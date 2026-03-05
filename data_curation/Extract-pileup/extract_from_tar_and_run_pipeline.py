#!/usr/bin/env python3
"""Extract selected BAM/BAI files from tar archives and run consensus pipeline.

Given metadata with columns:
  - SampleID
  - Reference_Taxon_ID

the script locates matching files inside one or more tar archives under:
  <tar_root>/viral_pipeline_run/results/<SampleID>/<Reference_Taxon_ID>/

It extracts BAM + index into `data/` and then calls `pipeline.consensus_and_pileup`
using references from:
  defaults/library_references/<Reference_Taxon_ID>_reference.fasta
"""

from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
import tarfile
from dataclasses import dataclass
from pathlib import Path
import re
import sys

import pandas as pd


PATH_REGEX = re.compile(
    r"(?:.*/)?viral_pipeline_run/results/(?P<sample>[^/]+)/(?P<taxon>[^/]+)/(?P<filename>[^/]+)$"
)


@dataclass(frozen=True)
class Target:
    sample_id: str
    taxon_id: str


@dataclass
class IndexedEntry:
    tar_path: Path
    bam_member: str | None = None
    bai_member: str | None = None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Extract selected BAM/BAI from tar and run consensus+pileup pipeline."
    )
    parser.add_argument(
        "--metadata",
        required=True,
        help="Path to metadata file containing SampleID and Reference_Taxon_ID",
    )
    parser.add_argument(
        "--tar",
        nargs="+",
        required=True,
        help="One or more tar archive paths to scan for BAM/BAI files",
    )
    parser.add_argument(
        "--data-dir",
        default="data",
        help="Destination directory for extracted BAM/BAI files (default: data)",
    )
    parser.add_argument(
        "--reference-dir",
        default="defaults/library_references",
        help="Directory containing <Reference_Taxon_ID>_reference.fasta (default: defaults/library_references)",
    )
    parser.add_argument(
        "--outdir",
        default="output",
        help="Output directory for consensus/pileup files (default: output)",
    )
    parser.add_argument(
        "--min-depth",
        type=int,
        default=10,
        help="Minimum depth required to call a consensus base (default: 10)",
    )
    parser.add_argument(
        "--samtools",
        default="samtools",
        help="Path to samtools binary (default: samtools)",
    )
    parser.add_argument(
        "--bcftools",
        default="bcftools",
        help="Path to bcftools binary (default: bcftools)",
    )
    parser.add_argument(
        "--overwrite-extracted",
        action="store_true",
        help="Overwrite existing extracted BAM/BAI in --data-dir",
    )
    parser.add_argument(
        "--skip-existing-output",
        action="store_true",
        help="Skip pipeline step when consensus output already exists",
    )
    parser.add_argument(
        "--jobs",
        type=int,
        default=1,
        help="Number of samples to process in parallel in batch mode (default: 1)",
    )
    parser.add_argument(
        "--mpileup-min-mapq",
        type=int,
        default=None,
        help="Optional bcftools mpileup minimum mapping quality (-q)",
    )
    parser.add_argument(
        "--mpileup-min-baseq",
        type=int,
        default=None,
        help="Optional bcftools mpileup minimum base quality (-Q)",
    )
    parser.add_argument(
        "--mpileup-max-depth",
        type=int,
        default=None,
        help="Optional bcftools mpileup max depth per input file (-d)",
    )
    parser.add_argument(
        "--pileup-only",
        action="store_true",
        help="Only generate pileup and mask BED; skip variant calling and consensus",
    )
    parser.add_argument(
        "--no-resume",
        action="store_true",
        help="Disable resume behavior for batch pipeline BAM mode",
    )
    return parser.parse_args()


def load_targets(metadata_path: Path) -> list[Target]:
    df = pd.read_csv(metadata_path, sep=None, engine="python")

    missing_cols = [c for c in ("SampleID", "Reference_Taxon_ID") if c not in df.columns]
    if missing_cols:
        raise ValueError(f"Metadata missing required columns: {', '.join(missing_cols)}")

    df = df[["SampleID", "Reference_Taxon_ID"]].dropna()
    df["SampleID"] = df["SampleID"].astype(str).str.strip()
    df["Reference_Taxon_ID"] = df["Reference_Taxon_ID"].astype(str).str.strip()
    df = df[(df["SampleID"] != "") & (df["Reference_Taxon_ID"] != "")]

    seen: set[tuple[str, str]] = set()
    targets: list[Target] = []
    for row in df.itertuples(index=False):
        key = (row.SampleID, row.Reference_Taxon_ID)
        if key in seen:
            continue
        seen.add(key)
        targets.append(Target(sample_id=row.SampleID, taxon_id=row.Reference_Taxon_ID))
    return targets


def classify_member(filename: str) -> str | None:
    lower = filename.lower()
    if lower.endswith(".bam") and not lower.endswith(".bam.bai"):
        return "bam"
    if lower.endswith(".bam.bai") or lower.endswith(".bai"):
        return "bai"
    return None


def index_tar_archives(tar_paths: list[Path]) -> dict[tuple[str, str], IndexedEntry]:
    index: dict[tuple[str, str], IndexedEntry] = {}

    for tar_path in tar_paths:
        if not tar_path.exists():
            print(f"[extract] Warning: tar file not found, skipping: {tar_path}")
            continue

        print(f"[extract] Indexing tar: {tar_path}")
        with tarfile.open(tar_path, "r:*") as tf:
            for member in tf.getmembers():
                if not member.isfile():
                    continue
                match = PATH_REGEX.match(member.name)
                if not match:
                    continue

                sample_id = match.group("sample")
                taxon_id = match.group("taxon")
                member_type = classify_member(match.group("filename"))
                if member_type is None:
                    continue

                key = (sample_id, taxon_id)
                entry = index.get(key)
                if entry is None:
                    entry = IndexedEntry(tar_path=tar_path)
                    index[key] = entry

                # Keep first match encountered (stable by tar and member order).
                if member_type == "bam" and entry.bam_member is None:
                    entry.bam_member = member.name
                elif member_type == "bai" and entry.bai_member is None:
                    entry.bai_member = member.name

    return index


def extract_member(tf: tarfile.TarFile, member_name: str, destination: Path) -> None:
    extracted = tf.extractfile(member_name)
    if extracted is None:
        raise RuntimeError(f"Failed to extract member: {member_name}")
    with extracted, destination.open("wb") as out_handle:
        shutil.copyfileobj(extracted, out_handle)


def ensure_extracted_files(
    target: Target,
    entry: IndexedEntry,
    data_dir: Path,
    overwrite: bool,
) -> tuple[Path, Path]:
    sample_dir = data_dir / target.sample_id / target.taxon_id
    sample_dir.mkdir(parents=True, exist_ok=True)

    bam_out = sample_dir / f"{target.sample_id}_{target.taxon_id}.bam"
    bai_out = sample_dir / f"{target.sample_id}_{target.taxon_id}.bam.bai"

    if bam_out.exists() and bai_out.exists() and not overwrite:
        print(f"[extract] Reusing existing BAM/BAI: {bam_out}")
        return bam_out, bai_out

    if entry.bam_member is None or entry.bai_member is None:
        raise RuntimeError(
            f"Incomplete archive entry for {target.sample_id}/{target.taxon_id}: "
            f"bam={entry.bam_member is not None}, bai={entry.bai_member is not None}"
        )

    print(
        f"[extract] Extracting {target.sample_id}/{target.taxon_id} "
        f"from {entry.tar_path.name}"
    )
    with tarfile.open(entry.tar_path, "r:*") as tf:
        extract_member(tf, entry.bam_member, bam_out)
        extract_member(tf, entry.bai_member, bai_out)

    return bam_out, bai_out


def _write_manifest(rows: list[dict], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["sample_id", "taxon_id", "bam_path", "reference_path"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)


def _run_batch_pipeline(
    script_dir: Path,
    manifest_path: Path,
    outdir: Path,
    min_depth: int,
    samtools: str,
    bcftools: str,
    jobs: int,
    mpileup_min_mapq: int | None,
    mpileup_min_baseq: int | None,
    mpileup_max_depth: int | None,
    pileup_only: bool,
    no_resume: bool,
) -> None:
    batch_path = script_dir / "batch_pipeline.py"
    cmd = [
        sys.executable,
        str(batch_path),
        "--bam-manifest",
        str(manifest_path),
        "--outdir",
        str(outdir),
        "--min-depth",
        str(min_depth),
        "--samtools",
        samtools,
        "--bcftools",
        bcftools,
        "--jobs",
        str(jobs),
    ]
    if mpileup_min_mapq is not None:
        cmd.extend(["--mpileup-min-mapq", str(mpileup_min_mapq)])
    if mpileup_min_baseq is not None:
        cmd.extend(["--mpileup-min-baseq", str(mpileup_min_baseq)])
    if mpileup_max_depth is not None:
        cmd.extend(["--mpileup-max-depth", str(mpileup_max_depth)])
    if pileup_only:
        cmd.append("--pileup-only")
    if no_resume:
        cmd.append("--no-resume")
    print(f"[extract] Running batch pipeline: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)


def _rewrite_pileup_first_col_to_sample(pileup_path: Path, sample_id: str) -> None:
    if not pileup_path.exists():
        return
    tmp_path = pileup_path.with_suffix(pileup_path.suffix + ".tmp")
    with pileup_path.open("r", encoding="utf-8") as src, tmp_path.open("w", encoding="utf-8") as dst:
        for raw_line in src:
            line = raw_line.rstrip("\n")
            if not line:
                dst.write(raw_line)
                continue
            parts = line.split("\t")
            if len(parts) >= 3:
                parts[0] = sample_id
                dst.write("\t".join(parts) + "\n")
            else:
                dst.write(raw_line)
    tmp_path.replace(pileup_path)


def main() -> int:
    args = parse_args()

    metadata_path = Path(args.metadata)
    tar_paths = [Path(p) for p in args.tar]
    data_dir = Path(args.data_dir)
    script_dir = Path(__file__).resolve().parent
    reference_dir = Path(args.reference_dir)
    if not reference_dir.is_absolute():
        reference_dir = script_dir / reference_dir
    outdir = Path(args.outdir)

    if not metadata_path.exists():
        print(f"[extract] Error: metadata file not found: {metadata_path}")
        return 2

    targets = load_targets(metadata_path)
    if not targets:
        print("[extract] No valid (SampleID, Reference_Taxon_ID) pairs found in metadata.")
        return 1

    print(f"[extract] Loaded {len(targets)} unique targets from metadata.")
    index = index_tar_archives(tar_paths)
    if not index:
        print("[extract] No matching BAM/BAI members found in provided tar files.")
        return 1

    processed = 0
    missing = 0
    errors = 0
    skipped_missing_reference = 0
    manifest_rows: list[dict] = []
    rewrite_targets: list[tuple[Path, str]] = []

    for target in targets:
        key = (target.sample_id, target.taxon_id)
        entry = index.get(key)
        if entry is None:
            print(f"[extract] Missing in tar: {target.sample_id}/{target.taxon_id}")
            missing += 1
            continue
        try:
            bam_path, _bai_path = ensure_extracted_files(
                target=target,
                entry=entry,
                data_dir=data_dir,
                overwrite=args.overwrite_extracted,
            )
            reference_path = reference_dir / f"{target.taxon_id}_reference.fasta"
            if not reference_path.exists():
                print(
                    f"[extract] Warning: reference not found for {target.taxon_id}: {reference_path}; skipping"
                )
                skipped_missing_reference += 1
                continue

            sample_outdir = outdir / target.sample_id
            consensus_path = sample_outdir / f"consensus_{target.sample_id}_{target.taxon_id}.fasta"
            pileup_path = sample_outdir / f"pileup_{target.sample_id}_{target.taxon_id}.txt"
            existing_target = pileup_path if args.pileup_only else consensus_path
            if args.skip_existing_output and existing_target.exists():
                print(f"[extract] Skipping existing output for {target.sample_id}/{target.taxon_id}")
                processed += 1
                rewrite_targets.append(
                    (pileup_path, target.sample_id)
                )
                continue

            manifest_rows.append(
                {
                    "sample_id": target.sample_id,
                    "taxon_id": target.taxon_id,
                    "bam_path": str(bam_path.resolve()),
                    "reference_path": str(reference_path.resolve()),
                }
            )
            rewrite_targets.append(
                (pileup_path, target.sample_id)
            )
            processed += 1
        except Exception as exc:
            errors += 1
            print(f"[extract] Error for {target.sample_id}/{target.taxon_id}: {exc}")

    if manifest_rows:
        try:
            manifest_path = outdir / "batch_bam_manifest.tsv"
            _write_manifest(manifest_rows, manifest_path)
            _run_batch_pipeline(
                script_dir=script_dir,
                manifest_path=manifest_path,
                outdir=outdir,
                min_depth=args.min_depth,
                samtools=args.samtools,
                bcftools=args.bcftools,
                jobs=args.jobs,
                mpileup_min_mapq=args.mpileup_min_mapq,
                mpileup_min_baseq=args.mpileup_min_baseq,
                mpileup_max_depth=args.mpileup_max_depth,
                pileup_only=args.pileup_only,
                no_resume=args.no_resume,
            )
        except Exception as exc:
            errors += 1
            print(f"[extract] Error while running batch pipeline: {exc}")

    for pileup_path, sample_id in rewrite_targets:
        try:
            _rewrite_pileup_first_col_to_sample(pileup_path, sample_id)
        except Exception as exc:
            errors += 1
            print(f"[extract] Error rewriting pileup header for {sample_id}: {exc}")

    print(
        f"[extract] Done. processed={processed} "
        f"missing_in_tar={missing} "
        f"skipped_missing_reference={skipped_missing_reference} "
        f"errors={errors} total={len(targets)}"
    )
    return 0 if errors == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
