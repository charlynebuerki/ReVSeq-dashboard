#!/usr/bin/env python3
import argparse
import csv
from concurrent.futures import ThreadPoolExecutor, as_completed
import re
import subprocess
import sys
import threading
from pathlib import Path

import pipeline as single_pipeline


READ_REGEX = re.compile(r"^(?P<base>.+?)(?:[_-](?:R)?)(?P<read>[12])$")


def find_pairs(fastq_dir):
    fastq_dir = Path(fastq_dir)
    fastq_files = sorted(fastq_dir.rglob("*.fastq.gz")) + sorted(
        fastq_dir.rglob("*.fq.gz")
    )

    pairs = {}
    skipped = []

    for fq in fastq_files:
        stem = fq.name
        for suffix in (".fastq.gz", ".fq.gz"):
            if stem.endswith(suffix):
                stem = stem[: -len(suffix)]
                break

        match = READ_REGEX.match(stem)
        if not match:
            skipped.append(fq)
            continue

        base = match.group("base")
        read = match.group("read")
        pairs.setdefault(base, {})[read] = fq

    return pairs, skipped


def derive_sample_id(fastq_path):
    name = Path(fastq_path).name
    for suffix in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break
    return name.replace("#", "_")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Batch runner for the viral consensus pipeline."
    )
    parser.add_argument(
        "--fastq-dir",
        default=None,
        help="Directory containing paired FASTQ.gz files",
    )
    parser.add_argument(
        "--bam-manifest",
        default=None,
        help=(
            "Optional TSV/CSV manifest for BAM-based batch mode with columns: "
            "sample_id, taxon_id, bam_path, reference_path"
        ),
    )
    parser.add_argument(
        "--reference",
        default=None,
        help="Path to reference FASTA (required for non-segmented mode)",
    )
    parser.add_argument(
        "--references",
        nargs="+",
        default=None,
        help="Segment reference FASTA paths (required for --segmented mode)",
    )
    parser.add_argument(
        "--segment-names",
        nargs="+",
        default=None,
        help="Optional segment labels aligned with --references order",
    )
    parser.add_argument(
        "--segmented",
        action="store_true",
        help="Use segmented_pipeline.py instead of pipeline.py",
    )
    parser.add_argument(
        "--taxon-id",
        required=False,
        default=None,
        help="Taxon ID to embed in headers (required in FASTQ mode)",
    )
    parser.add_argument("--outdir", default="output", help="Output directory")
    parser.add_argument(
        "--min-depth",
        type=int,
        default=10,
        help="Minimum depth required to call a consensus base (default: 10)",
    )
    parser.add_argument(
        "--minimap2",
        default="minimap2",
        help="Path to minimap2 binary (default: minimap2)",
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
        "--skip-mapping",
        action="store_true",
        help="Skip mapping step and use existing BAM if present",
    )
    parser.add_argument(
        "--progress-file",
        default=None,
        help="Path to resume progress file (default: <outdir>/batch_progress.txt)",
    )
    parser.add_argument(
        "--no-resume",
        action="store_true",
        help="Disable resume behavior and process all pairs",
    )
    parser.add_argument(
        "--delete-fastq",
        action="store_true",
        help="Delete FASTQ files after successful processing",
    )
    parser.add_argument(
        "--jobs",
        type=int,
        default=1,
        help="Number of samples to process in parallel (default: 1)",
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
    return parser.parse_args()


def main():
    args = parse_args()
    if args.jobs < 1:
        print("[batch] --jobs must be >= 1")
        return 2

    if args.bam_manifest:
        return run_bam_manifest_mode(args)

    if args.segmented:
        if not args.references:
            print("[batch] --references is required when --segmented is set")
            sys.exit(2)
        if args.segment_names and len(args.segment_names) != len(args.references):
            print("[batch] --segment-names must have same length as --references")
            sys.exit(2)
    else:
        if not args.reference:
            print("[batch] --reference is required when --segmented is not set")
            sys.exit(2)
    if not args.taxon_id:
        print("[batch] --taxon-id is required in FASTQ mode")
        sys.exit(2)
    if not args.fastq_dir:
        print("[batch] --fastq-dir is required in FASTQ mode")
        sys.exit(2)

    pairs, skipped = find_pairs(args.fastq_dir)

    if skipped:
        print("[batch] Skipped files with no read indicator:")
        for fq in skipped:
            print(f"[batch]   {fq}")

    to_process = []
    for base, reads in sorted(pairs.items()):
        if "1" in reads and "2" in reads:
            to_process.append((base, reads["1"], reads["2"]))
        else:
            missing = "2" if "1" in reads else "1"
            print(f"[batch] Missing read {missing} for base: {base}")

    if not to_process:
        print("[batch] No complete pairs found.")
        sys.exit(1)

    script_name = "segmented_pipeline.py" if args.segmented else "pipeline.py"
    pipeline_path = Path(__file__).parent / script_name

    progress_path = (
        Path(args.progress_file)
        if args.progress_file
        else Path(args.outdir) / "batch_progress.txt"
    )
    progress_path.parent.mkdir(parents=True, exist_ok=True)

    completed = set()
    if not args.no_resume and progress_path.exists():
        with progress_path.open("r") as handle:
            completed = {line.strip() for line in handle if line.strip()}
        if completed:
            print(f"[batch] Resume enabled; {len(completed)} samples already complete")

    progress_lock = threading.Lock()
    failures = []

    def process_fastq_pair(base, fq1, fq2):
        sample_id = derive_sample_id(fq1)
        if not args.no_resume and sample_id in completed:
            print(f"[batch] Skipping completed sample: {sample_id}")
            return True

        print(f"[batch] Processing pair: {fq1.name} / {fq2.name}")
        cmd = [
            sys.executable,
            str(pipeline_path),
            "--fastq1",
            str(fq1),
            "--fastq2",
            str(fq2),
            "--taxon-id",
            args.taxon_id,
            "--outdir",
            args.outdir,
            "--min-depth",
            str(args.min_depth),
            "--minimap2",
            args.minimap2,
            "--samtools",
            args.samtools,
            "--bcftools",
            args.bcftools,
        ]
        if args.segmented:
            cmd.extend(["--references", *args.references])
            if args.segment_names:
                cmd.extend(["--segment-names", *args.segment_names])
        else:
            cmd.extend(["--reference", args.reference])
            if args.mpileup_min_mapq is not None:
                cmd.extend(["--mpileup-min-mapq", str(args.mpileup_min_mapq)])
            if args.mpileup_min_baseq is not None:
                cmd.extend(["--mpileup-min-baseq", str(args.mpileup_min_baseq)])
            if args.mpileup_max_depth is not None:
                cmd.extend(["--mpileup-max-depth", str(args.mpileup_max_depth)])
        if args.skip_mapping:
            cmd.append("--skip-mapping")
        if args.pileup_only:
            cmd.append("--pileup-only")

        result = subprocess.run(cmd)
        if result.returncode != 0:
            print(f"[batch] Pipeline failed for pair: {fq1.name}")
            return False

        if not args.no_resume:
            with progress_lock:
                with progress_path.open("a") as handle:
                    handle.write(sample_id + "\n")
                completed.add(sample_id)

        if args.delete_fastq:
            try:
                fq1.unlink()
                fq2.unlink()
                print(f"[batch] Deleted FASTQ files for sample: {sample_id}")
            except Exception as exc:
                print(f"[batch] Warning: failed to delete FASTQ files: {exc}")
        return True

    if args.jobs == 1:
        for base, fq1, fq2 in to_process:
            ok = process_fastq_pair(base, fq1, fq2)
            if not ok:
                failures.append((base, fq1, fq2))
                break
    else:
        with ThreadPoolExecutor(max_workers=args.jobs) as executor:
            future_map = {
                executor.submit(process_fastq_pair, base, fq1, fq2): (base, fq1, fq2)
                for base, fq1, fq2 in to_process
            }
            for future in as_completed(future_map):
                base, fq1, fq2 = future_map[future]
                try:
                    ok = future.result()
                except Exception as exc:
                    print(f"[batch] Pipeline failed for pair: {fq1.name}: {exc}")
                    ok = False
                if not ok:
                    failures.append((base, fq1, fq2))

    if failures:
        print(f"[batch] {len(failures)} sample(s) failed.")
        return 1

    print("[batch] All samples completed.")
    return 0


def _read_manifest_rows(path: Path):
    text = path.read_text(encoding="utf-8")
    if not text.strip():
        return []

    # Prefer deterministic delimiter by file suffix, then content sniffing fallback.
    delimiters = []
    suffix = path.suffix.lower()
    if suffix in {".tsv", ".tab"}:
        delimiters = ["\t", ","]
    elif suffix == ".csv":
        delimiters = [",", "\t"]
    else:
        delimiters = ["\t", ","]

    # If one delimiter clearly dominates in header, prioritize it.
    header = text.splitlines()[0] if text.splitlines() else ""
    if header.count("\t") > header.count(","):
        delimiters = ["\t", ","]
    elif header.count(",") > header.count("\t"):
        delimiters = [",", "\t"]

    last_error = None
    for delimiter in delimiters:
        try:
            with path.open("r", encoding="utf-8", newline="") as handle:
                reader = csv.DictReader(handle, delimiter=delimiter)
                rows = list(reader)
            if rows:
                return rows
        except Exception as exc:  # pragma: no cover
            last_error = exc

    if last_error:
        raise last_error
    return []


def run_bam_manifest_mode(args):
    manifest_path = Path(args.bam_manifest)
    if not manifest_path.exists():
        print(f"[batch] Manifest not found: {manifest_path}")
        return 2

    rows = _read_manifest_rows(manifest_path)
    required = {"sample_id", "taxon_id", "bam_path", "reference_path"}
    if not rows:
        print(f"[batch] Empty manifest: {manifest_path}")
        return 1
    if not required.issubset(set(rows[0].keys())):
        print(
            "[batch] Manifest missing required columns: "
            f"{', '.join(sorted(required - set(rows[0].keys())))}"
        )
        return 2

    progress_path = (
        Path(args.progress_file)
        if args.progress_file
        else Path(args.outdir) / "batch_progress_bam.txt"
    )
    progress_path.parent.mkdir(parents=True, exist_ok=True)

    completed = set()
    if not args.no_resume and progress_path.exists():
        with progress_path.open("r", encoding="utf-8") as handle:
            completed = {line.strip() for line in handle if line.strip()}
        if completed:
            print(f"[batch] Resume enabled; {len(completed)} BAM targets already complete")

    progress_lock = threading.Lock()
    failures = []

    def process_row(row):
        sample_id = str(row["sample_id"]).strip()
        taxon_id = str(row["taxon_id"]).strip()
        bam_path = Path(str(row["bam_path"]).strip())
        reference_path = Path(str(row["reference_path"]).strip())
        key = f"{sample_id}|{taxon_id}"

        if not args.no_resume and key in completed:
            print(f"[batch] Skipping completed BAM target: {key}")
            return True
        if not bam_path.exists():
            print(f"[batch] Missing BAM: {bam_path}")
            return False
        if not reference_path.exists():
            print(f"[batch] Missing reference: {reference_path}")
            return False

        sample_outdir = Path(args.outdir) / sample_id
        sample_outdir.mkdir(parents=True, exist_ok=True)
        print(f"[batch] Processing BAM target: {key}")
        single_pipeline.consensus_and_pileup(
            reference=str(reference_path),
            bam_path=bam_path,
            taxon_id=taxon_id,
            sample_id=sample_id,
            outdir=sample_outdir,
            min_depth=args.min_depth,
            samtools_path=args.samtools,
            bcftools_path=args.bcftools,
            mpileup_min_mapq=args.mpileup_min_mapq,
            mpileup_min_baseq=args.mpileup_min_baseq,
            mpileup_max_depth=args.mpileup_max_depth,
            pileup_only=args.pileup_only,
        )

        if not args.no_resume:
            with progress_lock:
                with progress_path.open("a", encoding="utf-8") as handle:
                    handle.write(key + "\n")
                completed.add(key)
        return True

    if args.jobs == 1:
        for row in rows:
            ok = process_row(row)
            if not ok:
                failures.append(row)
                break
    else:
        with ThreadPoolExecutor(max_workers=args.jobs) as executor:
            future_map = {executor.submit(process_row, row): row for row in rows}
            for future in as_completed(future_map):
                row = future_map[future]
                try:
                    ok = future.result()
                except Exception as exc:
                    print(
                        "[batch] Pipeline failed for BAM target "
                        f"{row.get('sample_id')}|{row.get('taxon_id')}: {exc}"
                    )
                    ok = False
                if not ok:
                    failures.append(row)

    if failures:
        print(f"[batch] {len(failures)} BAM manifest target(s) failed.")
        return 1

    print("[batch] All BAM manifest targets completed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
