#!/usr/bin/env python3
"""Generate segmented consensus + pileup from paired FASTQ files.

This script is a segmented-virus counterpart to pipeline.py. It accepts
multiple segment references (for example 8 Influenza A segments), maps reads
once against a combined reference, and outputs:
  - one consensus FASTA containing one sequence per segment
  - one pileup file containing depth for all segments
"""

import argparse
import subprocess
import sys
import time
from pathlib import Path


def run_command(cmd, cwd=None):
    print(f"[segmented_pipeline] Running: {' '.join(cmd)}")
    result = subprocess.run(
        cmd, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"Command failed: {' '.join(cmd)}\nstdout:\n{result.stdout}\nstderr:\n{result.stderr}"
        )
    return result


def run_piped_command(cmd1, cmd2):
    print(f"[segmented_pipeline] Running: {' '.join(cmd1)} | {' '.join(cmd2)}")
    p1 = subprocess.Popen(
        cmd1, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=False
    )
    p2 = subprocess.Popen(
        cmd2, stdin=p1.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=False
    )
    p1.stdout.close()
    _out2, err2 = p2.communicate()
    err1 = p1.stderr.read()
    p1.stderr.close()

    if p1.wait() != 0:
        raise RuntimeError(
            f"Command failed: {' '.join(cmd1)}\nstderr:\n{err1.decode(errors='replace')}"
        )
    if p2.returncode != 0:
        raise RuntimeError(
            f"Command failed: {' '.join(cmd2)}\nstderr:\n{err2.decode(errors='replace')}"
        )


def _parse_fasta(path: Path):
    """Yield (header, sequence) from FASTA."""
    header = None
    seq_chunks = []
    with path.open("r", encoding="utf-8") as fh:
        for raw_line in fh:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_chunks)
                header = line[1:].split()[0]
                seq_chunks = []
            else:
                seq_chunks.append(line)
    if header is not None:
        yield header, "".join(seq_chunks)


def _write_fasta(records, path: Path):
    with path.open("w", encoding="utf-8") as fh:
        for header, sequence in records:
            fh.write(f">{header}\n")
            for i in range(0, len(sequence), 80):
                fh.write(sequence[i : i + 80] + "\n")


def build_segment_reference(
    references, output_path: Path, segment_names=None
):
    """Create a combined segmented reference FASTA and return segment ids."""
    records = []
    segment_ids = []

    if segment_names is not None and len(segment_names) != len(references):
        raise ValueError("--segment-names must have same length as --references")

    for idx, ref in enumerate(references):
        ref_path = Path(ref)
        if not ref_path.exists():
            raise FileNotFoundError(f"Reference not found: {ref_path}")

        parsed = list(_parse_fasta(ref_path))
        if not parsed:
            raise ValueError(f"No FASTA sequence found in: {ref_path}")
        if len(parsed) > 1:
            raise ValueError(
                f"Reference must contain one segment sequence per file: {ref_path}"
            )

        _orig_header, sequence = parsed[0]
        segment_id = (
            segment_names[idx]
            if segment_names is not None
            else ref_path.stem
        )
        records.append((segment_id, sequence))
        segment_ids.append(segment_id)

    _write_fasta(records, output_path)
    return segment_ids


def build_bam(minimap2_path, samtools_path, reference, fastq1, fastq2, bam_path):
    print("[segmented_pipeline] Mapping reads and building BAM")
    print(f"[segmented_pipeline] Combined reference: {reference}")
    print(f"[segmented_pipeline] FASTQ1: {fastq1}")
    print(f"[segmented_pipeline] FASTQ2: {fastq2}")

    cmd1 = [minimap2_path, "-ax", "sr", reference, fastq1, fastq2]
    cmd2 = [samtools_path, "sort", "-o", str(bam_path)]

    run_piped_command(cmd1, cmd2)
    print(f"[segmented_pipeline] BAM written: {bam_path}")

    run_command([samtools_path, "index", str(bam_path)])
    print(f"[segmented_pipeline] BAM index written: {bam_path}.bai")


def _rewrite_segmented_consensus_headers(consensus_path: Path, taxon_id: str):
    records = list(_parse_fasta(consensus_path))
    rewritten = []
    for segment_id, seq in records:
        rewritten.append((f"{taxon_id}|{segment_id}", seq))
    _write_fasta(rewritten, consensus_path)


def consensus_and_pileup_segmented(
    reference,
    bam_path,
    taxon_id,
    sample_id,
    outdir,
    min_depth=10,
    samtools_path="samtools",
    bcftools_path="bcftools",
):
    print("[segmented_pipeline] Building segmented consensus and pileup")
    print(f"[segmented_pipeline] Reference FASTA: {reference}")
    print(f"[segmented_pipeline] BAM: {bam_path}")

    t_start = time.perf_counter()

    outdir_path = Path(outdir)
    pileup_path = outdir_path / f"pileup_{sample_id}_{taxon_id}.txt"
    mask_path = outdir_path / f"mask_{sample_id}_{taxon_id}.bed"
    vcf_path = outdir_path / f"calls_{sample_id}_{taxon_id}.vcf.gz"
    consensus_path = outdir_path / f"consensus_{sample_id}_{taxon_id}.fasta"

    print(f"[segmented_pipeline] Calculating depth to {pileup_path}...")
    with open(pileup_path, "w") as f_out:
        subprocess.run(
            [samtools_path, "depth", "-a", str(bam_path)], stdout=f_out, check=True
        )

    print("[segmented_pipeline] Generating BED mask for low-coverage regions...")
    with open(pileup_path, "r") as f_in, open(mask_path, "w") as f_out:
        for line in f_in:
            chrom, pos, depth = line.strip().split("\t")
            if int(depth) < min_depth:
                pos_int = int(pos)
                f_out.write(f"{chrom}\t{pos_int - 1}\t{pos_int}\n")

    print("[segmented_pipeline] Calling variants (mpileup -> call -> norm)...")
    cmd_mpileup = [bcftools_path, "mpileup", "-Ou", "-f", reference, str(bam_path)]
    cmd_call = [bcftools_path, "call", "-Ou", "-mv"]
    cmd_norm = [bcftools_path, "norm", "-f", reference, "-Oz", "-o", str(vcf_path)]

    p1 = subprocess.Popen(cmd_mpileup, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p2 = subprocess.Popen(
        cmd_call, stdin=p1.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    p3 = subprocess.Popen(
        cmd_norm, stdin=p2.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )

    p1.stdout.close()
    p2.stdout.close()
    _out3, err3 = p3.communicate()

    if p3.returncode != 0:
        raise RuntimeError(
            f"Variant calling failed. bcftools norm stderr:\n{err3.decode(errors='replace')}"
        )

    print("[segmented_pipeline] Indexing VCF...")
    subprocess.run([bcftools_path, "index", str(vcf_path)], check=True)

    print("[segmented_pipeline] Generating segmented consensus FASTA...")
    cmd_consensus = [
        bcftools_path,
        "consensus",
        "-f",
        reference,
        "-m",
        str(mask_path),
        "-o",
        str(consensus_path),
        str(vcf_path),
    ]
    subprocess.run(cmd_consensus, check=True)

    _rewrite_segmented_consensus_headers(consensus_path, taxon_id)

    t_end = time.perf_counter()
    print(f"[segmented_pipeline] Completed in {t_end - t_start:.1f}s")

    return consensus_path, pileup_path, vcf_path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate segmented consensus genomes and combined pileup from paired FASTQ.gz files."
    )
    parser.add_argument("--fastq1", required=True, help="Path to R1 FASTQ.gz")
    parser.add_argument("--fastq2", required=True, help="Path to R2 FASTQ.gz")
    parser.add_argument(
        "--references",
        nargs="+",
        required=True,
        help="Paths to segment reference FASTA files (one per segment)",
    )
    parser.add_argument(
        "--segment-names",
        nargs="+",
        default=None,
        help="Optional segment labels aligned with --references order",
    )
    parser.add_argument(
        "--taxon-id", required=True, help="Taxon ID to embed in consensus headers"
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
    return parser.parse_args()


def derive_sample_id(fastq_path):
    name = Path(fastq_path).name
    for suffix in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break
    return name.replace("#", "_")


class Tee:
    def __init__(self, *streams):
        self.streams = streams

    def write(self, data):
        for stream in self.streams:
            stream.write(data)
            stream.flush()

    def flush(self):
        for stream in self.streams:
            stream.flush()

    def isatty(self):
        return False


def main():
    args = parse_args()
    sample_id = derive_sample_id(args.fastq1)

    logs_dir = Path("logs")
    logs_dir.mkdir(parents=True, exist_ok=True)
    log_path = logs_dir / f"{sample_id}.log"
    log_handle = open(log_path, "w")
    sys.stdout = Tee(sys.__stdout__, log_handle)
    sys.stderr = Tee(sys.__stderr__, log_handle)

    print("[segmented_pipeline] Starting pipeline")
    print(f"[segmented_pipeline] Sample ID: {sample_id}")
    print(f"[segmented_pipeline] Log file: {log_path}")

    outdir = Path(args.outdir) / sample_id
    outdir.mkdir(parents=True, exist_ok=True)

    combined_ref = outdir / f"combined_reference_{sample_id}_{args.taxon_id}.fasta"
    segment_ids = build_segment_reference(
        references=args.references,
        output_path=combined_ref,
        segment_names=args.segment_names,
    )
    print(f"[segmented_pipeline] Segments: {', '.join(segment_ids)}")

    bam_path = outdir / f"aligned_{sample_id}_{args.taxon_id}.bam"

    t_map_start = time.perf_counter()
    if bam_path.exists() and args.skip_mapping:
        print(f"[segmented_pipeline] BAM already exists; resuming after mapping: {bam_path}")
        t_map_end = time.perf_counter()
        print(f"[segmented_pipeline] Mapping skipped in {t_map_end - t_map_start:.1f}s")
    elif not args.skip_mapping:
        build_bam(
            args.minimap2,
            args.samtools,
            str(combined_ref),
            args.fastq1,
            args.fastq2,
            bam_path,
        )
        t_map_end = time.perf_counter()
        print(f"[segmented_pipeline] Mapping completed in {t_map_end - t_map_start:.1f}s")
    else:
        print(f"[segmented_pipeline] Skipping mapping; assuming BAM exists: {bam_path}")
        t_map_end = time.perf_counter()
        print(f"[segmented_pipeline] Mapping skipped in {t_map_end - t_map_start:.1f}s")

    consensus_path, pileup_path, vcf_path = consensus_and_pileup_segmented(
        reference=str(combined_ref),
        bam_path=bam_path,
        taxon_id=args.taxon_id,
        sample_id=sample_id,
        outdir=outdir,
        min_depth=args.min_depth,
        samtools_path=args.samtools,
        bcftools_path=args.bcftools,
    )

    print("\n--- Segmented Pipeline Summary ---")
    print(f"Combined reference: {combined_ref}")
    print(f"Alignments (BAM): {bam_path}")
    print(f"Variant Calls (VCF): {vcf_path}")
    print(f"Depth Profile (all segments): {pileup_path}")
    print(f"Consensus FASTA (all segments): {consensus_path}")
    print("[segmented_pipeline] Done")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        sys.exit(1)
