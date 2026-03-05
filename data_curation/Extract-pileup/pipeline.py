#!/usr/bin/env python3
import argparse
import os
import subprocess
import sys
import time
from pathlib import Path


def run_command(cmd, cwd=None):
    print(f"[pipeline] Running: {' '.join(cmd)}")
    result = subprocess.run(
        cmd, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"Command failed: {' '.join(cmd)}\nstdout:\n{result.stdout}\nstderr:\n{result.stderr}"
        )
    return result


def run_piped_command(cmd1, cmd2):
    print(f"[pipeline] Running: {' '.join(cmd1)} | {' '.join(cmd2)}")
    p1 = subprocess.Popen(
        cmd1, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=False
    )
    p2 = subprocess.Popen(
        cmd2,
        stdin=p1.stdout,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=False,
    )
    p1.stdout.close()
    out2, err2 = p2.communicate()
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


def build_bam(minimap2_path, samtools_path, reference, fastq1, fastq2, bam_path):
    print("[pipeline] Mapping reads and building BAM")
    print(f"[pipeline] Reference: {reference}")
    print(f"[pipeline] FASTQ1: {fastq1}")
    print(f"[pipeline] FASTQ2: {fastq2}")

    # -ax sr configures minimap2 specifically for short-read alignment
    cmd1 = [minimap2_path, "-ax", "sr", reference, fastq1, fastq2]
    cmd2 = [samtools_path, "sort", "-o", str(bam_path)]

    run_piped_command(cmd1, cmd2)
    print(f"[pipeline] BAM written: {bam_path}")

    run_command([samtools_path, "index", str(bam_path)])
    print(f"[pipeline] BAM index written: {bam_path}.bai")


def consensus_and_pileup(
    reference,
    bam_path,
    taxon_id,
    sample_id,
    outdir,
    min_depth=1,
    samtools_path="samtools",
    bcftools_path="bcftools",
    mpileup_min_mapq=None,
    mpileup_min_baseq=None,
    mpileup_max_depth=None,
    pileup_only=False,
):
    print("[pipeline] Building consensus and pileup using bcftools")
    print(f"[pipeline] Reference FASTA: {reference}")
    print(f"[pipeline] BAM: {bam_path}")

    t_start = time.perf_counter()

    outdir_path = Path(outdir)
    pileup_path = outdir_path / f"pileup_{sample_id}_{taxon_id}.txt"
    mask_path = outdir_path / f"mask_{sample_id}_{taxon_id}.bed"
    vcf_path = outdir_path / f"calls_{sample_id}_{taxon_id}.vcf.gz"
    consensus_path = outdir_path / f"consensus_{sample_id}_{taxon_id}.fasta"

    # 1. Generate Pileup (Depth) and Low-Coverage Mask BED
    print(f"[pipeline] Calculating depth to {pileup_path}...")
    with open(pileup_path, "w") as f_out:
        subprocess.run(
            [samtools_path, "depth", "-a", str(bam_path)], stdout=f_out, check=True
        )

    print("[pipeline] Generating BED mask for low-coverage regions...")
    with open(pileup_path, "r") as f_in, open(mask_path, "w") as f_out:
        for line in f_in:
            chrom, pos, depth = line.strip().split("\t")
            if int(depth) < min_depth:
                pos_int = int(pos)
                # BED format is 0-indexed, half-open
                f_out.write(f"{chrom}\t{pos_int - 1}\t{pos_int}\n")

    if pileup_only:
        t_end = time.perf_counter()
        print(f"[pipeline] Pileup-only mode completed in {t_end - t_start:.1f}s")
        return None, pileup_path, None

    # 2. Call Variants probabilistically using bcftools
    print("[pipeline] Calling variants (mpileup -> call -> norm)...")
    cmd_mpileup = [bcftools_path, "mpileup", "-Ou", "-f", reference]
    if mpileup_min_mapq is not None:
        cmd_mpileup.extend(["-q", str(mpileup_min_mapq)])
    if mpileup_min_baseq is not None:
        cmd_mpileup.extend(["-Q", str(mpileup_min_baseq)])
    if mpileup_max_depth is not None:
        cmd_mpileup.extend(["-d", str(mpileup_max_depth)])
    cmd_mpileup.append(str(bam_path))
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
    _, err3 = p3.communicate()

    if p3.returncode != 0:
        raise RuntimeError(
            f"Variant calling failed. bcftools norm stderr:\n{err3.decode(errors='replace')}"
        )

    # 3. Index the VCF
    print("[pipeline] Indexing VCF...")
    subprocess.run([bcftools_path, "index", str(vcf_path)], check=True)

    # 4. Generate the Final Consensus
    print("[pipeline] Generating masked consensus FASTA...")
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

    # Replace all FASTA headers in the output with the taxon_id
    with open(consensus_path, "r") as f:
        lines = f.readlines()
    with open(consensus_path, "w") as f:
        for line in lines:
            if line.startswith(">"):
                f.write(f">{taxon_id}\n")
            else:
                f.write(line)

    t_end = time.perf_counter()
    print(f"[pipeline] Consensus/pileup completed in {t_end - t_start:.1f}s")

    return consensus_path, pileup_path, vcf_path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate a robust viral consensus genome and depth map from paired FASTQ.gz files."
    )
    parser.add_argument("--fastq1", required=True, help="Path to R1 FASTQ.gz")
    parser.add_argument("--fastq2", required=True, help="Path to R2 FASTQ.gz")
    parser.add_argument("--reference", required=True, help="Path to reference FASTA")
    parser.add_argument(
        "--taxon-id", required=True, help="Taxon ID to embed in headers and filenames"
    )
    parser.add_argument("--outdir", default="output", help="Output directory")
    parser.add_argument(
        "--min-depth",
        type=int,
        default=10,
        help="Minimum depth required to call a consensus base (default: 1)",
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

    print("[pipeline] Starting pipeline")
    print(f"[pipeline] Sample ID: {sample_id}")
    print(f"[pipeline] Log file: {log_path}")

    outdir = Path(args.outdir) / sample_id
    outdir.mkdir(parents=True, exist_ok=True)

    bam_path = outdir / f"aligned_{sample_id}_{args.taxon_id}.bam"

    t_map_start = time.perf_counter()
    if bam_path.exists() and args.skip_mapping:
        print(f"[pipeline] BAM already exists; resuming after mapping: {bam_path}")
        t_map_end = time.perf_counter()
        print(f"[pipeline] Mapping skipped in {t_map_end - t_map_start:.1f}s")
    elif not args.skip_mapping:
        build_bam(
            args.minimap2,
            args.samtools,
            args.reference,
            args.fastq1,
            args.fastq2,
            bam_path,
        )
        t_map_end = time.perf_counter()
        print(f"[pipeline] Mapping completed in {t_map_end - t_map_start:.1f}s")
    else:
        print(f"[pipeline] Skipping mapping; assuming BAM exists: {bam_path}")
        t_map_end = time.perf_counter()
        print(f"[pipeline] Mapping skipped in {t_map_end - t_map_start:.1f}s")

    consensus_path, pileup_path, vcf_path = consensus_and_pileup(
        reference=args.reference,
        bam_path=bam_path,
        taxon_id=args.taxon_id,
        sample_id=sample_id,
        outdir=outdir,
        min_depth=args.min_depth,
        samtools_path=args.samtools,
        bcftools_path=args.bcftools,
        mpileup_min_mapq=args.mpileup_min_mapq,
        mpileup_min_baseq=args.mpileup_min_baseq,
        mpileup_max_depth=args.mpileup_max_depth,
        pileup_only=args.pileup_only,
    )

    print("\n--- Pipeline Summary ---")
    print(f"Alignments (BAM): {bam_path}")
    print(f"Depth Profile: {pileup_path}")
    if args.pileup_only:
        print("Variant Calls (VCF): skipped (pileup-only)")
        print("Consensus FASTA: skipped (pileup-only)")
    else:
        print(f"Variant Calls (VCF): {vcf_path}")
        print(f"Consensus FASTA: {consensus_path}")
    print("[pipeline] Done")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        sys.exit(1)
