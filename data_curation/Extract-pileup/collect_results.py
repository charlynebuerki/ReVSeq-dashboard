#!/usr/bin/env python3
import argparse
import shutil
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Collect consensus and pileup files from pipeline output into a single directory."
    )
    parser.add_argument(
        "--outdir",
        required=True,
        help="Pipeline output directory (e.g., output)",
    )
    parser.add_argument(
        "--dest",
        default="processed_pileup",
        help="Destination directory for collected files (default: processed_pileup)",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    outdir = Path(args.outdir)
    dest = Path(args.dest)

    if not outdir.exists():
        print(f"Error: output directory does not exist: {outdir}")
        return 1

    dest.mkdir(parents=True, exist_ok=True)

    consensus_files = list(outdir.rglob("consensus_*.fasta"))
    pileup_files = list(outdir.rglob("pileup_*.txt"))

    if not consensus_files and not pileup_files:
        print(f"No consensus or pileup files found in {outdir}")
        return 1

    print(f"[collect] Copying {len(consensus_files)} consensus file(s)...")
    for fasta in consensus_files:
        dest_path = dest / fasta.name
        shutil.copy2(fasta, dest_path)
        print(f"[collect]   {fasta.name}")

    print(f"[collect] Copying {len(pileup_files)} pileup file(s)...")
    for pileup in pileup_files:
        dest_path = dest / pileup.name
        shutil.copy2(pileup, dest_path)
        print(f"[collect]   {pileup.name}")

    total = len(consensus_files) + len(pileup_files)
    print(f"[collect] Collected {total} file(s) to {dest}")

    return 0


if __name__ == "__main__":
    import sys

    sys.exit(main())
