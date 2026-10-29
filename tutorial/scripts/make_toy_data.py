#!/usr/bin/env python3
"""Simulate a small, fully synthetic sequencing run for the ReVSeq tutorial.

For each toy sample we mutate a public reference genome (clade + private SNPs),
simulate paired-end reads with uneven coverage, and map them with minimap2.

Outputs (in --outdir):
    references/<virus>.fasta      reference used for mapping (one per virus)
    bams/<sample>_<virus>.bam     sorted + indexed BAM, one per sample and virus
    manifest.tsv                  BAM manifest for data_curation/Extract-pileup/batch_pipeline.py
    samples.tsv                   sample metadata (date, location, sequencing + PCR results)

Usage:
    python tutorial/scripts/make_toy_data.py --outdir tutorial/workdir/toy_data
"""

from __future__ import annotations

import argparse
import csv
import gzip
import random
import subprocess
from datetime import date, timedelta
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
LIBRARY = REPO / "data_curation/Extract-pileup/defaults/library_references"
FLU_SEGMENTS = ["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"]

# virus key (as used by the Nextstrain pipeline) -> reference, taxon and reported names
VIRUSES = {
    "rsv-a": {"taxon_id": "208893", "refs": ["208893_reference.fasta"], "seq": "rsv a", "pcr": "RSV - A/B"},
    "rsv-b": {"taxon_id": "208895", "refs": ["208895 _reference.fasta"], "seq": "rsv b", "pcr": "RSV - A/B"},
    "flu-a-h3n2": {"taxon_id": "41857", "refs": [f"h3n2/{s}.fa" for s in FLU_SEGMENTS], "seq": "influenza a h3n2", "pcr": "Influenza A"},
    "hpiv-3": {"taxon_id": "11216", "refs": ["11216_reference.fasta"], "seq": "parainfluenza 3", "pcr": "Parainfluenza 3"},
}
N_SAMPLES = {"rsv-a": 6, "rsv-b": 4, "flu-a-h3n2": 5, "hpiv-3": 4}
# Location spellings are deliberately mixed to show location harmonization.
LOCATIONS = ["ZH", "Zürich", "BS", "Basel-Stadt", "GE", "Genève", "BE", "VD", "TI", "Luzern"]
READ_LEN, FRAG_MEAN, ERROR_RATE = 150, 320, 0.002
# A shared zero-coverage region in the G gene of the first RSV-A samples, as a failed amplicon would
# leave. Part 2 of the tutorial shows how to find it in the dashboard.
DROPOUT = {"virus": "rsv-a", "n_samples": 3, "start": 4900, "end": 5250}


def read_fasta(path: Path) -> list[tuple[str, str]]:
    records, name, seq = [], None, []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if name:
                records.append((name, "".join(seq).upper()))
            name, seq = line[1:].split()[0], []
        elif line.strip():
            seq.append(line.strip())
    records.append((name, "".join(seq).upper()))
    return records


def mutate(seq: str, positions: list[int], rng: random.Random) -> str:
    seq = list(seq)
    for pos in positions:
        if seq[pos] in "ACGT":
            seq[pos] = rng.choice([b for b in "ACGT" if b != seq[pos]])
    return "".join(seq)


def revcomp(seq: str) -> str:
    return seq.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def add_errors(read: str, rng: random.Random) -> str:
    # At most one substitution per read (~ERROR_RATE per base): fast and good enough for a toy run.
    if rng.random() >= ERROR_RATE * len(read):
        return read
    i = rng.randrange(len(read))
    return read[:i] + rng.choice("ACGT") + read[i + 1:]


def simulate_reads(records, mean_depth, rng, fq1, fq2, prefix, dropout=None):
    """Write read pairs; per-300bp coverage weights mimic amplicon unevenness and dropouts.

    dropout: optional (start, end) interval that receives no reads at all.
    """
    n = 0
    for name, seq in records:
        n_bins = max(1, len(seq) // 300)
        weights = [rng.lognormvariate(0, 0.5) for _ in range(n_bins)]
        if rng.random() < 0.4:  # some samples have a poorly covered region
            start = rng.randrange(n_bins)
            for i in range(start, min(n_bins, start + rng.randint(1, 3))):
                weights[i] = 0.003
        n_pairs = int(mean_depth * len(seq) / (2 * READ_LEN))
        for bin_idx in rng.choices(range(n_bins), weights=weights, k=n_pairs):
            frag = max(READ_LEN, int(rng.gauss(FRAG_MEAN, 30)))
            start = min(max(0, bin_idx * 300 + rng.randrange(300) - frag // 2), max(0, len(seq) - frag))
            if dropout and start < dropout[1] and start + frag > dropout[0]:
                continue
            fragment = seq[start:start + frag]
            r1, r2 = fragment[:READ_LEN], revcomp(fragment)[:READ_LEN]
            qual = "I" * READ_LEN
            fq1.write(f"@{prefix}_{n}/1\n{add_errors(r1, rng)}\n+\n{qual[:len(r1)]}\n")
            fq2.write(f"@{prefix}_{n}/2\n{add_errors(r2, rng)}\n+\n{qual[:len(r2)]}\n")
            n += 1


def map_reads(reference: Path, fq1: Path, fq2: Path, bam: Path, threads: int):
    minimap = subprocess.Popen(
        ["minimap2", "-ax", "sr", "-t", str(threads), str(reference), str(fq1), str(fq2)],
        stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
    )
    subprocess.run(["samtools", "sort", "-@", str(threads), "-o", str(bam), "-"],
                   stdin=minimap.stdout, stderr=subprocess.DEVNULL, check=True)
    minimap.wait()
    subprocess.run(["samtools", "index", str(bam)], check=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--outdir", default="tutorial/workdir/toy_data")
    parser.add_argument("--seed", type=int, default=2025)
    parser.add_argument("--threads", type=int, default=2)
    args = parser.parse_args()
    rng = random.Random(args.seed)

    outdir = Path(args.outdir)
    (outdir / "references").mkdir(parents=True, exist_ok=True)
    (outdir / "bams").mkdir(exist_ok=True)

    # One reference FASTA per virus (flu: all 8 segments in one file).
    references = {}
    for virus, info in VIRUSES.items():
        records = [rec for ref in info["refs"] for rec in read_fasta(LIBRARY / ref)]
        references[virus] = records
        path = outdir / "references" / f"{virus}.fasta"
        path.write_text("".join(f">{name}\n{seq}\n" for name, seq in records))
        subprocess.run(["samtools", "faidx", str(path)], check=True)

    # Sample sheet: one row per (sample, virus); the last sample is an RSV-A + HPIV-3 co-infection.
    samples, start = [], date(2024, 10, 7)
    for virus, n in N_SAMPLES.items():
        for _ in range(n):
            samples.append({"virus": virus, "date": start + timedelta(days=rng.randrange(170))})
    samples.sort(key=lambda s: s["date"])
    for i, sample in enumerate(samples, start=1):
        sample.update(SampleID=f"TOY{i:03d}", location=rng.choice(LOCATIONS))
    coinfection = {"SampleID": f"TOY{len(samples) + 1:03d}", "date": date(2025, 1, 20), "location": "Kanton Bern"}
    samples += [{**coinfection, "virus": "rsv-a"}, {**coinfection, "virus": "hpiv-3", "pcr_missed": True}]

    # Two "clades" per virus so the phylogenies show some structure.
    clades = {
        virus: [[(r, rng.randrange(len(seq))) for r, (_, seq) in enumerate(recs) for _ in range(max(1, len(seq) // 2000))]
                for _ in range(2)]
        for virus, recs in references.items()
    }

    dropout_samples = [s["SampleID"] for s in samples if s["virus"] == DROPOUT["virus"]][:DROPOUT["n_samples"]]
    print(f"[toy] {DROPOUT['virus']} samples without coverage at {DROPOUT['start']}-{DROPOUT['end']}: {', '.join(dropout_samples)}")

    manifest, sheet = [], []
    for sample in samples:
        virus, info = sample["virus"], VIRUSES[sample["virus"]]
        snps = rng.choice(clades[virus]) + [(r, rng.randrange(len(s))) for r, (_, s) in enumerate(references[virus]) for _ in range(2)]
        genome = [(name, mutate(seq, [p for r, p in snps if r == idx], rng)) for idx, (name, seq) in enumerate(references[virus])]

        prefix = f"{sample['SampleID']}_{virus}"
        fq1, fq2, bam = (outdir / "bams" / f"{prefix}{suffix}" for suffix in ("_R1.fq.gz", "_R2.fq.gz", ".bam"))
        with gzip.open(fq1, "wt") as h1, gzip.open(fq2, "wt") as h2:
            dropout = (DROPOUT["start"], DROPOUT["end"]) if sample["SampleID"] in dropout_samples and virus == DROPOUT["virus"] else None
            simulate_reads(genome, rng.lognormvariate(5.5, 0.6), rng, h1, h2, prefix, dropout)
        map_reads(outdir / "references" / f"{virus}.fasta", fq1, fq2, bam, args.threads)
        fq1.unlink(), fq2.unlink()
        print(f"[toy] {bam.name}")

        manifest.append({"sample_id": sample["SampleID"], "taxon_id": info["taxon_id"],
                         "bam_path": str(bam.resolve()), "reference_path": str((outdir / "references" / f"{virus}.fasta").resolve())})
        sheet.append({"SampleID": sample["SampleID"], "virus": virus, "date": sample["date"].isoformat(),
                      "location": sample["location"], "virus_identified": info["seq"],
                      "Virus_PCR": "" if sample.get("pcr_missed") else info["pcr"]})

    for name, rows in (("manifest.tsv", manifest), ("samples.tsv", sheet)):
        with open(outdir / name, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=rows[0].keys(), delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
    print(f"[toy] {len(manifest)} BAMs, manifest and sample sheet written to {outdir}")


if __name__ == "__main__":
    main()
