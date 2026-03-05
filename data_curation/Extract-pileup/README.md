# Extract-pileup (BAM-first workflow)

This module generates:

- pileup depth files (`pileup_<sample>_<taxon>.txt`)
- consensus FASTA (`consensus_<sample>_<taxon>.fasta`) unless `--pileup-only` is used
- low-depth mask BED (`mask_<sample>_<taxon>.bed`)

The default workflow starts from existing BAM files.

## Requirements

- Python 3.9+
- `samtools`
- `bcftools`
- (only for FASTQ mode) `minimap2`

Install Python deps:

```bash
pip install -r requirements.txt
```

## 1. Preferred workflow: run from BAM manifest

Create a TSV/CSV manifest with columns:

- `sample_id`
- `taxon_id`
- `bam_path`
- `reference_path`

Example (`manifest.tsv`):

```tsv
sample_id taxon_id bam_path reference_path
S1 31631 /abs/path/S1.bam /abs/path/defaults/library_references/31631_reference.fasta
S2 31631 /abs/path/S2.bam /abs/path/defaults/library_references/31631_reference.fasta
```

Run batch BAM mode:

```bash
python batch_pipeline.py \
  --bam-manifest manifest.tsv \
  --outdir output \
  --min-depth 10 \
  --samtools samtools \
  --bcftools bcftools \
  --jobs 4
```

Pileup-only mode (skip consensus/VCF):

```bash
python batch_pipeline.py \
  --bam-manifest manifest.tsv \
  --outdir output \
  --pileup-only \
  --jobs 4
```

Notes:

- Resume file: `output/batch_progress_bam.txt`
- Disable resume with `--no-resume`
- Optional mpileup tuning: `--mpileup-min-mapq`, `--mpileup-min-baseq`, `--mpileup-max-depth`

## 2. Optional workflow: extract BAM/BAI from tar archives, then run batch

Use `extract_from_tar_and_run_pipeline.py` when your BAM/BAI files are inside `.tar` archives.

Expected metadata columns:

- `SampleID`
- `Reference_Taxon_ID`

Expected archive path pattern:

`<tar_root>/viral_pipeline_run/results/<SampleID>/<Reference_Taxon_ID>/`

Run:

```bash
python extract_from_tar_and_run_pipeline.py \
  --metadata data/oc43_metadata.tsv \
  --tar /path/to/run_1.tar /path/to/run_2.tar \
  --data-dir data \
  --reference-dir defaults/library_references \
  --outdir output \
  --jobs 4
```

This script:

1. extracts BAM+BAI into `data/<sample>/<taxon>/`
2. builds a BAM manifest
3. calls `batch_pipeline.py --bam-manifest ...`
4. rewrites pileup first column to sample name

## 3. FASTQ workflows (secondary)

Of note this is a minimal workflow, we encourage users to seek out a more robust consensus genome assembly pipeline tailored to their workflow.

Single sample from FASTQ:

```bash
python pipeline.py \
  --fastq1 data/sample_1.fastq.gz \
  --fastq2 data/sample_2.fastq.gz \
  --reference defaults/references/11216_reference.fasta \
  --taxon-id 11216 \
  --outdir output
```

Batch FASTQ mode:

```bash
python batch_pipeline.py \
  --fastq-dir data \
  --reference defaults/references/11216_reference.fasta \
  --taxon-id 11216 \
  --outdir output
```

Segmented FASTQ mode:

```bash
python batch_pipeline.py \
  --segmented \
  --fastq-dir data \
  --references refs/PB2.fasta refs/PB1.fasta refs/PA.fasta refs/HA.fasta refs/NP.fasta refs/NA.fasta refs/MP.fasta refs/NS.fasta \
  --segment-names PB2 PB1 PA HA NP NA MP NS \
  --taxon-id 11320 \
  --outdir output
```

Important: segmented processing currently applies to FASTQ mode, not BAM-manifest mode.

## 4. Collect all outputs

```bash
python collect_results.py --outdir output --dest processed_pileup
```

Copies all `consensus_*.fasta` and `pileup_*.txt` into `processed_pileup/`.
