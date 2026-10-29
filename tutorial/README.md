# ReVSeq Visuals tutorial

Build a dashboard for your own viral sequencing data, end to end, on a small synthetic dataset:
3 viruses — RSV (A and B), Influenza A (H3N2) and human parainfluenza virus 3 — 20 samples, one co-infection.

| Part | You will | Time |
|---|---|---|
| [1. Bioinformatics](01_bioinformatics.ipynb) | turn BAM files into sequencing-depth, phylogeny and metadata files | ~1 min |
| [2. Dashboard](02_dashboard.ipynb) | configure the dashboard, run it locally, and deploy it with Docker | ~1 min |

Every step is a notebook cell you can run, inspect and edit. All outputs go to `tutorial/workdir/` (git-ignored) and `dashboard/static/` (where the dashboard reads its data), so the tutorial must run in a **fresh clone** of the repository.

## Open the tutorial

**Option A — GitHub Codespaces (nothing to install).**
[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/charlynebuerki/ReVSeq-dashboard?quickstart=1)

Click the button (or, on the repository page, *Code → Codespaces → Create codespace*). The environment is built automatically from `.devcontainer/`. Open `tutorial/01_bioinformatics.ipynb` and pick the `revseq-tutorial` kernel. When Part 2 starts the dashboard, it makes port 4001 (the tree viewer) *Public* so the embedded trees can load; if the trees stay empty, set it by hand: *Ports* tab → right-click 4001 → *Port Visibility* → *Public*.

**Option B — your computer** (macOS or Linux, [conda](https://conda-forge.org/download/) required):

```bash
git clone https://github.com/charlynebuerki/ReVSeq-dashboard.git
```

```bash
cd ReVSeq-dashboard && conda env create -f tutorial/environment.yml
```

```bash
conda activate revseq-tutorial && jupyter lab tutorial/
```

The Docker steps at the end of Part 2 also need [Docker](https://docs.docker.com/get-docker/).

## What the toy data are
`tutorial/scripts/make_toy_data.py` mutates public reference genomes shipped with the repository, simulates paired-end reads with uneven coverage (including a planted coverage gap in the G gene of three RSV-A samples) and maps them with minimap2 (fixed random seed, so every run is identical). No patient data are involved. The HPIV-3 phylogeny uses a Nextclade dataset built from public GenBank sequences (`tutorial/data/nextclade_dataset_hpiv-3.zip`).

## Build the tutorial website
The notebooks are stored with their outputs, so the website can be built without re-running them:

```bash
jupyter-book build tutorial/
```

Open `tutorial/_build/html/index.html`.
