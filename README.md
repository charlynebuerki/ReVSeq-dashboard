# ReVSeq Visuals: A Rapid, Modular Environment for Viral Sequence Visualization

Web dashboard for respiratory virus surveillance, built with Django and Auspice.

## Configuration (`dashboard/config.py`)

The main dashboard behavior is controlled in `dashboard/config.py`.

### Global settings

- `DEFAULT_MODULES`
  - Default module toggles on strain pages:
    - `barplot_pcr`
    - `barplot_sequencing`
    - `map`
    - `pileup`
    - `tree`
- `DEFAULT_PILEUP_LEVELS`
  - Enabled pileup display modes (for example `all`, `substrain`, `individual`)
- `DEFAULT_PILEUP_MAX_INDIVIDUAL_TRACES`
  - Max number of traces rendered in individual pileup mode
- `MIXED_PAGE_ENABLED`
  - Enables/disables the mixed (co-infection) page

### Strain settings (`STRAIN_CONFIG`)

Each strain entry (for example `Influenza_A`, `RSV`, `SARS-CoV-2`) can define:

- `label`
  - UI display name
- `data_name`
  - Canonical strain name used for metadata harmonization
- `modules`
  - Per-strain module overrides
- `pileup_data_prefix`
  - Prefix used to load pileup JSONs from `dashboard/static/data/pileup`
- `pileup_segments`
  - Segment definitions for segmented viruses
- `pileup_default_segment`
  - Default selected segment in segmented pileup views
- `trees`
  - Auspice datasets to embed on the strain page (`title`, `dataset`)

### Supporting mappings

- `STRAIN_ORDER`
  - Controls strain ordering in UI/navigation
- `COLOR_BY_STRAIN`
  - Plot color mapping by canonical strain
- Harmonization helpers
  - Strain label normalization
  - Location normalization (for canton mapping)

## Project Overview

### `dashboard/`

- Django app serving:
  - home page (all strains)
  - strain-specific pages
  - mixed/co-infection page
- Modules include:
  - barplots (sequencing and PCR)
  - maps
  - pileup visualizations
  - embedded phylogenetic trees

### `data_curation/`

- Pipeline code used to prepare curated metadata, tree inputs, and pileup inputs.
- Includes:
  - `Nextstrain-pipelines/` for tree/metadata workflows
  - `Extract-pileup/` for consensus/pileup extraction workflows

## Runtime Configuration

Main environment variables:

- `NEXTSTRAIN_BASE_URL`
  - local example: `http://127.0.0.1:4001`
  - reverse-proxy deployment example: `/nextstrain`
- `DJANGO_DEBUG` (`1` for local dev, `0` for deployment)
- `DJANGO_ALLOWED_HOSTS` (comma-separated hostnames)
- `DJANGO_SECRET_KEY` (required when `DJANGO_DEBUG=0`)

## Run Locally (No Docker)

Terminal 1 (Auspice):

```bash
HOST=0.0.0.0 PORT=4001 npx auspice-revseq view --datasetDir dashboard/static/auspice
```

Terminal 2 (Django):

```bash
export NEXTSTRAIN_BASE_URL=http://127.0.0.1:4001
python manage.py runserver
```

Open:

- `http://127.0.0.1:8000/home/`

## Run Locally (Docker Compose)

```bash
docker compose up --build
```

Open:

- Dashboard: `http://localhost:8000/home/`
- Auspice: `http://localhost:4000/`

## Deploy with Docker Images

### 1. Build and push versioned images

```bash
export REVSEQ_TAG=2026.03.05-1
export DOCKERHUB_USER=<dockerhub_user>

docker build -f Dockerfile -t ${DOCKERHUB_USER}/revseq-dashboard:${REVSEQ_TAG} .
docker build -f Dockerfile_nextstrain -t ${DOCKERHUB_USER}/revseq-nextstrain:${REVSEQ_TAG} .

docker push ${DOCKERHUB_USER}/revseq-dashboard:${REVSEQ_TAG}
docker push ${DOCKERHUB_USER}/revseq-nextstrain:${REVSEQ_TAG}
```

### 2. Set deploy environment on host

In the same directory as `docker-compose.deploy.yml`, create `.env`:

```env
REVSEQ_TAG=2026.03.05-1
DJANGO_SECRET_KEY=<strong_random_secret>
```

Ensure `docker-compose.deploy.yml` references:

- `${REVSEQ_TAG}` in image tags
- `${DJANGO_SECRET_KEY}` for dashboard environment

### 3. Pull and recreate containers

```bash
docker compose -f docker-compose.deploy.yml pull
docker compose -f docker-compose.deploy.yml up -d --force-recreate
```

### 4. Verify deployment

```bash
docker compose -f docker-compose.deploy.yml ps
docker logs revseq-dashboard --tail 80
docker logs revseq-nextstrain --tail 80
```

Optional local checks on host:

```bash
curl -I http://127.0.0.1:8000/home/
curl -I http://127.0.0.1:4000
```

If deployed behind a reverse proxy, route:

- `/` to Django (`:8000`)
- `/nextstrain/` to Auspice (`:4000`)
- `/dist/` and `/charon/` to Auspice (`:4000`)
