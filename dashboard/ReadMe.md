# Dashboard App (`dashboard/`)

This folder contains the Django app that powers the ReVSeq web interface.

## High-Level Structure

- `views.py`
  - Entry points for page routes (`home`, `strain`, `mixed`).
  - Delegates data loading and plotting logic to service modules.
- `urls.py`
  - URL routing for dashboard pages.
- `config.py`
  - Central configuration for strains, module toggles, colors, tree datasets, pileup behavior, and label harmonization.
- `services/`
  - logic split by concern:
    - `data.py`: metadata loading, normalization, harmonization, caching
    - `home.py`: home page context and top-level plot generation
    - `strain.py`: strain page context and per-strain plot generation
    - `mixed.py`: mixed/co-infection matrix and mixed-page context
    - `pileup.py`: pileup figure generation and segmented/substrain handling
    - `plots.py`: reusable plotting helpers
    - `assets.py`: static asset/version helpers
- `templates/`
  - HTML templates for `home.html`, `strain.html`, and `mixed.html`.
- `static/`
  - Input data and visualization assets:
    - `data/`: metadata + pileup JSON inputs
    - `annotations/`: GenBank annotations used in pileup views
    - `auspice/`: Nextstrain/Auspice dataset JSON files
    - generated HTML artifacts (barplots/maps/pileups)

## Core Functionalities

- Home page (`/home/`)
  - Global sequencing and PCR barplots with map support.
- Strain page (`/strain/<strain_name>/`)
  - Module-based layout controlled by config (`barplot`, `map`, `pileup`, `tree`).
  - Supports sequencing/PCR toggles and strain-specific trees.
  - Pileup supports:
    - non-segmented and segmented viruses
    - all/substrain/individual modes (config-driven)
    - capped individual traces for performance
- Mixed page (`/mixed/`)
  - Co-infection matrix visualization and mixed-sample pileup access (if enabled).

## Request/Render Flow

1. A view receives the request.
2. Service layer loads normalized metadata (cached by source file mtime).
3. Plot functions generate or refresh HTML plot assets in `dashboard/static/`.
4. Template renders iframes and UI controls for source/mode switching.

## Notes

- Plot HTML files in `dashboard/static/barplots/` and `dashboard/static/pileup_html/` are generated artifacts.
- Input metadata and tree JSON files are expected under `dashboard/static/data/` and `dashboard/static/auspice/`.
- Local/hosted tree embedding behavior is controlled by `NEXTSTRAIN_BASE_URL`.
