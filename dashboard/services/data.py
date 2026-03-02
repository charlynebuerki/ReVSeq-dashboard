from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from dashboard.config import harmonize_canton, harmonize_detected_strain

DATA_DIR = Path("dashboard/static/data")
PCR_METADATA_PATH = DATA_DIR / "metadata_pcr.tsv"
SEQ_METADATA_PATH = DATA_DIR / "metadata_sequencing.tsv"
MIXED_PCR_PATH = DATA_DIR / "mixed_pcr.tsv"
MIXED_SEQ_PATH = DATA_DIR / "mixed_sequencing.tsv"
CONFIG_PATH = Path("dashboard/config.py")


@dataclass
class MetadataBundle:
    # Separate module-specific sources to avoid accidental cross-use.
    pcr: pd.DataFrame
    sequencing: pd.DataFrame


# In-process metadata caches, invalidated by file mtime changes.
_DASHBOARD_CACHE_KEY = None
_DASHBOARD_CACHE_BUNDLE: MetadataBundle | None = None
_MIXED_CACHE_KEY = None
_MIXED_CACHE_BUNDLE: MetadataBundle | None = None


def _mtime(path: Path):
    return path.stat().st_mtime if path.exists() else None


def _read_tsv(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t")


def _add_week(df: pd.DataFrame, date_col: str = "date") -> pd.DataFrame:
    if df.empty:
        return df
    df = df.copy()
    df[date_col] = pd.to_datetime(df[date_col], errors="coerce")
    df = df.dropna(subset=[date_col])
    # Keep weekly bucketing stable across all modules.
    df["week"] = df[date_col].dt.strftime("%Y/%W")
    return df


def _harmonize_column(df: pd.DataFrame, detected_col: str) -> pd.DataFrame:
    if df.empty:
        return df
    df = df.copy()
    source = "sequencing" if detected_col == "virus_identified_sequencing" else "pcr"
    harmonized = df[detected_col].map(lambda value: harmonize_detected_strain(value, source=source))
    df[["canonical_strain", "display_label"]] = pd.DataFrame(harmonized.tolist(), index=df.index)
    # Rows that cannot be mapped to configured strains are excluded.
    df = df.dropna(subset=["canonical_strain"])
    return df


def _normalize_canton_column(df: pd.DataFrame, canton_col: str = "canton") -> pd.DataFrame:
    if df.empty:
        return df
    df = df.copy()
    df[canton_col] = df[canton_col].map(harmonize_canton)
    # Keep only rows mappable to Swiss canton codes.
    return df.dropna(subset=[canton_col])


def _normalize_match_flag(series: pd.Series) -> pd.Series:
    """Normalize match flags to integer 0/1 for aggregations."""
    lowered = series.astype(str).str.strip().str.lower()
    truthy = {"1", "true", "t", "yes", "y"}
    return lowered.isin(truthy).astype(int)


def load_sequencing_metadata(path: Path | None = None) -> pd.DataFrame:
    path = path or SEQ_METADATA_PATH
    df = _read_tsv(path)
    if df.empty:
        return df

    df = df.rename(
        columns={
            "strain": "sample_id",
            "virus_identified": "virus_identified_sequencing",
            "location": "canton",
        }
    )
    df = _add_week(df)
    df = _harmonize_column(df, "virus_identified_sequencing")
    df = _normalize_canton_column(df, "canton")
    if "Match_PCR" in df.columns:
        df["match_pcr"] = _normalize_match_flag(df["Match_PCR"])
    else:
        df["match_pcr"] = 0
    df["substrain"] = df["display_label"]
    return df


def load_pcr_metadata(path: Path | None = None) -> pd.DataFrame:
    path = path or PCR_METADATA_PATH
    df = _read_tsv(path)
    if df.empty:
        return df

    df = df.rename(
        columns={
            "strain": "sample_id",
            "virus_identified_pcr": "virus_identified_pcr",
            "location": "canton",
        }
    )
    df = _add_week(df)
    df = _harmonize_column(df, "virus_identified_pcr")
    df = _normalize_canton_column(df, "canton")
    if "Match_Sequencing" in df.columns:
        df["match_sequencing"] = _normalize_match_flag(df["Match_Sequencing"])
    else:
        df["match_sequencing"] = 0
    # PCR plots should stay at canonical strain level.
    df["substrain"] = df["display_label"]
    return df


def load_dashboard_metadata() -> MetadataBundle:
    global _DASHBOARD_CACHE_KEY, _DASHBOARD_CACHE_BUNDLE
    cache_key = (_mtime(PCR_METADATA_PATH), _mtime(SEQ_METADATA_PATH), _mtime(CONFIG_PATH))
    if _DASHBOARD_CACHE_BUNDLE is not None and _DASHBOARD_CACHE_KEY == cache_key:
        return MetadataBundle(
            pcr=_DASHBOARD_CACHE_BUNDLE.pcr.copy(),
            sequencing=_DASHBOARD_CACHE_BUNDLE.sequencing.copy(),
        )

    sequencing = load_sequencing_metadata()
    pcr = load_pcr_metadata()

    if not pcr.empty:
        pcr["display_label"] = pcr["display_label"].fillna(pcr["canonical_strain"])
        pcr["substrain"] = pcr["substrain"].fillna(pcr["display_label"])
    if not sequencing.empty:
        sequencing["display_label"] = sequencing["display_label"].fillna(sequencing["canonical_strain"])
        sequencing["substrain"] = sequencing["substrain"].fillna(sequencing["display_label"])

    bundle = MetadataBundle(pcr=pcr, sequencing=sequencing)
    _DASHBOARD_CACHE_KEY = cache_key
    _DASHBOARD_CACHE_BUNDLE = bundle
    return MetadataBundle(pcr=bundle.pcr.copy(), sequencing=bundle.sequencing.copy())


def load_mixed_metadata(
    pcr_path: Path | None = None, sequencing_path: Path | None = None
) -> MetadataBundle:
    global _MIXED_CACHE_KEY, _MIXED_CACHE_BUNDLE
    pcr_path = pcr_path or MIXED_PCR_PATH
    sequencing_path = sequencing_path or MIXED_SEQ_PATH
    cache_key = (_mtime(pcr_path), _mtime(sequencing_path), str(pcr_path), str(sequencing_path))
    if _MIXED_CACHE_BUNDLE is not None and _MIXED_CACHE_KEY == cache_key:
        return MetadataBundle(
            pcr=_MIXED_CACHE_BUNDLE.pcr.copy(),
            sequencing=_MIXED_CACHE_BUNDLE.sequencing.copy(),
        )

    pcr = _read_tsv(pcr_path)
    sequencing = _read_tsv(sequencing_path)

    if not pcr.empty:
        # Mixed inputs are already harmonized by design; only normalize date/location.
        pcr = pcr.rename(columns={"location": "canton"})
        pcr = _add_week(pcr)
        pcr = _normalize_canton_column(pcr, "canton")
    if not sequencing.empty:
        sequencing = sequencing.rename(columns={"location": "canton"})
        sequencing = _add_week(sequencing)
        sequencing = _normalize_canton_column(sequencing, "canton")

    bundle = MetadataBundle(pcr=pcr, sequencing=sequencing)
    _MIXED_CACHE_KEY = cache_key
    _MIXED_CACHE_BUNDLE = bundle
    return MetadataBundle(pcr=bundle.pcr.copy(), sequencing=bundle.sequencing.copy())
