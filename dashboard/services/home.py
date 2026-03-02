from __future__ import annotations

from pathlib import Path
import hashlib
import pandas as pd

from .data import PCR_METADATA_PATH, SEQ_METADATA_PATH, load_dashboard_metadata
from .plots import make_weekly_strain_figure, build_weekly_canton_map


def _asset_is_stale(output_path: str | Path, input_paths: list[str | Path]) -> bool:
    output_path = Path(output_path)
    if not output_path.exists():
        return True

    output_mtime = output_path.stat().st_mtime
    # Rebuild when any input artifact (metadata or plotting code) is newer.
    for input_path in input_paths:
        input_path = Path(input_path)
        if input_path.exists() and input_path.stat().st_mtime > output_mtime:
            return True
    return False


def _asset_version(input_paths: list[str | Path]) -> str:
    mtimes = []
    for input_path in input_paths:
        p = Path(input_path)
        mtimes.append(str(p.stat().st_mtime) if p.exists() else "missing")
    return hashlib.md5("|".join(mtimes).encode("utf-8")).hexdigest()[:10]



def _aggregate_by_week(
    df: pd.DataFrame, strain_col: str = "canonical_strain", match_col: str | None = None
) -> pd.DataFrame:
    if df.empty:
        return pd.DataFrame(columns=["week", strain_col, "count", "week_start", "match_count", "match_over_total"])

    grouped = df.groupby(["week", strain_col]).size().reset_index(name="count")
    if match_col and match_col in df.columns:
        matches = df.groupby(["week", strain_col])[match_col].sum().reset_index(name="match_count")
        grouped = grouped.merge(matches, on=["week", strain_col], how="left")
    else:
        grouped["match_count"] = 0
    grouped["match_count"] = grouped["match_count"].fillna(0).astype(int)
    grouped["count"] = grouped["count"].fillna(0).astype(int)
    grouped["match_over_total"] = grouped["match_count"].astype(str) + "/" + grouped["count"].astype(str)
    grouped["week_start"] = pd.to_datetime(grouped["week"] + "/1", format="%Y/%W/%w")
    return grouped



def build_home_assets(bundle):
    plot_code_path = "dashboard/services/plots.py"
    config_path = "dashboard/config.py"
    seq_output = "dashboard/static/barplot_seq.html"
    pcr_output = "dashboard/static/barplot_pcr.html"
    map_seq_output = "dashboard/static/map_seq.html"
    map_pcr_output = "dashboard/static/map_pcr.html"

    # Pre-render static HTML plots used by iframe embeds.
    if _asset_is_stale(seq_output, [SEQ_METADATA_PATH, plot_code_path, config_path]):
        seq_grouped = _aggregate_by_week(bundle.sequencing, match_col="match_pcr")
        fig_seq = make_weekly_strain_figure(
            seq_grouped, "canonical_strain", "Sequencing", match_label="Match PCR"
        )
        fig_seq.write_html(seq_output, include_plotlyjs="cdn")

    if _asset_is_stale(pcr_output, [PCR_METADATA_PATH, plot_code_path, config_path]):
        pcr_grouped = _aggregate_by_week(bundle.pcr, match_col="match_sequencing")
        fig_pcr = make_weekly_strain_figure(
            pcr_grouped, "canonical_strain", "PCR", match_label="Match Sequencing"
        )
        fig_pcr.write_html(pcr_output, include_plotlyjs="cdn")

    if _asset_is_stale(
        map_seq_output,
        [SEQ_METADATA_PATH, "dashboard/static/swiss_cantons.geojson", plot_code_path, config_path],
    ):
        fig_map_seq = build_weekly_canton_map(bundle.sequencing)
        fig_map_seq.write_html(map_seq_output, include_plotlyjs="cdn")

    if _asset_is_stale(
        map_pcr_output,
        [PCR_METADATA_PATH, "dashboard/static/swiss_cantons.geojson", plot_code_path, config_path],
    ):
        fig_map_pcr = build_weekly_canton_map(bundle.pcr)
        fig_map_pcr.write_html(map_pcr_output, include_plotlyjs="cdn")



def get_home_context():
    bundle = load_dashboard_metadata()
    build_home_assets(bundle)

    # Requested behavior: sample count is sequencing-derived only.
    no_samples = bundle.sequencing["sample_id"].nunique() if not bundle.sequencing.empty else 0
    default_map = "seq" if not bundle.sequencing.empty else "pcr"
    version = _asset_version(
        [
            PCR_METADATA_PATH,
            SEQ_METADATA_PATH,
            "dashboard/services/plots.py",
            "dashboard/config.py",
            "dashboard/static/map_seq.html",
            "dashboard/static/map_pcr.html",
        ]
    )
    return {"no_samples": int(no_samples), "asset_version": version, "default_map": default_map}
