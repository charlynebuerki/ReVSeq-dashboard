from __future__ import annotations

import hashlib
from pathlib import Path

import pandas as pd

from dashboard.config import SLUG_BY_DATA_NAME

from .data import PCR_METADATA_PATH, SEQ_METADATA_PATH, load_dashboard_metadata
from .pileup import ensure_pileup_assets
from .plots import build_weekly_canton_map, make_weekly_substrain_figure


def _asset_is_stale(output_path: str | Path, input_paths: list[str | Path]) -> bool:
    output_path = Path(output_path)
    if not output_path.exists():
        return True

    output_mtime = output_path.stat().st_mtime
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



def _aggregate_strain_by_week(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return pd.DataFrame(columns=["week", "substrain", "count", "week_start"])

    grouped = df.groupby(["week", "substrain"]).size().reset_index(name="count")
    grouped["week_start"] = pd.to_datetime(grouped["week"] + "/1", format="%Y/%W/%w")
    return grouped



def _ensure_strain_histogram(df: pd.DataFrame, output_path: str, title: str, source_files: list[str | Path]):
    if not _asset_is_stale(output_path, source_files):
        return
    # Weekly bars with optional substrain split.
    grouped = _aggregate_strain_by_week(df)
    fig = make_weekly_substrain_figure(grouped, "substrain", title=title)
    fig.write_html(output_path, include_plotlyjs="cdn")



def _select_substrains(seq_df: pd.DataFrame, pcr_df: pd.DataFrame, strain_name: str):
    substrains = []
    if not seq_df.empty:
        substrains.extend(seq_df["substrain"].dropna().unique().tolist())
    if not pcr_df.empty:
        substrains.extend(pcr_df["substrain"].dropna().unique().tolist())

    # Preserve first-seen order while deduplicating across both sources.
    ordered = []
    for value in substrains:
        if value not in ordered and value != strain_name:
            ordered.append(value)
    return ordered



def get_strain_context(strain_slug, strain_config, request_host):
    modules = strain_config["modules"]
    strain_name = strain_config["data_name"]
    plot_code_path = "dashboard/services/plots.py"
    config_path = "dashboard/config.py"

    bundle = load_dashboard_metadata()
    seq_df = (
        bundle.sequencing[bundle.sequencing["canonical_strain"] == strain_name]
        if not bundle.sequencing.empty
        else bundle.sequencing
    )
    pcr_df = (
        bundle.pcr[bundle.pcr["canonical_strain"] == strain_name]
        if not bundle.pcr.empty
        else bundle.pcr
    )

    if modules["histogram_pcr"]:
        _ensure_strain_histogram(
            pcr_df,
            f"dashboard/static/barplots/{strain_slug}_pcr.html",
            title="PCR",
            source_files=[PCR_METADATA_PATH, plot_code_path, config_path],
        )
    if modules["histogram_sequencing"]:
        _ensure_strain_histogram(
            seq_df,
            f"dashboard/static/barplots/{strain_slug}_seq.html",
            title="Sequencing",
            source_files=[SEQ_METADATA_PATH, plot_code_path, config_path],
        )

    no_samples = int(seq_df["sample_id"].nunique()) if not seq_df.empty else 0
    substrains = _select_substrains(seq_df, pcr_df, strain_name)

    if modules["map"]:
        map_seq_output = f"dashboard/static/barplots/{strain_slug}_map_seq.html"
        map_pcr_output = f"dashboard/static/barplots/{strain_slug}_map_pcr.html"
        if _asset_is_stale(
            map_seq_output,
            [SEQ_METADATA_PATH, "dashboard/static/swiss_cantons.geojson", plot_code_path, config_path],
        ):
            fig_map_seq = build_weekly_canton_map(seq_df)
            fig_map_seq.write_html(map_seq_output, include_plotlyjs="cdn")
        if _asset_is_stale(
            map_pcr_output,
            [PCR_METADATA_PATH, "dashboard/static/swiss_cantons.geojson", plot_code_path, config_path],
        ):
            fig_map_pcr = build_weekly_canton_map(pcr_df)
            fig_map_pcr.write_html(map_pcr_output, include_plotlyjs="cdn")

    if modules["pileup"]:
        pileup_slug = SLUG_BY_DATA_NAME.get(strain_name, strain_slug)
        ensure_pileup_assets(pileup_slug, strain_name, substrains, no_samples)

    hostname = request_host.split(":")[0]
    nextstrain_host = f"{hostname}:4000"
    default_histogram = "seq" if modules["histogram_sequencing"] else "pcr"
    default_map = "seq" if not seq_df.empty else "pcr"
    default_source = default_histogram if (modules["histogram_sequencing"] or modules["histogram_pcr"]) else default_map
    version = _asset_version(
        [
            PCR_METADATA_PATH,
            SEQ_METADATA_PATH,
            plot_code_path,
            config_path,
            f"dashboard/static/barplots/{strain_slug}_seq.html",
            f"dashboard/static/barplots/{strain_slug}_pcr.html",
            f"dashboard/static/barplots/{strain_slug}_map_seq.html",
            f"dashboard/static/barplots/{strain_slug}_map_pcr.html",
        ]
    )

    return {
        "strain": strain_slug,
        "strain_name": strain_name,
        "no_samples": no_samples,
        "substrains": substrains,
        "nextstrain_host": nextstrain_host,
        "tree_sources": strain_config.get("trees", []),
        "modules": modules,
        "default_histogram": default_histogram,
        "default_map": default_map,
        "default_source": default_source,
        "asset_version": version,
    }
