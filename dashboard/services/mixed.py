from __future__ import annotations

from itertools import combinations

import numpy as np
import pandas as pd

from .data import load_dashboard_metadata, load_mixed_metadata
from .pileup import build_mixed_sample_pileup_assets
from .plots import build_coinfection_composite_figure

HOVER_SAMPLE_LIMIT = 2


def _as_bool(value) -> bool:
    return str(value).strip().lower() in {"1", "true", "t", "yes", "y"}


def _normalize_events(df: pd.DataFrame, source: str) -> pd.DataFrame:
    """Normalize row-wise mixed metadata into sample-level event rows."""
    if df.empty or "sample_id" not in df.columns or "canonical_strain" not in df.columns:
        return pd.DataFrame(
            columns=[
                "sample_id",
                "canonical_strain",
                "display_label",
                "date",
                "location",
                "match_value",
                "source",
            ]
        )

    out = df.copy()
    out["sample_id"] = out["sample_id"].astype(str).str.strip()
    out["canonical_strain"] = out["canonical_strain"].astype(str).str.strip()
    out = out[(out["sample_id"] != "") & (out["canonical_strain"] != "")]

    out["date"] = pd.to_datetime(out.get("date"), errors="coerce").dt.strftime("%Y-%m-%d")
    out["date"] = out["date"].fillna("NA")
    out["location"] = out.get("location", pd.Series(index=out.index, dtype=str)).fillna("NA").astype(str)

    if source == "sequencing":
        match_series = out["match_pcr"] if "match_pcr" in out.columns else pd.Series(0, index=out.index)
        out["match_value"] = match_series.map(_as_bool)
    else:
        match_series = (
            out["match_sequencing"]
            if "match_sequencing" in out.columns
            else pd.Series(0, index=out.index)
        )
        out["match_value"] = match_series.map(_as_bool)

    out["source"] = source
    out["display_label"] = out.get("display_label", out["canonical_strain"]).fillna(out["canonical_strain"])
    return out[
        [
            "sample_id",
            "canonical_strain",
            "display_label",
            "date",
            "location",
            "match_value",
            "source",
        ]
    ]


def _pair_records(events: pd.DataFrame) -> dict[tuple[str, str], list[dict]]:
    """Build pair-to-sample-details map, supporting 2+ strains per sample."""
    if events.empty:
        return {}

    pairs: dict[tuple[str, str], list[dict]] = {}
    for sample_id, sample_df in events.groupby("sample_id"):
        strains = sorted(sample_df["canonical_strain"].dropna().unique().tolist())
        if len(strains) < 2:
            continue

        date = sample_df["date"].iloc[0] if not sample_df["date"].empty else "NA"
        location = sample_df["location"].iloc[0] if not sample_df["location"].empty else "NA"
        match_value = bool(sample_df["match_value"].any())
        match_label = "Match PCR" if sample_df["source"].iloc[0] == "sequencing" else "Match Sequencing"
        detail = {
            "sample_id": sample_id,
            "date": date,
            "location": location,
            "match_label": match_label,
            "match_value": "Yes" if match_value else "No",
        }

        for a, b in combinations(strains, 2):
            key = tuple(sorted((a, b)))
            pairs.setdefault(key, []).append(detail)

    return pairs


def _build_source_matrix(events: pd.DataFrame, source_title: str, denominator_df: pd.DataFrame):
    pair_details = _pair_records(events)
    if not pair_details:
        empty = np.empty((0, 0))
        return {
            "html": build_coinfection_composite_figure(
                empty, [], np.empty((0, 0), dtype=object), [], [], title=f"{source_title} Co-infections"
            ),
            "count": 0,
            "strains": [],
        }

    strains = sorted({s for pair in pair_details.keys() for s in pair})
    index = {strain: i for i, strain in enumerate(strains)}
    z = np.full((len(strains), len(strains)), np.nan, dtype=float)
    hover = np.full((len(strains), len(strains)), "", dtype=object)

    total = 0
    for (a, b), details in pair_details.items():
        i, j = index[a], index[b]
        lo, hi = min(i, j), max(i, j)  # upper triangle only
        count = len(details)
        total += count
        z[lo, hi] = count

        sample_lines = [
            f"{d['sample_id']} | {d['date']} | {d['location']} | {d['match_label']}: {d['match_value']}"
            for d in details[:HOVER_SAMPLE_LIMIT]
        ]
        if len(details) > HOVER_SAMPLE_LIMIT:
            sample_lines.append(f"... (+{len(details) - HOVER_SAMPLE_LIMIT} more)")
        hover[lo, hi] = (
            f"{a} & {b}<br>"
            f"Pair count: {count}<br>"
            f"Samples:<br>{'<br>'.join(sample_lines)}"
        )

    # Frequency bar uses unique co-infected samples per virus as numerator,
    # denominator from all detections in the source metadata.
    sample_strain = events[["sample_id", "canonical_strain", "display_label"]].drop_duplicates()
    co_sample_count = sample_strain.groupby("sample_id")["canonical_strain"].nunique()
    coinfected_samples = set(co_sample_count[co_sample_count >= 2].index.tolist())
    coinf_events = sample_strain[sample_strain["sample_id"].isin(coinfected_samples)]

    denom = (
        denominator_df.groupby("canonical_strain").size().to_dict()
        if not denominator_df.empty and "canonical_strain" in denominator_df.columns
        else {}
    )
    denom_sub = (
        denominator_df.groupby(["canonical_strain", "display_label"]).size().to_dict()
        if not denominator_df.empty and {"canonical_strain", "display_label"}.issubset(denominator_df.columns)
        else {}
    )
    num_by_strain = (
        coinf_events.groupby("canonical_strain")["sample_id"].nunique().to_dict()
        if not coinf_events.empty
        else {}
    )
    num_by_sub = (
        coinf_events.groupby(["canonical_strain", "display_label"])["sample_id"].nunique().to_dict()
        if not coinf_events.empty
        else {}
    )

    freq = []
    freq_hover = []
    for strain in strains:
        numerator = int(num_by_strain.get(strain, 0))
        denominator = int(denom.get(strain, 0))
        pct = (100.0 * numerator / denominator) if denominator > 0 else 0.0
        freq.append(pct)

        hover_lines = [
            f"{strain}",
            f"Co-infected samples: {numerator}",
            f"Total detections: {denominator}",
            f"Frequency: {pct:.1f}%",
        ]
        sub_keys = [k for k in num_by_sub.keys() if k[0] == strain]
        if sub_keys:
            hover_lines.append("Substrain frequencies:")
            for _, sub in sorted(sub_keys, key=lambda x: x[1]):
                sub_num = int(num_by_sub.get((strain, sub), 0))
                sub_den = int(denom_sub.get((strain, sub), 0))
                # Some sources use mixed display labels (e.g. canonical vs "Human ...").
                # Fall back to canonical denominator so we never report x/0 when strain totals exist.
                if sub_den == 0:
                    sub_den = denominator
                sub_pct = (100.0 * sub_num / sub_den) if sub_den > 0 else 0.0
                hover_lines.append(f"- {sub}: {sub_pct:.1f}% ({sub_num}/{sub_den})")
        freq_hover.append("<br>".join(hover_lines))

    html = build_coinfection_composite_figure(
        z,
        strains,
        hover,
        freq,
        freq_hover,
        title=f"{source_title} Co-infections",
    )
    return {"html": html, "count": int(total), "strains": strains}


def get_mixed_pileup_context(selected_sample: str | None = None) -> dict:
    """Return sample options and pileup items for mixed-page async updates."""
    bundle = load_mixed_metadata()
    seq_events = _normalize_events(bundle.sequencing, source="sequencing")
    seq_per_sample = (
        seq_events.groupby("sample_id")["canonical_strain"].nunique()
        if not seq_events.empty
        else pd.Series(dtype=int)
    )
    coinfected_samples = sorted(seq_per_sample[seq_per_sample >= 2].index.tolist())
    active_sample = selected_sample if selected_sample in coinfected_samples else None
    if active_sample is None and coinfected_samples:
        active_sample = coinfected_samples[0]
    sample_rows = (
        seq_events[seq_events["sample_id"] == active_sample]
        if active_sample is not None
        else pd.DataFrame(columns=["canonical_strain", "display_label"])
    )
    mixed_pileup_items = build_mixed_sample_pileup_assets(active_sample or "", sample_rows)
    return {
        "mixed_pileup_sample_options": coinfected_samples,
        "mixed_pileup_selected_sample": active_sample,
        "mixed_pileup_items": mixed_pileup_items,
    }


def get_mixed_context(selected_sample: str | None = None):
    bundle = load_mixed_metadata()
    dashboard_bundle = load_dashboard_metadata()
    seq_events = _normalize_events(bundle.sequencing, source="sequencing")
    pcr_events = _normalize_events(bundle.pcr, source="pcr")

    seq_data = _build_source_matrix(
        seq_events, source_title="Sequencing", denominator_df=dashboard_bundle.sequencing
    )
    pcr_data = _build_source_matrix(
        pcr_events, source_title="PCR", denominator_df=dashboard_bundle.pcr
    )
    default_source = "sequencing" if seq_data["count"] > 0 else "pcr"

    pileup_ctx = get_mixed_pileup_context(selected_sample=selected_sample)

    return {
        "co_inf_mat_seq": seq_data["html"],
        "co_inf_mat_pcr": pcr_data["html"],
        "no_co_infections_seq": seq_data["count"],
        "no_co_infections_pcr": pcr_data["count"],
        "default_source": default_source,
        "strains": sorted(set(seq_data["strains"]) | set(pcr_data["strains"])),
        **pileup_ctx,
    }
