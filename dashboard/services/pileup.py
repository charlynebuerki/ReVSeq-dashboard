from __future__ import annotations

"""Build and cache interactive pileup visualizations.

This module reads precomputed pileup JSON assets and renders Plotly HTML views
for three resolutions:
- all samples for a strain
- per-substrain summaries
- individual sample traces
"""

import hashlib
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from Bio import SeqIO
from plotly.subplots import make_subplots
from dashboard.config import SLUG_BY_DATA_NAME, get_strain_config

PILEUP_DATA_DIR = Path("dashboard/static/data/pileup")
ANNOTATIONS_DIR = Path("dashboard/static/annotations")
PILEUP_OUTPUT_DIR = Path("dashboard/static/pileup_html")


def _asset_is_stale(output_path: Path, input_paths: list[Path]) -> bool:
    if not output_path.exists():
        return True
    output_mtime = output_path.stat().st_mtime
    for input_path in input_paths:
        if input_path.exists() and input_path.stat().st_mtime > output_mtime:
            return True
    return False


def _slugify(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", str(value).strip())


def _resolve_avg_path(prefixes: list[str], substrain: str | None = None) -> Path | None:
    """Resolve an avg JSON path from supported naming conventions."""
    for prefix in prefixes:
        if substrain is None:
            candidate = PILEUP_DATA_DIR / f"{prefix}_avg.json"
            if candidate.exists():
                return candidate
            continue

        clean = _slugify(substrain)
        lower_clean = clean.lower()
        kebab_clean = lower_clean.replace("_", "-")
        candidates = [
            PILEUP_DATA_DIR / f"{prefix}_{clean}_avg.json",
            PILEUP_DATA_DIR / f"{prefix}_{lower_clean}_avg.json",
            # Allow direct substrain filenames like rsv-a_avg.json.
            PILEUP_DATA_DIR / f"{clean}_avg.json",
            PILEUP_DATA_DIR / f"{lower_clean}_avg.json",
            PILEUP_DATA_DIR / f"{kebab_clean}_avg.json",
        ]
        for candidate in candidates:
            if candidate.exists():
                return candidate
    return None


def _resolve_individual_path(prefixes: list[str]) -> Path | None:
    """Resolve a strain-level individual-trace JSON path."""
    for prefix in prefixes:
        for suffix in ("indiv", "individual"):
            candidate = PILEUP_DATA_DIR / f"{prefix}_{suffix}.json"
            if candidate.exists():
                return candidate
    return None


def _resolve_individual_path_for_substrain(prefixes: list[str], substrain: str) -> Path | None:
    """Resolve a substrain-level individual-trace JSON path."""
    clean = _slugify(substrain)
    lower_clean = clean.lower()
    kebab_clean = lower_clean.replace("_", "-")
    for prefix in prefixes:
        for suffix in ("indiv", "individual"):
            candidates = [
                PILEUP_DATA_DIR / f"{prefix}_{clean}_{suffix}.json",
                PILEUP_DATA_DIR / f"{prefix}_{lower_clean}_{suffix}.json",
                PILEUP_DATA_DIR / f"{prefix}-{lower_clean}_{suffix}.json",
                PILEUP_DATA_DIR / f"{prefix}-{kebab_clean}_{suffix}.json",
                PILEUP_DATA_DIR / f"{clean}_{suffix}.json",
                PILEUP_DATA_DIR / f"{lower_clean}_{suffix}.json",
                PILEUP_DATA_DIR / f"{kebab_clean}_{suffix}.json",
            ]
            for candidate in candidates:
                if candidate.exists():
                    return candidate
    return None


def _resolve_segmented_individual_path(prefix: str) -> Path | None:
    for suffix in ("indiv", "individual"):
        candidate = PILEUP_DATA_DIR / f"{prefix}_{suffix}.json"
        if candidate.exists():
            return candidate
    return None


def _discover_segmented_subtypes(prefix_candidates: list[str]) -> list[dict]:
    """Discover subtype data prefixes from segmented avg files on disk."""
    found: dict[str, dict] = {}
    for prefix in prefix_candidates:
        # Single segmented file with no explicit subtype (e.g. flu-b_avg.json).
        base_candidate = PILEUP_DATA_DIR / f"{prefix}_avg.json"
        if base_candidate.exists():
            value = prefix.lower()
            found[value] = {"value": value, "label": prefix.upper(), "data_prefix": prefix}

        for path in sorted(PILEUP_DATA_DIR.glob(f"{prefix}-*_avg.json")) + sorted(
            PILEUP_DATA_DIR.glob(f"{prefix}_*_avg.json")
        ):
            stem = path.stem
            if not stem.endswith("_avg"):
                continue
            data_prefix = stem[: -len("_avg")]
            if data_prefix == prefix:
                # Already handled as base candidate; this is non-subtyped segmented data.
                continue
            subtype_value = data_prefix
            for sep in ("-", "_"):
                marker = f"{prefix}{sep}"
                if data_prefix.startswith(marker):
                    subtype_value = data_prefix[len(marker) :]
                    break
            subtype_key = subtype_value.lower()
            if subtype_key.startswith("h") and "n" in subtype_key and subtype_key[1].isdigit():
                label = subtype_key.upper()
            else:
                label = subtype_value.replace("-", " ").replace("_", " ").title()
            found[subtype_key] = {"value": subtype_key, "label": label, "data_prefix": data_prefix}
    return [found[k] for k in sorted(found.keys())]


def _extract_segmented_segment_names(avg_path: Path) -> list[str]:
    with open(avg_path, "r", encoding="utf-8") as fh:
        payload = json.load(fh)
    if isinstance(payload, dict) and isinstance(payload.get("segments"), dict):
        return sorted(payload["segments"].keys())
    if isinstance(payload, dict) and payload.get("type") == "segmented":
        for key in ("segments", "data", "values"):
            container = payload.get(key)
            if isinstance(container, dict):
                return sorted(container.keys())
    return []


def _coerce_avg_payload_to_df(payload: object, path: Path) -> pd.DataFrame:
    """Normalize avg/IQR payload content to a standard dataframe."""
    if isinstance(payload, dict) and "positions" in payload:
        df = pd.DataFrame(
            {
                "position": payload["positions"],
                "median": payload.get("median") or payload.get("avg") or payload.get("mean"),
                "q1": payload.get("q1") or payload.get("iqr_low") or payload.get("p25"),
                "q3": payload.get("q3") or payload.get("iqr_high") or payload.get("p75"),
            }
        )
    else:
        records = payload.get("records", payload) if isinstance(payload, dict) else payload
        df = pd.DataFrame(records)
        rename_map = {
            "pos": "position",
            "base": "position",
            "avg": "median",
            "mean": "median",
            "iqr_low": "q1",
            "iqr_high": "q3",
            "p25": "q1",
            "p75": "q3",
        }
        df = df.rename(columns={k: v for k, v in rename_map.items() if k in df.columns})

    required = {"position", "median", "q1", "q3"}
    if not required.issubset(df.columns):
        missing = required - set(df.columns)
        raise ValueError(f"Missing avg pileup columns in {path}: {missing}")

    df = df[["position", "median", "q1", "q3"]].copy()
    df["position"] = pd.to_numeric(df["position"], errors="coerce")
    df["median"] = pd.to_numeric(df["median"], errors="coerce")
    df["q1"] = pd.to_numeric(df["q1"], errors="coerce")
    df["q3"] = pd.to_numeric(df["q3"], errors="coerce")
    df = df.dropna(subset=["position", "median", "q1", "q3"]).sort_values("position")
    return df


def _extract_segment_payload(payload: object, segment_name: str) -> object | None:
    """Extract one segment payload from supported segmented JSON layouts."""
    if not isinstance(payload, dict):
        return None
    if isinstance(payload.get("segments"), dict):
        return payload["segments"].get(segment_name)
    if payload.get("type") == "segmented":
        for key in ("segments", "data", "values"):
            container = payload.get(key)
            if isinstance(container, dict) and segment_name in container:
                return container.get(segment_name)
    if segment_name in payload and isinstance(payload.get(segment_name), (dict, list)):
        return payload.get(segment_name)
    return None


def _load_avg_depth_json(path: Path, segment_name: str | None = None) -> pd.DataFrame:
    """Load avg/IQR depth JSON into normalized columns.

    When `segment_name` is provided, segmented payloads are resolved to that segment.
    """
    with open(path, "r", encoding="utf-8") as fh:
        payload = json.load(fh)

    if segment_name is not None:
        segment_payload = _extract_segment_payload(payload, segment_name)
        if segment_payload is None:
            raise ValueError(f"Missing segment '{segment_name}' in {path}")
        return _coerce_avg_payload_to_df(segment_payload, path)

    return _coerce_avg_payload_to_df(payload, path)


def _load_individual_depth_json(path: Path) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    """Load individual depth JSON into aligned arrays."""
    with open(path, "r", encoding="utf-8") as fh:
        payload = json.load(fh)

    if isinstance(payload, dict) and "positions" in payload and "samples" in payload:
        positions = np.asarray(payload["positions"], dtype=float)
        samples = {
            str(sample_id): np.asarray(depths, dtype=float)
            for sample_id, depths in payload["samples"].items()
        }
        return positions, samples

    records = payload.get("records", payload) if isinstance(payload, dict) else payload
    df = pd.DataFrame(records)
    rename_map = {
        "pos": "position",
        "base": "position",
        "sample": "sample_id",
        "sample_name": "sample_id",
        "value": "depth",
        "coverage": "depth",
    }
    df = df.rename(columns={k: v for k, v in rename_map.items() if k in df.columns})

    required = {"position", "sample_id", "depth"}
    if not required.issubset(df.columns):
        missing = required - set(df.columns)
        raise ValueError(f"Missing individual pileup columns in {path}: {missing}")

    df = df[["position", "sample_id", "depth"]].copy()
    df["position"] = pd.to_numeric(df["position"], errors="coerce")
    df["depth"] = pd.to_numeric(df["depth"], errors="coerce")
    df = df.dropna(subset=["position", "sample_id", "depth"])

    positions = np.sort(df["position"].unique())
    samples: dict[str, np.ndarray] = {}
    for sample_id, group in df.groupby("sample_id"):
        group = group.sort_values("position")
        aligned = pd.Series(group["depth"].values, index=group["position"].values).reindex(positions, fill_value=np.nan)
        samples[str(sample_id)] = aligned.values.astype(float)

    return positions.astype(float), samples


def _load_individual_depth_json_for_segment(
    path: Path, segment_name: str
) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    with open(path, "r", encoding="utf-8") as fh:
        payload = json.load(fh)
    if isinstance(payload, dict):
        segment_payload = _extract_segment_payload(payload, segment_name)
        if segment_payload is not None:
            return _load_individual_payload(segment_payload, path)
    return _load_individual_depth_json(path)


def _load_individual_payload(payload: object, path: Path) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    if isinstance(payload, dict) and "positions" in payload and "samples" in payload:
        positions = np.asarray(payload["positions"], dtype=float)
        samples = {
            str(sample_id): np.asarray(depths, dtype=float)
            for sample_id, depths in payload["samples"].items()
        }
        return positions, samples

    records = payload.get("records", payload) if isinstance(payload, dict) else payload
    df = pd.DataFrame(records)
    rename_map = {
        "pos": "position",
        "base": "position",
        "sample": "sample_id",
        "sample_name": "sample_id",
        "value": "depth",
        "coverage": "depth",
    }
    df = df.rename(columns={k: v for k, v in rename_map.items() if k in df.columns})

    required = {"position", "sample_id", "depth"}
    if not required.issubset(df.columns):
        missing = required - set(df.columns)
        raise ValueError(f"Missing individual pileup columns in {path}: {missing}")

    df = df[["position", "sample_id", "depth"]].copy()
    df["position"] = pd.to_numeric(df["position"], errors="coerce")
    df["depth"] = pd.to_numeric(df["depth"], errors="coerce")
    df = df.dropna(subset=["position", "sample_id", "depth"])

    positions = np.sort(df["position"].unique())
    samples: dict[str, np.ndarray] = {}
    for sample_id, group in df.groupby("sample_id"):
        group = group.sort_values("position")
        aligned = pd.Series(group["depth"].values, index=group["position"].values).reindex(
            positions, fill_value=np.nan
        )
        samples[str(sample_id)] = aligned.values.astype(float)

    return positions.astype(float), samples


def _load_annotation_features(annotation_path: Path) -> pd.DataFrame:
    """Load CDS/gene intervals from a GenBank annotation file."""
    if not annotation_path.exists():
        return pd.DataFrame(columns=["start", "end", "label"])

    record = SeqIO.read(str(annotation_path), "genbank")
    rows = []
    for feature in record.features:
        if feature.type not in {"CDS", "gene"}:
            continue
        start = int(feature.location.start) + 1
        end = int(feature.location.end)
        qualifiers = feature.qualifiers
        label = (
            qualifiers.get("gene", [None])[0]
            or qualifiers.get("product", [None])[0]
            or qualifiers.get("locus_tag", [None])[0]
            or feature.type
        )
        rows.append({"start": start, "end": end, "label": str(label)})

    if not rows:
        return pd.DataFrame(columns=["start", "end", "label"])

    df = pd.DataFrame(rows).drop_duplicates(subset=["start", "end", "label"]).sort_values("start")
    return df


def _load_annotation_length(annotation_path: Path) -> int | None:
    if not annotation_path.exists():
        return None
    record = SeqIO.read(str(annotation_path), "genbank")
    seq = getattr(record, "seq", None)
    if seq is None:
        return None
    length = len(seq)
    return length if length > 0 else None


def _resolve_annotation_max_len(annotation_path: Path, features: pd.DataFrame) -> int | None:
    if not features.empty:
        try:
            return int(features["end"].max())
        except Exception:
            pass
    return _load_annotation_length(annotation_path)


def _clamp_avg_df(avg_df: pd.DataFrame, max_len: int | None) -> pd.DataFrame:
    if max_len is None or avg_df.empty:
        return avg_df
    return avg_df[avg_df["position"] <= max_len].copy()


def _clamp_individual_series(
    positions: np.ndarray, sample_depths: dict[str, np.ndarray], max_len: int | None
) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    if max_len is None or positions.size == 0:
        return positions, sample_depths
    mask = positions <= max_len
    if not mask.any():
        return positions[:0], {k: v[:0] for k, v in sample_depths.items()}
    clamped_positions = positions[mask]
    clamped_samples = {k: v[mask] for k, v in sample_depths.items()}
    return clamped_positions, clamped_samples


def _add_annotation_track(fig: go.Figure, features: pd.DataFrame) -> None:
    if features.empty:
        return

    centers = []
    labels = []
    colors = ["rgba(15,122,167,0.35)", "rgba(29,63,114,0.35)"]

    for idx, row in features.reset_index(drop=True).iterrows():
        fig.add_shape(
            type="rect",
            x0=float(row["start"]),
            x1=float(row["end"]),
            y0=0,
            y1=1,
            line={"width": 0},
            fillcolor=colors[idx % 2],
            row=2,
            col=1,
        )
        centers.append((float(row["start"]) + float(row["end"])) / 2)
        labels.append(str(row["label"]))

    fig.add_trace(
        go.Scatter(
            x=centers,
            y=[0.5] * len(centers),
            mode="text",
            text=labels,
            hoverinfo="skip",
            showlegend=False,
            textfont={"size": 10, "color": "#1d2a3a"},
        ),
        row=2,
        col=1,
    )


def build_avg_pileup_figure(
    avg_df: pd.DataFrame, features: pd.DataFrame, title: str, max_len: int | None = None
) -> go.Figure:
    """Build the summary pileup view: median line + IQR + genome track."""
    avg_df = _clamp_avg_df(avg_df, max_len)
    fig = make_subplots(
        rows=2,
        cols=1,
        shared_xaxes=True,
        row_heights=[0.78, 0.22],
        vertical_spacing=0.03,
    )

    x = avg_df["position"].to_numpy()
    q1 = avg_df["q1"].clip(lower=0.1).to_numpy()
    q3 = avg_df["q3"].clip(lower=0.1).to_numpy()
    median = avg_df["median"].clip(lower=0.1).to_numpy()

    fig.add_trace(
        go.Scatter(x=x, y=q3, mode="lines", line={"width": 0}, showlegend=False, hoverinfo="skip"),
        row=1,
        col=1,
    )
    fig.add_trace(
        go.Scatter(
            x=x,
            y=q1,
            mode="lines",
            line={"width": 0},
            fill="tonexty",
            fillcolor="rgba(15,122,167,0.25)",
            name="IQR",
            hoverinfo="skip",
        ),
        row=1,
        col=1,
    )
    fig.add_trace(
        go.Scatter(
            x=x,
            y=median,
            mode="lines",
            line={"color": "#123b6d", "width": 1.5},
            name="Median",
            hovertemplate="Position=%{x}<br>Median depth=%{y:.1f}<extra></extra>",
        ),
        row=1,
        col=1,
    )

    _add_annotation_track(fig, features)

    fig.update_yaxes(type="log", title_text="Depth", row=1, col=1)
    fig.update_yaxes(visible=False, range=[0, 1], row=2, col=1)
    fig.update_xaxes(title_text="Genome position", row=2, col=1)
    fig.update_layout(
        template="simple_white",
        title=title,
        hovermode="x unified",
        margin={"l": 50, "r": 20, "t": 50, "b": 40},
        legend={"orientation": "h", "x": 0, "y": 1.02},
    )
    if max_len:
        fig.update_xaxes(range=[1, max_len])
    return fig


def build_multi_avg_pileup_figure(
    avg_by_substrain: dict[str, pd.DataFrame],
    features: pd.DataFrame,
    title: str,
    max_len: int | None = None,
) -> go.Figure:
    """Build a combined summary view with one median trace per substrain."""
    fig = make_subplots(
        rows=2,
        cols=1,
        shared_xaxes=True,
        row_heights=[0.78, 0.22],
        vertical_spacing=0.03,
    )
    colors = ["#123b6d", "#0f7aa7", "#2a9d8f", "#bc6c25", "#6a4c93"]

    for idx, (substrain, avg_df) in enumerate(sorted(avg_by_substrain.items(), key=lambda x: x[0])):
        avg_df = _clamp_avg_df(avg_df, max_len)
        x = avg_df["position"].to_numpy()
        q1 = avg_df["q1"].clip(lower=0.1).to_numpy()
        q3 = avg_df["q3"].clip(lower=0.1).to_numpy()
        median = avg_df["median"].clip(lower=0.1).to_numpy()
        color = colors[idx % len(colors)]

        fig.add_trace(
            go.Scatter(
                x=x,
                y=q3,
                mode="lines",
                line={"width": 0},
                showlegend=False,
                hoverinfo="skip",
                legendgroup=substrain,
            ),
            row=1,
            col=1,
        )
        fig.add_trace(
            go.Scatter(
                x=x,
                y=q1,
                mode="lines",
                line={"width": 0},
                fill="tonexty",
                fillcolor="rgba(15,122,167,0.12)",
                name=f"{substrain} IQR",
                hoverinfo="skip",
                legendgroup=substrain,
                showlegend=False,
            ),
            row=1,
            col=1,
        )
        fig.add_trace(
            go.Scatter(
                x=x,
                y=median,
                mode="lines",
                line={"color": color, "width": 1.8},
                name=substrain,
                hovertemplate=f"Substrain={substrain}<br>Position=%{{x}}<br>Median depth=%{{y:.1f}}<extra></extra>",
                legendgroup=substrain,
            ),
            row=1,
            col=1,
        )

    _add_annotation_track(fig, features)
    fig.update_yaxes(type="log", title_text="Depth", row=1, col=1)
    fig.update_yaxes(visible=False, range=[0, 1], row=2, col=1)
    fig.update_xaxes(title_text="Genome position", row=2, col=1)
    fig.update_layout(
        template="simple_white",
        title=title,
        hovermode="x unified",
        margin={"l": 50, "r": 20, "t": 50, "b": 40},
        legend={"orientation": "h", "x": 0, "y": 1.02},
    )
    if max_len:
        fig.update_xaxes(range=[1, max_len])
    return fig


def build_individual_pileup_figure(
    positions: np.ndarray,
    sample_depths: dict[str, np.ndarray],
    features: pd.DataFrame,
    title: str,
    max_traces: int,
    max_len: int | None = None,
) -> go.Figure:
    """Build the individual sample view with trace capping for responsiveness."""
    positions, sample_depths = _clamp_individual_series(positions, sample_depths, max_len)
    fig = make_subplots(
        rows=2,
        cols=1,
        shared_xaxes=True,
        row_heights=[0.78, 0.22],
        vertical_spacing=0.03,
    )

    sample_ids = sorted(sample_depths.keys())
    truncated = len(sample_ids) > max_traces
    if truncated:
        sample_ids = sample_ids[:max_traces]

    for sample_id in sample_ids:
        y = np.clip(sample_depths[sample_id], a_min=0.1, a_max=None)
        fig.add_trace(
            go.Scattergl(
                x=positions,
                y=y,
                mode="lines",
                line={"width": 1},
                name=sample_id,
                hovertemplate=f"Sample={sample_id}<br>Position=%{{x}}<br>Depth=%{{y:.1f}}<extra></extra>",
            ),
            row=1,
            col=1,
        )

    _add_annotation_track(fig, features)

    fig.update_yaxes(type="log", title_text="Depth", row=1, col=1)
    fig.update_yaxes(visible=False, range=[0, 1], row=2, col=1)
    fig.update_xaxes(title_text="Genome position", row=2, col=1)
    fig.update_layout(
        template="simple_white",
        title=title,
        hovermode="closest",
        margin={"l": 50, "r": 20, "t": 50, "b": 40},
        showlegend=False,
    )
    if max_len:
        fig.update_xaxes(range=[1, max_len])

    if truncated:
        fig.add_annotation(
            xref="paper",
            yref="paper",
            x=1,
            y=1.08,
            showarrow=False,
            text=f"Showing first {max_traces} of {len(sample_depths)} samples",
        )

    return fig


def _write_figure(fig: go.Figure, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.write_html(str(output_path), include_plotlyjs="cdn")


def _asset_version(paths: list[Path]) -> str:
    mtimes = [str(path.stat().st_mtime) if path.exists() else "missing" for path in paths]
    return hashlib.md5("|".join(mtimes).encode("utf-8")).hexdigest()[:10]


def _merge_individual_traces(
    per_substrain: dict[str, tuple[np.ndarray, dict[str, np.ndarray]]]
) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    """Merge multiple substrain individual datasets on a shared position axis."""
    all_positions = set()
    for positions, _samples in per_substrain.values():
        all_positions.update(int(p) for p in positions.tolist())
    merged_positions = np.array(sorted(all_positions), dtype=float)

    merged_samples: dict[str, np.ndarray] = {}
    for substrain, (positions, samples) in per_substrain.items():
        pos_index = pd.Index(positions.astype(float))
        for sample_id, values in samples.items():
            aligned = (
                pd.Series(values, index=pos_index)
                .reindex(merged_positions, fill_value=np.nan)
                .values.astype(float)
            )
            merged_samples[f"{substrain}:{sample_id}"] = aligned

    return merged_positions, merged_samples


def _build_segmented_subtype_context(
    strain_slug: str,
    data_prefix: str,
    segmented_subtypes: list[dict],
    segmented_segments: list[dict],
    levels: list[str],
    max_individual_traces: int,
    default_segment: str | None = None,
) -> dict:
    """Build segmented pileup assets with all/substrain/individual views."""
    assets: dict[str, object] = {
        "available_views": [],
        "default_view": None,
        "all_file": None,
        "individual_file": None,
        "substrain_files": {},
        "mode": "segmented",
        "segmented_all_files": {},
        "segmented_substrain_files": {},
        "segmented_individual_files": {},
        "segmented_subtypes": [],
        "segmented_segments": [],
        "segmented_default_subtype": None,
        "segmented_default_segment": None,
    }
    cache_inputs: list[Path] = [Path(__file__)]
    prefix_candidates = [data_prefix] if data_prefix else []

    auto_subtypes = segmented_subtypes or _discover_segmented_subtypes(prefix_candidates)
    subtype_options: list[dict] = []
    subtype_avg_paths: dict[str, Path] = {}
    subtype_indiv_paths: dict[str, Path] = {}
    for subtype in auto_subtypes:
        subtype_value = str(subtype.get("value") or subtype.get("key") or subtype.get("label") or "").strip()
        if not subtype_value:
            continue
        subtype_label = str(subtype.get("label") or subtype_value)
        subtype_prefix = str(subtype.get("data_prefix") or subtype_value)
        avg_path = _resolve_avg_path([subtype_prefix, subtype_value])
        if avg_path is None:
            continue
        subtype_options.append({"value": subtype_value, "label": subtype_label, "data_prefix": subtype_prefix})
        subtype_avg_paths[subtype_value] = avg_path
        indiv_path = _resolve_segmented_individual_path(subtype_prefix)
        if indiv_path is not None:
            subtype_indiv_paths[subtype_value] = indiv_path

    # Auto-discover segments from the first valid subtype avg file if not configured.
    if segmented_segments:
        segment_options = []
        for segment in segmented_segments:
            segment_value = str(segment.get("value") or segment.get("key") or segment.get("label") or "").strip()
            if not segment_value:
                continue
            segment_options.append(
                {
                    "value": segment_value,
                    "label": str(segment.get("label") or segment_value),
                    "annotation": str(segment.get("annotation") or ""),
                }
            )
    else:
        segment_options = []
        if subtype_options:
            first_sub = subtype_options[0]["value"]
            names = _extract_segmented_segment_names(subtype_avg_paths[first_sub])
            for name in names:
                segment_options.append(
                    {
                        "value": name,
                        "label": name,
                        "annotation": f"Influenza_A_{name}.gb",
                    }
                )

    all_files: dict[str, str] = {}
    substrain_files: dict[str, dict[str, str]] = {}
    individual_files: dict[str, dict[str, str]] = {}

    for segment in segment_options:
        segment_value = segment["value"]
        segment_label = segment["label"]
        annotation_name = segment.get("annotation") or ""
        if not annotation_name:
            continue
        annotation_path = ANNOTATIONS_DIR / annotation_name
        features = _load_annotation_features(annotation_path)
        max_len = _resolve_annotation_max_len(annotation_path, features)

        if "all" in levels:
            avg_by_subtype: dict[str, pd.DataFrame] = {}
            for subtype in subtype_options:
                subtype_value = subtype["value"]
                avg_path = subtype_avg_paths[subtype_value]
                try:
                    avg_by_subtype[subtype["label"]] = _load_avg_depth_json(
                        avg_path, segment_name=segment_value
                    )
                except ValueError:
                    continue
            if avg_by_subtype:
                output_path = PILEUP_OUTPUT_DIR / f"{strain_slug}_seg_all_{_slugify(segment_value)}.html"
                input_paths = [annotation_path, Path(__file__)] + [
                    subtype_avg_paths[s["value"]] for s in subtype_options if s["value"] in subtype_avg_paths
                ]
                if _asset_is_stale(output_path, input_paths):
                    fig = build_multi_avg_pileup_figure(
                        avg_by_subtype,
                        features,
                        title=f"Pileup - All samples - {segment_label}",
                        max_len=max_len,
                    )
                    _write_figure(fig, output_path)
                all_files[segment_value] = f"pileup_html/{output_path.name}"
                cache_inputs.extend(input_paths + [output_path])

        if "substrain" in levels:
            per_segment_subtypes: dict[str, str] = {}
            for subtype in subtype_options:
                subtype_value = subtype["value"]
                subtype_label = subtype["label"]
                avg_path = subtype_avg_paths[subtype_value]
                try:
                    avg_df = _load_avg_depth_json(avg_path, segment_name=segment_value)
                except ValueError:
                    continue
                output_path = PILEUP_OUTPUT_DIR / (
                    f"{strain_slug}_seg_sub_{_slugify(subtype_value)}_{_slugify(segment_value)}.html"
                )
                if _asset_is_stale(output_path, [avg_path, annotation_path, Path(__file__)]):
                    fig = build_avg_pileup_figure(
                        avg_df,
                        features,
                        title=f"Pileup - {subtype_label} - {segment_label}",
                        max_len=max_len,
                    )
                    _write_figure(fig, output_path)
                per_segment_subtypes[subtype_value] = f"pileup_html/{output_path.name}"
                cache_inputs.extend([avg_path, annotation_path, output_path])
            if per_segment_subtypes:
                substrain_files[segment_value] = per_segment_subtypes

        if "individual" in levels:
            per_segment_individual: dict[str, str] = {}
            for subtype in subtype_options:
                subtype_value = subtype["value"]
                subtype_label = subtype["label"]
                indiv_path = subtype_indiv_paths.get(subtype_value)
                if indiv_path is None:
                    continue
                try:
                    positions, sample_depths = _load_individual_depth_json_for_segment(
                        indiv_path, segment_name=segment_value
                    )
                except ValueError:
                    continue
                output_path = PILEUP_OUTPUT_DIR / (
                    f"{strain_slug}_seg_ind_{_slugify(subtype_value)}_{_slugify(segment_value)}_{max_individual_traces}.html"
                )
                if _asset_is_stale(output_path, [indiv_path, annotation_path, Path(__file__)]):
                    fig = build_individual_pileup_figure(
                        positions,
                        sample_depths,
                        features,
                        title=f"Pileup - Individuals - {subtype_label} - {segment_label}",
                        max_traces=max_individual_traces,
                        max_len=max_len,
                    )
                    _write_figure(fig, output_path)
                per_segment_individual[subtype_value] = f"pileup_html/{output_path.name}"
                cache_inputs.extend([indiv_path, annotation_path, output_path])
            if per_segment_individual:
                individual_files[segment_value] = per_segment_individual

    assets["segmented_subtypes"] = [{"value": s["value"], "label": s["label"]} for s in subtype_options]
    assets["segmented_segments"] = [{"value": s["value"], "label": s["label"]} for s in segment_options]
    assets["segmented_all_files"] = all_files
    assets["segmented_substrain_files"] = substrain_files
    assets["segmented_individual_files"] = individual_files

    if all_files:
        assets["available_views"].append("all")
    if substrain_files and len(subtype_options) > 1:
        assets["available_views"].append("substrain")
    if individual_files:
        assets["available_views"].append("individual")

    if assets["available_views"]:
        assets["default_view"] = "all" if "all" in assets["available_views"] else assets["available_views"][0]
        assets["enabled"] = True
    else:
        assets["enabled"] = False

    if assets["segmented_subtypes"]:
        assets["segmented_default_subtype"] = assets["segmented_subtypes"][0]["value"]

    segment_values = [s["value"] for s in assets["segmented_segments"]]
    if default_segment and default_segment in segment_values:
        assets["segmented_default_segment"] = default_segment
    elif "HA" in segment_values:
        assets["segmented_default_segment"] = "HA"
    elif segment_values:
        assets["segmented_default_segment"] = segment_values[0]

    assets["version"] = _asset_version(cache_inputs)
    return assets


def _resolve_annotation_for_segment(strain_config: dict, strain_slug: str, segment: str) -> str:
    for item in strain_config.get("pileup_segments", []):
        if str(item.get("value")) == str(segment) and item.get("annotation"):
            return str(item["annotation"])
    return strain_config.get("pileup_annotation", f"{strain_slug}.gb")


def _extract_segment_names_from_individual_json(path: Path) -> list[str]:
    with open(path, "r", encoding="utf-8") as fh:
        payload = json.load(fh)
    segments = payload.get("segments") if isinstance(payload, dict) else None
    if isinstance(segments, dict):
        return sorted(str(k) for k in segments.keys())
    return []


def _pick_sample_trace(sample_depths: dict[str, np.ndarray], sample_id: str) -> tuple[str | None, np.ndarray | None]:
    def _normalize_sample_token(value: str) -> str:
        token = str(value).strip()
        if token.lower().startswith("m2-"):
            return token[3:]
        return token

    if sample_id in sample_depths:
        return sample_id, sample_depths[sample_id]

    # Fallback: mixed metadata sample IDs can be "sample|date|location" while
    # pileup JSON sample IDs are just "sample".
    sample_id_raw = _normalize_sample_token(str(sample_id).strip())
    sample_id_base = _normalize_sample_token(sample_id_raw.split("|", 1)[0].strip())

    lowered = {
        _normalize_sample_token(str(k).strip()).lower(): k
        for k in sample_depths.keys()
    }
    key = lowered.get(sample_id_raw.lower())
    if key is not None:
        return str(key), sample_depths[str(key)]

    key = lowered.get(sample_id_base.lower())
    if key is not None:
        return str(key), sample_depths[str(key)]

    # Last fallback: match JSON key prefix before first "|" to metadata base.
    for original_key in sample_depths.keys():
        key_base = _normalize_sample_token(str(original_key).split("|", 1)[0].strip()).lower()
        if key_base and key_base == sample_id_base.lower():
            return str(original_key), sample_depths[str(original_key)]
    return None, None


def build_mixed_sample_pileup_assets(sample_id: str, sample_rows: pd.DataFrame) -> list[dict]:
    """Build one pileup iframe asset per virus detected in a mixed sample."""
    if not sample_id or sample_rows.empty:
        return []

    unique_rows = (
        sample_rows[["canonical_strain", "display_label"]]
        .dropna(subset=["canonical_strain"])
        .drop_duplicates()
    )
    assets: list[dict] = []

    for row in unique_rows.itertuples(index=False):
        canonical = str(row.canonical_strain)
        display = str(row.display_label) if pd.notna(row.display_label) else canonical
        label = display if display != canonical else canonical
        slug = SLUG_BY_DATA_NAME.get(canonical)
        if slug is None:
            assets.append(
                {
                    "label": label,
                    "available": False,
                    "message": "not available",
                }
            )
            continue

        strain_config = get_strain_config(slug) or {}
        if not strain_config.get("modules", {}).get("pileup", False):
            assets.append(
                {
                    "label": label,
                    "available": False,
                    "message": "not available",
                }
            )
            continue

        prefix = strain_config.get("pileup_data_prefix", slug)
        prefix_candidates: list[str] = []
        for candidate in (prefix, slug):
            if candidate and candidate not in prefix_candidates:
                prefix_candidates.append(candidate)

        indiv_path = None
        if display != canonical:
            indiv_path = _resolve_individual_path_for_substrain(prefix_candidates, display)
        if indiv_path is None:
            indiv_path = _resolve_individual_path(prefix_candidates)
        if indiv_path is None:
            assets.append({"label": label, "available": False, "message": "not available"})
            continue

        segment = None
        segment_names = _extract_segment_names_from_individual_json(indiv_path)
        if segment_names:
            configured_default = str(strain_config.get("pileup_default_segment", "")).strip()
            if configured_default and configured_default in segment_names:
                segment = configured_default
            elif "HA" in segment_names:
                segment = "HA"
            else:
                segment = segment_names[0]

        if segment is not None:
            annotation_name = _resolve_annotation_for_segment(strain_config, slug, segment)
            annotation_path = ANNOTATIONS_DIR / annotation_name
            positions, sample_depths = _load_individual_depth_json_for_segment(indiv_path, segment)
            title = f"Pileup - {label} - {sample_id} - {segment}"
        else:
            annotation_name = strain_config.get("pileup_annotation", f"{slug}.gb")
            annotation_path = ANNOTATIONS_DIR / annotation_name
            positions, sample_depths = _load_individual_depth_json(indiv_path)
            title = f"Pileup - {label} - {sample_id}"

        trace_name, trace = _pick_sample_trace(sample_depths, sample_id)
        if trace_name is None or trace is None:
            assets.append(
                {
                    "label": label,
                    "available": False,
                    "message": "not available",
                }
            )
            continue

        output_name = f"mixed_{_slugify(sample_id)}_{_slugify(label)}"
        if segment is not None:
            output_name = f"{output_name}_{_slugify(segment)}"
        output_path = PILEUP_OUTPUT_DIR / f"{output_name}.html"

        input_paths = [indiv_path, annotation_path, Path(__file__)]
        if _asset_is_stale(output_path, input_paths):
            features = _load_annotation_features(annotation_path)
            max_len = _resolve_annotation_max_len(annotation_path, features)
            fig = build_individual_pileup_figure(
                positions=positions,
                sample_depths={trace_name: trace},
                features=features,
                title=title,
                max_traces=1,
                max_len=max_len,
            )
            _write_figure(fig, output_path)

        assets.append(
            {
                "label": label,
                "available": True,
                "file": f"pileup_html/{output_path.name}",
                "version": _asset_version(input_paths + [output_path]),
            }
        )

    return assets


def build_pileup_context(
    strain_slug: str,
    strain_name: str,
    data_prefix: str,
    substrains: list[str],
    annotation_name: str,
    levels: list[str],
    max_individual_traces: int,
    segmented_subtypes: list[dict] | None = None,
    segmented_segments: list[dict] | None = None,
    segmented_default_segment: str | None = None,
) -> dict:
    """Build pileup asset metadata and render missing/stale HTML outputs.

    The function prefers direct strain files and falls back to combining
    substrain files when needed (for both `all` and `individual` views).
    """
    if segmented_segments:
        return _build_segmented_subtype_context(
            strain_slug=strain_slug,
            data_prefix=data_prefix,
            segmented_subtypes=segmented_subtypes or [],
            segmented_segments=segmented_segments,
            levels=levels,
            max_individual_traces=max_individual_traces,
            default_segment=segmented_default_segment,
        )

    annotation_path = ANNOTATIONS_DIR / annotation_name
    features = _load_annotation_features(annotation_path)
    max_len = _resolve_annotation_max_len(annotation_path, features)
    prefix_candidates = []
    for candidate in [data_prefix, strain_slug]:
        if candidate and candidate not in prefix_candidates:
            prefix_candidates.append(candidate)

    assets: dict[str, object] = {
        "available_views": [],
        "default_view": None,
        "all_file": None,
        "individual_file": None,
        "substrain_files": {},
    }
    cache_inputs: list[Path] = [Path(__file__), annotation_path]

    if "all" in levels:
        avg_path = _resolve_avg_path(prefix_candidates)
        if avg_path is not None:
            all_output = PILEUP_OUTPUT_DIR / f"{strain_slug}_all.html"
            if _asset_is_stale(all_output, [avg_path, annotation_path, Path(__file__)]):
                avg_df = _load_avg_depth_json(avg_path)
                fig = build_avg_pileup_figure(
                    avg_df, features, title="Pileup - All samples", max_len=max_len
                )
                _write_figure(fig, all_output)
            assets["all_file"] = f"pileup_html/{all_output.name}"
            assets["available_views"].append("all")
            cache_inputs.extend([avg_path, all_output])
        elif substrains:
            avg_by_substrain = {}
            input_paths: list[Path] = [annotation_path, Path(__file__)]
            for substrain in substrains:
                substrain_path = _resolve_avg_path(prefix_candidates, substrain)
                if substrain_path is None:
                    continue
                avg_by_substrain[substrain] = _load_avg_depth_json(substrain_path)
                input_paths.append(substrain_path)
            if avg_by_substrain:
                all_output = PILEUP_OUTPUT_DIR / f"{strain_slug}_all.html"
                if _asset_is_stale(all_output, input_paths):
                    fig = build_multi_avg_pileup_figure(
                        avg_by_substrain,
                        features,
                        title="Pileup - All samples",
                        max_len=max_len,
                    )
                    _write_figure(fig, all_output)
                assets["all_file"] = f"pileup_html/{all_output.name}"
                assets["available_views"].append("all")
                cache_inputs.extend(input_paths + [all_output])

    if "substrain" in levels:
        substrain_assets = {}
        for substrain in substrains:
            substrain_path = _resolve_avg_path(prefix_candidates, substrain)
            if substrain_path is None:
                continue
            clean = _slugify(substrain)
            output_path = PILEUP_OUTPUT_DIR / f"{strain_slug}_substrain_{clean}.html"
            if _asset_is_stale(output_path, [substrain_path, annotation_path, Path(__file__)]):
                avg_df = _load_avg_depth_json(substrain_path)
                fig = build_avg_pileup_figure(
                    avg_df, features, title=f"Pileup - {substrain}", max_len=max_len
                )
                _write_figure(fig, output_path)
            substrain_assets[substrain] = f"pileup_html/{output_path.name}"
            cache_inputs.extend([substrain_path, output_path])

        if substrain_assets:
            assets["substrain_files"] = substrain_assets
            assets["available_views"].append("substrain")

    if "individual" in levels:
        indiv_path = _resolve_individual_path(prefix_candidates)
        if indiv_path is not None:
            output_path = PILEUP_OUTPUT_DIR / f"{strain_slug}_individual_{max_individual_traces}.html"
            if _asset_is_stale(output_path, [indiv_path, annotation_path, Path(__file__)]):
                positions, sample_depths = _load_individual_depth_json(indiv_path)
                fig = build_individual_pileup_figure(
                    positions,
                    sample_depths,
                    features,
                    title="Pileup - Individual samples",
                    max_traces=max_individual_traces,
                    max_len=max_len,
                )
                _write_figure(fig, output_path)
            assets["individual_file"] = f"pileup_html/{output_path.name}"
            assets["available_views"].append("individual")
            cache_inputs.extend([indiv_path, output_path])
        elif substrains:
            per_substrain: dict[str, tuple[np.ndarray, dict[str, np.ndarray]]] = {}
            input_paths: list[Path] = [annotation_path, Path(__file__)]
            for substrain in substrains:
                substrain_indiv = _resolve_individual_path_for_substrain(prefix_candidates, substrain)
                if substrain_indiv is None:
                    continue
                per_substrain[substrain] = _load_individual_depth_json(substrain_indiv)
                input_paths.append(substrain_indiv)
            if per_substrain:
                output_path = PILEUP_OUTPUT_DIR / f"{strain_slug}_individual_{max_individual_traces}.html"
                if _asset_is_stale(output_path, input_paths):
                    positions, sample_depths = _merge_individual_traces(per_substrain)
                    fig = build_individual_pileup_figure(
                        positions,
                        sample_depths,
                        features,
                        title="Pileup - Individual samples",
                        max_traces=max_individual_traces,
                        max_len=max_len,
                    )
                    _write_figure(fig, output_path)
                assets["individual_file"] = f"pileup_html/{output_path.name}"
                assets["available_views"].append("individual")
                cache_inputs.extend(input_paths + [output_path])

    if assets["available_views"]:
        assets["default_view"] = "all" if "all" in assets["available_views"] else assets["available_views"][0]

    assets["enabled"] = bool(assets["available_views"])
    assets["version"] = _asset_version(cache_inputs)
    return assets
