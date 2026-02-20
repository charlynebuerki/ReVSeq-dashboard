from __future__ import annotations

import json
from functools import lru_cache
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from dna_features_viewer import BiopythonTranslator
from matplotlib import pyplot as plt

from dashboard.config import COLOR_BY_STRAIN


def pileup_plot(coverage_file, annotation_file, out_file, figsize=(20, 5), height_ratios=(4, 1.5)):
    coverage = pd.read_csv(coverage_file)
    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=figsize, sharex=True, gridspec_kw={"height_ratios": height_ratios}
    )
    ax1.plot(coverage["idx"], coverage["mean"])
    ax1.fill_between(coverage["idx"], coverage["ci_lower"], coverage["ci_upper"], color="b", alpha=0.15)
    ax1.set_ylabel("Sequencing Depth", fontsize=10)
    ax1.set_title("Sequencing Depth", fontsize=15, loc="left", pad=20)
    ax1.set_yscale("log")
    ax1.set_ylim(ymin=1)
    ax1.get_xaxis().set_visible(False)
    ax1.axhline(y=10, color="b", linestyle="-", label="DP10")
    ax1.legend(bbox_to_anchor=(1, 1), loc="upper left")

    graphic_record = BiopythonTranslator().translate_record(annotation_file)
    graphic_record.plot(ax=ax2, strand_in_label_threshold=4, with_ruler=True)
    ax2.set_xlabel("Position")
    ax2.get_yaxis().set_visible(False)
    fig.savefig(out_file, bbox_inches="tight")



def make_weekly_strain_figure(df: pd.DataFrame, strain_col: str, title: str):
    fig = px.bar(
        df,
        x="week_start",
        y="count",
        color=strain_col,
        barmode="stack",
        labels={"week_start": "Week", "count": "No. infections", strain_col: "Strain"},
        color_discrete_map=COLOR_BY_STRAIN,
        template="simple_white",
        title=title,
    )
    # Weekly data; month ticks/grid for readability on long timelines.
    fig.update_xaxes(
        tickformat="%m-%y",
        dtick="M1",
        ticklabelmode="period",
        showgrid=True,
        minor=dict(showgrid=False),
    )
    fig.update_yaxes(showgrid=True)
    return fig



def make_weekly_substrain_figure(df: pd.DataFrame, substrain_col: str, title: str):
    if df[substrain_col].nunique() > 1:
        fig = px.bar(
            df,
            x="week_start",
            y="count",
            color=substrain_col,
            barmode="stack",
            labels={"week_start": "Week", "count": "No. infections", substrain_col: "Substrain"},
            template="simple_white",
            title=title,
        )
    else:
        fig = px.bar(
            df,
            x="week_start",
            y="count",
            labels={"week_start": "Week", "count": "No. infections"},
            template="simple_white",
            title=title,
        )

    # Keep the same axis semantics for strain and substrain views.
    fig.update_xaxes(
        tickformat="%m-%y",
        dtick="M1",
        ticklabelmode="period",
        showgrid=True,
        minor=dict(showgrid=False),
    )
    fig.update_yaxes(showgrid=True)
    fig.update_layout(legend=dict(xanchor="left", x=0, yanchor="bottom", y=1))
    return fig



@lru_cache(maxsize=1)
def _load_geojson(geojson_path: str):
    with open(geojson_path, "r", encoding="utf-8") as fh:
        geojson = json.load(fh)
    canton_names = {f["properties"]["canton"]: f["properties"]["NAME"] for f in geojson["features"]}
    canton_codes = sorted(canton_names.keys())
    return geojson, canton_codes, canton_names


def build_weekly_canton_map(df: pd.DataFrame, geojson_path: str = "dashboard/static/swiss_cantons.geojson"):
    # Return an explicit empty-state chart so templates can always render an iframe.
    if df.empty:
        fig = go.Figure()
        fig.update_layout(
            template="simple_white",
            title="No data available",
            xaxis={"visible": False},
            yaxis={"visible": False},
            annotations=[{"text": "No mapped samples for this selection", "xref": "paper", "yref": "paper", "x": 0.5, "y": 0.5, "showarrow": False}],
        )
        return fig

    geojson, canton_codes, canton_names = _load_geojson(geojson_path)
    data = df[["week", "canton"]].dropna().copy()
    data = data[data["canton"].isin(canton_codes)]
    if data.empty:
        return build_weekly_canton_map(pd.DataFrame())

    counts = data.groupby(["week", "canton"]).size().reset_index(name="No. infections")
    weeks = sorted(counts["week"].unique().tolist())
    full_index = pd.MultiIndex.from_product([weeks, canton_codes], names=["week", "canton"])
    counts = counts.set_index(["week", "canton"]).reindex(full_index, fill_value=0).reset_index()
    counts["NAME"] = counts["canton"].map(canton_names)

    max_count = int(counts["No. infections"].max()) if not counts.empty else 0
    fig = px.choropleth_mapbox(
        counts,
        geojson=geojson,
        locations="canton",
        featureidkey="properties.canton",
        hover_name="NAME",
        zoom=5.5,
        center={"lat": 46.8, "lon": 8.3},
        color="No. infections",
        color_continuous_scale="viridis",
        animation_frame="week",
        mapbox_style="carto-positron",
        range_color=[0, max(max_count, 1)],
    )
    fig.update_layout(
        margin={"l": 0, "r": 0, "t": 30, "b": 0},
        coloraxis_colorbar={"title": "Count"},
    )
    return fig



def build_coinfection_heatmap(co_infections_np, labels):
    fig = go.Figure(
        go.Heatmap(
            z=co_infections_np[1:, :-1][::-1],
            x=labels[:-1],
            y=labels[1:][::-1],
            colorscale="Viridis",
        )
    )
    fig.update_layout({"paper_bgcolor": "rgba(0, 0, 0, 0)", "plot_bgcolor": "rgba(0, 0, 0, 0)"})
    fig = fig.update_traces(
        text=np.nan_to_num(co_infections_np[1:, :-1][::-1], nan=0).astype(int).astype(str),
        texttemplate="%{text}",
        hovertemplate=None,
        showscale=False,
    )
    from plotly.offline import plot

    return plot(fig, output_type="div")
