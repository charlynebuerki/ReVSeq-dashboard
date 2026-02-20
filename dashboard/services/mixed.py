from __future__ import annotations

import numpy as np
import pandas as pd

from .data import load_mixed_metadata
from .plots import build_coinfection_heatmap



def _coerce_pair_frame(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return pd.DataFrame(columns=["strain_1", "strain_2", "week", "canton"])
    missing = {"strain_1", "strain_2", "week"} - set(df.columns)
    if missing:
        return pd.DataFrame(columns=["strain_1", "strain_2", "week", "canton"])
    return df[["strain_1", "strain_2", "week", "canton"]].dropna(subset=["strain_1", "strain_2"])



def _build_coinfection_matrix(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return pd.DataFrame()

    strains = sorted(set(df["strain_1"]).union(set(df["strain_2"])))
    matrix = pd.DataFrame(0, index=strains, columns=strains, dtype=int)

    # Symmetric pair counts: (A,B) increments both [A,B] and [B,A].
    for row in df.itertuples(index=False):
        s1, s2 = row.strain_1, row.strain_2
        if s1 == s2:
            continue
        matrix.loc[s1, s2] += 1
        matrix.loc[s2, s1] += 1
    return matrix



def get_mixed_context():
    bundle = load_mixed_metadata()
    # Combine PCR and sequencing mixed calls into one co-infection matrix.
    mixed_df = pd.concat([_coerce_pair_frame(bundle.pcr), _coerce_pair_frame(bundle.sequencing)], ignore_index=True)

    co_infections = _build_coinfection_matrix(mixed_df)
    if co_infections.empty:
        return {
            "co_inf_mat": "",
            "strains": [],
            "pairs": [],
            "no_co_infections": 0,
            "show_pair_plots": False,
        }

    # Lower triangle sum gives unique pair totals (upper triangle masked).
    co_infections_np = co_infections.values.astype(float)
    co_infections_np[np.triu_indices(len(co_infections_np))] = np.nan
    no_co_infections = int(np.nansum(co_infections_np))

    co_inf_mat = build_coinfection_heatmap(co_infections_np, co_infections.columns.tolist())

    pairs = [
        {
            "first": co_infections.index[i],
            "second": co_infections.index[j],
            "first_clean": co_infections.index[i].replace(" ", "_").replace("/", "_"),
            "second_clean": co_infections.index[j].replace(" ", "_").replace("/", "_"),
        }
        for i, j in np.argwhere(np.array(co_infections) > 0)
        if j > i
    ]

    return {
        "co_inf_mat": co_inf_mat,
        "strains": co_infections.index.tolist(),
        "pairs": pairs,
        "no_co_infections": no_co_infections,
        "show_pair_plots": False,
    }
