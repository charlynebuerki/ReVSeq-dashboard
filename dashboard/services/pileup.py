from __future__ import annotations

import os
from pathlib import Path

import pandas as pd
from dna_features_viewer import BiopythonTranslator
from matplotlib import pyplot as plt

from .plots import pileup_plot


def _resolve_substrain_csv(strain_slug: str, substrain: str) -> str | None:
    clean_substrain = substrain.replace(" ", "_")
    base = Path("dashboard/static/pileup")

    candidates = [base / f"{strain_slug}_{clean_substrain}.csv"]

    # Backward-compatible source names for RSV subtype files.
    if strain_slug == "RSV":
        if clean_substrain == "RSV-A":
            candidates.append(base / "RSV_Respiratory_syncytial_virus_(type_A).csv")
        elif clean_substrain == "RSV-B":
            candidates.append(base / "RSV_Human_Respiratory_syncytial_virus_9320_(type_B).csv")

    for candidate in candidates:
        if candidate.is_file():
            return str(candidate)
    return None


def ensure_pileup_assets(strain_slug, strain_name_long, substrains, no_samples):
    match strain_name_long:
        case "Influenza A":
            if not os.path.isfile("dashboard/static/pileup/Influenza_A_PB1_all.png"):
                segments = ["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"]
                for seg in segments:
                    for substrain in ["all", "H1N1", "H3N2"]:
                        pileup_plot(
                            f"dashboard/static/pileup/Influenza_A_{seg}_{substrain}.csv",
                            f"dashboard/static/annotations/Influenza_A_{seg}.gb",
                            f"dashboard/static/pileup/Influenza_A_{seg}_{substrain}.png",
                        )
        case "Influenza B":
            if not os.path.isfile("dashboard/static/pileup/Influenza_B_PB1.png"):
                segments = ["PB2", "PB1", "HA", "NP", "NA", "MP", "NS"]
                for seg in segments:
                    pileup_plot(
                        f"dashboard/static/pileup/Influenza_B_{seg}_all.csv",
                        f"dashboard/static/annotations/Influenza_B_{seg}.gb",
                        f"dashboard/static/pileup/Influenza_B_{seg}_all.png",
                    )
        case _:
            if strain_slug == "Adenovirus":
                figsize = (20, 6.5)
                height_ratios = [4, 2.5]
            else:
                figsize = (20, 5)
                height_ratios = [4, 1]

            pileup_plot(
                f"dashboard/static/pileup/{strain_slug}_all.csv",
                f"dashboard/static/annotations/{strain_slug}.gb",
                f"dashboard/static/pileup/{strain_slug}_all.png",
                figsize=figsize,
                height_ratios=height_ratios,
            )

            if len(substrains) > 1:
                for substrain in substrains:
                    clean_substrain = substrain.replace(" ", "_")
                    source_csv = _resolve_substrain_csv(strain_slug, substrain)
                    if source_csv is None:
                        continue
                    pileup_plot(
                        source_csv,
                        f"dashboard/static/annotations/{strain_slug}.gb",
                        f"dashboard/static/pileup/{strain_slug}_{clean_substrain}.png",
                        figsize=figsize,
                        height_ratios=height_ratios,
                    )
            elif no_samples <= 25:
                coverage = pd.read_csv(f"dashboard/static/pileup/{strain_slug}_all_indiv.csv", index_col=0)
                fig, (ax1, ax2) = plt.subplots(
                    2,
                    1,
                    figsize=figsize,
                    sharex=True,
                    gridspec_kw={"height_ratios": height_ratios},
                )
                coverage.plot(ax=ax1)
                ax1.set_ylabel("Sequencing Depth", fontsize=10)
                ax1.set_title("Sequencing Depth", fontsize=15, loc="left", pad=20)
                ax1.set_yscale("log")
                ax1.set_ylim(ymin=1)
                ax1.get_xaxis().set_visible(False)
                ax1.axhline(y=10, color="b", linestyle="-", label="DP10")
                ax1.legend(title="Sample", bbox_to_anchor=(1, 1), loc="upper left")

                graphic_record = BiopythonTranslator().translate_record(
                    f"dashboard/static/annotations/{strain_slug}.gb"
                )
                graphic_record.plot(ax=ax2, strand_in_label_threshold=4, with_ruler=True)
                ax2.set_xlabel("Position")
                ax2.get_yaxis().set_visible(False)
                fig.savefig(f"dashboard/static/pileup/{strain_slug}_all_indiv.png", bbox_inches="tight")
