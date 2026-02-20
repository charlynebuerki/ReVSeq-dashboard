from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import MagicMock, patch

import pandas as pd
from django.test import TestCase

from dashboard.config import harmonize_canton, harmonize_detected_strain
from dashboard.services import data as data_service
from dashboard.services.home import _aggregate_by_week, get_home_context
from dashboard.services.mixed import get_mixed_context
from dashboard.services.strain import get_strain_context


class HarmonizationTests(TestCase):
    def test_harmonize_detected_strain_aliases(self):
        self.assertEqual(
            harmonize_detected_strain("Respiratory syncytial virus (type A)", source="sequencing"),
            ("RSV - A/B", "RSV-A"),
        )
        self.assertEqual(
            harmonize_detected_strain("Respiratory syncytial virus (type A)", source="pcr"),
            ("RSV - A/B", "RSV - A/B"),
        )
        self.assertEqual(
            harmonize_detected_strain("Influenza A virus (A/Texas/50/2012(H3N2))", source="sequencing")[1],
            "H3N2",
        )
        self.assertEqual(
            harmonize_detected_strain("Influenza A", source="pcr")[1],
            "Influenza A",
        )
        self.assertEqual(
            harmonize_detected_strain("Human coronavirus OC43", source="pcr"),
            ("coronavirus OC43", "coronavirus OC43"),
        )
        self.assertEqual(harmonize_detected_strain("Unknown Virus", source="sequencing"), (None, None))

    def test_harmonize_canton_aliases(self):
        self.assertEqual(harmonize_canton("ZH"), "ZH")
        self.assertEqual(harmonize_canton("Zürich"), "ZH")
        self.assertEqual(harmonize_canton("Geneve"), "GE")
        self.assertEqual(harmonize_canton("Unknown"), None)


class MetadataPipelineTests(TestCase):
    def test_load_dashboard_metadata_keeps_pcr_substrain_canonical(self):
        with TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            pcr_file = tmp / "metadata_pcr.tsv"
            seq_file = tmp / "metadata_sequencing.tsv"

            pd.DataFrame(
                [
                    {
                        "strain": "sample-1",
                        "date": "2024-01-05",
                        "location": "ZH",
                        "Match_Sequencing": True,
                        "virus_identified_pcr": "RSV - A/B",
                    },
                    {
                        "strain": "sample-2",
                        "date": "2024-01-08",
                        "location": "BE",
                        "Match_Sequencing": False,
                        "virus_identified_pcr": "SARS-CoV-2",
                    },
                ]
            ).to_csv(pcr_file, sep="\t", index=False)

            pd.DataFrame(
                [
                    {
                        "strain": "sample-1",
                        "virus_identified": "Human Respiratory syncytial virus 9320 (type B)",
                        "date": "2024-01-05",
                        "location": "ZH",
                        "Match_PCR": True,
                    },
                    {
                        "strain": "sample-2",
                        "virus_identified": "SARS-CoV-2",
                        "date": "2024-01-08",
                        "location": "BE",
                        "Match_PCR": True,
                    },
                ]
            ).to_csv(seq_file, sep="\t", index=False)

            with patch.object(data_service, "PCR_METADATA_PATH", pcr_file), patch.object(
                data_service, "SEQ_METADATA_PATH", seq_file
            ):
                bundle = data_service.load_dashboard_metadata()

            self.assertEqual(set(bundle.pcr["canonical_strain"].unique()), {"RSV - A/B", "SARS-CoV-2"})
            rsv_row = bundle.pcr[bundle.pcr["canonical_strain"] == "RSV - A/B"].iloc[0]
            self.assertEqual(rsv_row["substrain"], "RSV - A/B")
            self.assertEqual(rsv_row["display_label"], "RSV - A/B")

            seq_rsv_row = bundle.sequencing[bundle.sequencing["canonical_strain"] == "RSV - A/B"].iloc[0]
            self.assertEqual(seq_rsv_row["display_label"], "RSV-B")


class AggregationTests(TestCase):
    def test_home_aggregate_by_week(self):
        df = pd.DataFrame(
            [
                {"week": "2024/01", "canonical_strain": "RSV - A/B"},
                {"week": "2024/01", "canonical_strain": "RSV - A/B"},
                {"week": "2024/01", "canonical_strain": "Influenza A"},
                {"week": "2024/02", "canonical_strain": "Influenza A"},
            ]
        )
        grouped = _aggregate_by_week(df)

        rsv_week = grouped[(grouped["week"] == "2024/01") & (grouped["canonical_strain"] == "RSV - A/B")]
        inf_week = grouped[(grouped["week"] == "2024/02") & (grouped["canonical_strain"] == "Influenza A")]
        self.assertEqual(int(rsv_week.iloc[0]["count"]), 2)
        self.assertEqual(int(inf_week.iloc[0]["count"]), 1)


class StrainContextTests(TestCase):
    def test_strain_context_uses_separate_inputs_for_pcr_and_sequencing(self):
        with TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            pcr_file = tmp / "metadata_pcr.tsv"
            seq_file = tmp / "metadata_sequencing.tsv"

            pd.DataFrame(
                [
                    {
                        "strain": "sample-a",
                        "date": "2024-01-01",
                        "location": "ZH",
                        "Match_Sequencing": True,
                        "virus_identified_pcr": "RSV - A/B",
                    },
                    {
                        "strain": "sample-b",
                        "date": "2024-01-02",
                        "location": "BE",
                        "Match_Sequencing": True,
                        "virus_identified_pcr": "RSV - A/B",
                    },
                ]
            ).to_csv(pcr_file, sep="\t", index=False)

            pd.DataFrame(
                [
                    {
                        "strain": "sample-a",
                        "virus_identified": "Respiratory syncytial virus (type A)",
                        "date": "2024-01-01",
                        "location": "ZH",
                        "Match_PCR": True,
                    },
                    {
                        "strain": "sample-b",
                        "virus_identified": "Human Respiratory syncytial virus 9320 (type B)",
                        "date": "2024-01-02",
                        "location": "BE",
                        "Match_PCR": True,
                    },
                ]
            ).to_csv(seq_file, sep="\t", index=False)

            strain_config = {
                "data_name": "RSV - A/B",
                "modules": {
                    "histogram_pcr": True,
                    "histogram_sequencing": True,
                    "map": True,
                    "pileup": False,
                    "tree": False,
                },
                "trees": [],
            }

            with patch.object(data_service, "PCR_METADATA_PATH", pcr_file), patch.object(
                data_service, "SEQ_METADATA_PATH", seq_file
            ), patch("dashboard.services.strain._ensure_strain_histogram") as histogram_mock, patch(
                "dashboard.services.strain.build_weekly_canton_map",
                return_value=MagicMock(write_html=lambda *_args, **_kwargs: None),
            ):
                context = get_strain_context("RSV", strain_config, "localhost:8000")

            self.assertEqual(context["strain_name"], "RSV - A/B")
            self.assertEqual(context["no_samples"], 2)
            self.assertEqual(histogram_mock.call_count, 2)

            pcr_df = histogram_mock.call_args_list[0].args[0]
            seq_df = histogram_mock.call_args_list[1].args[0]
            self.assertTrue((pcr_df["canonical_strain"] == "RSV - A/B").all())
            self.assertTrue((seq_df["canonical_strain"] == "RSV - A/B").all())


class SampleCountTests(TestCase):
    def test_home_sample_count_uses_sequencing_only(self):
        with TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            pcr_file = tmp / "metadata_pcr.tsv"
            seq_file = tmp / "metadata_sequencing.tsv"

            pd.DataFrame(
                [
                    {
                        "strain": "sample-seq",
                        "date": "2024-01-01",
                        "location": "ZH",
                        "Match_Sequencing": True,
                        "virus_identified_pcr": "RSV - A/B",
                    },
                    {
                        "strain": "sample-pcr-only",
                        "date": "2024-01-02",
                        "location": "BE",
                        "Match_Sequencing": False,
                        "virus_identified_pcr": "RSV - A/B",
                    },
                ]
            ).to_csv(pcr_file, sep="\t", index=False)

            pd.DataFrame(
                [
                    {
                        "strain": "sample-seq",
                        "virus_identified": "Respiratory syncytial virus (type A)",
                        "date": "2024-01-01",
                        "location": "ZH",
                        "Match_PCR": True,
                    }
                ]
            ).to_csv(seq_file, sep="\t", index=False)

            with patch.object(data_service, "PCR_METADATA_PATH", pcr_file), patch.object(
                data_service, "SEQ_METADATA_PATH", seq_file
            ), patch("dashboard.services.home.build_home_assets"):
                context = get_home_context()

            self.assertEqual(context["no_samples"], 1)


class MixedContextTests(TestCase):
    def test_mixed_context_reads_new_mixed_files(self):
        with TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            mixed_pcr = tmp / "mixed_pcr.tsv"
            mixed_seq = tmp / "mixed_sequencing.tsv"

            pd.DataFrame(
                [
                    {
                        "strain_1": "RSV - A/B",
                        "strain_2": "SARS-CoV-2",
                        "date": "2024-01-01",
                        "location": "ZH",
                        "Match_Sequencing": True,
                    }
                ]
            ).to_csv(mixed_pcr, sep="\t", index=False)

            pd.DataFrame(
                [
                    {
                        "strain_1": "RSV - A/B",
                        "strain_2": "Influenza A",
                        "date": "2024-01-01",
                        "location": "ZH",
                        "Match_PCR": True,
                    }
                ]
            ).to_csv(mixed_seq, sep="\t", index=False)

            with patch.object(data_service, "MIXED_PCR_PATH", mixed_pcr), patch.object(
                data_service, "MIXED_SEQ_PATH", mixed_seq
            ), patch("dashboard.services.mixed.build_coinfection_heatmap", return_value="<div>heatmap</div>"):
                context = get_mixed_context()

            self.assertEqual(context["co_inf_mat"], "<div>heatmap</div>")
            self.assertEqual(context["no_co_infections"], 2)
            self.assertGreaterEqual(len(context["pairs"]), 2)
