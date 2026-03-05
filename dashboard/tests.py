from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import MagicMock, patch

import pandas as pd
from django.test import TestCase

from dashboard import config as dashboard_config
from dashboard.config import harmonize_canton, harmonize_detected_strain
from dashboard.services import data as data_service
from dashboard.services.home import _aggregate_by_week, get_home_context
from dashboard.services.mixed import get_mixed_context
from dashboard.services import pileup as pileup_service
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
                {"week": "2024/01", "canonical_strain": "RSV - A/B", "match_pcr": 1},
                {"week": "2024/01", "canonical_strain": "RSV - A/B", "match_pcr": 0},
                {"week": "2024/01", "canonical_strain": "Influenza A", "match_pcr": 1},
                {"week": "2024/02", "canonical_strain": "Influenza A", "match_pcr": 1},
            ]
        )
        grouped = _aggregate_by_week(df, match_col="match_pcr")

        rsv_week = grouped[(grouped["week"] == "2024/01") & (grouped["canonical_strain"] == "RSV - A/B")]
        inf_week = grouped[(grouped["week"] == "2024/02") & (grouped["canonical_strain"] == "Influenza A")]
        self.assertEqual(int(rsv_week.iloc[0]["count"]), 2)
        self.assertEqual(int(inf_week.iloc[0]["count"]), 1)
        self.assertEqual(int(rsv_week.iloc[0]["match_count"]), 1)
        self.assertEqual(int(inf_week.iloc[0]["match_count"]), 1)
        self.assertEqual(rsv_week.iloc[0]["match_over_total"], "1/2")
        self.assertEqual(inf_week.iloc[0]["match_over_total"], "1/1")


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
                    {
                        "strain": "sample-a",
                        "date": "2024-01-03",
                        "location": "ZH",
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
                    "barplot_pcr": True,
                    "barplot_sequencing": True,
                    "map": True,
                    "pileup": False,
                    "tree": False,
                },
                "trees": [],
            }

            with patch.object(data_service, "PCR_METADATA_PATH", pcr_file), patch.object(
                data_service, "SEQ_METADATA_PATH", seq_file
            ), patch("dashboard.services.strain._ensure_strain_barplot") as histogram_mock, patch(
                "dashboard.services.strain.build_weekly_canton_map",
                return_value=MagicMock(write_html=lambda *_args, **_kwargs: None),
            ):
                context = get_strain_context("RSV", strain_config)

            self.assertEqual(context["strain_name"], "RSV - A/B")
            self.assertEqual(context["no_sequences"], 2)
            self.assertEqual(context["no_detections"], 3)
            self.assertEqual(histogram_mock.call_count, 2)

            pcr_df = histogram_mock.call_args_list[0].args[0]
            seq_df = histogram_mock.call_args_list[1].args[0]
            self.assertTrue((pcr_df["canonical_strain"] == "RSV - A/B").all())
            self.assertTrue((seq_df["canonical_strain"] == "RSV - A/B").all())


class SampleCountTests(TestCase):
    def test_home_sample_counts_use_sequencing_unique_and_pcr_rows(self):
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
                        "strain": "sample-seq",
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

            self.assertEqual(context["no_sequences"], 1)
            self.assertEqual(context["no_detections"], 2)


class MixedContextTests(TestCase):
    def test_mixed_context_builds_source_specific_pair_counts(self):
        with TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            mixed_pcr = tmp / "metadata_co_infection_pcr.tsv"
            mixed_seq = tmp / "metadata_co_infection_sequencing.tsv"

            pd.DataFrame(
                [
                    {
                        "strain": "sample-pcr-1",
                        "virus_identified_pcr": "SARS-CoV-2",
                        "date": "2024-01-01",
                        "location": "ZH",
                        "Match_sequencing": True,
                    },
                    {
                        "strain": "sample-pcr-1",
                        "virus_identified_pcr": "Influenza A",
                        "date": "2024-01-01",
                        "location": "ZH",
                        "Match_sequencing": True,
                    },
                    {
                        "strain": "sample-pcr-2",
                        "virus_identified_pcr": "SARS-CoV-2",
                        "date": "2024-01-02",
                        "location": "BE",
                        "Match_sequencing": False,
                    },
                    {
                        "strain": "sample-pcr-2",
                        "virus_identified_pcr": "Influenza A",
                        "date": "2024-01-02",
                        "location": "BE",
                        "Match_sequencing": False,
                    },
                    {
                        "strain": "sample-pcr-2",
                        "virus_identified_pcr": "coronavirus OC43",
                        "date": "2024-01-02",
                        "location": "BE",
                        "Match_sequencing": False,
                    },
                ]
            ).to_csv(mixed_pcr, sep="\t", index=False)

            pd.DataFrame(
                [
                    {
                        "strain": "sample-seq-1",
                        "virus_identified": "SARS-CoV-2",
                        "date": "2024-01-01",
                        "location": "ZH",
                        "Match_PCR": True,
                    },
                    {
                        "strain": "sample-seq-1",
                        "virus_identified": "Influenza A",
                        "date": "2024-01-01",
                        "location": "ZH",
                        "Match_PCR": True,
                    },
                    {
                        "strain": "sample-seq-2",
                        "virus_identified": "SARS-CoV-2",
                        "date": "2024-01-03",
                        "location": "BS",
                        "Match_PCR": True,
                    },
                    {
                        "strain": "sample-seq-2",
                        "virus_identified": "Influenza A",
                        "date": "2024-01-03",
                        "location": "BS",
                        "Match_PCR": True,
                    },
                    {
                        "strain": "sample-seq-2",
                        "virus_identified": "coronavirus OC43",
                        "date": "2024-01-03",
                        "location": "BS",
                        "Match_PCR": True,
                    },
                ]
            ).to_csv(mixed_seq, sep="\t", index=False)

            with patch.object(data_service, "MIXED_PCR_PATH", mixed_pcr), patch.object(
                data_service, "MIXED_SEQ_PATH", mixed_seq
            ), patch("dashboard.services.mixed.build_coinfection_composite_figure", return_value="<div>heatmap</div>"), patch(
                "dashboard.services.mixed.build_mixed_sample_pileup_assets",
                return_value=[{"label": "Influenza A", "available": False, "message": "not available"}],
            ) as pileup_builder:
                context = get_mixed_context(selected_sample="sample-seq-2")

            self.assertEqual(context["co_inf_mat_seq"], "<div>heatmap</div>")
            self.assertEqual(context["no_co_infections_seq"], 4)
            self.assertEqual(context["no_co_infections_pcr"], 4)
            self.assertIn(context["default_source"], {"sequencing", "pcr"})
            self.assertEqual(context["mixed_pileup_selected_sample"], "sample-seq-2")
            self.assertIn("sample-seq-1", context["mixed_pileup_sample_options"])
            self.assertIn("sample-seq-2", context["mixed_pileup_sample_options"])
            pileup_builder.assert_called_once()
            self.assertEqual(pileup_builder.call_args[0][0], "sample-seq-2")


class MixedPageConfigTests(TestCase):
    def test_mixed_page_can_be_disabled_from_config(self):
        with patch.object(dashboard_config, "MIXED_PAGE_ENABLED", False):
            options = dashboard_config.get_strain_options()
            self.assertTrue(all(opt["value"] != "/mixed" for opt in options))

            response = self.client.get("/mixed/")
            self.assertEqual(response.status_code, 404)

    def test_mixed_pileup_data_endpoint_returns_json(self):
        payload = {
            "mixed_pileup_sample_options": ["sample-a"],
            "mixed_pileup_selected_sample": "sample-a",
            "mixed_pileup_items": [{"label": "SARS-CoV-2", "available": False, "message": "not available"}],
        }
        with patch("dashboard.views.get_mixed_pileup_context", return_value=payload):
            response = self.client.get("/mixed/pileup/?sample=sample-a")
        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.json(), payload)


class PileupTests(TestCase):
    def _write_genbank(self, path: Path):
        content = """LOCUS       TESTSEQ                 20 bp    DNA     linear   SYN 01-JAN-1980
FEATURES             Location/Qualifiers
     gene            1..8
                     /gene="GeneA"
     CDS             9..16
                     /gene="GeneB"
ORIGIN
        1 atgcatgcat gcatgcatgc
//
"""
        path.write_text(content, encoding="utf-8")

    def test_build_pileup_context_all_substrain_individual(self):
        with TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            data_dir = tmp / "data"
            ann_dir = tmp / "annotations"
            out_dir = tmp / "out"
            data_dir.mkdir(parents=True)
            ann_dir.mkdir(parents=True)
            out_dir.mkdir(parents=True)

            (data_dir / "RSV_avg.json").write_text(
                '{"positions":[1,2,3],"median":[10,20,30],"q1":[5,10,20],"q3":[15,25,40]}',
                encoding="utf-8",
            )
            (data_dir / "RSV_RSV-A_avg.json").write_text(
                '{"positions":[1,2,3],"median":[8,15,25],"q1":[4,8,12],"q3":[12,20,30]}',
                encoding="utf-8",
            )
            (data_dir / "RSV_indiv.json").write_text(
                '{"positions":[1,2,3],"samples":{"sample-1":[10,20,30],"sample-2":[5,8,12]}}',
                encoding="utf-8",
            )
            self._write_genbank(ann_dir / "RSV.gb")

            with patch.object(pileup_service, "PILEUP_DATA_DIR", data_dir), patch.object(
                pileup_service, "ANNOTATIONS_DIR", ann_dir
            ), patch.object(pileup_service, "PILEUP_OUTPUT_DIR", out_dir):
                context = pileup_service.build_pileup_context(
                    strain_slug="RSV",
                    strain_name="RSV - A/B",
                    data_prefix="RSV",
                    substrains=["RSV-A", "RSV-B"],
                    annotation_name="RSV.gb",
                    levels=["all", "substrain", "individual"],
                    max_individual_traces=10,
                )

            self.assertTrue(context["enabled"])
            self.assertIn("all", context["available_views"])
            self.assertIn("substrain", context["available_views"])
            self.assertIn("individual", context["available_views"])
            self.assertIn("RSV-A", context["substrain_files"])
            self.assertTrue((out_dir / "RSV_all.html").exists())
            self.assertTrue((out_dir / "RSV_substrain_RSV-A.html").exists())
            self.assertTrue((out_dir / "RSV_individual_10.html").exists())

            individual_html = (out_dir / "RSV_individual_10.html").read_text(encoding="utf-8")
            self.assertIn("sample-1", individual_html)

    def test_build_pileup_context_respects_configured_levels(self):
        with TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            data_dir = tmp / "data"
            ann_dir = tmp / "annotations"
            out_dir = tmp / "out"
            data_dir.mkdir(parents=True)
            ann_dir.mkdir(parents=True)
            out_dir.mkdir(parents=True)

            (data_dir / "SARS-CoV-2_avg.json").write_text(
                '{"positions":[1,2,3],"median":[10,20,30],"q1":[5,10,20],"q3":[15,25,40]}',
                encoding="utf-8",
            )
            self._write_genbank(ann_dir / "SARS-CoV-2.gb")

            with patch.object(pileup_service, "PILEUP_DATA_DIR", data_dir), patch.object(
                pileup_service, "ANNOTATIONS_DIR", ann_dir
            ), patch.object(pileup_service, "PILEUP_OUTPUT_DIR", out_dir):
                context = pileup_service.build_pileup_context(
                    strain_slug="SARS-CoV-2",
                    strain_name="SARS-CoV-2",
                    data_prefix="SARS-CoV-2",
                    substrains=[],
                    annotation_name="SARS-CoV-2.gb",
                    levels=["all"],
                    max_individual_traces=10,
                )

            self.assertTrue(context["enabled"])
            self.assertEqual(context["available_views"], ["all"])
            self.assertEqual(context["default_view"], "all")
            self.assertIsNone(context["individual_file"])

    def test_build_pileup_context_resolves_direct_substrain_file_names(self):
        with TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            data_dir = tmp / "data"
            ann_dir = tmp / "annotations"
            out_dir = tmp / "out"
            data_dir.mkdir(parents=True)
            ann_dir.mkdir(parents=True)
            out_dir.mkdir(parents=True)

            # Direct substrain naming used in real RSV inputs.
            (data_dir / "rsv-a_avg.json").write_text(
                '{"positions":[1,2,3],"median":[8,15,25],"q1":[4,8,12],"q3":[12,20,30]}',
                encoding="utf-8",
            )
            self._write_genbank(ann_dir / "RSV.gb")

            with patch.object(pileup_service, "PILEUP_DATA_DIR", data_dir), patch.object(
                pileup_service, "ANNOTATIONS_DIR", ann_dir
            ), patch.object(pileup_service, "PILEUP_OUTPUT_DIR", out_dir):
                context = pileup_service.build_pileup_context(
                    strain_slug="RSV",
                    strain_name="RSV - A/B",
                    data_prefix="RSV",
                    substrains=["RSV-A"],
                    annotation_name="RSV.gb",
                    levels=["substrain"],
                    max_individual_traces=10,
                )

            self.assertTrue(context["enabled"])
            self.assertIn("substrain", context["available_views"])
            self.assertIn("RSV-A", context["substrain_files"])

    def test_build_pileup_context_rebuilds_individual_for_new_trace_cap(self):
        with TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            data_dir = tmp / "data"
            ann_dir = tmp / "annotations"
            out_dir = tmp / "out"
            data_dir.mkdir(parents=True)
            ann_dir.mkdir(parents=True)
            out_dir.mkdir(parents=True)

            (data_dir / "RSV_indiv.json").write_text(
                '{"positions":[1,2,3],"samples":{"sample-1":[10,20,30],"sample-2":[5,8,12],"sample-3":[7,9,11]}}',
                encoding="utf-8",
            )
            self._write_genbank(ann_dir / "RSV.gb")

            with patch.object(pileup_service, "PILEUP_DATA_DIR", data_dir), patch.object(
                pileup_service, "ANNOTATIONS_DIR", ann_dir
            ), patch.object(pileup_service, "PILEUP_OUTPUT_DIR", out_dir):
                ctx_3 = pileup_service.build_pileup_context(
                    strain_slug="RSV",
                    strain_name="RSV - A/B",
                    data_prefix="RSV",
                    substrains=[],
                    annotation_name="RSV.gb",
                    levels=["individual"],
                    max_individual_traces=3,
                )
                ctx_1 = pileup_service.build_pileup_context(
                    strain_slug="RSV",
                    strain_name="RSV - A/B",
                    data_prefix="RSV",
                    substrains=[],
                    annotation_name="RSV.gb",
                    levels=["individual"],
                    max_individual_traces=1,
                )

            self.assertNotEqual(ctx_3["individual_file"], ctx_1["individual_file"])
            html_cap_1 = (out_dir / "RSV_individual_1.html").read_text(encoding="utf-8")
            self.assertIn("Showing first 1 of 3 samples", html_cap_1)
            self.assertIn("Sample=sample-1", html_cap_1)
            self.assertNotIn("Sample=sample-2", html_cap_1)

    def test_build_pileup_context_builds_all_from_substrain_avg_files(self):
        with TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            data_dir = tmp / "data"
            ann_dir = tmp / "annotations"
            out_dir = tmp / "out"
            data_dir.mkdir(parents=True)
            ann_dir.mkdir(parents=True)
            out_dir.mkdir(parents=True)

            (data_dir / "rsv-a_avg.json").write_text(
                '{"positions":[1,2,3],"median":[8,15,25],"q1":[4,8,12],"q3":[12,20,30]}',
                encoding="utf-8",
            )
            (data_dir / "rsv-b_avg.json").write_text(
                '{"positions":[1,2,3],"median":[6,12,18],"q1":[3,6,9],"q3":[9,16,24]}',
                encoding="utf-8",
            )
            self._write_genbank(ann_dir / "RSV.gb")

            with patch.object(pileup_service, "PILEUP_DATA_DIR", data_dir), patch.object(
                pileup_service, "ANNOTATIONS_DIR", ann_dir
            ), patch.object(pileup_service, "PILEUP_OUTPUT_DIR", out_dir):
                context = pileup_service.build_pileup_context(
                    strain_slug="RSV",
                    strain_name="RSV - A/B",
                    data_prefix="RSV",
                    substrains=["RSV-A", "RSV-B"],
                    annotation_name="RSV.gb",
                    levels=["all"],
                    max_individual_traces=10,
                )

            self.assertTrue(context["enabled"])
            self.assertIn("all", context["available_views"])
            self.assertEqual(context["all_file"], "pileup_html/RSV_all.html")
            all_html = (out_dir / "RSV_all.html").read_text(encoding="utf-8")
            self.assertIn('"name":"RSV-A"', all_html)
            self.assertIn('"name":"RSV-B"', all_html)

    def test_build_pileup_context_builds_individual_from_substrain_indiv_files(self):
        with TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            data_dir = tmp / "data"
            ann_dir = tmp / "annotations"
            out_dir = tmp / "out"
            data_dir.mkdir(parents=True)
            ann_dir.mkdir(parents=True)
            out_dir.mkdir(parents=True)

            (data_dir / "rsv-a_indiv.json").write_text(
                '{"positions":[1,2,3],"samples":{"sample-a1":[10,20,30],"sample-a2":[8,12,14]}}',
                encoding="utf-8",
            )
            (data_dir / "rsv-b_indiv.json").write_text(
                '{"positions":[1,2,3],"samples":{"sample-b1":[6,8,11]}}',
                encoding="utf-8",
            )
            self._write_genbank(ann_dir / "RSV.gb")

            with patch.object(pileup_service, "PILEUP_DATA_DIR", data_dir), patch.object(
                pileup_service, "ANNOTATIONS_DIR", ann_dir
            ), patch.object(pileup_service, "PILEUP_OUTPUT_DIR", out_dir):
                context = pileup_service.build_pileup_context(
                    strain_slug="RSV",
                    strain_name="RSV - A/B",
                    data_prefix="RSV",
                    substrains=["RSV-A", "RSV-B"],
                    annotation_name="RSV.gb",
                    levels=["individual"],
                    max_individual_traces=10,
                )

            self.assertTrue(context["enabled"])
            self.assertIn("individual", context["available_views"])
            self.assertEqual(context["individual_file"], "pileup_html/RSV_individual_10.html")
            indiv_html = (out_dir / "RSV_individual_10.html").read_text(encoding="utf-8")
            self.assertIn("Sample=RSV-A:sample-a1", indiv_html)
            self.assertIn("Sample=RSV-B:sample-b1", indiv_html)

    def test_build_pileup_context_segmented_subtype_segment_files(self):
        with TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            data_dir = tmp / "data"
            ann_dir = tmp / "annotations"
            out_dir = tmp / "out"
            data_dir.mkdir(parents=True)
            ann_dir.mkdir(parents=True)
            out_dir.mkdir(parents=True)

            (data_dir / "flu-a-h1n1_avg.json").write_text(
                (
                    '{"virus":"flu-a-h1n1","is_segmented":true,"segments":'
                    '{"HA":{"positions":[1,2,3],"median":[8,15,25],"q1":[4,8,12],"q3":[12,20,30]},'
                    '"PB1":{"positions":[1,2,3],"median":[9,14,19],"q1":[5,8,11],"q3":[13,19,24]}}}'
                ),
                encoding="utf-8",
            )
            (data_dir / "flu-a-h3n2_avg.json").write_text(
                (
                    '{"virus":"flu-a-h3n2","is_segmented":true,"segments":'
                    '{"HA":{"positions":[1,2,3],"median":[7,13,21],"q1":[3,7,10],"q3":[11,18,28]},'
                    '"PB1":{"positions":[1,2,3],"median":[10,12,17],"q1":[6,7,10],"q3":[14,16,22]}}}'
                ),
                encoding="utf-8",
            )
            self._write_genbank(ann_dir / "Influenza_A_HA.gb")
            self._write_genbank(ann_dir / "Influenza_A_PB1.gb")

            with patch.object(pileup_service, "PILEUP_DATA_DIR", data_dir), patch.object(
                pileup_service, "ANNOTATIONS_DIR", ann_dir
            ), patch.object(pileup_service, "PILEUP_OUTPUT_DIR", out_dir):
                context = pileup_service.build_pileup_context(
                    strain_slug="Influenza_A",
                    strain_name="Influenza A",
                    data_prefix="Influenza_A",
                    substrains=[],
                    annotation_name="Influenza_A_HA.gb",
                    levels=["all", "substrain", "individual"],
                    max_individual_traces=10,
                    segmented_subtypes=[
                        {"value": "h1n1", "label": "H1N1", "data_prefix": "flu-a-h1n1"},
                        {"value": "h3n2", "label": "H3N2", "data_prefix": "flu-a-h3n2"},
                    ],
                    segmented_segments=[
                        {"value": "HA", "label": "HA", "annotation": "Influenza_A_HA.gb"},
                        {"value": "PB1", "label": "PB1", "annotation": "Influenza_A_PB1.gb"},
                    ],
                )

            self.assertTrue(context["enabled"])
            self.assertEqual(context.get("mode"), "segmented")
            self.assertIn("all", context.get("available_views", []))
            self.assertIn("substrain", context.get("available_views", []))
            self.assertEqual(context.get("default_view"), "all")
            self.assertIn("HA", context.get("segmented_all_files", {}))
            self.assertIn("HA", context.get("segmented_substrain_files", {}))
            self.assertIn("h1n1", context["segmented_substrain_files"]["HA"])
            self.assertTrue((out_dir / "Influenza_A_seg_all_HA.html").exists())
            self.assertTrue((out_dir / "Influenza_A_seg_sub_h1n1_HA.html").exists())

    def test_build_pileup_context_segmented_auto_detect_subtypes_and_default_segment(self):
        with TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            data_dir = tmp / "data"
            ann_dir = tmp / "annotations"
            out_dir = tmp / "out"
            data_dir.mkdir(parents=True)
            ann_dir.mkdir(parents=True)
            out_dir.mkdir(parents=True)

            segmented_payload = (
                '{"virus":"flu-a-h1n1","is_segmented":true,"segments":'
                '{"HA":{"positions":[1,2],"median":[8,9],"q1":[4,5],"q3":[12,13]},'
                '"PB1":{"positions":[1,2],"median":[6,7],"q1":[3,4],"q3":[9,10]}}}'
            )
            (data_dir / "flu-a-h1n1_avg.json").write_text(segmented_payload, encoding="utf-8")
            (data_dir / "flu-a-h3n2_avg.json").write_text(segmented_payload, encoding="utf-8")
            self._write_genbank(ann_dir / "Influenza_A_HA.gb")
            self._write_genbank(ann_dir / "Influenza_A_PB1.gb")

            with patch.object(pileup_service, "PILEUP_DATA_DIR", data_dir), patch.object(
                pileup_service, "ANNOTATIONS_DIR", ann_dir
            ), patch.object(pileup_service, "PILEUP_OUTPUT_DIR", out_dir):
                context = pileup_service.build_pileup_context(
                    strain_slug="Influenza_A",
                    strain_name="Influenza A",
                    data_prefix="flu-a",
                    substrains=[],
                    annotation_name="Influenza_A_HA.gb",
                    levels=["all", "substrain", "individual"],
                    max_individual_traces=10,
                    segmented_segments=[
                        {"value": "PB1", "label": "PB1", "annotation": "Influenza_A_PB1.gb"},
                        {"value": "HA", "label": "HA", "annotation": "Influenza_A_HA.gb"},
                    ],
                    segmented_default_segment="HA",
                )

            subtype_values = [x["value"] for x in context.get("segmented_subtypes", [])]
            self.assertIn("h1n1", subtype_values)
            self.assertIn("h3n2", subtype_values)
            self.assertEqual(context.get("segmented_default_segment"), "HA")

    def test_build_pileup_context_segmented_auto_detect_single_non_subtyped_source(self):
        with TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            data_dir = tmp / "data"
            ann_dir = tmp / "annotations"
            out_dir = tmp / "out"
            data_dir.mkdir(parents=True)
            ann_dir.mkdir(parents=True)
            out_dir.mkdir(parents=True)

            segmented_payload = (
                '{"virus":"flu-b","is_segmented":true,"segments":'
                '{"HA":{"positions":[1,2],"median":[8,9],"q1":[4,5],"q3":[12,13]},'
                '"PB1":{"positions":[1,2],"median":[6,7],"q1":[3,4],"q3":[9,10]}}}'
            )
            (data_dir / "flu-b_avg.json").write_text(segmented_payload, encoding="utf-8")
            self._write_genbank(ann_dir / "Influenza_B_HA.gb")
            self._write_genbank(ann_dir / "Influenza_B_PB1.gb")

            with patch.object(pileup_service, "PILEUP_DATA_DIR", data_dir), patch.object(
                pileup_service, "ANNOTATIONS_DIR", ann_dir
            ), patch.object(pileup_service, "PILEUP_OUTPUT_DIR", out_dir):
                context = pileup_service.build_pileup_context(
                    strain_slug="Influenza_B",
                    strain_name="Influenza B",
                    data_prefix="flu-b",
                    substrains=[],
                    annotation_name="Influenza_B_HA.gb",
                    levels=["all", "substrain", "individual"],
                    max_individual_traces=10,
                    segmented_segments=[
                        {"value": "PB1", "label": "PB1", "annotation": "Influenza_B_PB1.gb"},
                        {"value": "HA", "label": "HA", "annotation": "Influenza_B_HA.gb"},
                    ],
                    segmented_default_segment="HA",
                )

            self.assertTrue(context.get("enabled"))
            subtype_values = [x["value"] for x in context.get("segmented_subtypes", [])]
            self.assertIn("flu-b", subtype_values)
            self.assertIn("all", context.get("available_views", []))
            self.assertNotIn("substrain", context.get("available_views", []))
            self.assertIn("HA", context.get("segmented_all_files", {}))

    def test_build_mixed_sample_pileup_assets_non_segmented(self):
        with TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            data_dir = tmp / "data"
            ann_dir = tmp / "annotations"
            out_dir = tmp / "out"
            data_dir.mkdir(parents=True)
            ann_dir.mkdir(parents=True)
            out_dir.mkdir(parents=True)

            (data_dir / "sars-cov-2_indiv.json").write_text(
                '{"positions":[1,2,3],"samples":{"sample-x":[10,20,30]}}',
                encoding="utf-8",
            )
            self._write_genbank(ann_dir / "SARS-CoV-2.gb")

            sample_rows = pd.DataFrame(
                [{"canonical_strain": "SARS-CoV-2", "display_label": "SARS-CoV-2"}]
            )
            with patch.object(pileup_service, "PILEUP_DATA_DIR", data_dir), patch.object(
                pileup_service, "ANNOTATIONS_DIR", ann_dir
            ), patch.object(pileup_service, "PILEUP_OUTPUT_DIR", out_dir):
                assets = pileup_service.build_mixed_sample_pileup_assets("sample-x", sample_rows)

            self.assertEqual(len(assets), 1)
            self.assertTrue(assets[0]["available"])
            self.assertTrue((out_dir / "mixed_sample-x_SARS-CoV-2.html").exists())
            html = (out_dir / "mixed_sample-x_SARS-CoV-2.html").read_text(encoding="utf-8")
            self.assertIn("sample-x", html)

    def test_build_mixed_sample_pileup_assets_segmented_uses_default_segment(self):
        with TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            data_dir = tmp / "data"
            ann_dir = tmp / "annotations"
            out_dir = tmp / "out"
            data_dir.mkdir(parents=True)
            ann_dir.mkdir(parents=True)
            out_dir.mkdir(parents=True)

            (data_dir / "flu-b_indiv.json").write_text(
                (
                    '{"virus":"flu-b","is_segmented":true,"segments":'
                    '{"HA":{"positions":[1,2],"samples":{"sample-y":[8,9]}},'
                    '"PB1":{"positions":[1,2],"samples":{"sample-y":[4,5]}}}}'
                ),
                encoding="utf-8",
            )
            self._write_genbank(ann_dir / "Influenza_B_HA.gb")
            self._write_genbank(ann_dir / "Influenza_B_PB1.gb")

            sample_rows = pd.DataFrame(
                [{"canonical_strain": "Influenza B", "display_label": "Influenza B"}]
            )
            with patch.object(pileup_service, "PILEUP_DATA_DIR", data_dir), patch.object(
                pileup_service, "ANNOTATIONS_DIR", ann_dir
            ), patch.object(pileup_service, "PILEUP_OUTPUT_DIR", out_dir):
                assets = pileup_service.build_mixed_sample_pileup_assets("sample-y", sample_rows)

            self.assertEqual(len(assets), 1)
            self.assertTrue(assets[0]["available"])
            self.assertTrue((out_dir / "mixed_sample-y_Influenza_B_HA.html").exists())
            html = (out_dir / "mixed_sample-y_Influenza_B_HA.html").read_text(encoding="utf-8")
            self.assertIn("Influenza B - sample-y - HA", html)

    def test_build_mixed_sample_pileup_assets_matches_sample_prefix_from_metadata(self):
        with TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            data_dir = tmp / "data"
            ann_dir = tmp / "annotations"
            out_dir = tmp / "out"
            data_dir.mkdir(parents=True)
            ann_dir.mkdir(parents=True)
            out_dir.mkdir(parents=True)

            (data_dir / "sars-cov-2_indiv.json").write_text(
                '{"positions":[1,2,3],"samples":{"sample-z":[11,21,31]}}',
                encoding="utf-8",
            )
            self._write_genbank(ann_dir / "SARS-CoV-2.gb")

            sample_rows = pd.DataFrame(
                [{"canonical_strain": "SARS-CoV-2", "display_label": "SARS-CoV-2"}]
            )
            with patch.object(pileup_service, "PILEUP_DATA_DIR", data_dir), patch.object(
                pileup_service, "ANNOTATIONS_DIR", ann_dir
            ), patch.object(pileup_service, "PILEUP_OUTPUT_DIR", out_dir):
                assets = pileup_service.build_mixed_sample_pileup_assets(
                    "sample-z|2026-03-05|ZH",
                    sample_rows,
                )

            self.assertEqual(len(assets), 1)
            self.assertTrue(assets[0]["available"])
            self.assertTrue((out_dir / "mixed_sample-z_2026-03-05_ZH_SARS-CoV-2.html").exists())

    def test_build_mixed_sample_pileup_assets_matches_sample_prefix_with_m2_tag(self):
        with TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            data_dir = tmp / "data"
            ann_dir = tmp / "annotations"
            out_dir = tmp / "out"
            data_dir.mkdir(parents=True)
            ann_dir.mkdir(parents=True)
            out_dir.mkdir(parents=True)

            (data_dir / "sars-cov-2_indiv.json").write_text(
                '{"positions":[1,2,3],"samples":{"sample-z":[11,21,31]}}',
                encoding="utf-8",
            )
            self._write_genbank(ann_dir / "SARS-CoV-2.gb")

            sample_rows = pd.DataFrame(
                [{"canonical_strain": "SARS-CoV-2", "display_label": "SARS-CoV-2"}]
            )
            with patch.object(pileup_service, "PILEUP_DATA_DIR", data_dir), patch.object(
                pileup_service, "ANNOTATIONS_DIR", ann_dir
            ), patch.object(pileup_service, "PILEUP_OUTPUT_DIR", out_dir):
                assets = pileup_service.build_mixed_sample_pileup_assets(
                    "m2-sample-z|2026-03-05|ZH",
                    sample_rows,
                )

            self.assertEqual(len(assets), 1)
            self.assertTrue(assets[0]["available"])
            self.assertTrue((out_dir / "mixed_m2-sample-z_2026-03-05_ZH_SARS-CoV-2.html").exists())
