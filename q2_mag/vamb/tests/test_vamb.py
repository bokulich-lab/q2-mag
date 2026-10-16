# ----------------------------------------------------------------------------
# Copyright (c) 2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import glob
import os
import shutil
import tempfile
import pandas as pd
import unittest
from pathlib import Path
from unittest.mock import ANY, call, patch

from q2_types.per_sample_sequences import (
    BAMDirFmt,
    ContigSequencesDirFmt,
    MultiFASTADirectoryFormat,
)
from rachis.plugin.testing import TestPluginBase

from q2_mag.vamb.vamb import (
    _assert_reference_integrity,
    _assert_samples,
    _assert_taxids,
    _bin_contigs_vamb,
    _convert_taxonomy,
    _process_sample,
    _run_vamb,
    _process_vamb_arg,
    bin_contigs_vamb,
)


class TestVAMB(TestPluginBase):
    package = "q2_mag.vamb.tests"

    def test_assert_reference_integrity_ok(self):
        contigs_path = self.get_data_path("contigs")
        maps_path = self.get_data_path("maps")

        sample_set = {
            "samp1": {
                "contigs": str(Path(contigs_path) / "samp1_contigs.fa"),
                "map": str(Path(maps_path) / "samp1_alignment.bam"),
            },
            "samp2": {
                "contigs": str(Path(contigs_path) / "samp2_contigs.fa"),
                "map": str(Path(maps_path) / "samp2_alignment.bam"),
            },
        }

        _assert_reference_integrity(sample_set, "contigs", "map")

    def test_assert_reference_integrity_mismatch(self):
        contigs_path = self.get_data_path("contigs")
        maps_path = self.get_data_path("maps")

        sample_set = {
            "samp1": {
                "contigs": str(Path(contigs_path) / "samp1_contigs.fa"),
                "map": str(Path(maps_path) / "samp2_alignment.bam"),
            },
            "samp2": {
                "contigs": str(Path(contigs_path) / "samp2_contigs.fa"),
                "map": str(Path(maps_path) / "samp1_alignment.bam"),
            },
        }

        with self.assertRaisesRegex(
            Exception,
            "Alignment maps do not match the corresponding contigs in at least one "
            "sample. The following samples had a mismatch in count, name, order, and "
            "or length:",
        ):
            _assert_reference_integrity(sample_set, "contigs", "map")

    def test_assert_samples_ok(self):
        contigs_path = self.get_data_path("contigs")
        maps_path = self.get_data_path("maps")

        contigs = ContigSequencesDirFmt(contigs_path, mode="r")
        maps = BAMDirFmt(maps_path, mode="r")

        obs_samples = _assert_samples(contigs, maps)

        exp_samples = {
            "samp1": {
                "contigs": str(Path(contigs_path) / "samp1_contigs.fa"),
                "map": str(Path(maps_path) / "samp1_alignment.bam"),
            },
            "samp2": {
                "contigs": str(Path(contigs_path) / "samp2_contigs.fa"),
                "map": str(Path(maps_path) / "samp2_alignment.bam"),
            },
        }
        self.assertDictEqual(exp_samples, obs_samples)

    def test_assert_samples_uneven(self):
        contigs_path = self.get_data_path("contigs")
        with tempfile.TemporaryDirectory() as maps_path:
            map_path = Path(self.get_data_path("maps")) / "samp1_alignment.bam"
            shutil.copy(map_path, maps_path)

            contigs = ContigSequencesDirFmt(contigs_path, mode="r")
            maps = BAMDirFmt(maps_path, mode="r")

            with self.assertRaisesRegex(
                Exception,
                "Contigs and alignment maps should belong to the same sample "
                "set. You provided contigs for samples: samp1,samp2 but maps "
                "for samples: samp1. Please check your inputs and try again.",
            ):
                _assert_samples(contigs, maps)

    def test_assert_samples_non_matching(self):
        with tempfile.TemporaryDirectory() as tempdir:
            contigs_path = Path(tempdir) / "contigs-path"
            maps_path = Path(tempdir) / "maps-path"
            os.makedirs(contigs_path)
            os.makedirs(maps_path)

            contig_path = Path(self.get_data_path("contigs")) / "samp1_contigs.fa"
            map_path = Path(self.get_data_path("maps")) / "samp2_alignment.bam"

            shutil.copy(contig_path, contigs_path)
            shutil.copy(map_path, maps_path)

            contigs = ContigSequencesDirFmt(contigs_path, mode="r")
            maps = BAMDirFmt(maps_path, mode="r")

            with self.assertRaisesRegex(
                Exception,
                "Contigs and alignment maps should belong to the same sample"
                " set. You provided contigs for samples: samp1 but maps "
                "for samples: samp2. Please check your inputs and try again.",
            ):
                _assert_samples(contigs, maps)

    @patch("subprocess.run")
    def test_run_vamb_ok(self, p1):
        fake_props = {"map": "/some/where/map.bam", "contigs": "/a/b/co.fa"}
        fake_args = ["--verbose", "--minfasta", "1000", "-m", "2000"]

        with tempfile.TemporaryDirectory() as fake_loc:
            obs_fp = _run_vamb(
                "default", "samp1", fake_props, fake_loc, fake_args, "", None
            )
            exp_fp = os.path.join(fake_loc, "samp1")

            self.assertEqual(exp_fp, obs_fp)

            exp_cmd = [
                "vamb",
                "bin",
                "default",
                "--fasta",
                fake_props["contigs"],
                "--bamdir",
                ANY,
                "--outdir",
                os.path.join(fake_loc, "samp1", "bin"),
                "-o",
                "",
            ]
            exp_cmd.extend(fake_args)
            p1.assert_called_once_with(exp_cmd, check=True)

    @patch("subprocess.run")
    def test_run_taxvamb_ok(self, p1):
        fake_props = {"map": "/some/where/map.bam", "contigs": "/a/b/co.fa"}
        fake_args = ["--verbose", "--minfasta", "1000", "-m", "2000"]
        taxonomy_fp = self.get_data_path("/taxonomy/taxonomy.tsv")
        taxonomy_df = pd.read_csv(taxonomy_fp, sep="\t")

        with tempfile.TemporaryDirectory() as fake_loc:
            obs_fp = _run_vamb(
                "default", "samp1", fake_props, fake_loc, fake_args, "", taxonomy_df
            )
            exp_fp = os.path.join(fake_loc, "samp1")

            self.assertEqual(exp_fp, obs_fp)

            exp_cmd = [
                "vamb",
                "bin",
                "default",
                "--fasta",
                fake_props["contigs"],
                "--bamdir",
                ANY,
                "--outdir",
                os.path.join(fake_loc, "samp1", "bin"),
                "-o",
                "",
            ]
            exp_cmd.extend(fake_args)
            exp_cmd.extend(["--taxonomy", os.path.join(fake_loc, "taxonomy.tsv")])
            p1.assert_called_once_with(exp_cmd, check=True)

    @patch("tempfile.TemporaryDirectory")
    @patch("q2_mag.vamb.vamb.uuid4")
    @patch("q2_mag.vamb.vamb._run_vamb")
    def test_process_sample(self, p1, p2, p3):
        fake_props = {
            "map": "some/where/samp1_alignment.bam",
            "contigs": "some/where/samp1_contigs.fasta",
        }
        fake_args = ["--verbose", "--minfasta", "1000", "-m", "2000"]

        p2.side_effect = [
            "522775d4-b1c6-4ee3-8b47-cd990f17eb8b",
            "684db670-6304-4f33-a0ea-7f570532e178",
            "37356c23-b8db-4bbe-b4c9-d35e1cef615b",
            "51c19113-31f0-4e4c-bbb3-b9df26b949f3",
        ]
        fake_temp_dir = tempfile.mkdtemp()
        p3.return_value.__enter__.return_value = fake_temp_dir

        with tempfile.TemporaryDirectory() as fake_loc:
            p1.return_value = os.path.join(fake_loc, "samp1")

            # copy two expected bins to the new location
            samp1_bins_fp = self.get_data_path("bins-no-uuid/samp1")
            shutil.copytree(
                samp1_bins_fp,
                os.path.join(fake_loc, "samp1", "bin", "bins"),
                dirs_exist_ok=True,
            )

            _process_sample(
                "samp1", fake_props, "default", False, fake_args, fake_loc, None
            )

            # find the newly formed bins
            obs_bins = set(
                [
                    x.split("/")[-1]
                    for x in glob.glob(os.path.join(fake_loc, "samp1", "*.fa"))
                ]
            )
            exp_bins = {
                "522775d4-b1c6-4ee3-8b47-cd990f17eb8b.fa",
                "684db670-6304-4f33-a0ea-7f570532e178.fa",
            }
            self.assertSetEqual(exp_bins, obs_bins)

            p1.assert_called_once_with(
                "default", "samp1", fake_props, ANY, fake_args, "", None
            )

    @patch("q2_mag.vamb.vamb.ContigSequencesDirFmt")
    @patch("q2_mag.vamb.vamb.MultiFASTADirectoryFormat")
    @patch("q2_mag.vamb.vamb._process_sample")
    def test_bin_contigs_vamb(self, p1, p2, p3):
        input_contigs = self.get_data_path("contigs")
        input_maps = self.get_data_path("maps")
        contigs = ContigSequencesDirFmt(input_contigs, mode="r")
        maps = BAMDirFmt(input_maps, mode="r")

        args = ["--verbose", "--minfasta", "1000", "-m", "2000"]

        mock_bins = MultiFASTADirectoryFormat(self.get_data_path("bins"), "r")
        p2.return_value = mock_bins

        mock_unbinned = ContigSequencesDirFmt(
            self.get_data_path("contigs/samp1_contigs.fa"), "r"
        )
        p3.return_value = mock_unbinned

        obs_bins, obs_map = _bin_contigs_vamb(contigs, maps, None, False, args)

        self.assertIsInstance(obs_bins, MultiFASTADirectoryFormat)
        p1.assert_has_calls(
            [
                call(
                    "samp1",
                    {
                        "contigs": self.get_data_path("/contigs/samp1_contigs.fa"),
                        "map": self.get_data_path("/maps/samp1_alignment.bam"),
                    },
                    "default",
                    False,
                    args,
                    str(mock_bins),
                    None,
                ),
                call(
                    "samp2",
                    {
                        "contigs": self.get_data_path("/contigs/samp2_contigs.fa"),
                        "map": self.get_data_path("/maps/samp2_alignment.bam"),
                    },
                    "default",
                    False,
                    args,
                    str(mock_bins),
                    None,
                ),
            ]
        )

        # find the newly formed bins
        obs_bins = sorted(
            [
                sorted(
                    [
                        "/".join(x.split("/")[-2:])
                        for x in glob.glob(
                            os.path.join(str(obs_bins), f"samp{y}", "*.fa")
                        )
                    ]
                )
                for y in (1, 2)
            ]
        )
        exp_bins = [
            [
                "samp1/522775d4-b1c6-4ee3-8b47-cd990f17eb8b.fa",
                "samp1/684db670-6304-4f33-a0ea-7f570532e178.fa",
            ],
            [
                "samp2/37356c23-b8db-4bbe-b4c9-d35e1cef615b.fa",
                "samp2/51c19113-31f0-4e4c-bbb3-b9df26b949f3.fa",
            ],
        ]
        self.assertListEqual(exp_bins, obs_bins)

    @patch("q2_mag.vamb.vamb.MultiFASTADirectoryFormat")
    @patch("q2_mag.vamb.vamb._process_sample")
    def test_bin_contigs_vamb_no_mags(self, p1, p2):
        input_contigs = self.get_data_path("contigs")
        input_maps = self.get_data_path("maps")
        contigs = ContigSequencesDirFmt(input_contigs, mode="r")
        maps = BAMDirFmt(input_maps, mode="r")

        args = ["--verbose", "--minfasta", "1000", "-m", "2000"]

        mock_bins = MultiFASTADirectoryFormat()
        p2.return_value = mock_bins

        with self.assertRaisesRegex(ValueError, "No MAGs were formed"):
            _bin_contigs_vamb(contigs, maps, None, False, args)

    @patch("q2_mag.vamb.vamb._bin_contigs_vamb")
    @patch("q2_mag.vamb.vamb._process_common_input_params")
    def test_bin_contigs_vamb_wrapper(self, p1, p2):
        input_contigs = self.get_data_path("contigs")
        input_maps = self.get_data_path("maps")
        contigs = ContigSequencesDirFmt(input_contigs, mode="r")
        maps = BAMDirFmt(input_maps, mode="r")

        p1.return_value = ["--minfasta", "1000", "-m", "2000"]
        p2.return_value = ("bins", {"contigA": "bin1"})

        obs_bins = bin_contigs_vamb(
            contigs=contigs,
            alignment_maps=maps,
            min_contig_len=1000,
            minfasta=2000,
            threads=8,
            seed="12345",
        )
        exp_bins = ("bins", {"contigA": "bin1"})

        p1.assert_called_once_with(
            processing_func=_process_vamb_arg,
            params={
                "min_contig_len": 1000,
                "minfasta": 2000,
                "threads": 8,
                "seed": "12345",
            },
        )
        p2.assert_called_once_with(
            contigs=contigs,
            alignment_maps=maps,
            taxonomy=None,
            multi_split=False,
            common_args=p1.return_value,
        )
        self.assertTupleEqual(exp_bins, obs_bins)

    def test_bin_contigs_vamb_wrapper_taxid_map_without_taxonomy(self):
        taxid_to_contigs_fp = self.get_data_path(
            "/taxid-to-contigids-maps/feature_map.json"
        )
        taxid_to_contigs_df = pd.read_json(
            taxid_to_contigs_fp, lines=True, dtype={"name": str}
        )
        taxid_to_contigs_dict = taxid_to_contigs_df.set_index("name")[
            "members"
        ].to_dict()

        with self.assertRaisesRegex(
            ValueError, "Received `taxid_map`, but no `taxonomy` input"
        ):
            bin_contigs_vamb(
                contigs=None,
                alignment_maps=None,
                taxid_map=taxid_to_contigs_dict,
                taxonomy=None,
            )

    def test_bin_contigs_vamb_wrapper_taxonomy_without_taxid_map(self):
        taxid_to_full_string_fp = self.get_data_path(
            "/taxid-to-full-string-maps/taxonomy.tsv"
        )
        taxid_to_full_string_df = pd.read_csv(
            taxid_to_full_string_fp, sep="\t", index_col="Feature ID"
        )
        taxid_to_full_string_df.index = taxid_to_full_string_df.index.astype(str)

        with self.assertRaisesRegex(
            ValueError, "Received `taxonomy`, but no `taxid_map` input"
        ):
            bin_contigs_vamb(
                contigs=None,
                alignment_maps=None,
                taxid_map=None,
                taxonomy=taxid_to_full_string_df,
            )

    def test_bin_contigs_vamb_wrapper_no_predictor_without_taxonomic_info(self):
        with self.assertRaisesRegex(
            ValueError, "`--p-no-predictor` flag is set, but no taxonomic information"
        ):
            bin_contigs_vamb(
                contigs=None,
                alignment_maps=None,
                taxid_map=None,
                taxonomy=None,
                no_predictor=True,
            )

    @patch("q2_mag.vamb.vamb._bin_contigs_vamb")
    @patch("q2_mag.vamb.vamb._process_common_input_params")
    @patch("q2_mag.vamb.vamb._convert_taxonomy")
    def test_bin_contigs_vamb_wrapper_with_taxonomic_information(self, p1, p2, p3):
        taxid_to_contigs_fp = self.get_data_path(
            "/taxid-to-contigids-maps/feature_map.json"
        )
        taxid_to_contigs_df = pd.read_json(
            taxid_to_contigs_fp, lines=True, dtype={"name": str}
        )
        taxid_to_contigs_dict = taxid_to_contigs_df.set_index("name")[
            "members"
        ].to_dict()

        taxid_to_full_string_fp = self.get_data_path(
            "/taxid-to-full-string-maps/taxonomy.tsv"
        )
        taxid_to_full_string_df = pd.read_csv(
            taxid_to_full_string_fp, sep="\t", index_col="Feature ID"
        )
        taxid_to_full_string_df.index = taxid_to_full_string_df.index.astype(str)

        input_contigs = self.get_data_path("contigs")
        input_maps = self.get_data_path("maps")
        contigs = ContigSequencesDirFmt(input_contigs, mode="r")
        maps = BAMDirFmt(input_maps, mode="r")

        p1.return_value = taxid_to_full_string_df
        p2.return_value = ["--minfasta", "1000"]
        p3.return_value = ("bins", {"contigA": "bin1"})

        obs_bins = bin_contigs_vamb(
            contigs=contigs,
            alignment_maps=maps,
            taxid_map=taxid_to_contigs_dict,
            taxonomy=taxid_to_full_string_df,
        )

        p1.assert_called_once_with(
            taxid_map=taxid_to_contigs_dict,
            taxonomy=taxid_to_full_string_df,
        )

        p2.assert_called_once_with(
            processing_func=_process_vamb_arg,
            params={
                "min_contig_len": 2000,
                "minfasta": 2000,
                "threads": 8,
                "seed": None,
                "no_predictor": False,
            },
        )

        p3.assert_called_once_with(
            contigs=contigs,
            alignment_maps=maps,
            taxonomy=taxid_to_full_string_df,
            multi_split=False,
            common_args=p2.return_value,
        )
        self.assertTupleEqual(("bins", {"contigA": "bin1"}), obs_bins)

    def test_assert_taxids_ok(self):
        taxid_to_contigs_fp = self.get_data_path(
            "/taxid-to-contigids-maps/feature_map.json"
        )
        taxid_to_full_string_fp = self.get_data_path(
            "/taxid-to-full-string-maps/taxonomy.tsv"
        )
        taxid_to_contigs_df = pd.read_json(
            taxid_to_contigs_fp, lines=True, dtype={"name": str}
        )
        taxid_to_contigs_dict = taxid_to_contigs_df.set_index("name")[
            "members"
        ].to_dict()
        taxid_to_full_string_df = pd.read_csv(
            taxid_to_full_string_fp, sep="\t", index_col="Feature ID"
        )
        taxid_to_full_string_df.index = taxid_to_full_string_df.index.astype(str)

        _assert_taxids(
            taxid_map=taxid_to_contigs_dict, taxonomy=taxid_to_full_string_df
        )

    def test_assert_taxids_missing_taxid(self):
        taxid_to_contigs_fp = self.get_data_path(
            "/taxid-to-contigids-maps/feature_map.json"
        )
        taxid_to_full_string_fp = self.get_data_path(
            "/taxid-to-full-string-maps/taxonomy.tsv"
        )
        taxid_to_contigs_df = pd.read_json(
            taxid_to_contigs_fp, lines=True, dtype={"name": str}
        )
        taxid_to_contigs_dict = taxid_to_contigs_df.set_index("name")[
            "members"
        ].to_dict()
        taxid_to_full_string_df = pd.read_csv(
            taxid_to_full_string_fp, sep="\t", index_col="Feature ID"
        )
        taxid_to_full_string_df.index = taxid_to_full_string_df.index.astype(str)

        taxid_to_full_string_df = taxid_to_full_string_df.drop(index="1485")

        with self.assertRaisesRegex(
            AssertionError,
            "The following taxonomy IDs are missing: 1485",
        ):
            _assert_taxids(
                taxid_map=taxid_to_contigs_dict, taxonomy=taxid_to_full_string_df
            )

    def test_convert_taxonomy(self):
        taxid_to_contigs_fp = self.get_data_path(
            "/taxid-to-contigids-maps/feature_map.json"
        )
        taxid_to_full_string_fp = self.get_data_path(
            "/taxid-to-full-string-maps/taxonomy.tsv"
        )
        taxid_to_contigs_df = pd.read_json(
            taxid_to_contigs_fp, lines=True, dtype={"name": str}
        )
        taxid_to_contigs_dict = taxid_to_contigs_df.set_index("name")[
            "members"
        ].to_dict()
        taxid_to_full_string_df = pd.read_csv(
            taxid_to_full_string_fp, sep="\t", index_col="Feature ID"
        )
        taxid_to_full_string_df.index = taxid_to_full_string_df.index.astype(str)

        obs = (
            _convert_taxonomy(
                taxid_map=taxid_to_contigs_dict, taxonomy=taxid_to_full_string_df
            )
            .sort_values("contigs")
            .reset_index(drop=True)
        )

        exp_taxonomy_fp = self.get_data_path("/taxonomy/taxonomy.tsv")
        exp = (
            pd.read_csv(exp_taxonomy_fp, sep="\t")
            .sort_values("contigs")
            .reset_index(drop=True)
        )

        pd.testing.assert_frame_equal(obs, exp, check_dtype=False)


if __name__ == "__main__":
    unittest.main()
