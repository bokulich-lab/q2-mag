# ----------------------------------------------------------------------------
# Copyright (c) 2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import glob

import shutil

import os
import tempfile
import unittest
import uuid
from pathlib import Path
from q2_types.per_sample_sequences import ContigSequencesDirFmt, BAMDirFmt
from q2_types.per_sample_sequences import MultiFASTADirectoryFormat
from unittest.mock import MagicMock, patch, call

from qiime2.plugin.testing import TestPluginBase

from q2_mag.semibin2.semibin2 import (
    bin_contigs_semibin2,
    _bin_contigs_semibin2,
    _bin_partition_contigs_semibin2,
    _filter_alignment_maps,
    _run_semibin2,
    _process_sample,
)


class TestSemibin2(TestPluginBase):
    package = "q2_mag.semibin2.tests"

    def test_filter_alignment_maps(self):
        contigs_path = self.get_data_path("contigs")
        maps = BAMDirFmt(self.get_data_path("maps"), mode="r")

        contigs_one_sample = ContigSequencesDirFmt()
        shutil.copy(
            os.path.join(contigs_path, "samp1_contigs.fa"), contigs_one_sample.path
        )

        obs = _filter_alignment_maps(contigs_one_sample, maps)

        obs_maps = obs.file_dict()
        exp_map = maps.path / "samp1_alignment.bam"

        self.assertListEqual(["samp1_alignment"], list(obs_maps))
        self.assertEqual(
            exp_map.read_bytes(),
            Path(obs_maps["samp1_alignment"]).read_bytes(),
        )

    @patch("subprocess.run")
    def test_run_semibin2_ok(self, p1):
        fake_props = {"map": "/some/where/map.bam", "contigs": "/a/b/co.fa"}
        fake_args = [
            "--epochs",
            "16",
            "--min-len",
            "1000",
        ]

        with tempfile.TemporaryDirectory() as fake_loc:
            obs_fp = _run_semibin2(
                "samp1", fake_props, fake_loc, "single_easy_bin", fake_args
            )
            exp_fp = os.path.join(fake_loc, "samp1")

            self.assertEqual(exp_fp, obs_fp)

            exp_cmd = [
                "SemiBin2",
                "single_easy_bin",
                "--input-fasta",
                fake_props["contigs"],
                "--input-bam",
                fake_props["map"],
                "--output",
                os.path.join(fake_loc, "samp1", "bin"),
                "--compression",
                "none",
                "--verbose",
            ]
            exp_cmd.extend(fake_args)
            p1.assert_called_once_with(exp_cmd, check=True)

    @patch("tempfile.TemporaryDirectory")
    @patch("q2_mag.semibin2.semibin2.uuid4")
    @patch("q2_mag.semibin2.semibin2._run_semibin2")
    def test_process_sample(self, p1, p2, p3):
        fake_props = {
            "map": "some/where/samp1_alignment.bam",
            "contigs": "some/where/samp1_contigs.fasta",
        }
        fake_args = [
            "--compression",
            "none",
            "--epochs",
            "16",
            "--min-len",
            "1000",
        ]

        fake_temp_dir = tempfile.mkdtemp()
        p2.side_effect = [
            uuid.UUID("93fa393d-504a-4257-8a38-931af0cc4b1f"),
            uuid.UUID("a69ab2c1-b23a-400f-be45-7dfc35bf343c"),
        ]
        p3.return_value.__enter__.return_value = fake_temp_dir

        with tempfile.TemporaryDirectory() as fake_loc:
            p1.return_value = os.path.join(fake_loc, "bins", "samp1")

            # copy two expected bins to the new location
            samp1_bins_fp = self.get_data_path("bins-no-uuid/samp1")
            shutil.copytree(
                samp1_bins_fp,
                os.path.join(fake_loc, "bins", "samp1", "bin", "output_bins"),
                dirs_exist_ok=True,
            )

            _process_sample("samp1", fake_props, "single_easy_bin", fake_args, fake_loc)

            # find the newly formed bins
            obs_bins = set(
                [
                    x.split("/")[-1]
                    for x in glob.glob(os.path.join(fake_loc, "samp1", "*.fa"))
                ]
            )
            exp_bins = {
                "93fa393d-504a-4257-8a38-931af0cc4b1f.fa",
                "a69ab2c1-b23a-400f-be45-7dfc35bf343c.fa",
            }
            self.assertSetEqual(exp_bins, obs_bins)

    @patch("q2_mag.semibin2.semibin2._filter_alignment_maps")
    @patch("q2_mag.semibin2.semibin2.ContigSequencesDirFmt")
    @patch("q2_mag.semibin2.semibin2.MultiFASTADirectoryFormat")
    @patch("q2_mag.semibin2.semibin2._process_sample")
    def test_bin_partition_contigs_semibin2(self, p1, p2, _, filter_alignment_maps):
        input_contigs = self.get_data_path("contigs")
        input_maps = self.get_data_path("maps")
        contigs = ContigSequencesDirFmt(input_contigs, mode="r")
        maps = BAMDirFmt(input_maps, mode="r")
        filter_alignment_maps.return_value = maps

        args = [
            "--compression",
            "none",
            "--epochs",
            "16",
            "--min-len",
            "1000",
        ]

        mock_bins = MultiFASTADirectoryFormat(self.get_data_path("bins"), "r")
        p2.return_value = mock_bins

        obs_bins, obs_map = _bin_partition_contigs_semibin2(
            contigs, maps, "single_easy_bin", args
        )

        filter_alignment_maps.assert_called_once_with(contigs, maps)
        self.assertIsInstance(obs_bins, MultiFASTADirectoryFormat)
        p1.assert_has_calls(
            [
                call(
                    "samp1",
                    {
                        "contigs": self.get_data_path("/contigs/samp1_contigs.fa"),
                        "map": self.get_data_path("/maps/samp1_alignment.bam"),
                    },
                    "single_easy_bin",
                    args,
                    str(mock_bins),
                ),
                call(
                    "samp2",
                    {
                        "contigs": self.get_data_path("/contigs/samp2_contigs.fa"),
                        "map": self.get_data_path("/maps/samp2_alignment.bam"),
                    },
                    "single_easy_bin",
                    args,
                    str(mock_bins),
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

    @patch("q2_mag.semibin2.semibin2.MultiFASTADirectoryFormat")
    @patch("q2_mag.semibin2.semibin2._process_sample")
    def test_bin_partition_contigs_semibin2_no_mags(self, p1, p2):
        input_contigs = self.get_data_path("contigs")
        input_maps = self.get_data_path("maps")
        contigs = ContigSequencesDirFmt(input_contigs, mode="r")
        maps = BAMDirFmt(input_maps, mode="r")

        args = [
            "--compression",
            "none",
            "--epochs",
            "16",
            "--min-len",
            "1000",
        ]

        mock_bins = MultiFASTADirectoryFormat()
        p2.return_value = mock_bins

        with self.assertRaisesRegex(ValueError, "No MAGs were formed"):
            _bin_partition_contigs_semibin2(contigs, maps, "single_easy_bin", args)

    @patch("q2_mag.semibin2.semibin2._bin_partition_contigs_semibin2")
    @patch("q2_mag.semibin2.semibin2._process_common_input_params")
    def test_bin_contigs_semibin2_action(self, p1, p2):
        p1.return_value = ["--epochs", "16"]
        p2.return_value = ("bins", {"contigA": "bin1"})

        contigs = MagicMock()
        bams = MagicMock()

        result = _bin_contigs_semibin2(
            contigs=contigs,
            alignment_maps=bams,
            epochs=10,
        )

        p1.assert_called_once()

        p2.assert_called_once_with(
            contigs=contigs,
            alignment_maps=bams,
            mode="single_easy_bin",
            common_args=["--epochs", "16"],
        )

        assert result == ("bins", {"contigA": "bin1"})

    def test_bin_contigs_semibin2_pipeline(self):
        contigs = ContigSequencesDirFmt(self.get_data_path("contigs"), mode="r")
        alignment_maps = BAMDirFmt(self.get_data_path("maps"), mode="r")
        partitions = MagicMock()
        partitions.values.return_value = ["partition1", "partition2"]

        partition_contigs = MagicMock(return_value=(partitions,))
        bin_partition = MagicMock(
            side_effect=[
                ("mags1", "contig_map1"),
                ("mags2", "contig_map2"),
            ]
        )
        collate_mags = MagicMock(return_value=("collated_mags",))
        collate_contig_maps = MagicMock(return_value=("collated_contig_map",))

        actions = {
            ("types", "partition_contigs"): partition_contigs,
            ("mag", "_bin_contigs_semibin2"): bin_partition,
            ("types", "collate_sample_data_mags"): collate_mags,
            ("mag", "collate_contig_maps"): collate_contig_maps,
        }
        ctx = MagicMock()
        ctx.get_action.side_effect = lambda plugin, action: actions[(plugin, action)]

        obs = bin_contigs_semibin2(
            ctx=ctx,
            contigs=contigs,
            alignment_maps=alignment_maps,
            epochs=10,
            num_partitions=2,
        )

        partition_contigs.assert_called_once_with(contigs, 2)
        worker_params = {
            "training_type": None,
            "orf_finder": "fast-naive",
            "environment": "global",
            "engine": "auto",
            "sequencing_type": "short_read",
            "minfasta_kbs": 200,
            "no_recluster": False,
            "epochs": 10,
            "batch_size": 2048,
            "max_node": 1,
            "max_edges": 200,
            "ratio": 0.05,
            "threads": 1,
            "min_len": None,
            "ml_threshold": None,
            "random_seed": None,
            "debug": False,
        }
        self.assertListEqual(
            [
                call("partition1", alignment_maps, **worker_params),
                call("partition2", alignment_maps, **worker_params),
            ],
            bin_partition.call_args_list,
        )
        collate_mags.assert_called_once_with(["mags1", "mags2"])
        collate_contig_maps.assert_called_once_with(["contig_map1", "contig_map2"])
        self.assertTupleEqual(("collated_mags", "collated_contig_map"), obs)


if __name__ == "__main__":
    unittest.main()
