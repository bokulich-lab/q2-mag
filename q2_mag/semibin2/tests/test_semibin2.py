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
import unittest
import uuid
from pathlib import Path
from unittest.mock import MagicMock, call, patch

from q2_types.per_sample_sequences import (
    BAMDirFmt,
    ContigSequencesDirFmt,
    MultiFASTADirectoryFormat,
)
from qiime2 import Artifact
from qiime2.plugin.testing import TestPluginBase
from rachis.sdk import Context

from q2_mag.semibin2.semibin2 import (
    _bin_contigs_semibin2,
    _bin_partition_contigs_semibin2,
    _filter_alignment_maps,
    _process_sample,
    _run_semibin2,
    bin_contigs_semibin2,
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

        self.assertListEqual(["samp1"], list(obs_maps))
        self.assertEqual(
            exp_map.read_bytes(),
            Path(obs_maps["samp1"]).read_bytes(),
        )

    @patch("subprocess.run")
    def test_run_semibin2_ok(self, p1):
        fake_props = {"map": ["/some/where/map.bam"], "contigs": "/a/b/co.fa"}
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
                *fake_props["map"],
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
            contigs, maps, "single-sample", args
        )

        filter_alignment_maps.assert_called_once_with(contigs, maps)
        self.assertIsInstance(obs_bins, MultiFASTADirectoryFormat)
        p1.assert_has_calls(
            [
                call(
                    "samp1",
                    {
                        "contigs": self.get_data_path("/contigs/samp1_contigs.fa"),
                        "map": [self.get_data_path("/maps/samp1_alignment.bam")],
                    },
                    "single_easy_bin",
                    args,
                    str(mock_bins),
                    separator="",
                ),
                call(
                    "samp2",
                    {
                        "contigs": self.get_data_path("/contigs/samp2_contigs.fa"),
                        "map": [self.get_data_path("/maps/samp2_alignment.bam")],
                    },
                    "single_easy_bin",
                    args,
                    str(mock_bins),
                    separator="",
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
            _bin_partition_contigs_semibin2(contigs, maps, "single-sample", args)

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
            binning_workflow="single-sample",
        )

        p1.assert_called_once()

        p2.assert_called_once_with(
            contigs=contigs,
            alignment_maps=bams,
            binning_workflow="single-sample",
            common_args=["--epochs", "16"],
        )

        assert result == ("bins", {"contigA": "bin1"})

    def test_bin_contigs_semibin2_pipeline(self):
        contigs = Artifact.import_data(
            "SampleData[Contigs]", self.get_data_path("contigs")
        )
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
        collated_contig_map = MagicMock()
        collated_contig_map.view.return_value = {
            "mag1": ["contig2", "contig1"],
            "mag2": ["contig3"],
        }
        collate_contig_maps = MagicMock(return_value=(collated_contig_map,))
        filter_contigs = MagicMock(return_value=("unbinned_contigs",))

        actions = {
            ("types", "partition_contigs"): partition_contigs,
            ("mag", "_bin_contigs_semibin2"): bin_partition,
            ("types", "collate_sample_data_mags"): collate_mags,
            ("types", "collate_contig_maps"): collate_contig_maps,
            ("assembly", "filter_contigs"): filter_contigs,
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
                call("partition1", alignment_maps, "single-sample", **worker_params),
                call("partition2", alignment_maps, "single-sample", **worker_params),
            ],
            bin_partition.call_args_list,
        )
        collate_mags.assert_called_once_with(["mags1", "mags2"])
        collate_contig_maps.assert_called_once_with(["contig_map1", "contig_map2"])
        filter_contigs.assert_called_once_with(
            contigs,
            on="contig",
            ids=["contig1", "contig2", "contig3"],
            exclude_ids=True,
        )
        self.assertTupleEqual(
            ("collated_mags", collated_contig_map, "unbinned_contigs"), obs
        )

    def test_bin_contigs_semibin2_coassembly_invalid_maps(self):
        self.package = "q2_mag.tests"
        contigs = Artifact.import_data(
            "FeatureData[Contig % Properties('coassembled')]",
            self.get_data_path("coassembled-contigs"),
        )
        alignment_maps = Artifact.import_data(
            "SampleData[AlignmentMap % Properties('sorted')]",
            self.get_data_path("coassembled-maps"),
        )

        with self.assertRaisesRegex(ValueError, "Expected maps of the type"):
            bin_contigs_semibin2(
                ctx=None,
                contigs=contigs,
                alignment_maps=alignment_maps,
            )

    def test_bin_contigs_semibin2_multi_sample_invalid_maps(self):
        self.package = "q2_mag.tests"
        contigs = Artifact.import_data(
            "FeatureData[Contig % Properties('pooled')]",
            self.get_data_path("pooled-contigs/colon-separator"),
        )
        alignment_maps = Artifact.import_data(
            "SampleData[AlignmentMap % Properties('sorted')]",
            self.get_data_path("pooled-maps/colon-separator"),
        )

        with self.assertRaisesRegex(ValueError, "Expected maps of the type"):
            bin_contigs_semibin2(
                ctx=None,
                contigs=contigs,
                alignment_maps=alignment_maps,
            )

    def test_bin_contigs_semibin2_unsupported_contigs(self):
        self.package = "q2_mag.tests"
        contigs = Artifact.import_data(
            "FeatureData[Sequence]",
            self.get_data_path("sequences/dna-sequences.fasta"),
        )

        with self.assertRaises(NotImplementedError):
            bin_contigs_semibin2(
                ctx=None,
                contigs=contigs,
                alignment_maps=None,
            )

    @patch("q2_mag.semibin2.semibin2._bin_partition_contigs_semibin2")
    @patch("q2_mag.semibin2.semibin2._process_common_input_params")
    def test_bin_contigs_semibin2_multi_sample_no_environment(self, p1, p2):
        contigs = ContigSequencesDirFmt(self.get_data_path("contigs"), mode="r")
        alignment_maps = BAMDirFmt(self.get_data_path("maps"), mode="r")
        p1.return_value = ["--epochs", "15"]

        _bin_contigs_semibin2(
            contigs=contigs,
            alignment_maps=alignment_maps,
            binning_workflow="multi-sample",
            environment="global",
        )

        p1.assert_called_once()
        self.assertNotIn("environment", p1.call_args.kwargs["params"])
        p2.assert_called_once_with(
            contigs=contigs,
            alignment_maps=alignment_maps,
            binning_workflow="multi-sample",
            common_args=["--epochs", "15"],
        )

    @patch("q2_mag.semibin2.semibin2.MultiFASTADirectoryFormat")
    @patch("q2_mag.semibin2.semibin2._process_sample")
    def test_bin_partition_contigs_semibin2_multi_sample(self, p1, p2):
        self.package = "q2_mag.tests"
        contigs = ContigSequencesDirFmt(
            self.get_data_path("pooled-contigs/semicolon-separator"), mode="r"
        )
        maps = BAMDirFmt(
            self.get_data_path("pooled-maps/semicolon-separator"), mode="r"
        )
        args = ["--epochs", "15"]

        mock_bins = MultiFASTADirectoryFormat()
        bin_id = "522775d4-b1c6-4ee3-8b47-cd990f17eb8b"
        os.makedirs(mock_bins.path / "pooled")
        shutil.copy(
            contigs.path / "pooled_contigs.fa",
            mock_bins.path / "pooled" / f"{bin_id}.fa",
        )
        p2.return_value = mock_bins

        obs_bins, obs_map = _bin_partition_contigs_semibin2(
            contigs, maps, "multi-sample", args
        )

        p1.assert_called_once_with(
            "pooled",
            {
                "contigs": self.get_data_path(
                    "pooled-contigs/semicolon-separator/pooled_contigs.fa"
                ),
                "map": [
                    self.get_data_path("pooled-maps/semicolon-separator/sample1.bam"),
                    self.get_data_path("pooled-maps/semicolon-separator/sample2.bam"),
                ],
            },
            "multi_easy_bin",
            args + ["--separator", ";"],
            str(mock_bins),
            separator=";",
        )
        self.assertEqual(obs_bins, mock_bins)
        self.assertDictEqual(obs_map, {bin_id: ["sample1;contig1", "sample2;contig1"]})

    @patch("q2_assembly.filter.filter_contigs")
    @patch("q2_mag.semibin2.semibin2._bin_contigs_semibin2")
    @patch("rachis.sdk.Context.get_action")
    def test_bin_contigs_semibin2_pipeline_multi_sample(self, p1, p2, p3):
        self.package = "q2_mag.tests"
        contigs = Artifact.import_data(
            "FeatureData[Contig % Properties('pooled')]",
            self.get_data_path("pooled-contigs/colon-separator"),
        )
        maps = Artifact.import_data(
            "SampleData[AlignmentMap % Properties('sorted', 'pooled')]",
            self.get_data_path("pooled-maps/colon-separator"),
        )
        contig_map = Artifact.import_data(
            "FeatureMap[MAGtoContigs]",
            {
                "522775d4-b1c6-4ee3-8b47-cd990f17eb8b": [
                    "sample2:contig1",
                    "sample1:contig1",
                ],
            },
        )
        p1.side_effect = [p2, p3]
        p2.return_value = ("mags", contig_map)
        p3.return_value = ("unbinned_contigs",)

        obs = bin_contigs_semibin2(
            Context(), contigs, maps, epochs=10, num_partitions=2
        )

        p1.assert_has_calls(
            [
                call("mag", "_bin_contigs_semibin2"),
                call("assembly", "filter_contigs"),
            ]
        )
        p2.assert_called_once_with(
            contigs,
            maps,
            "multi-sample",
            training_type=None,
            orf_finder="fast-naive",
            environment="global",
            engine="auto",
            sequencing_type="short_read",
            minfasta_kbs=200,
            no_recluster=False,
            epochs=10,
            batch_size=2048,
            max_node=1,
            max_edges=200,
            ratio=0.05,
            threads=1,
            min_len=None,
            ml_threshold=None,
            random_seed=None,
            debug=False,
        )
        p3.assert_called_once_with(
            contigs,
            on="contig",
            ids=["sample1:contig1", "sample2:contig1"],
            exclude_ids=True,
        )
        self.assertTupleEqual(("mags", contig_map, "unbinned_contigs"), obs)

    @patch("q2_mag.semibin2.semibin2.MultiFASTADirectoryFormat")
    @patch("q2_mag.semibin2.semibin2._process_sample")
    def test_bin_partition_contigs_semibin2_coassembly(self, p1, p2):
        self.package = "q2_mag.tests"
        contigs = ContigSequencesDirFmt(
            self.get_data_path("coassembled-contigs"), mode="r"
        )
        maps = BAMDirFmt(self.get_data_path("coassembled-maps"), mode="r")
        args = ["--epochs", "15"]

        mock_bins = MultiFASTADirectoryFormat()
        bin_id = "522775d4-b1c6-4ee3-8b47-cd990f17eb8b"
        os.makedirs(mock_bins.path / "pooled")
        shutil.copy(
            contigs.path / "pooled.fasta",
            mock_bins.path / "pooled" / f"{bin_id}.fa",
        )
        p2.return_value = mock_bins

        obs_bins, obs_map = _bin_partition_contigs_semibin2(
            contigs, maps, "co-assembly", args
        )

        p1.assert_called_once_with(
            "pooled",
            {
                "contigs": self.get_data_path("coassembled-contigs/pooled.fasta"),
                "map": [
                    self.get_data_path("coassembled-maps/sample1.bam"),
                    self.get_data_path("coassembled-maps/sample2.bam"),
                ],
            },
            "single_easy_bin",
            args,
            str(mock_bins),
            separator="",
        )
        self.assertEqual(obs_bins, mock_bins)
        self.assertDictEqual(obs_map, {bin_id: ["pooled:contig1", "pooled:contig2"]})

    @patch("q2_assembly.filter.filter_contigs")
    @patch("q2_mag.semibin2.semibin2._bin_contigs_semibin2")
    @patch("rachis.sdk.Context.get_action")
    def test_bin_contigs_semibin2_pipeline_coassembly(self, p1, p2, p3):
        self.package = "q2_mag.tests"
        contigs = Artifact.import_data(
            "FeatureData[Contig % Properties('coassembled')]",
            self.get_data_path("coassembled-contigs"),
        )
        maps = Artifact.import_data(
            "SampleData[AlignmentMap % Properties('sorted', 'pooled')]",
            self.get_data_path("coassembled-maps"),
        )
        contig_map = Artifact.import_data(
            "FeatureMap[MAGtoContigs]",
            {
                "522775d4-b1c6-4ee3-8b47-cd990f17eb8b": [
                    "pooled:contig2",
                    "pooled:contig1",
                ],
            },
        )
        p1.side_effect = [p2, p3]
        p2.return_value = ("mags", contig_map)
        p3.return_value = ("unbinned_contigs",)

        obs = bin_contigs_semibin2(
            Context(), contigs, maps, epochs=10, num_partitions=2
        )

        p1.assert_has_calls(
            [
                call("mag", "_bin_contigs_semibin2"),
                call("assembly", "filter_contigs"),
            ]
        )
        p2.assert_called_once_with(
            contigs,
            maps,
            "co-assembly",
            training_type=None,
            orf_finder="fast-naive",
            environment="global",
            engine="auto",
            sequencing_type="short_read",
            minfasta_kbs=200,
            no_recluster=False,
            epochs=10,
            batch_size=2048,
            max_node=1,
            max_edges=200,
            ratio=0.05,
            threads=1,
            min_len=None,
            ml_threshold=None,
            random_seed=None,
            debug=False,
        )
        p3.assert_called_once_with(
            contigs,
            on="contig",
            ids=["pooled:contig1", "pooled:contig2"],
            exclude_ids=True,
        )
        self.assertTupleEqual(("mags", contig_map, "unbinned_contigs"), obs)


if __name__ == "__main__":
    unittest.main()
