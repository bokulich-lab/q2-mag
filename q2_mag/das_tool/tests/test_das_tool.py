# ----------------------------------------------------------------------------
# Copyright (c) 2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import glob
import io
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest
import warnings
from unittest.mock import patch, ANY, call

import pandas as pd
from q2_types.genome_data import ProteinsDirectoryFormat
from q2_types.per_sample_sequences import (
    ContigSequencesDirFmt,
    MultiFASTADirectoryFormat,
)
from qiime2.plugin.testing import TestPluginBase

from q2_mag.das_tool.das_tool import (
    _append_summary,
    _get_sample_ids,
    _generate_labels,
    _parse_labels,
    _process_das_tool_arg,
    _run_das_tool,
    _write_contig2bin_map,
    _refine_bins_das_tool,
    refine_bins_das_tool,
)


class TestDASTool(TestPluginBase):
    package = "q2_mag.das_tool.tests"

    def _get_concatenated_outputs(self):
        summary = pd.read_csv(
            self.get_data_path("summaries/summary_concatenated.tsv"),
            sep="\t",
            index_col="id",
        )
        evaluation = pd.read_csv(
            self.get_data_path("summaries/input_bins_evaluation_concatenated.eval"),
            sep="\t",
            index_col="id",
        )
        return summary, evaluation

    def test_process_das_tool_arg(self):
        self.assertEqual(
            _process_das_tool_arg("score_threshold", 0.6),
            ["--score_threshold", "0.6"],
        )
        self.assertEqual(_process_das_tool_arg("debug", True), ["--debug"])

    def test_parse_labels(self):
        self.assertEqual(_parse_labels("metabat,semibin", 2), ["metabat", "semibin"])

    def test_parse_labels_wrong_count(self):
        with self.assertRaisesRegex(
            ValueError, "number of labels provided is different"
        ):
            _parse_labels("metabat,semibin", 3)

    def test_parse_labels_duplicate_labels(self):
        with self.assertRaisesRegex(ValueError, "Duplicate labels detected"):
            _parse_labels("metabat,metabat", 2)

    def test_generate_labels(self):
        self.assertEqual(_generate_labels(3), ["binning_1", "binning_2", "binning_3"])

    def test_get_sample_ids_inconsistent_bin_samples(self):
        input_contigs = self.get_data_path("contigs")
        input_bins_1 = self.get_data_path("bins-single-sample")
        input_bins_2 = self.get_data_path("bins")
        contigs = ContigSequencesDirFmt(input_contigs, mode="r")
        bins_1 = MultiFASTADirectoryFormat(input_bins_1, mode="r")
        bins_2 = MultiFASTADirectoryFormat(input_bins_2, mode="r")

        with self.assertRaisesRegex(
            ValueError,
            r"Sample IDs must stay consistent across all bins\.\n"
            r"Missing sample IDs by bin:\n"
            r"- metabat: samp2",
        ):
            _get_sample_ids(
                ["metabat", "semibin"],
                contigs,
                None,
                bins_1,
                bins_2,
            )

    def test_get_sample_ids_missing_contig_sample(self):
        input_contigs = self.get_data_path("contigs-single-sample")
        input_bins = self.get_data_path("bins")
        contigs = ContigSequencesDirFmt(input_contigs, mode="r")
        bins = MultiFASTADirectoryFormat(input_bins, mode="r")

        with self.assertRaisesRegex(
            ValueError,
            "Missing from contigs: samp2",
        ):
            _get_sample_ids(
                ["metabat"],
                contigs,
                None,
                bins,
            )

    def test_get_sample_ids_missing_protein_sample(self):
        input_contigs = self.get_data_path("contigs")
        input_bins = self.get_data_path("bins")
        input_proteins = self.get_data_path("proteins-single-sample")
        contigs = ContigSequencesDirFmt(input_contigs, mode="r")
        bins = MultiFASTADirectoryFormat(input_bins, mode="r")
        proteins = ProteinsDirectoryFormat(input_proteins, mode="r")

        with self.assertRaisesRegex(
            ValueError,
            "Missing from proteins: samp2",
        ):
            _get_sample_ids(
                ["metabat"],
                contigs,
                proteins,
                bins,
            )

    def test_get_sample_ids_allows_extra_contig_and_protein_samples(self):
        input_contigs = self.get_data_path("contigs")
        input_bins = self.get_data_path("bins-single-sample")
        input_proteins = self.get_data_path("proteins")
        contigs = ContigSequencesDirFmt(input_contigs, mode="r")
        bins = MultiFASTADirectoryFormat(input_bins, mode="r")
        proteins = ProteinsDirectoryFormat(input_proteins, mode="r")

        observed = _get_sample_ids(
            ["metabat"],
            contigs,
            proteins,
            bins,
        )

        self.assertEqual(observed, ["samp1"])

    def test_write_contig2bin_map(self):
        input_bins = self.get_data_path("bins")
        bins = MultiFASTADirectoryFormat(input_bins, mode="r")

        with tempfile.TemporaryDirectory() as fake_loc:
            obs = _write_contig2bin_map(bins, "samp1", "metabat", fake_loc)
            with open(obs) as fh:
                lines = sorted(line.strip().split("\t") for line in fh)

        self.assertIn(
            [
                "NODE_12_length_1973_cov_1.555266",
                "522775d4-b1c6-4ee3-8b47-cd990f17eb8b",
            ],
            lines,
        )

    def test_append_summary(self):
        summary1_path = self.get_data_path("summaries/run1.tsv")
        summary2_path = self.get_data_path("summaries/run2.tsv")
        runs_concatenated_path = self.get_data_path("summaries/runs_concatenated.tsv")

        summaries = _append_summary(
            sample_id="sample1", summary=summary1_path, summaries=None
        )
        summaries = _append_summary(
            sample_id="sample2", summary=summary2_path, summaries=summaries
        )

        obs = summaries
        exp = pd.read_csv(
            runs_concatenated_path,
            sep="\t",
            index_col="id",
        )

        pd.testing.assert_frame_equal(obs, exp)

    @patch("q2_mag.das_tool.das_tool._collect_refined_bins")
    @patch("q2_mag.das_tool.das_tool._append_summary")
    @patch("q2_mag.das_tool.das_tool._run_das_tool")
    @patch("q2_mag.das_tool.das_tool._get_sample_ids")
    @patch("q2_mag.das_tool.das_tool._parse_labels")
    def test_refine_bins_das_tool(self, p1, p2, p3, p4, p5):
        bins1_dirpath = self.get_data_path("mags/binning_1")
        bins2_dirpath = self.get_data_path("mags/binning_2")
        contigs_dirpath = self.get_data_path("contigs")
        bins1 = MultiFASTADirectoryFormat(bins1_dirpath, mode="r")
        bins2 = MultiFASTADirectoryFormat(bins2_dirpath, mode="r")
        contigs = ContigSequencesDirFmt(contigs_dirpath, mode="r")
        cat_summary, cat_evaluation = self._get_concatenated_outputs()

        p1.return_value = ["binning_1", "binning_2"]
        p2.return_value = ["samp2"]
        p3.return_value = (
            "some/where/refined_bins_dir",
            self.get_data_path("summaries/summary.tsv"),
            self.get_data_path("summaries/input_bins_evaluation.eval"),
        )
        p4.side_effect = [cat_summary, cat_evaluation]
        p5.return_value = 1

        args = [
            "--search_engine",
            "diamond",
            "--score_threshold",
            "0.1",
            "--threads",
            "2",
            "--debug",
        ]

        obs_bins, obs_summary, obs_evaluation = _refine_bins_das_tool(
            bins=[bins1, bins2],
            contigs=contigs,
            proteins=None,
            labels="binning_1,binning_2",
            common_args=args,
        )

        self.assertIsInstance(obs_bins, MultiFASTADirectoryFormat)

        self.assertEqual(
            set(obs_summary.to_dataframe()["sample_id"]),
            {"samp2"},
        )
        self.assertEqual(
            set(obs_evaluation.to_dataframe()["sample_id"]),
            {"samp2"},
        )

        p1.assert_called_once_with("binning_1,binning_2", 2)
        p2.assert_called_once_with(
            ["binning_1", "binning_2"],
            contigs,
            None,
            bins1,
            bins2,
        )
        p3.assert_called_once_with(
            sample_id="samp2",
            bins=[bins1, bins2],
            contigs_fp=self.get_data_path("contigs/samp2_contigs.fa"),
            proteins_fp=None,
            labels=["binning_1", "binning_2"],
            common_args=args,
            output_dir=ANY,
        )
        p4.assert_has_calls(
            [
                call(
                    "samp2",
                    self.get_data_path("summaries/summary.tsv"),
                    None,
                ),
                call(
                    "samp2",
                    self.get_data_path("summaries/input_bins_evaluation.eval"),
                    None,
                ),
            ]
        )
        p5.assert_called_once_with("samp2", ANY, obs_bins)

    @patch("q2_mag.das_tool.das_tool._collect_refined_bins")
    @patch("q2_mag.das_tool.das_tool._append_summary")
    @patch("q2_mag.das_tool.das_tool._run_das_tool")
    @patch("q2_mag.das_tool.das_tool._get_sample_ids")
    @patch("q2_mag.das_tool.das_tool._parse_labels")
    def test_refine_bins_das_tool_warns_and_continues_when_samples_fail(
        self, p1, p2, p3, p4, p5
    ):
        bins1_dirpath = self.get_data_path("mags/binning_1")
        bins2_dirpath = self.get_data_path("mags/binning_2")
        contigs_dirpath = self.get_data_path("contigs")
        bins1 = MultiFASTADirectoryFormat(bins1_dirpath, mode="r")
        bins2 = MultiFASTADirectoryFormat(bins2_dirpath, mode="r")
        contigs = ContigSequencesDirFmt(contigs_dirpath, mode="r")
        cat_summary, cat_evaluation = self._get_concatenated_outputs()

        p1.return_value = ["binning_1", "binning_2"]
        p2.return_value = ["samp1", "samp2"]
        p3.side_effect = [
            subprocess.CalledProcessError(
                1,
                ["DAS_Tool"],
                output="Error:  No bins with bin-score >0.1 found.",
            ),
            (
                "some/where/refined_bins_dir",
                self.get_data_path("summaries/summary.tsv"),
                self.get_data_path("summaries/input_bins_evaluation.eval"),
            ),
        ]
        p4.side_effect = [cat_summary, cat_evaluation]
        p5.side_effect = [0, 1]

        args = [
            "--search_engine",
            "diamond",
            "--score_threshold",
            "0.1",
            "--threads",
            "2",
            "--debug",
        ]

        with self.assertWarnsRegex(UserWarning, r"No bins produced for sample samp1\."):
            obs_bins, _, _ = _refine_bins_das_tool(
                bins=[bins1, bins2],
                contigs=contigs,
                proteins=None,
                labels="binning_1,binning_2",
                common_args=args,
            )

        self.assertIsInstance(obs_bins, MultiFASTADirectoryFormat)

        p1.assert_called_once_with("binning_1,binning_2", 2)
        p2.assert_called_once_with(
            ["binning_1", "binning_2"],
            contigs,
            None,
            bins1,
            bins2,
        )
        self.assertEqual(
            p3.call_args_list,
            [
                call(
                    sample_id="samp1",
                    bins=[bins1, bins2],
                    contigs_fp=self.get_data_path("contigs/samp1_contigs.fa"),
                    proteins_fp=None,
                    labels=["binning_1", "binning_2"],
                    common_args=args,
                    output_dir=ANY,
                ),
                call(
                    sample_id="samp2",
                    bins=[bins1, bins2],
                    contigs_fp=self.get_data_path("contigs/samp2_contigs.fa"),
                    proteins_fp=None,
                    labels=["binning_1", "binning_2"],
                    common_args=args,
                    output_dir=ANY,
                ),
            ],
        )
        self.assertEqual(
            p4.call_args_list,
            [
                call(
                    "samp2",
                    self.get_data_path("summaries/summary.tsv"),
                    None,
                ),
                call(
                    "samp2",
                    self.get_data_path("summaries/input_bins_evaluation.eval"),
                    None,
                ),
            ],
        )
        self.assertEqual(
            p5.call_args_list,
            [
                call("samp1", ANY, obs_bins),
                call("samp2", ANY, obs_bins),
            ],
        )

    @patch("q2_mag.das_tool.das_tool._collect_refined_bins")
    @patch("q2_mag.das_tool.das_tool._append_summary")
    @patch("q2_mag.das_tool.das_tool._run_das_tool")
    @patch("q2_mag.das_tool.das_tool._get_sample_ids")
    @patch("q2_mag.das_tool.das_tool._parse_labels")
    def test_run_das_tool_with_proteins(self, p1, p2, p3, p4, p5):
        bins1_dirpath = self.get_data_path("mags/binning_1")
        bins2_dirpath = self.get_data_path("mags/binning_2")
        contigs_dirpath = self.get_data_path("contigs")
        bins1 = MultiFASTADirectoryFormat(bins1_dirpath, mode="r")
        bins2 = MultiFASTADirectoryFormat(bins2_dirpath, mode="r")
        contigs = ContigSequencesDirFmt(contigs_dirpath, mode="r")
        proteins_dirpath = self.get_data_path("proteins")
        proteins = ProteinsDirectoryFormat(proteins_dirpath, mode="r")
        cat_summary, cat_evaluation = self._get_concatenated_outputs()

        p1.return_value = ["binning_1", "binning_2"]
        p2.return_value = ["samp2"]
        p3.return_value = (
            "some/where/refined_bins_dir",
            self.get_data_path("summaries/summary.tsv"),
            self.get_data_path("summaries/input_bins_evaluation.eval"),
        )
        p4.side_effect = [cat_summary, cat_evaluation]
        p5.return_value = 1

        args = [
            "--search_engine",
            "diamond",
            "--score_threshold",
            "0.1",
            "--threads",
            "2",
            "--debug",
        ]

        obs_bins, obs_summary, obs_evaluation = _refine_bins_das_tool(
            bins=[bins1, bins2],
            contigs=contigs,
            proteins=proteins,
            labels="binning_1,binning_2",
            common_args=args,
        )

        self.assertIsInstance(obs_bins, MultiFASTADirectoryFormat)

        self.assertEqual(
            set(obs_summary.to_dataframe()["sample_id"]),
            {"samp2"},
        )
        self.assertEqual(
            set(obs_evaluation.to_dataframe()["sample_id"]),
            {"samp2"},
        )

        p1.assert_called_once_with("binning_1,binning_2", 2)
        p2.assert_called_once_with(
            ["binning_1", "binning_2"],
            contigs,
            proteins,
            bins1,
            bins2,
        )
        p3.assert_called_once_with(
            sample_id="samp2",
            bins=[bins1, bins2],
            contigs_fp=self.get_data_path("contigs/samp2_contigs.fa"),
            proteins_fp=self.get_data_path("proteins/samp2.fasta"),
            labels=["binning_1", "binning_2"],
            common_args=args,
            output_dir=ANY,
        )
        p4.assert_has_calls(
            [
                call(
                    "samp2",
                    self.get_data_path("summaries/summary.tsv"),
                    None,
                ),
                call(
                    "samp2",
                    self.get_data_path("summaries/input_bins_evaluation.eval"),
                    None,
                ),
            ]
        )
        p5.assert_called_once_with("samp2", ANY, obs_bins)

    @patch("sys.stderr", new_callable=io.StringIO)
    @patch("sys.stdout", new_callable=io.StringIO)
    @patch("q2_mag.das_tool.das_tool.run_command")
    def test_run_das_tool_prints_streams_on_error(self, p1, p2, p3):
        bins1_dirpath = self.get_data_path("mags/binning_1")
        bins2_dirpath = self.get_data_path("mags/binning_2")
        bins1 = MultiFASTADirectoryFormat(bins1_dirpath, mode="r")
        bins2 = MultiFASTADirectoryFormat(bins2_dirpath, mode="r")
        contigs_fp = self.get_data_path("contigs/samp2_contigs.fa")

        p1.side_effect = subprocess.CalledProcessError(
            1,
            ["DAS_Tool"],
            output="DAS Tool output\n",
            stderr="DAS Tool error\n",
        )

        with self.assertRaises(subprocess.CalledProcessError):
            _run_das_tool(
                sample_id="samp1",
                bins=[bins1, bins2],
                contigs_fp=contigs_fp,
                proteins_fp=None,
                labels=["binning_1", "binning_2"],
                common_args=["--threads", "1"],
                output_dir=self.temp_dir.name,
            )

        self.assertEqual(p2.getvalue(), "DAS Tool output\n")
        self.assertEqual(p3.getvalue(), "DAS Tool error\n")

    @patch("subprocess.run")
    def test_refine_bins_das_tool_fails_when_all_samples_fail(self, p1):
        bins1_dirpath = self.get_data_path("mags/binning_1")
        bins2_dirpath = self.get_data_path("mags/binning_2")
        contigs_dirpath = self.get_data_path("contigs")
        bins1 = MultiFASTADirectoryFormat(bins1_dirpath, mode="r")
        bins2 = MultiFASTADirectoryFormat(bins2_dirpath, mode="r")
        contigs = ContigSequencesDirFmt(contigs_dirpath, mode="r")

        p1.side_effect = subprocess.CalledProcessError(
            1,
            ["DAS_Tool"],
            output=(
                "Error:  No bins with bin-score >0.6 found. Adjust score_threshold "
                "to report bins with lower quality.\n"
            ),
        )

        with warnings.catch_warnings(record=True) as warning_records:
            with self.assertRaisesRegex(ValueError, "No refined MAGs were formed"):
                refine_bins_das_tool(
                    bins=[bins1, bins2],
                    contigs=contigs,
                    score_threshold=0.6,
                )

        warning_messages = [str(record.message) for record in warning_records]
        self.assertIn("No bins produced for sample samp1.", warning_messages)
        self.assertIn("No bins produced for sample(s): samp1.", warning_messages)

    @patch("subprocess.run")
    def test_refine_bins_das_tool_raises_on_other_errors(self, p1):
        bins1_dirpath = self.get_data_path("mags/binning_1")
        bins2_dirpath = self.get_data_path("mags/binning_2")
        contigs_dirpath = self.get_data_path("contigs")
        bins1 = MultiFASTADirectoryFormat(bins1_dirpath, mode="r")
        bins2 = MultiFASTADirectoryFormat(bins2_dirpath, mode="r")
        contigs = ContigSequencesDirFmt(contigs_dirpath, mode="r")

        p1.side_effect = subprocess.CalledProcessError(
            1,
            ["DAS_Tool"],
            output="Error: Another DAS Tool error.\n",
        )

        with self.assertRaisesRegex(Exception, "error was encountered"):
            refine_bins_das_tool(
                bins=[bins1, bins2],
                contigs=contigs,
                score_threshold=0.6,
            )

    def test_refine_bins_das_tool_requires_two_binnings(self):
        input_contigs = self.get_data_path("contigs")
        input_bins = self.get_data_path("bins")
        contigs = ContigSequencesDirFmt(input_contigs, mode="r")
        bins = MultiFASTADirectoryFormat(input_bins, mode="r")

        with self.assertRaisesRegex(ValueError, "at least two binning methods"):
            refine_bins_das_tool(bins=[bins], contigs=contigs)

    def test_run_das_tool(self):
        # TODO: Rearrange the tests according to the order they occur.
        # TODO: See if more processes can be patched.
        # TODO: Check coverage for current tests.
        bins1_dirpath = self.get_data_path("mags/binning_1")
        bins2_dirpath = self.get_data_path("mags/binning_2")
        bins1 = MultiFASTADirectoryFormat(bins1_dirpath, mode="r")
        bins2 = MultiFASTADirectoryFormat(bins2_dirpath, mode="r")
        contigs_fp = self.get_data_path("contigs/samp2_contigs.fa")

        args = [
            "--search_engine",
            "diamond",
            "--score_threshold",
            "0.1",
            "--threads",
            "2",
            "--debug",
        ]

        with tempfile.TemporaryDirectory() as fake_loc:
            obs_bins, obs_summary, obs_evaluation = _run_das_tool(
                sample_id="samp1",
                bins=[
                    bins1,
                    bins2,
                ],
                contigs_fp=contigs_fp,
                proteins_fp=None,
                labels=["binning_1", "binning_2"],
                common_args=args,
                output_dir=fake_loc,
            )

            obs_bin_fps = [
                fp.name
                for fp in Path(obs_bins).iterdir()
                if fp.suffix in {".fa", ".fasta", ".fna"}
            ]
            self.assertEqual(len(obs_bin_fps), 1)
            obs_bin_fp = Path(obs_bins) / obs_bin_fps[0]
            obs_summary_fp = Path(obs_summary)
            obs_evaluation_fp = Path(obs_evaluation)

            exp_bin_fp = Path(
                self.get_data_path(
                    "refined_bins/samp2/" "05684760-dd45-4e0d-a5d4-63483f4b92f9.fasta"
                )
            )
            exp_summary_fp = Path(self.get_data_path("summaries/summary.tsv"))
            exp_evaluation_fp = Path(
                self.get_data_path("summaries/input_bins_evaluation.eval")
            )

            self.assertMultiLineEqual(
                exp_bin_fp.read_text(),
                obs_bin_fp.read_text(),
            )
            self.assertMultiLineEqual(
                exp_summary_fp.read_text(),
                obs_summary_fp.read_text(),
            )
            self.assertMultiLineEqual(
                exp_evaluation_fp.read_text(),
                obs_evaluation_fp.read_text(),
            )


if __name__ == "__main__":
    unittest.main()
