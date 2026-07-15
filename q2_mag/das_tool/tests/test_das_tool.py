# ----------------------------------------------------------------------------
# Copyright (c) 2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import glob
import os
import shutil
import subprocess
import tempfile
import unittest
import warnings
from unittest.mock import patch

import pandas as pd
from q2_types.per_sample_sequences import (
    ContigSequencesDirFmt,
    MultiFASTADirectoryFormat,
)
from qiime2.plugin.testing import TestPluginBase

from q2_mag.das_tool.das_tool import (
    _append_summary,
    _get_sample_ids,
    _get_sample_proteins,
    _generate_labels,
    _parse_labels,
    _process_das_tool_arg,
    _run_das_tool,
    _write_contig2bin_map,
    refine_bins_das_tool,
)


class _Bins:
    def __init__(self, bins, sample_id="samp1"):
        if isinstance(bins, dict) and all(isinstance(v, dict) for v in bins.values()):
            self._sample_dict = bins
        else:
            self._sample_dict = {sample_id: bins}

    def sample_dict(self):
        return self._sample_dict


class _Contigs:
    def __init__(self, contig_fp, sample_id="samp1"):
        if isinstance(contig_fp, dict):
            self._sample_dict = contig_fp
        else:
            self._sample_dict = {sample_id: contig_fp}

    def sample_dict(self):
        return self._sample_dict


class TestDASTool(TestPluginBase):
    package = "q2_mag.das_tool.tests"

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

    def test_get_sample_ids_mismatched_samples(self):
        with self.assertRaisesRegex(ValueError, "same sample IDs"):
            _get_sample_ids(
                _Contigs({"samp1": "samp1_contigs.fa"}),
                _Bins({"samp2": {"bin_1": "bin_1.fa"}}),
            )

    def test_get_sample_proteins(self):
        with tempfile.TemporaryDirectory() as tempdir:
            expected = {
                "samp1": os.path.join(tempdir, "samp1.fa"),
                "samp2": os.path.join(tempdir, "samp2.faa"),
                "samp3": os.path.join(tempdir, "samp3.fasta"),
            }

            for fp in expected.values():
                with open(fp, "w") as fh:
                    fh.write(">protein\nM\n")

            with open(os.path.join(tempdir, "samp4.txt"), "w") as fh:
                fh.write("ignored\n")

            obs = _get_sample_proteins(tempdir)

        self.assertEqual(obs, expected)

    def test_write_contig2bin_map(self):
        bins = MultiFASTADirectoryFormat(self.get_data_path("sample_data_mags"), "r")

        with tempfile.TemporaryDirectory() as tempdir:
            obs = _write_contig2bin_map(bins, "sample1", "metabat", tempdir)
            with open(obs) as fh:
                lines = sorted(line.strip().split("\t") for line in fh)

        self.assertIn(
            ["NZ_00000000.1_contig1", "24dee6fe-9b84-45bb-8145-de7b092533a1"],
            lines,
        )

    def test_append_summary(self):
        with tempfile.TemporaryDirectory() as tempdir:
            summary1 = os.path.join(tempdir, "sample1_summary.tsv")
            summary2 = os.path.join(tempdir, "sample2_summary.tsv")

            pd.DataFrame(
                {
                    "bin": ["metabat.1", "metabat.2"],
                    "bin_set": ["metabat", "metabat"],
                    "bin_score": ["1", "0.8"],
                }
            ).to_csv(summary1, sep="\t", index=False)
            pd.DataFrame(
                {"bin": ["semibin.1"], "bin_set": ["semibin"], "bin_score": ["0.9"]}
            ).to_csv(summary2, sep="\t", index=False)

            summaries = _append_summary("sample1", summary1)
            summaries = _append_summary("sample2", summary2, summaries)

        self.assertEqual(summaries.index.name, "id")
        self.assertEqual(list(summaries.index), ["metabat.1", "metabat.2", "semibin.1"])
        self.assertEqual(
            summaries.to_dict(orient="records"),
            [
                {
                    "sample_id": "sample1",
                    "bin_set": "metabat",
                    "bin_score": 1.0,
                },
                {
                    "sample_id": "sample1",
                    "bin_set": "metabat",
                    "bin_score": 0.8,
                },
                {
                    "sample_id": "sample2",
                    "bin_set": "semibin",
                    "bin_score": 0.9,
                },
            ],
        )

    @patch("subprocess.run")
    def test_refine_bins_das_tool(self, subp_run):
        bins_path = self.get_data_path("bins")
        contig_fp = os.path.join(self.get_data_path("contigs"), "samp1_contigs.fa")

        def _mock_das_tool(cmd, check):
            output_prefix = cmd[cmd.index("--outputbasename") + 1]
            output_dir = f"{output_prefix}_DASTool_bins"
            os.makedirs(output_dir)
            with open(os.path.join(output_dir, "refined.fa"), "w") as fh:
                fh.write(">NZ_00000000.1_contig1\nACGT\n")
            pd.DataFrame(
                {"bin": ["refined"], "bin_set": ["DASTool"], "bin_score": ["1"]}
            ).to_csv(f"{output_prefix}_DASTool_summary.tsv", sep="\t", index=False)
            pd.DataFrame(
                {"bin": ["input"], "bin_set": ["metabat"], "bin_score": ["0.9"]}
            ).to_csv(f"{output_prefix}_allBins.eval", sep="\t", index=False)

        subp_run.side_effect = _mock_das_tool

        obs, summary, input_bins_evaluation = refine_bins_das_tool(
            bins=[
                _Bins({"bin_1": os.path.join(bins_path, "bin_1_samp1.fa")}),
                _Bins({"bin_2": os.path.join(bins_path, "bin_2_samp1.fa")}),
            ],
            contigs=_Contigs(contig_fp),
            search_engine="diamond",
            score_threshold=0.6,
            threads=2,
            debug=True,
            labels="metabat,semibin",
        )

        self.assertIsInstance(obs, MultiFASTADirectoryFormat)
        self.assertEqual(list(summary.to_dataframe().index), ["refined"])
        self.assertEqual(list(summary.to_dataframe()["sample_id"]), ["samp1"])
        self.assertEqual(list(input_bins_evaluation.to_dataframe().index), ["input"])
        self.assertEqual(
            list(input_bins_evaluation.to_dataframe()["sample_id"]),
            ["samp1"],
        )
        self.assertEqual(len(subp_run.call_args_list), 1)

        first_cmd = subp_run.call_args_list[0].args[0]
        self.assertEqual(first_cmd[0], "DAS_Tool")
        self.assertIn("--write_bins", first_cmd)
        self.assertIn("--write_bin_evals", first_cmd)
        self.assertIn("--bins", first_cmd)
        self.assertIn("--labels", first_cmd)
        self.assertEqual(first_cmd[first_cmd.index("--labels") + 1], "metabat,semibin")
        self.assertIn("--search_engine", first_cmd)
        self.assertIn("diamond", first_cmd)
        self.assertIn("--score_threshold", first_cmd)
        self.assertIn("0.6", first_cmd)
        self.assertIn("--threads", first_cmd)
        self.assertIn("2", first_cmd)
        self.assertIn("--debug", first_cmd)

        obs_bins = sorted(
            "/".join(fp.split("/")[-2:])
            for fp in glob.glob(os.path.join(str(obs), "*", "*.fa"))
        )
        self.assertEqual(len(obs_bins), 1)
        self.assertTrue(obs_bins[0].startswith("samp1/"))

    @patch("subprocess.run")
    def test_refine_bins_das_tool_warns_and_continues_when_samples_fail(self, subp_run):
        bins_path = self.get_data_path("bins")
        contig_fp = os.path.join(self.get_data_path("contigs"), "samp1_contigs.fa")
        bin_1 = os.path.join(bins_path, "bin_1_samp1.fa")
        bin_2 = os.path.join(bins_path, "bin_2_samp1.fa")

        def _mock_das_tool(cmd, check):
            output_prefix = cmd[cmd.index("--outputbasename") + 1]
            if output_prefix.endswith("samp1"):
                output_dir = f"{output_prefix}_DASTool_bins"
                os.makedirs(output_dir)
                with open(os.path.join(output_dir, "refined.fa"), "w") as fh:
                    fh.write(">NZ_00000000.1_contig1\nACGT\n")
                pd.DataFrame(
                    {"bin": ["refined"], "bin_set": ["DASTool"], "bin_score": ["1"]}
                ).to_csv(f"{output_prefix}_DASTool_summary.tsv", sep="\t", index=False)
                pd.DataFrame(
                    {"bin": ["input"], "bin_set": ["metabat"], "bin_score": ["0.9"]}
                ).to_csv(f"{output_prefix}_allBins.eval", sep="\t", index=False)
            else:
                raise subprocess.CalledProcessError(1, cmd)

        subp_run.side_effect = _mock_das_tool

        with warnings.catch_warnings(record=True) as warning_records:
            warnings.simplefilter("always")
            obs, summary, input_bins_evaluation = refine_bins_das_tool(
                bins=[
                    _Bins(
                        {
                            "samp1": {"bin_1": bin_1},
                            "samp2": {"bin_1": bin_1},
                            "samp3": {"bin_1": bin_1},
                        }
                    ),
                    _Bins(
                        {
                            "samp1": {"bin_2": bin_2},
                            "samp2": {"bin_2": bin_2},
                            "samp3": {"bin_2": bin_2},
                        }
                    ),
                ],
                contigs=_Contigs(
                    {"samp1": contig_fp, "samp2": contig_fp, "samp3": contig_fp}
                ),
                score_threshold=0.6,
            )

        warning_messages = [str(record.message) for record in warning_records]
        self.assertTrue(
            any("lowering the score threshold" in msg for msg in warning_messages)
        )
        self.assertIn("No bins produced for sample(s): samp2, samp3", warning_messages)
        self.assertIsInstance(obs, MultiFASTADirectoryFormat)
        self.assertEqual(set(obs.sample_dict()), {"samp1", "samp2", "samp3"})
        self.assertEqual(obs.sample_dict()["samp2"], {})
        self.assertEqual(obs.sample_dict()["samp3"], {})
        self.assertEqual(list(summary.to_dataframe()["sample_id"]), ["samp1"])
        self.assertEqual(
            list(input_bins_evaluation.to_dataframe()["sample_id"]), ["samp1"]
        )
        self.assertEqual(len(subp_run.call_args_list), 3)

    @patch("q2_mag.das_tool.das_tool.run_command")
    def test_run_das_tool_with_proteins(self, run_command):
        contigs_path = self.get_data_path("contigs")
        bins_path = self.get_data_path("bins")

        contig_path = os.path.join(contigs_path, "samp1_contigs.fa")
        bin1_path = os.path.join(bins_path, "bin_1_samp1.fa")
        bin2_path = os.path.join(bins_path, "bin_2_samp1.fa")
        proteins_path = os.path.join(self.temp_dir.name, "samp1.faa")

        with open(proteins_path, "w") as fh:
            fh.write(">protein\nM\n")

        obs_bins, obs_summary, obs_evaluation = _run_das_tool(
            sample_id="samp1",
            bins=[
                _Bins({"bin_1": bin1_path}),
                _Bins({"bin_2": bin2_path}),
            ],
            contigs_fp=contig_path,
            proteins_fp=proteins_path,
            labels=["binning_1", "binning_2"],
            common_args=["--threads", "1"],
            output_dir=self.temp_dir.name,
        )

        cmd = run_command.call_args.args[0]
        self.assertIn("--proteins", cmd)
        self.assertEqual(cmd[cmd.index("--proteins") + 1], proteins_path)
        self.assertEqual(
            obs_bins, os.path.join(self.temp_dir.name, "samp1_DASTool_bins")
        )
        self.assertEqual(
            obs_summary, os.path.join(self.temp_dir.name, "samp1_DASTool_summary.tsv")
        )
        self.assertEqual(
            obs_evaluation, os.path.join(self.temp_dir.name, "samp1_allBins.eval")
        )

    @patch("subprocess.run")
    def test_refine_bins_das_tool_fails_when_all_samples_fail(self, subp_run):
        bins_path = self.get_data_path("bins")
        contig_fp = os.path.join(self.get_data_path("contigs"), "samp1_contigs.fa")

        subp_run.side_effect = subprocess.CalledProcessError(1, ["DAS_Tool"])

        with warnings.catch_warnings(record=True) as warning_records:
            warnings.simplefilter("always")
            with self.assertRaisesRegex(ValueError, "No refined MAGs were formed"):
                refine_bins_das_tool(
                    bins=[
                        _Bins({"bin_1": os.path.join(bins_path, "bin_1_samp1.fa")}),
                        _Bins({"bin_2": os.path.join(bins_path, "bin_2_samp1.fa")}),
                    ],
                    contigs=_Contigs(contig_fp),
                    score_threshold=0.6,
                )

        warning_messages = [str(record.message) for record in warning_records]
        self.assertTrue(
            any("lowering the score threshold" in msg for msg in warning_messages)
        )
        self.assertIn("No bins produced for sample(s): samp1", warning_messages)

    def test_refine_bins_das_tool_requires_two_binnings(self):
        bins = MultiFASTADirectoryFormat(self.get_data_path("sample_data_mags"), "r")
        contigs = ContigSequencesDirFmt(self.get_data_path("contigs"), "r")

        with self.assertRaisesRegex(ValueError, "at least two binning methods"):
            refine_bins_das_tool(bins=[bins], contigs=contigs)

    @unittest.skipUnless(shutil.which("DAS_Tool"), "DAS_Tool is not installed")
    def test_run_das_tool(self):
        contigs_path = self.get_data_path("contigs")
        bins_path = self.get_data_path("bins")

        contig_path = os.path.join(contigs_path, "samp1_contigs.fa")
        bin1_path = os.path.join(bins_path, "bin_1_samp1.fa")
        bin2_path = os.path.join(bins_path, "bin_2_samp1.fa")

        with tempfile.TemporaryDirectory() as tempdir:
            obs_bins, obs_summary, obs_evaluation = _run_das_tool(
                sample_id="samp1",
                bins=[
                    _Bins({"bin_1": bin1_path}),
                    _Bins({"bin_2": bin2_path}),
                ],
                contigs_fp=contig_path,
                proteins_fp=None,
                labels=["binning_1", "binning_2"],
                common_args=[
                    "--threads",
                    "1",
                    "--score_threshold",
                    "0.01",
                ],
                output_dir=tempdir,
            )

            self.assertEqual(obs_bins, os.path.join(tempdir, "samp1_DASTool_bins"))
            self.assertEqual(
                obs_summary, os.path.join(tempdir, "samp1_DASTool_summary.tsv")
            )
            self.assertEqual(
                obs_evaluation, os.path.join(tempdir, "samp1_allBins.eval")
            )
            self.assertTrue(
                os.path.exists(os.path.join(tempdir, "binning_1_contig2bin.tsv"))
            )
            self.assertTrue(
                os.path.exists(os.path.join(tempdir, "binning_2_contig2bin.tsv"))
            )
            self.assertTrue(os.path.exists(obs_summary))
            self.assertTrue(os.path.exists(obs_evaluation))


if __name__ == "__main__":
    unittest.main()
