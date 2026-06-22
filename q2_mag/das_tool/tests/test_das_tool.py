# ----------------------------------------------------------------------------
# Copyright (c) 2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import glob
import os
import tempfile
import unittest
from unittest.mock import patch

from q2_types.per_sample_sequences import (
    ContigSequencesDirFmt,
    MultiFASTADirectoryFormat,
)
from qiime2.plugin.testing import TestPluginBase

from q2_mag.das_tool.das_tool import (
    _process_das_tool_arg,
    _write_contig2bin_map,
    refine_bins_das_tool,
)


class TestDASTool(TestPluginBase):
    package = "q2_annotate.tests"

    def test_process_das_tool_arg(self):
        self.assertEqual(
            _process_das_tool_arg("score_threshold", 0.6),
            ["--score_threshold", "0.6"],
        )
        self.assertEqual(_process_das_tool_arg("debug", True), ["--debug"])

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

    @patch("subprocess.run")
    def test_refine_bins_das_tool(self, subp_run):
        bins = MultiFASTADirectoryFormat(self.get_data_path("sample_data_mags"), "r")
        contigs = ContigSequencesDirFmt(self.get_data_path("contigs"), "r")

        def _mock_das_tool(cmd, check):
            output_prefix = cmd[cmd.index("--outputbasename") + 1]
            output_dir = f"{output_prefix}_DASTool_bins"
            os.makedirs(output_dir)
            with open(os.path.join(output_dir, "refined.fa"), "w") as fh:
                fh.write(">NZ_00000000.1_contig1\nACGT\n")

        subp_run.side_effect = _mock_das_tool

        obs = refine_bins_das_tool(
            bins=[bins, bins],
            contigs=contigs,
            search_engine="diamond",
            score_threshold=0.6,
            threads=2,
            debug=True,
        )

        self.assertIsInstance(obs, MultiFASTADirectoryFormat)
        self.assertEqual(len(subp_run.call_args_list), 2)

        first_cmd = subp_run.call_args_list[0].args[0]
        self.assertEqual(first_cmd[0], "DAS_Tool")
        self.assertIn("--write_bins", first_cmd)
        self.assertIn("--bins", first_cmd)
        self.assertIn("--labels", first_cmd)
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
        self.assertEqual(len(obs_bins), 2)
        self.assertTrue(obs_bins[0].startswith("sample1/"))
        self.assertTrue(obs_bins[1].startswith("sample2/"))

    def test_refine_bins_das_tool_requires_two_binnings(self):
        bins = MultiFASTADirectoryFormat(self.get_data_path("sample_data_mags"), "r")
        contigs = ContigSequencesDirFmt(self.get_data_path("contigs"), "r")

        with self.assertRaisesRegex(ValueError, "at least two binning methods"):
            refine_bins_das_tool(bins=[bins], contigs=contigs)


if __name__ == "__main__":
    unittest.main()
