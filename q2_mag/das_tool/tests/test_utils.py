# ----------------------------------------------------------------------------
# Copyright (c) 2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ---------------------------------------------------------------------------

import io
import subprocess
from contextlib import redirect_stderr, redirect_stdout

from rachis.plugin.testing import TestPluginBase

from q2_mag.das_tool.utils import _print_streams


class TestDASToolUtils(TestPluginBase):
    package = "q2_mag.das_tool.tests"

    def setUp(self):
        super().setUp()
        self.result = subprocess.CompletedProcess(
            args=["DAS_Tool", "--threads", "4"],
            returncode=0,
            stdout="DAS Tool version 1.1.7\n",
            stderr="No bins with bin-score >0.4 found.\n",
        )

    def test_print_streams(self):
        obs_stdout = io.StringIO()
        obs_stderr = io.StringIO()

        with redirect_stdout(obs_stdout), redirect_stderr(obs_stderr):
            _print_streams(self.result)

        self.assertEqual(obs_stdout.getvalue(), self.result.stdout)
        self.assertEqual(obs_stderr.getvalue(), self.result.stderr)

    def test_print_streams_stderr_empty(self):
        self.result.stderr = ""

        obs_stdout = io.StringIO()
        obs_stderr = io.StringIO()

        with redirect_stdout(obs_stdout), redirect_stderr(obs_stderr):

            _print_streams(self.result)

        self.assertEqual(obs_stdout.getvalue(), self.result.stdout)
        self.assertEqual(obs_stderr.getvalue(), "")

    def test_print_streams_stdout_empty(self):
        self.result.stdout = ""

        obs_stdout = io.StringIO()
        obs_stderr = io.StringIO()

        with redirect_stdout(obs_stdout), redirect_stderr(obs_stderr):

            _print_streams(self.result)

        self.assertEqual(obs_stdout.getvalue(), "")
        self.assertEqual(obs_stderr.getvalue(), self.result.stderr)
