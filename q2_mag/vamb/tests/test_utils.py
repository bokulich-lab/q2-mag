# ----------------------------------------------------------------------------
# Copyright (c) 2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest
from rachis.plugin.testing import TestPluginBase

from q2_mag.vamb.utils import _process_vamb_arg


class TestVAMBUtils(TestPluginBase):
    package = "q2_mag.semibin2.tests"

    def test_process_vamb_arg_min_contig_len(self):
        obs = _process_vamb_arg("min_contig_len", 2000)
        exp = ["-m", "2000"]
        self.assertListEqual(obs, exp)

    def test_process_vamb_arg_threads(self):
        obs = _process_vamb_arg("threads", 8)
        exp = ["-p", "8"]
        self.assertListEqual(obs, exp)

    def test_process_vamb_arg_option(self):
        obs = _process_vamb_arg("minfasta", 1000)
        exp = ["--minfasta", "1000"]
        self.assertListEqual(obs, exp)

    def test_process_vamb_arg_bool(self):
        obs = _process_vamb_arg("verbose", True)
        exp = ["--verbose"]
        self.assertListEqual(obs, exp)


if __name__ == "__main__":
    unittest.main()
