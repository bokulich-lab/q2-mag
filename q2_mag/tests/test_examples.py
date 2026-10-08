# ----------------------------------------------------------------------------
# Copyright (c) 2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from qiime2.plugin.testing import TestPluginBase


class TestUsageExamples(TestPluginBase):
    package = "q2_mag.tests"

    def test_examples(self):
        self.execute_examples()
