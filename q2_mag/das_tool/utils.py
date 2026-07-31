# ----------------------------------------------------------------------------
# Copyright (c) 2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import sys


def _print_streams(obj):
    if obj.stdout:
        print(obj.stdout, end="", file=sys.stdout)
    if obj.stderr:
        print(obj.stderr, end="", file=sys.stderr)
