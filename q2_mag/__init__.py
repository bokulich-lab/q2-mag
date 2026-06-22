# flake8: noqa
# ----------------------------------------------------------------------------
# Copyright (c) 2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from . import abundance, busco, dereplication, filtering, metabat2, utils, semibin2, das_tool

try:
    from ._version import __version__
except ModuleNotFoundError:
    __version__ = "0.0.0+notfound"

__all__ = [
    "abundance",
    "busco",
    "das_tool",
    "dereplication",
    "filtering",
    "metabat2",
    "semibin2",
    "utils",
]
