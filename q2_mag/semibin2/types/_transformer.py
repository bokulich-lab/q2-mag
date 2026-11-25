# ----------------------------------------------------------------------------
# Copyright (c) 2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from q2_types.feature_data_mag import MAGSequencesDirFmt
from q2_types.per_sample_sequences import MultiFASTADirectoryFormat
from qiime2.util import duplicate

from ...plugin_setup import plugin


@plugin.register_transformer
def _1(source: MultiFASTADirectoryFormat) -> MAGSequencesDirFmt:
    result = MAGSequencesDirFmt()
    for fasta in source.path.glob("*/*.fa*"):
        duplicate(str(fasta), str(result.path / fasta.name))
    return result
