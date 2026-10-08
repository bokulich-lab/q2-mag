# ----------------------------------------------------------------------------
# Copyright (c) 2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import importlib.resources

import qiime2


def _data_path(name):
    return str(importlib.resources.files("q2_mag.semibin2.tests") / "data" / name)


def contigs_factory():
    return qiime2.Artifact.import_data("SampleData[Contigs]", _data_path("contigs"))


def alignment_maps_factory():
    return qiime2.Artifact.import_data(
        'SampleData[AlignmentMap % Properties("sorted")]', _data_path("maps")
    )


def bin_contigs_semibin2_example(use):
    contigs = use.init_artifact("contigs", contigs_factory)
    alignment_maps = use.init_artifact("alignment_maps", alignment_maps_factory)

    mags, contig_map = use.action(
        use.UsageAction(plugin_id="mag", action_id="bin_contigs_semibin2"),
        use.UsageInputs(
            contigs=contigs,
            alignment_maps=alignment_maps,
            epochs=1,
            batch_size=256,
            engine="cpu",
            threads=1,
            no_recluster=True,
            random_seed=1,
            num_partitions=2,
        ),
        use.UsageOutputNames(mags="mags", contig_map="contig_map"),
    )

    mags.assert_output_type("SampleData[MAGs]")
    contig_map.assert_output_type("FeatureMap[MAGtoContigs]")
