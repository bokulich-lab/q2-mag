# ----------------------------------------------------------------------------
# Copyright (c) 2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import glob
import os.path
import re
import shutil
import tempfile
from uuid import uuid4

from q2_types.per_sample_sequences import (
    BAMDirFmt,
    ContigSequencesDirFmt,
    MultiFASTADirectoryFormat,
)

from q2_mag.utils import _process_common_input_params, run_command
from q2_mag.metabat2.metabat2 import _assert_samples, _generate_contig_map
from q2_mag.semibin2.utils import _process_semibin2_arg


def _filter_alignment_maps(
    contigs: ContigSequencesDirFmt, alignment_maps: BAMDirFmt
) -> BAMDirFmt:
    samples = contigs.sample_dict()
    maps_by_sample = alignment_maps.file_dict()

    filtered_maps = BAMDirFmt()
    for sample in samples:
        shutil.copy(maps_by_sample[f"{sample}_alignment"], filtered_maps.path)

    return filtered_maps


def _run_semibin2(samp_name, samp_props, loc, mode, common_args):
    bins_dp = os.path.join(loc, samp_name)
    bins_prefix = os.path.join(bins_dp, "bin")
    os.makedirs(bins_dp)
    cmd = [
        "SemiBin2",
        mode,
        "--input-fasta",
        samp_props["contigs"],
        "--input-bam",
        samp_props["map"],
        "--output",
        bins_prefix,
        "--compression",
        "none",
        "--verbose",
    ]
    cmd.extend(common_args)
    run_command(cmd, verbose=True)
    return bins_dp


def _process_sample(samp_name, samp_props, mode, common_args, result_loc):
    with tempfile.TemporaryDirectory() as tmp:
        output_dp = _run_semibin2(samp_name, samp_props, tmp, mode, common_args)
        bins_dp = os.path.join(output_dp, "bin/output_bins/")

        all_outputs = glob.glob(os.path.join(bins_dp, "*.fa"))
        all_bins = [x for x in all_outputs if re.search(r"SemiBin\_[0-9]+\.fa$", x)]

        # rename using UUID v4
        bin_dest_dir = os.path.join(str(result_loc), samp_name)
        os.makedirs(bin_dest_dir, exist_ok=True)
        for old_bin in all_bins:
            new_bin = os.path.join(bin_dest_dir, f"{uuid4()}.fa")
            shutil.move(old_bin, new_bin)


def _bin_partition_contigs_semibin2(
    contigs: ContigSequencesDirFmt,
    alignment_maps: BAMDirFmt,
    mode: str,
    common_args: list,
) -> (MultiFASTADirectoryFormat, dict):
    filtered_alignment_maps = _filter_alignment_maps(contigs, alignment_maps)
    sample_set = _assert_samples(contigs, filtered_alignment_maps)

    bins = MultiFASTADirectoryFormat()
    for samp, props in sample_set.items():
        _process_sample(samp, props, mode, common_args, str(bins))

    if not glob.glob(os.path.join(str(bins), "*/*.fa")):
        raise ValueError(
            "No MAGs were formed during binning, please check your inputs."
        )

    contig_map = _generate_contig_map(bins)

    return bins, contig_map


def _bin_contigs_semibin2(
    contigs: ContigSequencesDirFmt,
    alignment_maps: BAMDirFmt,
    # mode: str,
    training_type: str | None = None,
    orf_finder: str = "fast-naive",
    environment: str = "global",
    engine: str = "auto",
    sequencing_type: str = "short_read",
    minfasta_kbs: int = 200,
    no_recluster: bool = False,
    epochs: int = 15,
    batch_size: int = 2048,
    max_node: int = 1,
    max_edges: int = 200,
    ratio: float = 0.05,
    threads: int | None = 1,
    min_len: int | None = None,
    ml_threshold: int | None = None,
    random_seed: int | None = None,
    debug: bool = False,
) -> (MultiFASTADirectoryFormat, dict):
    kwargs = {
        k: v
        for k, v in locals().items()
        if k not in ["contigs", "alignment_maps", "mode", "training_type"]
    }

    # mode = "single_easy_bin" if mode == "single" else "multi_easy_bin"
    mode = "single_easy_bin"

    common_args = _process_common_input_params(
        processing_func=_process_semibin2_arg, params=kwargs
    )

    return _bin_partition_contigs_semibin2(
        contigs=contigs,
        alignment_maps=alignment_maps,
        mode=mode,
        common_args=common_args,
    )


def bin_contigs_semibin2(
    ctx,
    contigs,
    alignment_maps,
    training_type=None,
    orf_finder="fast-naive",
    environment="global",
    engine="auto",
    sequencing_type="short_read",
    minfasta_kbs=200,
    no_recluster=False,
    epochs=15,
    batch_size=2048,
    max_node=1,
    max_edges=200,
    ratio=0.05,
    threads=1,
    min_len=None,
    ml_threshold=None,
    random_seed=None,
    debug=False,
    num_partitions=None,
):
    kwargs = {
        key: value
        for key, value in locals().items()
        if key not in {"ctx", "contigs", "alignment_maps", "num_partitions"}
    }

    partition_contigs = ctx.get_action("types", "partition_contigs")
    bin_partition = ctx.get_action("mag", "_bin_contigs_semibin2")
    collate_mags = ctx.get_action("types", "collate_sample_data_mags")
    collate_contig_maps = ctx.get_action("mag", "collate_contig_maps")

    (partitioned_contigs,) = partition_contigs(contigs, num_partitions)
    mags = []
    contig_maps = []
    for contig_partition in partitioned_contigs.values():
        partition_mags, partition_contig_map = bin_partition(
            contig_partition, alignment_maps, **kwargs
        )
        mags.append(partition_mags)
        contig_maps.append(partition_contig_map)

    (collated_mags,) = collate_mags(mags)
    (collated_contig_map,) = collate_contig_maps(contig_maps)

    return collated_mags, collated_contig_map
