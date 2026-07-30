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
from pathlib import Path
from uuid import uuid4

from q2_types.per_sample_sequences import (
    BAMDirFmt,
    ContigSequencesDirFmt,
    MultiFASTADirectoryFormat,
)
from q2_mag.metabat2.metabat2 import _generate_contig_map
from q2_mag.utils import _process_common_input_params
from q2_mag.vamb.utils import _process_vamb_arg


def _run_vamb(
    binner: str,
    samp_name: str,
):
    cmd = [
        "vamb",
        "bin",
        binner,
        "--fasta",
        samp_props["contigs"],
        "--bamdir",
        samp_props["map"],
        "--output",
        bins_prefix,
        ""
    ]

    return


def _process_sample(samp_name, samp_props, binner, common_args, result_loc):
    with tempfile.TemporaryDirectory() as tmp:
        output_dp = _run_vamb(binner, samp_name, samp_props, tmp, binner, common_args)
        bins_dp = os.path.join(output_dp, "bin/output_bins/")

        all_outputs = glob.glob(os.path.join(bins_dp, "*.fa"))
        all_bins = [x for x in all_outputs if re.search(r"SemiBin\_[0-9]+\.fa$", x)]

        # rename using UUID v4
        bin_dest_dir = os.path.join(str(result_loc), samp_name)
        os.makedirs(bin_dest_dir, exist_ok=True)
        for old_bin in all_bins:
            new_bin = os.path.join(bin_dest_dir, f"{uuid4()}.fa")
            shutil.move(old_bin, new_bin)

    return


def _assert_samples(
    fasta: ContigSequencesDirFmt,
    bamdir: BAMDirFmt,
    taxonomy: dict | None,
) -> dict:
    fasta_fps = fasta.sample_dict.values()
    bam_fps = glob.glob(os.path.join(str(bamdir), "*.bam"))
    fasta_fps, bam_fps = sorted(fasta_fps), sorted(bam_fps)

    fasta_samples = fasta.sample_dict.keys()
    bam_samples = [Path(fp).stem.rsplit("_alignment", 1)[0] for fp in bam_fps]

    if set(fasta_samples) != set(bam_samples):
        raise Exception(
            "Contigs and alignment maps should belong to the same sample set. "
            f'You provided contigs for samples: {",".join(fasta_samples)} '
            f'but maps for samples: {",".join(bam_samples)}. Please check '
            "your inputs and try again."
        )

    return {
        s: {"contigs": fasta_fps[i], "map": bam_fps[i]}
        for i, s in enumerate(fasta_fps)
    }


def _bin_contigs_vamb(
    fasta: ContigSequencesDirFmt,
    bamdir: BAMDirFmt,
    taxonomy: dict | None,
    common_args: list,
) -> (MultiFASTADirectoryFormat, dict):
    binner = "taxvamb" if taxonomy is not None else "default"
    sample_set = _assert_samples(fasta, bamdir)

    bins = MultiFASTADirectoryFormat()
    for samp, props in sample_set.items():
        _process_sample(samp, props)


def bin_contigs_vamb(
    fasta: ContigSequencesDirFmt,
    bamdir: BAMDirFmt,
    min_contig_len: int = 2000,
    minfasta: int = 2000,
    threads: int = 8,
    seed: int | None = None
) -> (MultiFASTADirectoryFormat, dict):
    kwargs = {
        k: v
        for k, v in locals().items()
        if k not in ["fasta", "bamdir"]
    }

    common_args = _process_common_input_params(
        processing_func=_process_vamb_arg, params=kwargs
    )

    return _bin_contigs_vamb(
        fasta=fasta,
        bamdir=bamdir,
        common_args=common_args
    )


def bin_contigs_taxvamb(
    fasta: ContigSequencesDirFmt,
    bamdir: BAMDirFmt,
    taxonomy: dict,
    min_contig_len: int = 2000,
    min_bin_size: int = 2000,
    no_predictor: bool = False,
    threads: int = 8,
    seed: int | None = None
) -> (MultiFASTADirectoryFormat, dict):
    kwargs = {
        k: v
        for k, v in locals().items()
        if k not in ["fasta", "bamdir"]
    }

    common_args = _process_common_input_params(
        processing_func=_process_vamb_arg, params=kwargs
    )

    return _bin_contigs_vamb(
        fasta=fasta,
        bamdir=bamdir,
        taxonomy=taxonomy
        common_args=common_args
    )
