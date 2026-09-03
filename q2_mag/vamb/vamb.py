# ----------------------------------------------------------------------------
# Copyright (c) 2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import glob
import os.path
import shutil
import tempfile
import pysam
from pathlib import Path
from uuid import uuid4

from q2_types.per_sample_sequences import (
    BAMDirFmt,
    ContigSequencesDirFmt,
    MultiFASTADirectoryFormat,
)
from q2_mag.metabat2.metabat2 import _generate_contig_map
from q2_mag.utils import _process_common_input_params, run_command
from q2_mag.vamb.utils import _process_vamb_arg


def _run_vamb(
    binner: str,
    samp_name: str,
    samp_props: dict[str],
    loc: str,
    common_args: list[str],
    binsplit_separator: str,
):
    bins_dp = os.path.join(loc, samp_name)
    bins_prefix = os.path.join(bins_dp, "bin")
    os.makedirs(bins_dp)

    # VAMB expects a directory of BAM files and not a single file
    bamdir = os.path.join(loc, f"{samp_name}_bams")
    os.makedirs(bamdir)
    os.symlink(
        samp_props["map"],
        os.path.join(bamdir, os.path.basename(samp_props["map"])),
    )

    cmd = [
        "vamb",
        "bin",
        binner,
        "--fasta",
        samp_props["contigs"],
        "--bamdir",
        bamdir,
        "--outdir",
        bins_prefix,
        "-o",
        binsplit_separator,
        "--norefcheck",  # References across maps, contigs, and taxonomy already checked
    ]
    cmd.extend(common_args)
    run_command(cmd, verbose=True)
    return bins_dp


def _process_sample(
    samp_name, samp_props, binner, multi_split, common_args, result_loc
):
    binsplit_separator = "C" if multi_split else ""

    with tempfile.TemporaryDirectory() as tmp:
        output_dp = _run_vamb(
            binner, samp_name, samp_props, tmp, common_args, binsplit_separator
        )
        bins_dp = os.path.join(output_dp, "bin", "bins")

        all_bins = glob.glob(os.path.join(bins_dp, "*.fna"))

        # rename using UUID v4
        bin_dest_dir = os.path.join(str(result_loc), samp_name)
        os.makedirs(bin_dest_dir, exist_ok=True)
        for old_bin in all_bins:
            new_bin = os.path.join(bin_dest_dir, f"{uuid4()}.fa")
            shutil.move(old_bin, new_bin)

    return


def _assert_reference_integrity(
    sample_set: dict[str], contig_key: str, map_key: str
) -> None:
    """Verify that each sample's FASTA and BAM references are identical.

    Parameters
    ----------
    sample_set : dict[str]
        Mapping of sample identifiers to their associated file paths.
    contig_key : str
        Key for a sample's contig FASTA path.
    map_key : str
        Key for a sample's alignment-map path.

    Raises
    ------
    ValueError
        If a sample's references differ in count, name, order, or length.
    """
    failed_samples = {}

    for samp, props in sample_set.items():
        with (
            pysam.FastaFile(props[contig_key]) as fasta,
            pysam.AlignmentFile(props[map_key], "r") as bam,
        ):
            fasta_records = tuple(zip(fasta.references, fasta.lengths))
            bam_records = tuple(zip(bam.references, bam.lengths))

            if fasta_records != bam_records:
                failed_samples[samp] = (len(fasta_records), len(bam_records))

    if failed_samples:
        failed_sample_details = ", ".join(
            (
                "\n  "
                f"Sample {sample}: {fasta_count} contigs, {bam_count} BAM references"
            )
            for sample, (fasta_count, bam_count) in failed_samples.items()
        )

        raise ValueError(
            "Alignment maps do not match the corresponding contigs in at least one "
            "sample. The following samples had a mismatch in count, name, order, and "
            f"or length: {failed_sample_details}."
        )

    return None


def _assert_samples(
    contigs: ContigSequencesDirFmt,
    alignment_maps: BAMDirFmt,
    taxonomy: dict | None = None,
) -> dict:
    fasta_fps = contigs.sample_dict().values()
    bam_fps = glob.glob(os.path.join(str(alignment_maps), "*.bam"))
    fasta_fps, bam_fps = sorted(fasta_fps), sorted(bam_fps)

    fasta_samples = contigs.sample_dict().keys()
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
        for i, s in enumerate(fasta_samples)
    }


def _bin_contigs_vamb(
    contigs: ContigSequencesDirFmt,
    alignment_maps: BAMDirFmt,
    taxonomy: dict | None,
    multi_split: bool,
    common_args: list,
) -> (MultiFASTADirectoryFormat, dict):
    binner = "taxvamb" if taxonomy is not None else "default"
    sample_set = _assert_samples(contigs, alignment_maps, taxonomy)
    _assert_reference_integrity(sample_set, "contigs", "map")

    bins = MultiFASTADirectoryFormat()
    for samp, props in sample_set.items():
        _process_sample(samp, props, binner, multi_split, common_args, str(bins))

    if not glob.glob(os.path.join(str(bins), "*/*.fa")):
        raise ValueError(
            "No MAGs were formed during binning, please check your inputs."
        )

    contig_map = _generate_contig_map(bins)

    return bins, contig_map


def bin_contigs_vamb(
    contigs: ContigSequencesDirFmt,
    alignment_maps: BAMDirFmt,
    # multi_split: bool = False,
    min_contig_len: int = 2000,
    minfasta: int = 2000,
    threads: int = 8,
    seed: int | None = None,
) -> (MultiFASTADirectoryFormat, dict):
    multi_split = False  # Placeholder until multi split is supported

    kwargs = {
        k: v
        for k, v in locals().items()
        if k not in ["contigs", "alignment_maps", "multi_split"]
    }

    common_args = _process_common_input_params(
        processing_func=_process_vamb_arg, params=kwargs
    )

    return _bin_contigs_vamb(
        contigs=contigs,
        alignment_maps=alignment_maps,
        taxonomy=None,
        multi_split=multi_split,
        common_args=common_args,
    )


# def bin_contigs_taxvamb(
#     fasta: ContigSequencesDirFmt,
#     bamdir: BAMDirFmt,
#     taxonomy: dict,
#     multi_split: bool = False,
#     min_contig_len: int = 2000,
#     minfasta: int = 2000,
#     threads: int = 8,
#     seed: int | None = None,
#     no_predictor: bool = False,
# ) -> (MultiFASTADirectoryFormat, dict):
#     kwargs = {
#         k: v for k, v in locals().items() if k not in ["fasta", "bamdir",
# "multi_split"]
#     }

#     common_args = _process_common_input_params(
#         processing_func=_process_vamb_arg, params=kwargs
#     )

#     return _bin_contigs_vamb(
#         fasta=fasta,
#         bamdir=bamdir,
#         taxonomy=taxonomy,
#         multi_split=multi_split,
#         common_args=common_args,
#     )
