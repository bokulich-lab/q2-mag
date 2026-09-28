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
import pandas as pd
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
    taxonomy: pd.DataFrame,
):
    bins_dp = os.path.join(loc, samp_name)
    bins_prefix = os.path.join(bins_dp, "bin")
    os.makedirs(bins_dp)

    if taxonomy != None:
        taxonomy_fp = os.path.join(loc, "taxonomy.tsv")
        taxonomy.to_csv(taxonomy_fp, sep='\t')

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
    ]
    cmd.extend(common_args)
    run_command(cmd, verbose=True)
    return bins_dp


def _process_sample(
    samp_name, samp_props, binner, multi_split, common_args, result_loc, taxonomy
):
    binsplit_separator = "C" if multi_split else ""

    with tempfile.TemporaryDirectory() as tmp:
        output_dp = _run_vamb(
            binner, samp_name, samp_props, tmp, common_args, binsplit_separator,
            taxonomy
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
) -> dict:
    fasta_fps = sorted(contigs.sample_dict().values())
    bam_fps = sorted(glob.glob(os.path.join(str(alignment_maps), "*.bam")))

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


def _assert_taxids(taxid_map: dict, taxonomy: pd.DataFrame) -> list[str]:
    """Asserts that all taxonomy IDs in the provided ``taxid_map`` are also found in
    ``taxonomy``. Returns a list of the taxonomy IDs found in both objects.

    Args:
        taxid_map (dict): Mapping between NCBI taxonomy IDs and contig IDs. 0 represents
            unclassified contigs.
        taxonomy (pd.DataFrame): Table mapping NCBI taxonomy IDs to full taxonomy
            strings.

    Raises:
        AssertionError: if 1+ taxonomy IDs are missing in ``taxonomy``.

    Returns:
        list[str]: List of taxonomy IDs found in both objects.
    """
    taxids = set(taxid_map.keys())
    feature_ids = set(taxonomy.index)
    missing_taxids = taxids - feature_ids

    if missing_taxids:
        raise AssertionError(
            "At least one taxonomy ID is missing in taxonomy. "
            f"The following taxonomy IDs are missing: {','.join(missing_taxids)}"
        )

    return list(taxids)


def _convert_taxonomy(taxid_map: dict, taxonomy: pd.DataFrame) -> pd.DataFrame:
    """Convert the provided taxonomy annotations into the format expected by TaxVAMB.

    Args:
        taxid_map (dict): Mapping between NCBI taxonomy IDs and contig IDs. 0 represents
            unclassified contigs.
        taxonomy (pd.DataFrame): Table mapping NCBI taxonomy IDs to full taxonomy
            strings.

    Returns:
        pd.DataFrame: Table mapping contigs to their predictions.

    Examples:
        The following is an example of a valid taxonomy file:

            contigs	predictions
            S18C13	Bacteria;Bacillota;Clostridia;Eubacteriales
            S18C25	Bacteria;Pseudomonadota
            S18C67	Bacteria;Bacillota;Bacilli;Bacillales;Staphylococcaceae
    """
    taxids = _assert_taxids(taxid_map=taxid_map, taxonomy=taxonomy)

    contigs = []
    predictions = []

    for taxid in taxids:
        contig_ids = taxid_map[taxid]
        contigs.extend(contig_ids)
        predictions.extend(
            [taxonomy.loc[taxid, "Taxon"]] * len(contig_ids)
        )

    return pd.DataFrame(
        data = zip(contigs, predictions),
        columns = ["contigs", "predictions"]
    )


def _bin_contigs_vamb(
    contigs: ContigSequencesDirFmt,
    alignment_maps: BAMDirFmt,
    taxonomy: dict | None,
    multi_split: bool,
    common_args: list,
) -> (MultiFASTADirectoryFormat, dict):
    binner = "taxvamb" if taxonomy is not None else "default"
    sample_set = _assert_samples(contigs, alignment_maps)
    _assert_reference_integrity(sample_set, "contigs", "map")

    bins = MultiFASTADirectoryFormat()
    for samp, props in sample_set.items():
        _process_sample(samp, props, binner, multi_split, common_args, str(bins),
                        taxonomy)

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


def bin_contigs_taxvamb(
    contigs: ContigSequencesDirFmt,
    alignment_maps: BAMDirFmt,
    taxid_map: dict,
    taxonomy: pd.DataFrame,
    # multi_split: bool = False,
    min_contig_len: int = 2000,
    minfasta: int = 2000,
    threads: int = 8,
    seed: int | None = None,
    no_predictor: bool = False,
) -> (MultiFASTADirectoryFormat, dict):
    multi_split = False  # Placeholder until multi split is supported

    kwargs = {
        k: v
        for k, v in locals().items()
        if k
        not in ["contigs", "alignment_maps", "taxid_map", "taxonomy", "multi_split"]
    }

    common_args = _process_common_input_params(
        processing_func=_process_vamb_arg, params=kwargs
    )

    converted_taxonomy = _convert_taxonomy(taxid_map=taxid_map, taxonomy=taxonomy)

    return _bin_contigs_vamb(
        contigs=contigs,
        alignment_maps=alignment_maps,
        taxonomy=converted_taxonomy,
        multi_split=multi_split,
        common_args=common_args,
    )
