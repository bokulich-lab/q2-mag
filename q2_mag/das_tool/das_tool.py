# ----------------------------------------------------------------------------
# Copyright (c) 2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import glob
import os
import shutil
import tempfile
from pathlib import Path
from uuid import uuid4

import skbio.io
from q2_types.genome_data import ProteinsDirectoryFormat
from q2_types.per_sample_sequences import (
    ContigSequencesDirFmt,
    MultiFASTADirectoryFormat,
)

from q2_annotate._utils import _process_common_input_params, run_command


def _process_das_tool_arg(arg_key, arg_val):
    if isinstance(arg_val, bool) and arg_val:
        return [f"--{arg_key}"]
    else:
        return [f"--{arg_key}", str(arg_val)]


def _get_sample_ids(*formats):
    sample_ids = None
    for fmt in formats:
        ids = set(fmt.sample_dict())
        if sample_ids is None:
            sample_ids = ids
        elif sample_ids != ids:
            raise ValueError(
                "Bins and contigs must contain the same sample IDs. "
                f"Observed sample sets: {sorted(sample_ids)} and {sorted(ids)}."
            )
    return sorted(sample_ids)


def _get_sample_proteins(proteins: ProteinsDirectoryFormat | None) -> dict:
    if proteins is None:
        return {}

    sample_proteins = {}
    fasta_paths = []
    for ext in ("*.fa", "*.faa", "*.fasta"):
        fasta_paths.extend(glob.glob(os.path.join(str(proteins), ext)))

    for fp in fasta_paths:
        sample_proteins[Path(fp).stem] = fp

    return sample_proteins


def _write_contig2bin_map(bins, sample_id, label, output_dir):
    sample_bins = bins.sample_dict()[sample_id]
    output_fp = os.path.join(output_dir, f"{label}_contig2bin.tsv")

    with open(output_fp, "w") as fh:
        for bin_id, bin_fp in sorted(sample_bins.items()):
            for seq in skbio.io.read(bin_fp, format="fasta", verify=False):
                fh.write(f"{seq.metadata['id']}\t{bin_id}\n")

    return output_fp


def _run_das_tool(
    sample_id,
    bins,
    contigs_fp,
    proteins_fp,
    labels,
    common_args,
    output_dir,
):
    contig2bin_fps = []
    for idx, binning in enumerate(bins, start=1):
        label = labels[idx - 1]
        contig2bin_fps.append(
            _write_contig2bin_map(binning, sample_id, label, output_dir)
        )

    output_prefix = os.path.join(output_dir, sample_id)
    cmd = [
        "DAS_Tool",
        "--bins",
        ",".join(contig2bin_fps),
        "--labels",
        ",".join(labels),
        "--contigs",
        contigs_fp,
        "--outputbasename",
        output_prefix,
        "--write_bins",
    ]

    if proteins_fp is not None:
        cmd.extend(["--proteins", proteins_fp])

    cmd.extend(common_args)
    run_command(cmd, verbose=True)

    return f"{output_prefix}_DASTool_bins"


def _collect_refined_bins(sample_id, das_tool_bins_dir, refined_bins):
    sample_output_dir = os.path.join(str(refined_bins), sample_id)
    os.makedirs(sample_output_dir, exist_ok=True)

    refined_bin_fps = []
    for ext in ("*.fa", "*.fasta", "*.fna"):
        refined_bin_fps.extend(glob.glob(os.path.join(das_tool_bins_dir, ext)))

    for src in sorted(refined_bin_fps):
        shutil.copy(src, os.path.join(sample_output_dir, f"{uuid4()}.fa"))

    return len(refined_bin_fps)


def _parse_labels(labels: str, n_bins: int) -> list[str]:
    split_labels = labels.split(",")
    if len(split_labels) != n_bins:
        raise ValueError(
            "The number of labels provided is different from the number of bins."
        )

    if len(split_labels) != len(set(split_labels)):
        duplicate_labels = set(
            [label for label in split_labels if split_labels.count(label) > 1]
        )
        raise ValueError(
            "Duplicate labels detected. Each label provided must be a unique string. "
            f"The following labels appear more than once: {",".join(duplicate_labels)}."
        )

    return split_labels


def _generate_labels(n_bins: int) -> list[str]:
    labels = []
    for idx in range(1, n_bins + 1):
        label = f"binning_{idx}"
        labels.append(label)
    return labels


def _refine_bins_das_tool(
    bins: MultiFASTADirectoryFormat,
    contigs: ContigSequencesDirFmt,
    proteins: ProteinsDirectoryFormat | None,
    labels: str | None,
    common_args: list,
) -> MultiFASTADirectoryFormat:
    if len(bins) < 2:
        raise ValueError("DAS Tool requires bins from at least two binning methods.")

    if labels is not None:
        labels = _parse_labels(labels, len(bins))
    else:
        labels = _generate_labels(len(bins))

    sample_ids = _get_sample_ids(contigs, *bins)
    sample_proteins = _get_sample_proteins(proteins)
    refined_bins = MultiFASTADirectoryFormat()
    num_refined_bins = 0

    with tempfile.TemporaryDirectory() as tmp:
        for sample_id in sample_ids:
            das_tool_bins_dir = _run_das_tool(
                sample_id=sample_id,
                bins=bins,
                contigs_fp=contigs.sample_dict()[sample_id],
                proteins_fp=sample_proteins.get(sample_id),
                labels=labels,
                common_args=common_args,
                output_dir=tmp,
            )
            num_refined_bins += _collect_refined_bins(
                sample_id, das_tool_bins_dir, refined_bins
            )

    if num_refined_bins == 0:
        raise ValueError(
            "No refined MAGs were formed by DAS Tool, please check your inputs."
        )

    return refined_bins


def refine_bins_das_tool(
    bins: MultiFASTADirectoryFormat,
    contigs: ContigSequencesDirFmt,
    proteins: ProteinsDirectoryFormat | None = None,
    labels: str | None = None,
    search_engine: str = "diamond",
    score_threshold: float = 0.5,
    duplicate_penalty: float = 0.6,
    megabin_penalty: float = 0.5,
    max_iter_post_threshold: int = 10,
    threads: int = 1,
    debug: bool | None = None,
) -> MultiFASTADirectoryFormat:
    kwargs = {
        k: v for k, v in locals().items() if k not in ["bins", "contigs", "proteins", "labels"]
    }

    common_args = _process_common_input_params(
        processing_func=_process_das_tool_arg, params=kwargs
    )

    return _refine_bins_das_tool(
        bins=bins,
        contigs=contigs,
        proteins=proteins,
        labels=labels,
        common_args=common_args,
    )
