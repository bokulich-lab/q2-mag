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
import subprocess
import tempfile
import warnings
from uuid import uuid4

import pandas as pd
import rachis
import skbio.io
from q2_mag.utils import _process_common_input_params, run_command
from q2_types.genome_data import ProteinsDirectoryFormat
from q2_types.per_sample_sequences import (
    ContigSequencesDirFmt,
    MultiFASTADirectoryFormat,
)
from q2_mag.das_tool.utils import _print_streams

SUMMARY_DTYPES = {
    "bin": str,
    "bin_set": str,
    "unique_SCGs": int,
    "redundant_SCGs": int,
    "SCG_set": str,
    "size": int,
    "contigs": int,
    "N50": int,
    "bin_score": float,
    "SCG_completeness": float,
    "SCG_redundancy": float,
}


def _process_das_tool_arg(arg_key, arg_val):
    if isinstance(arg_val, bool) and arg_val:
        return [f"--{arg_key}"]
    else:
        return [f"--{arg_key}", str(arg_val)]


def _get_sample_ids(
    labels: list[str],
    contigs: ContigSequencesDirFmt,
    proteins: ProteinsDirectoryFormat | None = None,
    *bins: MultiFASTADirectoryFormat,
) -> list[str]:
    # 1. Ensure that sample IDs are consistent across all bins.
    bin_sample_sets = [set(x.sample_dict()) for x in bins]
    bin_samples = set().union(*bin_sample_sets)

    inconsistent_bin_samples = {
        idx: missing
        for idx, bin_set in enumerate(bin_sample_sets)
        if (missing := bin_samples - bin_set)
    }

    if len(inconsistent_bin_samples) >= 1:
        missing_rows = []

        for idx, missing in inconsistent_bin_samples.items():
            missing_sample_ids = ", ".join(sorted(missing))
            missing_rows.append(f"- {labels[idx]}: {missing_sample_ids}")

        raise ValueError(
            "Sample IDs must stay consistent across all bins.\n"
            "Missing sample IDs by bin:\n" + "\n".join(missing_rows)
        )

    # 2. Assert that the contig sample IDs contain all sample IDs in bins.
    contig_samples = set(contigs.sample_dict())
    missing = bin_samples - contig_samples

    if len(missing) >= 1:
        raise ValueError(
            "Contigs must contain all sample IDs present in bins. "
            f"Missing from contigs: {', '.join(sorted(missing)) or 'none'}."
        )

    # 3. If proteins are provided, assert that their sample IDs match those of bins.
    if proteins is not None:
        protein_samples = set(proteins.file_dict().keys())
        missing = bin_samples - protein_samples

        if len(missing) >= 1:
            raise ValueError(
                "Proteins must contain all sample IDs present in bins. "
                f"Missing from proteins: {', '.join(sorted(missing)) or 'none'}."
            )

    return sorted(bin_samples)


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

    # Copy inputs into the temp workspace so DAS Tool cannot create extra files
    # inside QIIME's cached artifact directories.
    input_dir = os.path.join(output_dir, f"{sample_id}_inputs")
    os.makedirs(input_dir, exist_ok=True)
    staged_contigs_fp = shutil.copy(contigs_fp, input_dir)
    if proteins_fp is not None:
        proteins_fp = shutil.copy(proteins_fp, input_dir)

    # Avoid DAS_Tool locale warnings.
    env = os.environ.copy()
    env["LANG"] = "C"
    env["LC_ALL"] = "C"

    output_prefix = os.path.join(output_dir, sample_id)
    cmd = [
        "DAS_Tool",
        "--bins",
        ",".join(contig2bin_fps),
        "--labels",
        ",".join(labels),
        "--contigs",
        staged_contigs_fp,
        "--outputbasename",
        output_prefix,
        "--write_bins",
        "--write_bin_evals",
    ]

    if proteins_fp is not None:
        cmd.extend(["--proteins", proteins_fp])

    cmd.extend(common_args)

    try:
        _print_streams(run_command(cmd, verbose=True, env=env, pipe=True))
    except subprocess.CalledProcessError as error:
        _print_streams(error)
        raise

    return (
        f"{output_prefix}_DASTool_bins",
        f"{output_prefix}_DASTool_summary.tsv",
        f"{output_prefix}_allBins.eval",
    )


def _collect_refined_bins(sample_id, das_tool_bins_dir, refined_bins):
    sample_output_dir = os.path.join(str(refined_bins), sample_id)
    os.makedirs(sample_output_dir, exist_ok=True)

    refined_bin_fps = []
    for ext in ("*.fa", "*.fasta", "*.fna"):
        refined_bin_fps.extend(glob.glob(os.path.join(das_tool_bins_dir, ext)))

    for src in sorted(refined_bin_fps):
        shutil.copy(src, os.path.join(sample_output_dir, f"{uuid4()}.fa"))

    return len(refined_bin_fps)


def _append_summary(sample_id, summary, summaries=None):
    summary_df = pd.read_csv(summary, sep="\t", dtype=SUMMARY_DTYPES)

    summary_df.insert(0, "sample_id", str(sample_id))
    summary_df.set_index("bin", drop=True, inplace=True)
    summary_df.index.name = "id"

    if summaries is None:
        return summary_df

    return pd.concat([summaries, summary_df], sort=False)


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
            f"The following labels appear more than once: {','.join(duplicate_labels)}."
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
    proteins: ProteinsDirectoryFormat,
    labels: str | None,
    common_args: list,
) -> MultiFASTADirectoryFormat:
    if len(bins) < 2:
        raise ValueError("DAS Tool requires bins from at least two binning methods.")

    if labels is not None:
        labels = _parse_labels(labels, len(bins))
    else:
        labels = _generate_labels(len(bins))

    protein_records = proteins.file_dict() if proteins is not None else {}
    sample_ids = _get_sample_ids(labels, contigs, proteins, *bins)
    concatenated_summary = None
    concatenated_evaluation = None
    refined_bins = MultiFASTADirectoryFormat()
    num_refined_bins = 0
    failed_samples = []

    with tempfile.TemporaryDirectory() as tmp:
        for sample_id in sample_ids:
            das_tool_bins_dir = os.path.join(tmp, f"{sample_id}_DASTool_bins")

            # In cases where 1+ samples fail to produce any bins, we issue a warning
            # but continue the analysis for the other samples, given that at least one
            # bin is produced.
            try:
                _, summary, evaluation = _run_das_tool(
                    sample_id=sample_id,
                    bins=bins,
                    contigs_fp=contigs.sample_dict()[sample_id],
                    proteins_fp=protein_records.get(sample_id),
                    labels=labels,
                    common_args=common_args,
                    output_dir=tmp,
                )
            except subprocess.CalledProcessError as e:
                no_bins = any(
                    line.startswith("Error:  No bins with bin-score >")
                    for line in (e.stdout or "").splitlines()
                )
                if no_bins:
                    warnings.warn(
                        f"No bins produced for sample {sample_id}.",
                        UserWarning,
                    )
                    failed_samples.append(sample_id)
                else:
                    raise Exception(
                        "An error was encountered while running DAS Tool, "
                        f"(return code {e.returncode}), please inspect "
                        "stdout and stderr to learn more."
                    ) from e
            else:
                # Do not append summary files if no refined bins were recovered
                concatenated_summary = _append_summary(
                    sample_id, summary, concatenated_summary
                )
                concatenated_evaluation = _append_summary(
                    sample_id, evaluation, concatenated_evaluation
                )

            # Always record the sample ID even if no refined bins were recovered for
            # that sample.
            num_refined_bins += _collect_refined_bins(
                sample_id, das_tool_bins_dir, refined_bins
            )

    # Report the names of samples that failed to produce any bins
    if failed_samples:
        warnings.warn(
            f"No bins produced for sample(s): {', '.join(failed_samples)}.",
            UserWarning,
        )

    if num_refined_bins == 0:
        raise ValueError(
            "No refined MAGs were formed by DAS Tool, please check your inputs. "
            "If DAS Tool is filtering out all bins, try lowering the score threshold."
        )

    return (
        refined_bins,
        rachis.Metadata(concatenated_summary),
        rachis.Metadata(concatenated_evaluation),
    )


def refine_bins_das_tool(
    bins: MultiFASTADirectoryFormat,
    contigs: ContigSequencesDirFmt,
    proteins: ProteinsDirectoryFormat = None,
    labels: str | None = None,
    search_engine: str = "diamond",
    score_threshold: float = 0.5,
    duplicate_penalty: float = 0.6,
    megabin_penalty: float = 0.5,
    max_iter_post_threshold: int = 10,
    threads: int = 1,
    debug: bool | None = None,
) -> (MultiFASTADirectoryFormat, rachis.Metadata, rachis.Metadata):
    kwargs = {
        k: v
        for k, v in locals().items()
        if k not in ["bins", "contigs", "proteins", "labels"]
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
