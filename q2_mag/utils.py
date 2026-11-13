# ----------------------------------------------------------------------------
# Copyright (c) 2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import os
import subprocess
from typing import List, Union

import pandas as pd
import skbio
from q2_assembly._action_params import ALLOWED_SEPARATORS
from q2_assembly._utils import _split_contig_id
from q2_types.feature_data import DNASequencesDirectoryFormat

from q2_types.feature_data_mag import MAGSequencesDirFmt
from q2_types.per_sample_sequences import (
    BAMDirFmt,
    ContigSequencesDirFmt,
    MultiMAGSequencesDirFmt,
)

EXTERNAL_CMD_WARNING = (
    "Running external command line application(s). "
    "This may print messages to stdout and/or stderr.\n"
    "The command(s) being run are below. These commands "
    "cannot be manually re-run as they will depend on "
    "temporary files that no longer exist."
)


def run_command(cmd, env=None, verbose=True, pipe=False, **kwargs):
    if verbose:
        print(EXTERNAL_CMD_WARNING)
        print("\nCommand:", end=" ")
        print(" ".join(cmd), end="\n\n")

    if pipe:
        result = subprocess.run(
            cmd, env=env, check=True, capture_output=True, text=True
        )
        return result

    if env:
        subprocess.run(cmd, env=env, check=True, **kwargs)
    else:
        subprocess.run(cmd, check=True, **kwargs)


def _construct_param(arg_name):
    """Converts argument name into a command line parameter."""
    return f'--{arg_name.replace("_", "-")}'


def _process_common_input_params(processing_func, params: dict) -> List[str]:
    """Converts provided arguments and their values.

    Conversion is entirely dependent on the passed 'processing_func'
    that processes individual arguments. The output is a list of
    parameters with their values that can be directly passed to the
    respective command.

    Arguments without any value are skipped.
    Any other argument is processed using the 'processing_func' and
    appended to the final list.

    Args:
        processing_func: Function to be used for formatting a single argument.
        params (dict): Dictionary of parameter: value pairs to be processed.

    Returns:
        processed_args (list): List of processed arguments and their values.

    """
    processed_args = []
    for arg_key, arg_val in params.items():
        # This if condition excludes arguments which are falsy
        # (False, None, "", []), except for integers and floats.
        if type(arg_val) == int or type(arg_val) == float or arg_val:  # noqa: E721
            processed_args.extend(processing_func(arg_key, arg_val))

    return processed_args


def colorify(string: str):
    return "%s%s%s" % ("\033[1;32m", string, "\033[0m")


def _get_separator(
    header_id: str,
    sample_id: str,
    separator: str | None = None,
) -> str:
    """Extract the separator separating the sample ID from the rest of the header. The
    separator is a single character following the sample ID:

        <sample_id><separator><contig_id>

    Args:
        header_id (str): Full header ID.
        sample_id (str): Sample ID. Should be a subset of the full header ID.
        separator (str, optional): Expected separator character (if any).
            Defaults to "".

    Raises:
        AssertionError: If no separator follows the sample ID, the extracted
            character is not a valid separator, or it differs from the expected
            (provided) separator.

    Returns:
        str: The separator character following the sample ID.
    """
    if not separator:
        if len(header_id) < len(sample_id):
            raise AssertionError(
                "Found sample ID but no separator found in contig header. "
            )

        separator = header_id[len(sample_id)]

        if not separator in ALLOWED_SEPARATORS:
            raise AssertionError(
                f"{separator} is not a valid separator. Must be one of: "
                f"Valid separators are: {", ".join(ALLOWED_SEPARATORS)}"
            )
    elif header_id[len(sample_id)] != separator:
        raise AssertionError(
            "Found sample ID but no separator found in contig header. "
            f"Valid separators are: {", ".join(ALLOWED_SEPARATORS)}"
        )

    return separator


def _get_sample_ids(
    contigs: ContigSequencesDirFmt, expected_samples: list[str]
) -> tuple[set[str], str]:
    """Check that the expected sample IDs exists within the provided contig headers."""
    fasta_fps = contigs.sample_dict().values()
    headers = []
    n_samples = len(expected_samples)
    samples_found = set()
    separator = ""

    for fasta_fp in fasta_fps:
        for contig in skbio.io.read(fasta_fp, format="fasta"):
            header_id = contig.metadata["id"]
            headers.append(header_id)

            # We expected Sample IDs to be ordered. We therefore, use the last matching
            # sample ID instead of iterating over a fixed list.
            i = 0
            j = 0

            while i != n_samples:
                current_sample = expected_samples[j]

                if header_id.startswith(current_sample):
                    separator = _get_separator(header_id, current_sample, separator)
                    if current_sample not in samples_found:
                        samples_found.add(current_sample)
                else:
                    j += 1 if j <= n_samples else 0
                i += 1

    return (samples_found, separator)


def _align_samples_before_binning(
    contigs: ContigSequencesDirFmt,
    alignment_maps: BAMDirFmt,
    binning_workflow: str = "single-sample",  # Default provides backwards compatability
) -> dict[str, dict[str, str]]:
    """Validate that sample IDs match and pair each sample's FASTA and BAM paths."""
    if binning_workflow not in ("single-sample", "multi-sample", "co-assembly"):
        raise ValueError(
            "Binning workflow must be one of single-sample, multi-sample, or "
            f"co-assembly. You provided: {binning_workflow}."
        )

    contigs_by_sample = contigs.sample_dict()
    maps_by_sample = alignment_maps.file_dict()
    bam_samples = maps_by_sample.keys()

    # TODO: single-sample binning doesn't expect a separator
    if binning_workflow in ("multi-sample", "co-assembly"):
        fasta_samples, separator = _get_sample_ids(contigs, list(bam_samples))

        return (
            {
                "pooled": {
                    "contigs": next(iter(contigs_by_sample.values())),
                    "map": list(maps_by_sample.values()),
                }
            },
            separator,
        )
    else:
        fasta_samples = contigs_by_sample.keys()

        if fasta_samples != bam_samples:
            raise Exception(
                "Contigs and alignment maps should belong to the same sample set. "
                f'You provided contigs for samples: {", ".join(fasta_samples)} '
                f'but maps for samples: {", ".join(bam_samples)}. Please check '
                "your inputs and try again."
            )

        return (
            {
                sample: {"contigs": contig_path, "map": maps_by_sample[sample]}
                for sample, contig_path in contigs_by_sample.items()
            },
            separator,
        )


def get_feature_lengths(
    features: Union[
        MAGSequencesDirFmt,
        MultiMAGSequencesDirFmt,
        ContigSequencesDirFmt,
        DNASequencesDirectoryFormat,
    ],
) -> pd.DataFrame:
    """Calculate lengths of features in a feature data object."""
    ids, lengths = [], []
    per_sequence = False

    # for FeatureData[Sequence]
    if isinstance(features, DNASequencesDirectoryFormat):
        sample_dict = {"": {"": os.path.join(str(features), "dna-sequences.fasta")}}
        per_sequence = True

        # For SampleData[Contigs]
    elif isinstance(features, ContigSequencesDirFmt):
        sample_dict = {"": features.sample_dict()}
        per_sequence = True

    # For SampleData[MAGs]
    elif isinstance(features, MultiMAGSequencesDirFmt):
        sample_dict = features.sample_dict()

    # For FeatureData[MAG]
    else:
        sample_dict = {"": features.feature_dict()}

    for _, file_dict in sample_dict.items():
        for _id, fp in file_dict.items():
            sequences = skbio.io.read(fp, format="fasta", verify=False)

            if per_sequence:
                for seq in sequences:
                    ids.append(seq.metadata["id"])
                    lengths.append(len(seq))
            else:
                ids.append(_id)
                lengths.append(sum(len(seq) for seq in sequences))

    df = pd.DataFrame({"id": ids, "length": lengths})
    df.set_index("id", inplace=True)
    return df
