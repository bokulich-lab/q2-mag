# ----------------------------------------------------------------------------
# Copyright (c) 2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import os
import subprocess
from typing import List

from q2_mag.utils import colorify, run_command
from q2_mag.busco.types import BuscoDatabaseDirFmt

# OrthoDB dataset versions that can be requested from BUSCO
# Default is set explicitly so that any change to the datasets
# being fetched is intentional
DATASETS_VERSIONS = ["odb10", "odb12", "odb12.2"]
DEFAULT_DATASETS_VERSION = "odb12.2"

VERSIONED_SELECTIONS = {
    "all": "all_{version}",
    "prokaryota": "all_prokaryota_{version}",
    "eukaryota": "all_eukaryota_{version}",
}


def _process_lineages(lineages):
    """Process and validate BUSCO lineage selections.

    This function validates and normalizes lineage selections for BUSCO
    database downloads. It handles special cases like 'all' and domain-level
    lineages (prokaryota, eukaryota, virus).

    Parameters
    ----------
    lineages : list
        List of lineage names to process. Can include 'all' or domain names.

    Returns
    -------
    list
        Processed list of lineage names.

    Raises
    ------
    ValueError
        If no lineages are provided (empty list or None).
    """
    if not lineages:
        raise ValueError("No lineages provided.")

    if "all" in lineages:
        lineages = ["all"]
        return lineages

    domain_lineages = ["prokaryota", "eukaryota", "virus"]
    domain_matches = [lin for lin in lineages if lin in domain_lineages]
    if domain_matches:
        print(
            colorify(
                f"Domain lineages were provided ({', '.join(domain_matches)}) - "
                "other lineages, if any, will be ignored."
            )
        )
        lineages = domain_matches

    return lineages


def _resolve_download_targets(lineages, datasets_version):
    """Pin unversioned BUSCO lineage selections to a dataset version.

    BUSCO resolves 'all', 'prokaryota' and 'eukaryota' using a dataset
    version that is hardcoded per BUSCO release. To avoid the fetched datasets
    silently changing with the installed BUSCO version, those selections are
    replaced with their explicitly versioned equivalents (e.g. 'all_odb12').
    The 'virus' selection and individually named lineages are not affected.

    Parameters
    ----------
    lineages : list
        List of processed lineage names (see `_process_lineages`).
    datasets_version : str
        OrthoDB dataset version to be used, e.g. 'odb12'.

    Returns
    -------
    list
        List of download targets to be passed to `busco --download`.

    Raises
    ------
    ValueError
        If the dataset version is not supported.
    """
    if datasets_version not in DATASETS_VERSIONS:
        raise ValueError(
            f"Unsupported datasets version: '{datasets_version}'. "
            f"Supported versions are: {', '.join(DATASETS_VERSIONS)}."
        )

    return [
        (
            VERSIONED_SELECTIONS[lineage].format(version=datasets_version)
            if lineage in VERSIONED_SELECTIONS
            else lineage
        )
        for lineage in lineages
    ]


def delete_symlinks(root_dir):
    """Delete all symbolic links in a directory tree.

    This function walks through a directory tree and removes any symbolic
    links found in the file system. It does not follow symbolic links while
    traversing the directory structure.

    Parameters
    ----------
    root_dir : str
        Path to the root directory to start searching for symbolic links.
    """
    for dirpath, dirnames, filenames in os.walk(root_dir, followlinks=False):
        for _file in filenames:
            if os.path.islink(os.path.join(dirpath, _file)):
                os.unlink(os.path.join(dirpath, _file))


def fetch_busco_db(
    lineages: List[str] = None,
    datasets_version: str = DEFAULT_DATASETS_VERSION,
) -> BuscoDatabaseDirFmt:
    busco_db = BuscoDatabaseDirFmt(path=None, mode="w")

    lineages = _process_lineages(lineages)
    targets = _resolve_download_targets(lineages, datasets_version)

    print(colorify(f"Fetching lineages: {', '.join(targets)}."))
    cmd = ["busco", "--download_path", str(busco_db), "--download", *targets]
    try:
        run_command(cmd)
    except subprocess.CalledProcessError as e:
        raise Exception(f"Error during BUSCO database download: {e.returncode}")

    # There is a symlink in the BUSCO database that needs to be removed
    delete_symlinks(str(busco_db))

    print(
        colorify(
            "Download completed. \n"
            "Copying files from temporary directory to the final location..."
        )
    )

    return busco_db
