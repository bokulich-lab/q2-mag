# ----------------------------------------------------------------------------
# Copyright (c) 2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import os
from q2_mag.busco.database import (
    fetch_busco_db,
    delete_symlinks,
    _resolve_download_targets,
)
from unittest.mock import patch
from qiime2.plugin.testing import TestPluginBase


class TestFetchBUSCO(TestPluginBase):
    package = "q2_mag.busco.tests"

    @patch("subprocess.run")
    def test_fetch_busco_db_virus(self, subp_run):
        busco_db = fetch_busco_db(lineages=["virus"])
        cmd = ["busco", "--download_path", str(busco_db), "--download", "virus"]
        subp_run.assert_called_once_with(cmd, check=True)

    @patch("subprocess.run")
    def test_fetch_busco_db_prok_euk(self, subp_run):
        busco_db = fetch_busco_db(lineages=["prokaryota", "eukaryota"])
        cmd = [
            "busco",
            "--download_path",
            str(busco_db),
            "--download",
            "all_prokaryota_odb12.2",
            "all_eukaryota_odb12.2",
        ]
        subp_run.assert_called_once_with(cmd, check=True)

    @patch("subprocess.run")
    def test_fetch_busco_db_all(self, subp_run):
        busco_db = fetch_busco_db(lineages=["all"])
        cmd = ["busco", "--download_path", str(busco_db), "--download", "all_odb12.2"]
        subp_run.assert_called_once_with(cmd, check=True)

    @patch("subprocess.run")
    def test_fetch_busco_db_two_lineages(self, subp_run):
        busco_db = fetch_busco_db(lineages=["lineage1", "lineage2"])
        cmd = [
            "busco",
            "--download_path",
            str(busco_db),
            "--download",
            "lineage1",
            "lineage2",
        ]
        subp_run.assert_called_once_with(cmd, check=True)

    @patch("subprocess.run")
    def test_fetch_busco_db_lineages_with_all(self, subp_run):
        busco_db = fetch_busco_db(lineages=["lineage1", "lineage2", "all"])
        cmd = ["busco", "--download_path", str(busco_db), "--download", "all_odb12.2"]
        subp_run.assert_called_once_with(cmd, check=True)

    @patch("subprocess.run")
    def test_fetch_busco_db_lineages_with_domains(self, subp_run):
        busco_db = fetch_busco_db(
            lineages=["lineage1", "lineage2", "prokaryota", "virus"]
        )
        cmd = [
            "busco",
            "--download_path",
            str(busco_db),
            "--download",
            "all_prokaryota_odb12.2",
            "virus",
        ]
        subp_run.assert_called_once_with(cmd, check=True)

    @patch("subprocess.run")
    def test_fetch_busco_db_no_lineage(self, subp_run):
        with self.assertRaisesRegex(ValueError, "No lineages provided."):
            fetch_busco_db()

    @patch("subprocess.run")
    def test_fetch_busco_db_all_datasets_version(self, subp_run):
        for version in ["odb10", "odb12", "odb12.2"]:
            with self.subTest(version=version):
                subp_run.reset_mock()
                busco_db = fetch_busco_db(lineages=["all"], datasets_version=version)
                cmd = [
                    "busco",
                    "--download_path",
                    str(busco_db),
                    "--download",
                    f"all_{version}",
                ]
                subp_run.assert_called_once_with(cmd, check=True)

    @patch("subprocess.run")
    def test_fetch_busco_db_domains_datasets_version(self, subp_run):
        busco_db = fetch_busco_db(
            lineages=["prokaryota", "eukaryota", "virus"], datasets_version="odb12.2"
        )
        cmd = [
            "busco",
            "--download_path",
            str(busco_db),
            "--download",
            "all_prokaryota_odb12.2",
            "all_eukaryota_odb12.2",
            "virus",
        ]
        subp_run.assert_called_once_with(cmd, check=True)

    @patch("subprocess.run")
    def test_fetch_busco_db_named_lineages_ignore_datasets_version(self, subp_run):
        busco_db = fetch_busco_db(
            lineages=["bacteria_odb10", "archaea_odb12"], datasets_version="odb12.2"
        )
        cmd = [
            "busco",
            "--download_path",
            str(busco_db),
            "--download",
            "bacteria_odb10",
            "archaea_odb12",
        ]
        subp_run.assert_called_once_with(cmd, check=True)

    @patch("subprocess.run")
    def test_fetch_busco_db_invalid_datasets_version(self, subp_run):
        with self.assertRaisesRegex(ValueError, "Unsupported datasets version: 'odb9'"):
            fetch_busco_db(lineages=["all"], datasets_version="odb9")
        subp_run.assert_not_called()

    def test_resolve_download_targets_default_selections_are_versioned(self):
        # None of the unversioned selections should ever be passed to BUSCO
        # as-is, since BUSCO would resolve them using its own default version.
        obs = _resolve_download_targets(["all", "prokaryota", "eukaryota"], "odb12")
        exp = ["all_odb12", "all_prokaryota_odb12", "all_eukaryota_odb12"]
        self.assertListEqual(obs, exp)

    def test_delete_symlinks(self):
        temp_dir = self.temp_dir.name
        sub_dir = os.path.join(temp_dir, "subdir")
        os.makedirs(sub_dir)

        # Create symlinks
        symlink1 = os.path.join(temp_dir, "symlink1")
        symlink2 = os.path.join(sub_dir, "symlink2")

        # Run delete function
        delete_symlinks(temp_dir)

        # Check that symlinks were deleted
        self.assertFalse(os.path.exists(symlink1))
        self.assertFalse(os.path.exists(symlink2))
