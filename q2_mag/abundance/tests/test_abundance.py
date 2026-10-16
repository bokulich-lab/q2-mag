# ----------------------------------------------------------------------------
# Copyright (c) 2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import shutil

import biom
import numpy as np
import pandas as pd
import pysam
import qiime2
from parameterized import parameterized
from qiime2.plugin.testing import TestPluginBase

from q2_mag.abundance.abundance import rpkm, tpm, estimate_abundance
from q2_types.per_sample_sequences import BAMDirFmt


class TestAbundance(TestPluginBase):
    package = "q2_mag.abundance.tests"

    def setUp(self):
        super().setUp()
        self.df = pd.DataFrame(
            {
                "sample-id": [*["s1"] * 4, *["s2"] * 4, *["s3"] * 4],
                "mag-id": ["m1", "m2", "m3", "m4"] * 3,
                "length": [2000, 4000, 1000, 10000] * 3,
                "numreads": [10, 20, 5, 0, 12, 25, 8, 0, 30, 60, 15, 1],
            }
        )
        self.mags_derep = qiime2.Artifact.import_data(
            "FeatureData[MAG]", self.get_data_path("mag-sequences")
        )
        self.mag_length_df = pd.read_csv(
            self.get_data_path("mag-length.tsv"), sep="\t", index_col=0
        )
        self.mag_length = qiime2.Artifact.import_data(
            'FeatureData[SequenceCharacteristics % Properties("length")]',
            self.mag_length_df,
        )
        self.mapped_reads = qiime2.Artifact.import_data(
            "FeatureData[AlignmentMap]", self.get_data_path("reads-mapped")
        )

    def test_rpkm(self):
        obs = rpkm(self.df)
        exp = pd.Series(
            [
                142857.14,
                142857.14,
                142857.14,
                0.0,
                133333.33,
                138888.89,
                177777.78,
                0.0,
                141509.43,
                141509.43,
                141509.43,
                943.40,
            ]
        )
        pd.testing.assert_series_equal(obs, exp)

    def test_tpm(self):
        obs = tpm(self.df)
        exp = pd.Series(
            [
                333333.33,
                333333.33,
                333333.33,
                0.0,
                296296.3,
                308641.98,
                395061.73,
                0.0,
                332594.24,
                332594.24,
                332594.24,
                2217.29,
            ]
        )
        pd.testing.assert_series_equal(obs, exp)

    @parameterized.expand(["legacy", "plain", "mixed", "internal_alignment"])
    def test_estimate_mag_abundance(self, naming):
        maps = BAMDirFmt()
        sample_ids = []
        for i, source in enumerate(
            sorted(self.mapped_reads.view(BAMDirFmt).path.iterdir()), start=1
        ):
            sample_id = f"sample{i}"
            if naming == "internal_alignment" and i == 1:
                sample_id = "sample1_alignment_batch"
            sample_ids.append(sample_id)
            suffix = (
                "_alignment"
                if naming == "legacy" or (naming == "mixed" and i % 2)
                else ""
            )
            shutil.copyfile(source, maps.path / f"{sample_id}{suffix}.bam")

        obs = estimate_abundance(
            alignment_maps=maps,
            feature_lengths=self.mag_length_df,
            metric="rpkm",
        )
        exp = pd.DataFrame(
            {
                "2317a4bd-778d-46fa-901c-88428b2b863e": [
                    0.0,
                    285.7,
                    142.9,
                    271.4,
                    14.3,
                ],
                "890676ee-1310-477d-b713-457ef661b8f3": [
                    274.0,
                    0.0,
                    137.0,
                    13.7,
                    260.3,
                ],
            },
            index=pd.Index(sample_ids, name="sample-id"),
        )
        exp.columns.name = "feature-id"
        pd.testing.assert_frame_equal(obs, exp, check_exact=False, rtol=0.1)

        # just double-check that we have what we expect:
        # 2317: E. coli, 8906: M. tuberculosis
        #           |     2317    |   8906
        # sample1   |   0 reads   | 100 reads
        # sample2   |   100 reads | 0 reads
        # sample3   |   50 reads  | 50 reads
        # sample4   |   95 reads  | 5 reads
        # sample5   |   5 reads   | 95 reads
        obs_rel = obs.div(obs.sum(axis=1), axis=0)
        exp_rel = pd.DataFrame(
            {
                "2317a4bd-778d-46fa-901c-88428b2b863e": [0.0, 1.0, 0.5, 0.95, 0.05],
                "890676ee-1310-477d-b713-457ef661b8f3": [1.0, 0.0, 0.5, 0.05, 0.95],
            },
            index=pd.Index(sample_ids, name="sample-id"),
        )
        exp_rel.columns.name = "feature-id"
        pd.testing.assert_frame_equal(obs_rel, exp_rel, check_exact=False, rtol=0.1)


class TestContigAbundance(TestPluginBase):
    package = "q2_mag.abundance.tests"

    @parameterized.expand(
        [
            ("raw_contigs", ("k141_0", "k141_1")),
            (
                "sample_underscore",
                ("sample_A:dPj72SpWBxMXAVyPtYqJGA", "sample_A:2Ta5AWdPGsJiCUZJ4NeP92"),
            ),
            (
                "underscore_separator",
                ("pooled_dPj72SpWBxMXAVyPtYqJGA", "pooled_2Ta5AWdPGsJiCUZJ4NeP92"),
            ),
            ("matching_prefix", ("contig", "contig_part")),
        ]
    )
    def test_estimate_contig_abundance(self, name, contig_ids):
        maps = BAMDirFmt()
        sample_id = "sample_A"
        lengths = [1000, 2000]
        header = {
            "HD": {"VN": "1.6", "SO": "coordinate"},
            "SQ": [
                {"SN": identifier, "LN": length}
                for identifier, length in zip(contig_ids, lengths)
            ],
        }
        with pysam.AlignmentFile(
            str(maps.path / f"{sample_id}.bam"), "wb", header=header
        ) as bam:
            for reference_id, count in enumerate([2, 1]):
                for i in range(count):
                    read = pysam.AlignedSegment()
                    read.query_name = f"read{reference_id}_{i}"
                    read.query_sequence = "ACGT" * 25
                    read.flag = 0
                    read.reference_id = reference_id
                    read.reference_start = i * 10
                    read.mapping_quality = 60
                    read.cigar = [(0, 100)]
                    read.query_qualities = pysam.qualitystring_to_array("I" * 100)
                    bam.write(read)

        alignment = qiime2.Artifact.import_data("SampleData[AlignmentMap]", maps)
        feature_lengths = qiime2.Artifact.import_data(
            "FeatureData[SequenceCharacteristics % Properties('length')]",
            pd.DataFrame({"length": lengths}, index=pd.Index(contig_ids, name="id")),
        )
        action = self.plugin.methods["estimate_abundance"]
        for metric, expected in (
            ("rpkm", [2e6 / 3, 1e6 / 6]),
            ("tpm", [8e5, 2e5]),
        ):
            with self.subTest(metric=metric):
                (abundance,) = action(alignment, feature_lengths, metric=metric)
                abundance.validate("max")
                table = abundance.view(biom.Table)
                self.assertEqual(list(table.ids(axis="sample")), [sample_id])
                self.assertEqual(set(table.ids(axis="observation")), set(contig_ids))
                observed = table.to_dataframe(dense=True).loc[
                    list(contig_ids), sample_id
                ]
                np.testing.assert_allclose(observed.to_numpy(), expected)
