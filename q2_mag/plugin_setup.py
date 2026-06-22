# ----------------------------------------------------------------------------
# Copyright (c) 2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import importlib

from qiime2.plugin import Metadata
from q2_mag.busco.types import (
    BUSCOResultsFormat,
    BUSCOResultsDirectoryFormat,
    BuscoDatabaseDirFmt,
    BUSCOResults,
    BUSCO,
)
from q2_types.distance_matrix import DistanceMatrix
from q2_types.feature_data import (
    FeatureData,
    Sequence,
    SequenceCharacteristics,
)
from q2_types.feature_table import (
    FeatureTable,
    PresenceAbsence,
    Frequency,
)
from q2_types.per_sample_sequences import (
    MAGs,
    Contigs,
)
from q2_types.sample_data import SampleData
from q2_types.feature_map import FeatureMap, MAGtoContigs
from qiime2.core.type import (
    Bool,
    Range,
    Int,
    Str,
    Float,
    List,
    Choices,
    Visualization,
    Properties,
    TypeMap,
)
from qiime2.plugin import Plugin, Citations
import q2_mag
from q2_types.feature_data_mag import MAG
from q2_types.per_sample_sequences import AlignmentMap
from q2_types.reference_db import ReferenceDB
from q2_types.genome_data import GenomeData, Proteins

from q2_mag import __version__

citations = Citations.load("citations.bib", package="q2_mag")

plugin = Plugin(
    name="mag",
    version=__version__,
    website="https://github.com",
    package="q2_mag",
    description=(
        "Rachis plugin for contig binning, bin quality assessment "
        "and MAG dereplication."
    ),
    short_description=(
        "Rachis plugin for contig binning, bin quality assessment "
        "and MAG dereplication."
    ),
    citations=[],
)

importlib.import_module("q2_mag.metabat2")
importlib.import_module("q2_mag.das_tool")

partition_params = {"num_partitions": Int % Range(1, None)}
partition_param_descriptions = {
    "num_partitions": (
        "The number of partitions to split the data "
        "into. Defaults to partitioning into individual "
        "samples."
    )
}

plugin.methods.register_function(
    function=q2_mag.metabat2.bin_contigs_metabat,
    inputs={
        "contigs": SampleData[Contigs],
        "alignment_maps": SampleData[AlignmentMap % Properties("sorted")],
    },
    parameters={
        "min_contig": Int % Range(1500, None),
        "max_p": Int % Range(1, 100),
        "min_s": Int % Range(1, 100),
        "max_edges": Int % Range(1, None),
        "p_tnf": Int % Range(0, 100),
        "no_add": Bool,
        "min_cv": Int % Range(1, None),
        "min_cv_sum": Int % Range(1, None),
        "min_cls_size": Int % Range(1, None),
        "num_threads": Int % Range(0, None),
        "seed": Int % Range(0, None),
        "debug": Bool,
        "verbose": Bool,
    },
    outputs=[
        ("mags", SampleData[MAGs]),
        ("contig_map", FeatureMap[MAGtoContigs]),
        ("unbinned_contigs", SampleData[Contigs % Properties("unbinned")]),  # ??
    ],
    input_descriptions={
        "contigs": "Contigs to be binned.",
        "alignment_maps": "Reads-to-contig alignment maps.",
    },
    parameter_descriptions={
        "min_contig": "Minimum size of a contig for binning.",
        "max_p": (
            'Percentage of "good" contigs considered for binning '
            "decided by connection among contigs. The greater, the "
            "more sensitive."
        ),
        "min_s": (
            "Minimum score of a edge for binning. The greater, the more specific."
        ),
        "max_edges": (
            "Maximum number of edges per node. The greater, the more sensitive."
        ),
        "p_tnf": (
            "TNF probability cutoff for building TNF graph. Use it to "
            "skip the preparation step. (0: auto)"
        ),
        "no_add": "Turning off additional binning for lost or small contigs.",
        "min_cv": "Minimum mean coverage of a contig in each library for binning.",
        "min_cv_sum": (
            "Minimum total effective mean coverage of a contig "
            "(sum of depth over minCV) for binning."
        ),
        "min_cls_size": "Minimum size of a bin as the output.",
        "num_threads": "Number of threads to use (0: use all cores).",
        "seed": "For exact reproducibility. (0: use random seed)",
        "debug": "Debug output.",
        "verbose": "Verbose output.",
    },
    output_descriptions={
        "mags": "The resulting MAGs.",
        "contig_map": (
            "Mapping of MAG identifiers to the contig identifiers "
            "contained in each MAG."
        ),
        "unbinned_contigs": "Contigs that were not binned into any MAG.",
    },
    name="Bin contigs into MAGs using MetaBAT 2.",
    description=("This method uses MetaBAT 2 to bin provided contigs into MAGs."),
    citations=[
        citations["kang2019"],
        citations["heng2009samtools"],
        citations["scikit_bio_release"],
    ],
)

plugin.methods.register_function(
    function=q2_mag.dereplication.dereplicate_mags,
    inputs={
        "mags": SampleData[MAGs],
        "distance_matrix": DistanceMatrix,
    },
    parameters={
        "threshold": Float % Range(0, 1, inclusive_end=True),
        "metadata": Metadata,
        "metadata_column": Str,
        "find_max": Bool,
    },
    outputs=[
        ("dereplicated_mags", FeatureData[MAG]),
        ("table", FeatureTable[PresenceAbsence]),
    ],
    input_descriptions={
        "mags": "MAGs to be dereplicated.",
        "distance_matrix": "Matrix of distances between MAGs.",
    },
    parameter_descriptions={
        "threshold": ("Similarity threshold required to consider two bins identical."),
        "metadata": "Metadata table.",
        "metadata_column": (
            "Name of the metadata column used to choose the "
            "most representative bins."
        ),
        "find_max": (
            "Set to True to choose the bin with the highest value in "
            "the metadata column. Set to False to choose the bin "
            "with the lowest value."
        ),
    },
    output_descriptions={
        "dereplicated_mags": "Dereplicated MAGs.",
        "table": "Mapping between MAGs and samples.",
    },
    name="Dereplicate MAGs from multiple samples.",
    description=(
        "This method dereplicates MAGs from multiple samples "
        "using distances between them found in the provided "
        "distance matrix. For each cluster of similar MAGs, "
        "the longest one will be selected as the representative. If "
        "metadata is given as input, the MAG with the "
        "highest or lowest value in the specified metadata column "
        'is chosen, depending on the parameter "find-max". '
        "If there are MAGs with identical values, the longer one is "
        "chosen. For example an artifact of type BUSCOResults can be "
        "passed as metadata, and the dereplication can be done by "
        'highest "completeness".'
    ),
    citations=[],
)

busco_params = {
    "mode": Str % Choices(["genome"]),
    "lineage_dataset": Str,
    "augustus": Bool,
    "augustus_parameters": Str,
    "augustus_species": Str,
    "cpu": Int % Range(1, None),
    "contig_break": Int % Range(0, None),
    "evalue": Float % Range(0, None, inclusive_start=False),
    "limit": Int % Range(1, 20),
    "long": Bool,
    "metaeuk_parameters": Str,
    "metaeuk_rerun_parameters": Str,
    "miniprot": Bool,
    "additional_metrics": Bool,
}
busco_param_descriptions = {
    "mode": (
        "Specify which BUSCO analysis mode to run."
        "Currently only the 'genome' option is supported, "
        "for genome assemblies. In the future modes for transcriptome "
        "assemblies and for magd gene sets (proteins) will be made "
        "available."
    ),
    "lineage_dataset": (
        "Specify the name of the BUSCO lineage to be used. "
        "To see all possible options run `busco --list-datasets`."
    ),
    "augustus": "Use augustus gene predictor for eukaryote runs.",
    "augustus_parameters": (
        "Pass additional arguments to Augustus. All arguments should be contained "
        "within a single string with no white space, with each argument "
        "separated by a comma. Example: '--PARAM1=VALUE1,--PARAM2=VALUE2'."
    ),
    "augustus_species": "Specify a species for Augustus training.",
    "cpu": "Specify the number (N=integer) of threads/cores to use.",
    "contig_break": (
        "Number of contiguous Ns to signify a break between contigs. "
        "See https://gitlab.com/ezlab/busco/-/issues/691 for a "
        "more detailed explanation."
    ),
    "evalue": ("E-value cutoff for BLAST searches. Allowed formats, 0.001 or 1e-03."),
    "limit": (
        "How many candidate regions (contig or transcript) to consider per BUSCO."
    ),
    "long": (
        "Optimization Augustus self-training mode (Default: Off); "
        "adds considerably to the run time, "
        "but can improve results for some non-model organisms."
    ),
    "metaeuk_parameters": (
        "Pass additional arguments to Metaeuk for the first run. All arguments "
        "should be contained within a single string with no white space, with each "
        "argument separated by a comma. Example: `--PARAM1=VALUE1,--PARAM2=VALUE2`."
    ),
    "metaeuk_rerun_parameters": (
        "Pass additional arguments to Metaeuk for the second run. All arguments "
        "should be contained within a single string with no white space, with "
        "each argument separated by a comma. Example: "
        "`--PARAM1=VALUE1,--PARAM2=VALUE2`."
    ),
    "miniprot": "Use miniprot gene predictor for eukaryote runs.",
    "additional_metrics": (
        "Adds completeness and contamination values to the BUSCO "
        "report. Check here for documentation: https://github.com/"
        "metashot/busco?tab=readme-ov-file#documetation"
    ),
}


plugin.methods.register_function(
    function=q2_mag.busco.collate_busco_results,
    inputs={"results": List[BUSCOResults]},
    parameters={},
    outputs={"collated_results": BUSCOResults},
    name="Collate BUSCO results.",
    description="Collates BUSCO results.",
)

plugin.visualizers.register_function(
    function=q2_mag.busco._visualize_busco,
    inputs={
        "results": BUSCOResults,
        "unbinned_contigs": SampleData[Contigs],
    },
    parameters={},
    input_descriptions={
        "results": "BUSCO results table.",
        "unbinned_contigs": "Contigs which were not assigned to any bin.",
    },
    parameter_descriptions={},
    name="Visualize BUSCO results.",
    description=("This method generates a visualization from the BUSCO results table."),
    citations=[citations["manni_busco_2021"]],
)

plugin.methods.register_function(
    function=q2_mag.busco._evaluate_busco,
    inputs={
        "mags": SampleData[MAGs] | FeatureData[MAG],
        "unbinned_contigs": SampleData[Contigs],
        "db": ReferenceDB[BUSCO],
    },
    parameters=busco_params,
    outputs={"results": BUSCOResults},
    input_descriptions={
        "mags": "MAGs to be analyzed.",
        "db": "BUSCO database.",
        "unbinned_contigs": "Contigs which were not assigned to any bin.",
    },
    parameter_descriptions=busco_param_descriptions,
    output_descriptions={"results": "BUSCO result table."},
    name="Evaluate quality of the generated MAGs using BUSCO.",
    description=(
        "This method uses BUSCO to assess the quality of assembled MAGs "
        "and generates a table summarizing the results."
    ),
    citations=[citations["manni_busco_2021"]],
)

plugin.pipelines.register_function(
    function=q2_mag.busco.evaluate_busco,
    inputs={
        "mags": SampleData[MAGs] | FeatureData[MAG],
        "unbinned_contigs": SampleData[Contigs],
        "db": ReferenceDB[BUSCO],
    },
    parameters={**busco_params, **partition_params},
    outputs={"results": BUSCOResults, "visualization": Visualization},
    input_descriptions={
        "mags": "MAGs to be analyzed.",
        "db": "BUSCO database.",
        "unbinned_contigs": "Contigs which were not assigned to any bin.",
    },
    parameter_descriptions={**busco_param_descriptions, **partition_param_descriptions},
    output_descriptions={
        "results": "BUSCO result table.",
        "visualization": "Visualization of the BUSCO results.",
    },
    name="Evaluate quality of the generated MAGs using BUSCO.",
    description=(
        "This method uses BUSCO to assess the quality of assembled "
        "MAGs and generates a table summarizing the results."
    ),
    citations=[citations["manni_busco_2021"]],
)

plugin.methods.register_function(
    function=q2_mag.busco.fetch_busco_db,
    inputs={},
    outputs=[("db", ReferenceDB[BUSCO])],
    output_descriptions={"db": "BUSCO database for the specified lineages."},
    parameters={
        "lineages": List[Str],
    },
    parameter_descriptions={
        "lineages": (
            "Lineages to be included in the database. Can be any "
            "valid BUSCO lineage or any of the following: 'all' "
            "(for all lineages), 'prokaryota', 'eukaryota', 'virus'."
        ),
    },
    name="Download BUSCO database.",
    description=(
        "Downloads BUSCO database for the specified lineage. "
        "Output can be used to run BUSCO with the 'evaluate-busco' "
        "action."
    ),
    citations=[citations["manni_busco_2021"]],
)

filter_params = {
    "metadata": Metadata,
    "where": Str,
    "exclude_ids": Bool,
    "remove_empty": Bool,
}
filter_param_descriptions = {
    "metadata": (
        "Sample metadata indicating which MAG ids to filter. "
        "The optional `where` parameter may be used to filter ids "
        "based on specified conditions in the metadata. The "
        "optional `exclude_ids` parameter may be used to exclude "
        "the ids specified in the metadata from the filter."
    ),
    "where": (
        "Optional SQLite WHERE clause specifying MAG metadata "
        "criteria that must be met to be included in the filtered "
        "data. If not provided, all MAGs in `metadata` that are "
        "also in the MAG data will be retained."
    ),
    "exclude_ids": (
        "Defaults to False. If True, the MAGs selected by "
        "the `metadata` and optional `where` parameter will be "
        "excluded from the filtered data."
    ),
}

plugin.methods.register_function(
    function=q2_mag.filtering.filter_derep_mags,
    inputs={"mags": FeatureData[MAG]},
    parameters=filter_params,
    outputs={"filtered_mags": FeatureData[MAG]},
    input_descriptions={"mags": "Dereplicated MAGs to filter."},
    parameter_descriptions={
        **filter_param_descriptions,
        "remove_empty": "Remove empty MAGs.",
    },
    name="Filter dereplicated MAGs.",
    description="Filter dereplicated MAGs based on metadata.",
)

plugin.methods.register_function(
    function=q2_mag.filtering.filter_mags,
    inputs={"mags": SampleData[MAGs]},
    parameters={
        **filter_params,
        "on": Str % Choices(["sample", "mag"]),
    },
    outputs={"filtered_mags": SampleData[MAGs]},
    input_descriptions={"mags": "MAGs to filter."},
    parameter_descriptions={
        **filter_param_descriptions,
        "on": "Whether to filter based on sample or MAG metadata.",
        "remove_empty": "Remove empty MAGs.",
    },
    name="Filter MAGs.",
    description="Filter MAGs based on metadata.",
)

plugin.methods.register_function(
    function=q2_mag.utils.get_feature_lengths,
    inputs={
        "features": FeatureData[MAG | Sequence] | SampleData[MAGs | Contigs],
    },
    parameters={},
    outputs=[("lengths", FeatureData[SequenceCharacteristics % Properties("length")])],
    input_descriptions={"features": "Features to get lengths for."},
    parameter_descriptions={},
    output_descriptions={
        "lengths": "Feature lengths.",
    },
    name="Get feature lengths.",
    description="This method extract lengths for the provided feature set.",
    citations=[],
)

M_abundance_in, P_abundance_out = TypeMap(
    {
        Str % Choices(["rpkm"]): Properties("rpkm"),
        Str % Choices(["tpm"]): Properties("tpm"),
    }
)

plugin.methods.register_function(
    function=q2_mag.abundance.estimate_abundance,
    inputs={
        "alignment_maps": FeatureData[AlignmentMap] | SampleData[AlignmentMap],
        "feature_lengths": FeatureData[SequenceCharacteristics % Properties("length")],
    },
    parameters={
        "metric": M_abundance_in,
        "min_mapq": Int % Range(0, 255),
        "min_query_len": Int % Range(0, None),
        "min_base_quality": Int % Range(0, None),
        "min_read_len": Int % Range(0, None),
        "threads": Int % Range(1, None),
    },
    outputs=[
        ("abundances", FeatureTable[Frequency % P_abundance_out]),
    ],
    input_descriptions={
        "alignment_maps": (
            "Bowtie2 alignment maps between reads and features "
            "for which the abundance should be estimated."
        ),
        "feature_lengths": "Table containing length of every feature (MAG/contig).",
    },
    parameter_descriptions={
        "metric": "Metric to be used as a proxy of feature abundance.",
        "min_mapq": "Minimum mapping quality.",
        "min_query_len": "Minimum query length.",
        "min_base_quality": "Minimum base quality.",
        "min_read_len": "Minimum read length.",
        "threads": "Number of threads to pass to samtools.",
    },
    output_descriptions={
        "abundances": "Feature abundances.",
    },
    name="Estimate feature (MAG/contig) abundance.",
    description=(
        "This method estimates MAG/contig abundances by mapping the "
        "reads to them and calculating respective metric values"
        "which are then used as a proxy for the frequency."
    ),
    citations=[],
)

plugin.methods.register_function(
    function=q2_mag.semibin2.bin_contigs_semibin2,
    inputs={
        "contigs": SampleData[Contigs],
        "alignment_maps": SampleData[AlignmentMap % Properties("sorted")],
    },
    parameters={
        # "mode": Str % Choices("single", "multi"),
        "training_type": Str % Choices("semi", "self"),
        "orf_finder": Str % Choices("fast-naive", "prodigal", "fraggenescan"),
        "environment": Str
        % Choices(
            "human_gut",
            "dog_gut",
            "ocean",
            "soil",
            "cat_gut",
            "human_oral",
            "mouse_gut",
            "pig_gut",
            "built_environment",
            "wastewater",
            "chicken_caecum",
            "global",
        ),
        "engine": Str % Choices("auto", "gpu", "cpu"),
        "sequencing_type": Str % Choices("short_read", "long_read"),
        "minfasta_kbs": Int % Range(1, None),
        "no_recluster": Bool,
        "epochs": Int % Range(1, None),
        "batch_size": Int % Range(1, None),
        "max_node": Int % Range(1, None),
        "max_edges": Int % Range(1, None),
        "ratio": Float % Range(0.0, None),
        "min_len": Int % Range(1, None),
        "ml_threshold": Int % Range(1, None),
        "threads": Int % Range(0, None),
        "random_seed": Int % Range(0, None),
        "debug": Bool,
    },
    outputs=[
        ("mags", SampleData[MAGs]),
        ("contig_map", FeatureMap[MAGtoContigs]),
    ],
    input_descriptions={
        "contigs": "Contigs to be binned.",
        "alignment_maps": "Reads-to-contig alignment maps.",
    },
    parameter_descriptions={
        # "mode": "Binning mode controlling how coverage is used for embedding.",
        "training_type": "Training algorithm used to train the model.",
        "orf_finder": "Gene predictor used to estimate the number of bins.",
        "environment": "Which pre-trained model to use.",
        "engine": "Device used to train the model.",
        "sequencing_type": (
            "Specify whether your data consists of short- or long-reads. For hybrid "
            "data (long- and short-reads), it is recommended to use the long-reads "
            "pipeline."
        ),
        "minfasta_kbs": "Minimum bin size in kilo-basepairs.",
        "no_recluster": (
            "Do not recluster bins. This saves a small amount of time, but "
            "pre-reclustering bins are always output."
        ),
        "epochs": "Number of epochs used in the training process.",
        "batch_size": "Number of epochs used in the training process.",
        "max_node": "Percentage of contigs that considered to be binned.",
        "max_edges": "The maximum number of edges that can be connected to one contig.",
        "ratio": (
            "If the ratio of the number of base pairs of contigs between 1000-2500 bp "
            "smaller than this value, the minimal length will be set as 1000 bp, "
            "otherwise 2500 bp. Setting --p-min-length overrules this parameter."
        ),
        "min_len": (
            "Minimal contig length (bp) to include in binning. Contigs shorter than "
            "this length are excluded. This parameter overrules --p-ratio."
        ),
        "ml_threshold": (
            "Length threshold for generating must-link constraints. (By default, the "
            "threshold is calculated from the contig, and the default minimum value is "
            "4,000 bp)."
        ),
        "threads": "Number of threads to use (0: use all cores).",
        "random_seed": "For exact reproducibility. (0: use random seed)",
        "debug": "Debug output.",
    },
    output_descriptions={
        "mags": "The resulting MAGs.",
        "contig_map": (
            "Mapping of MAG identifiers to the contig identifiers "
            "contained in each MAG."
        ),
    },
    name="Bin contigs into MAGs using SemiBin2.",
    description=(
        "This method uses SemiBin2 to bin provided contigs into MAGs. Note that "
        "SemiBin2 does not output what it considers ‘unbinned’ contigs"
    ),
    citations=[
        citations["pan_deep_2022"],
        citations["pan_semibin2_2023"],
    ],
)


plugin.methods.register_function(
    function=q2_mag.das_tool.refine_bins_das_tool,
    inputs={
        "bins": List[SampleData[MAGs]],
        "contigs": SampleData[Contigs],
        "proteins": GenomeData[Proteins],
    },
    parameters={
        "labels": Str,
        # NOTE: usearch as a search engine is disabled until added as a dependency
        "search_engine": Str % Choices("diamond", "blastp"),  # TODO: add usearch
        "score_threshold": Float % Range(0.0, 1.0),
        "duplicate_penalty": Float % Range(0.0, 3.0),
        "megabin_penalty": Float % Range(0.0, 3.0),
        "max_iter_post_threshold": Int % Range(0, None),
        "threads": Int % Range(0, None),
        "debug": Bool,
    },
    outputs={
        "refined_bins": SampleData[MAGs],
    },
    input_descriptions={
        "bins": "Bins produced by a metagenomic binning tool.",
        "contigs": (
            "Contig sequences used for metagenomic binning. These sequences must "
            "correspond to the contigs that were used to generate the bins."
        ),
        "proteins": (
            "Predicted protein sequences derived from the input contigs using "
            "Prodigal."
        ),
    },
    parameter_descriptions={
        "labels": (
            "Comma-separated list of metagenomic binning methods names. "
            "The number of labels must match the number of bins. "
            "If unspecified, binning methods will be labeled binning_1, binning_2,..."
            "Duplicate label names are not allowed."
            "For exampe: `metabat,semibin`."
        ),
        "search_engine": "Engine used for single copy gene identification.",
        "score_threshold": (
            "Score threshold until selection algorithm will keep selecting bins."
        ),
        "duplicate_penalty": (
            "Penalty for duplicate single copy genes per bin (weight b)."
        ),
        "megabin_penalty": "Penalty for megabins (weight c).",
        "max_iter_post_threshold": (
            "Maximum number of iterations after reaching score threshold."
        ),
        "threads": "Number of threads to use (0: use all cores).",
        "debug": "Debug output.",
    },
    output_descriptions={
        "refined_bins": "The binned contigs created by DAS Tool.",
    },
    name="Refine bins produced by 2+ binning methods using DAS Tool.",
    description=(
        "This method uses DAS Tool to integrate multiple binning prediction to "
        "produce an optimized, non-redundant set of bins."),
    citations=[
        citations["sieber2018recovery"],
    ],
)

plugin.register_semantic_types(BUSCOResults, BUSCO)
plugin.register_formats(
    BUSCOResultsFormat, BUSCOResultsDirectoryFormat, BuscoDatabaseDirFmt
)
plugin.register_semantic_type_to_format(
    BUSCOResults, artifact_format=BUSCOResultsDirectoryFormat
)
plugin.register_semantic_type_to_format(ReferenceDB[BUSCO], BuscoDatabaseDirFmt)
importlib.import_module("q2_mag.busco.types._transformer")
