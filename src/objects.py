from collections import namedtuple

FastaDatabase = namedtuple("FastaDatabase", ["file_path", "proteins", "kmers"])
Precursor = namedtuple(
    "Precursor",
    ["id", "description", "mass", "charge", "retention_time", "abundance", "fragments"],
)
Fragment = namedtuple(
    "Fragment",
    [
        "id",
        "precursor_id",
        "precursor_mass",
        "precursor_charge",
        "mz_value",
        "abundance",
    ],
)
ExperimentParameters = namedtuple(
    "ExperimentParameters",
    [
        "id",
        "precursors",
        "sqllite_database",
        "max_peptide_length",
        "ppm_tolerance",
        "precursor_tolerance",
        "number_hybrids",
        "number_natives",
        "target_seq",
        "number_decimal_places",
    ],
)
KMer = namedtuple(
    "KMer",
    [
        "fragment_id",
        "precursor_mass",
        "precursor_charge",
        "protein_id",
        "kmer_mass",
        "location_start",
        "location_end",
        "ion",
        "charge",
        "subsequence",
        "kmer_type",
    ],
)
MatchedFragment = namedtuple("MatchedFragment", ["fragment", "kmers"])
Protein = namedtuple("Protein", ["id", "description", "sequence"])
MatchedProtein = namedtuple("MatchedProtein", ["protein", "kmers"])
Cluster = namedtuple("Cluster", ["protein", "ion", "longest_kmer", "score"])
ExtendedCluster = namedtuple("ExtendedCluster", ["cluster", "extended_sequence"])
Peptide = namedtuple(
    "Peptide", ["peptide_type", "b_extended_cluster", "y_extended_cluster", "score"]
)
AlignedPeptide = namedtuple(
    "AlignedPeptide",
    [
        "hybrid",
        "left_protein",
        "right_protein",
        "sequence",
        "b_score",
        "y_score",
        "total_score",
        "total_gaussian_score",
        "extensions",
        "precursor_mass",
        "precursor_charge",
        "total_mass_error",
        "total_count",
    ],
)
