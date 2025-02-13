from collections import defaultdict
from itertools import groupby
from operator import attrgetter, itemgetter
from typing import List

from src.erik import Peak, relative_ppm_tolerance_in_daltons
from src.lookups.constants import AMINO_ACIDS
from src.lookups.data_classes import Peak
from src.lookups.protein_product_ion_db import ProteinProductIonDb
from src.objects import (
    AlignedPeptide,
    Cluster,
    ExperimentParameters,
    ExtendedCluster,
    Fragment,
    KMer,
    MatchedFragment,
    MatchedProtein,
    Peptide,
    Precursor,
    Protein,
)


def get_filtered_fragments(
    precursors: List[Precursor], ppm_tolerance
) -> List[Fragment]:
    filtered_fragments = []
    for precursor in precursors:
        for fragment in precursor.fragments:
            if fragment.abundance >= ppm_tolerance:
                filtered_fragments.append(fragment)
    return filtered_fragments


def ppm_to_da(mass, ppm_tolerance):
    return abs((ppm_tolerance / 1000000) * mass)


def get_matching_ions_for_peak(
    peak: Peak,
):
    pass


def get_matched_fragment(
    fragment: Peak,
    db: ProteinProductIonDb,
    ppm_tolerance: float,
    number_decimal_places: int,
):
    fragment_id = fragment.id
    precursor_mass = fragment.precursor_mass
    precursor_charge = fragment.precursor_charge
    mz_value = fragment.mz_value
    adjusted_tolerance = ppm_to_da(mz_value, ppm_tolerance)
    b_rows, y_rows = db.query_mass_kmers(
        fragment_id,
        precursor_mass,
        precursor_charge,
        mz_value,
        adjusted_tolerance,
        number_decimal_places,
    )
    b_kmers = [KMer(*row) for row in b_rows]
    y_kmers = [KMer(*row) for row in y_rows]
    kmers = b_kmers + y_kmers
    kmers = [kmer._replace(ion="b" if kmer.ion == 0 else "y") for kmer in kmers]
    matched_fragment = MatchedFragment(fragment=fragment, kmers=kmers)
    return matched_fragment


def get_matched_fragments(
    fragments: List[Fragment],
    sqllite_database: ProteinProductIonDb,
    ppm_tolerance,
    number_decimal_places,
):
    matched_fragments = []
    for fragment in fragments:
        matched_fragment = get_matched_fragment(
            fragment, sqllite_database, ppm_tolerance, number_decimal_places
        )
        if len(matched_fragment.kmers) > 0:
            matched_fragments.append(matched_fragment)
    return matched_fragments


def get_matched_proteins(matched_fragments, sqllite_database: ProteinProductIonDb):
    all_kmers = []
    for matched_fragment in matched_fragments:
        all_kmers.extend(matched_fragment.kmers)
    sorted_kmers = sorted(all_kmers, key=lambda k: k.protein_id)
    grouped_kmers = {
        protein_id: list(group)
        for protein_id, group in groupby(sorted_kmers, key=lambda k: k.protein_id)
    }

    matched_proteins = []
    for protein_id, kmers in grouped_kmers.items():
        protein_row = sqllite_database.get_protein_by_id(protein_id)
        protein = Protein(*protein_row)
        matched_protein = MatchedProtein(protein=protein, kmers=kmers)
        matched_proteins.append(matched_protein)

    return matched_proteins


def get_clusters(matched_proteins):
    clusters = []
    for matched_protein in matched_proteins:
        protein = matched_protein.protein

        b_kmers = [k for k in matched_protein.kmers if k.ion == "b"]
        y_kmers = [k for k in matched_protein.kmers if k.ion == "y"]

        for key, group in groupby(
            sorted(b_kmers, key=attrgetter("location_start")),
            key=attrgetter("location_start"),
        ):
            kmers_in_group = list(group)
            longest_kmer = max(kmers_in_group, key=lambda kmer: len(kmer.subsequence))
            score = len(kmers_in_group)
            clusters.append(
                Cluster(
                    protein=matched_protein.protein,
                    ion="b",
                    longest_kmer=longest_kmer,
                    score=score,
                )
            )

        for key, group in groupby(
            sorted(y_kmers, key=attrgetter("location_end")),
            key=attrgetter("location_end"),
        ):
            kmers_in_group = list(group)
            longest_kmer = max(kmers_in_group, key=lambda kmer: len(kmer.subsequence))
            score = len(kmers_in_group)
            clusters.append(
                Cluster(
                    protein=matched_protein.protein,
                    ion="y",
                    longest_kmer=longest_kmer,
                    score=score,
                )
            )

    return clusters


def get_extended_clusters(clusters):
    extended_clusters = []
    bad_chars = ["B", "X", "U", "Z", "O", "J"]
    for cluster in clusters:
        protein = cluster.protein
        longest_kmer = cluster.longest_kmer
        precursor_mass = longest_kmer.precursor_mass
        cumulative_mass = longest_kmer.kmer_mass
        extended_sequence = ""
        if longest_kmer.ion == "b":
            position = longest_kmer.location_end
            while position < len(protein.sequence):
                next_aa = protein.sequence[position]
                if next_aa in bad_chars:
                    position += 1
                    continue
                next_aa_mass = AMINO_ACIDS.get(next_aa, 0)
                if cumulative_mass + next_aa_mass > precursor_mass:
                    break
                extended_sequence += next_aa
                cumulative_mass += next_aa_mass
                position += 1
        elif longest_kmer.ion == "y":
            position = longest_kmer.location_start - 1
            while position >= 0:
                next_aa = protein.sequence[position]
                if next_aa in bad_chars:
                    position -= 1
                    continue
                next_aa_mass = AMINO_ACIDS.get(next_aa, 0)
                if cumulative_mass + next_aa_mass > precursor_mass:
                    break
                extended_sequence = next_aa + extended_sequence
                cumulative_mass += next_aa_mass
                position -= 1
        extended_cluster = ExtendedCluster(
            cluster=cluster, extended_sequence=extended_sequence
        )

    extended_clusters.append(extended_cluster)

    return extended_clusters


def get_best_extended_cluster_match(
    target_extended_cluster, complimentary_extended_clusters
):
    target_precursor_mass = target_extended_cluster.cluster.longest_kmer.precursor_mass
    target_weight = target_extended_cluster.cluster.longest_kmer.kmer_mass
    best_extended_cluster = None
    max_peptide_weight = 0

    for complimentary_extended_cluster in complimentary_extended_clusters:
        cluster = complimentary_extended_cluster.cluster
        compliment_weight = cluster.longest_kmer.kmer_mass
        peptide_weight = target_weight + compliment_weight
        if (
            peptide_weight <= target_precursor_mass
            and peptide_weight > max_peptide_weight
        ):
            best_extended_cluster = complimentary_extended_cluster
            max_peptide_weight = peptide_weight
    return best_extended_cluster


def get_peptide_type(b_cluster, y_cluster, sqllite_database):
    if b_cluster is None or y_cluster is None:
        return "u"
    if b_cluster.protein.id != y_cluster.protein.id:
        return "h"
    else:
        return "n"


def get_kmer_overlap(b_seq, y_seq):
    min_len = min(len(b_seq), len(y_seq))
    overlap = 0

    for i in range(1, min_len + 1):
        if b_seq[-i:] == y_seq[:i]:
            overlap = i
    return overlap


def get_peptide_score(b_cluster, y_cluster, protein_length):
    combined_sequence_length = len(b_cluster.longest_kmer) + len(y_cluster.longest_kmer)
    return combined_sequence_length / protein_length


def get_peptide(
    target_extended_cluster, complimentary_extended_clusters, sqllite_database
):
    target_cluster = target_extended_cluster.cluster
    protein_length = len(target_cluster.protein.sequence)
    if target_cluster.ion == "b":
        y_cluster = get_best_extended_cluster_match(
            target_extended_cluster, complimentary_extended_clusters
        )
        if y_cluster is not None:
            peptide_type = get_peptide_type(target_cluster, y_cluster, sqllite_database)
            score = get_peptide_score(target_cluster, y_cluster, protein_length)
            peptide = Peptide(
                peptide_type=peptide_type,
                b_extended_cluster=target_extended_cluster,
                y_extended_cluster=y_cluster,
                score=score,
            )
            return peptide
    else:
        b_cluster = get_best_extended_cluster_match(
            target_extended_cluster, complimentary_extended_clusters
        )
        if b_cluster is not None:
            peptide_type = get_peptide_type(b_cluster, target_cluster, sqllite_database)
            score = get_peptide_score(b_cluster, target_cluster, protein_length)
            peptide = Peptide(
                peptide_type=peptide_type,
                b_extended_cluster=b_cluster,
                y_extended_cluster=target_extended_cluster,
                score=score,
            )
            return peptide


def get_peptides(extended_clusters, sqllite_database):
    peptides = []
    b_extended_clusters = [
        extended_cluster
        for extended_cluster in extended_clusters
        if extended_cluster.cluster.ion == "b"
    ]
    y_extended_clusters = [
        extended_cluster
        for extended_cluster in extended_clusters
        if extended_cluster.cluster.ion == "y"
    ]
    for b_extended_cluster in b_extended_clusters:
        peptide = get_peptide(b_extended_cluster, y_extended_clusters, sqllite_database)
        peptides.append(peptide)
    for y_extended_cluster in y_extended_clusters:
        peptide = get_peptide(y_extended_cluster, b_extended_clusters, sqllite_database)
        peptides.append(peptide)
    return peptides


def get_aligned_peptide(peptide):
    if peptide is not None:
        b_extended_cluster = peptide.b_extended_cluster
        y_extended_cluster = peptide.y_extended_cluster
        b_cluster = b_extended_cluster.cluster
        y_cluster = y_extended_cluster.cluster

        hybrid = peptide.peptide_type
        left_protein = b_cluster.protein.id
        right_protein = y_cluster.protein.id
        sequence = (
            b_cluster.longest_kmer
            + b_extended_cluster.extended_sequence
            + y_extended_cluster.extended_sequence
            + y_cluster.longest_kmer
        )
        b_score = b_cluster.score
        y_score = y_cluster.score
        total_score = b_score + y_score
        total_gaussian_score = peptide.score
        extensions = [
            b_extended_cluster.extended_sequence,
            y_extended_cluster.extended_sequence,
        ]
        precursor_mass = b_cluster.longest_kmer.precursor_mass
        precursor_charge = b_cluster.longest_kmer.precursor_charge
        total_mass_error = 0
        total_count = 0
        aligned_peptide = AlignedPeptide(
            hybrid=hybrid,
            left_protein=left_protein,
            right_protein=right_protein,
            sequence=sequence,
            b_score=b_score,
            y_score=y_score,
            total_score=total_score,
            total_gaussian_score=total_gaussian_score,
            extensions=extensions,
            precursor_mass=precursor_mass,
            precursor_charge=precursor_charge,
            total_mass_error=total_mass_error,
            total_count=total_count,
        )
        return aligned_peptide


def get_aligned_peptides(experiment_parameters):
    precursors = experiment_parameters.precursors
    ppm_tolerance = experiment_parameters.ppm_tolerance
    filtered_fragments = get_filtered_fragments(precursors, ppm_tolerance)
    sqllite_database = experiment_parameters.sqllite_database
    number_decimal_places = experiment_parameters.number_decimal_places
    matched_fragments = get_matched_fragments(
        filtered_fragments, sqllite_database, ppm_tolerance, number_decimal_places
    )
    matched_proteins = get_matched_proteins(matched_fragments, sqllite_database)
    clusters = get_clusters(matched_proteins)
    extended_clusters = get_extended_clusters(clusters)
    peptides = get_peptides(extended_clusters, sqllite_database)
    aligned_peptides = []
    for peptide in peptides:
        aligned_peptide = get_aligned_peptide(peptide)
        aligned_peptides.append(aligned_peptide)
    return aligned_peptides
