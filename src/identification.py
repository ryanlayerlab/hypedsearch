from objects import ExperimentParameters, Precursor, KMer, MatchedFragment, Protein, MatchedProtein, Cluster, Peptide, ExtendedCluster
from lookups.constants import AMINO_ACIDS
from itertools import groupby
from operator import itemgetter
from collections import defaultdict

def get_filtered_fragments(precursors, ppm_tolerance):
    filtered_fragments = []
    for precursor in precursors:
        for fragment in precursor.fragments:
            if fragment.abundance >= ppm_tolerance:
                filtered_fragments.append(fragment)
    return filtered_fragments

def ppm_to_da(mass, ppm_tolerance):
    return abs((ppm_tolerance / 1000000)*mass)

def get_matched_fragment(fragment, sqllite_database, ppm_tolerance, number_decimal_places):
    fragment_id = fragment.id
    precursor_mass = fragment.precursor_mass
    precursor_charge = fragment.precursor_charge
    mz_value = fragment.mz_value
    adjusted_tolerance =  ppm_to_da(mz_value, ppm_tolerance) 
    b_rows, y_rows = sqllite_database.query_mass_kmers(fragment_id, precursor_mass, precursor_charge ,mz_value, adjusted_tolerance, number_decimal_places)
    b_kmers = [KMer(*row) for row in b_rows]
    y_kmers = [KMer(*row) for row in y_rows]
    kmers = b_kmers + y_kmers
    kmers = [
        kmer._replace(ion='b' if kmer.ion == 0 else 'y') 
        for kmer in kmers]    
    matched_fragment = MatchedFragment(fragment=fragment, kmers=kmers)
    return matched_fragment    

def get_matched_fragments(fragments,sqllite_database,ppm_tolerance,number_decimal_places):
    matched_fragments = []
    for fragment in fragments:
        matched_fragment = get_matched_fragment(fragment,sqllite_database,ppm_tolerance,number_decimal_places)
        if len(matched_fragment.kmers) > 0:
            matched_fragments.append(matched_fragment)
    return matched_fragments

def get_matched_proteins(matched_fragments,sqllite_database):
    all_kmers = []
    for matched_fragment in matched_fragments:
        all_kmers.extend(matched_fragment.kmers)
    sorted_kmers = sorted(all_kmers, key=lambda k: k.protein_id)
    grouped_kmers = {protein_id: list(group) for protein_id, group in groupby(sorted_kmers, key=lambda k: k.protein_id)}

    matched_proteins = []
    for protein_id, kmers in grouped_kmers.items(): 
        protein_row = sqllite_database.get_protein(protein_id)
        protein = Protein(*protein_row)
        matched_protein = MatchedProtein(protein=protein, kmers=kmers)
        matched_proteins.append(matched_protein)

    return matched_proteins

def get_clusters(matched_proteins):
    clusters = []
    for matched_protein in matched_proteins:
        protein = matched_protein.protein
        kmers = matched_protein.kmers
        longest_kmer = max(kmers, key=lambda k: k.location_end - k.location_start)
        ion_of_longest = longest_kmer.ion
        score = sum(1 for kmer in kmers if kmer.ion == ion_of_longest)
        cluster = Cluster(protein=protein, ion=ion_of_longest, longest_kmer=longest_kmer, score=score)
        clusters.append(cluster)    
    return clusters

def get_extended_clusters(clusters):
    extended_clusters = []
    bad_chars = ['B', 'X', 'U', 'Z', 'O', 'J']
    for cluster in clusters:
        protein = cluster.protein
        longest_kmer = cluster.longest_kmer
        precursor_mass = longest_kmer.precursor_mass
        cumulative_mass = longest_kmer.kmer_mass
        extended_sequence = ''
        if longest_kmer.ion == 'b':
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
        elif longest_kmer.ion == 'y':
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
        extended_cluster = ExtendedCluster(cluster=cluster, extended_sequence=extended_sequence)
    
    extended_clusters.append(extended_cluster)

    return extended_clusters

def get_best_cluster_match(target_cluster, complimentary_clusters):
    target_precursor_maas = target_cluster.longest_kmer.precursor_mass
    target_weight = target_cluster.longest_kmer.kmer_mass
    best_cluster = None
    max_peptide_weight = 0

    for cluster in complimentary_clusters:
        compliment_weight = cluster.longest_kmer.kmer_mass
        peptide_weight = target_weight + compliment_weight
        if peptide_weight <= target_precursor_maas and peptide_weight > max_peptide_weight:
            best_cluster = cluster
            max_peptide_weight = peptide_weight

    return best_cluster

def get_peptide_type(b_cluster,y_cluster, sqllite_database):
    if b_cluster is None or y_cluster is None:
        return 'u'
    if b_cluster.protein.id != y_cluster.protein.id:
        return 'h'
    else:
        return 'n'

def get_peptide(target_cluster, complimentary_clusters,sqllite_database):    
    if target_cluster.cluster_type == 'b':
        y_cluster = get_best_cluster_match(target_cluster, complimentary_clusters)
        peptide_type = get_peptide_type(target_cluster,y_cluster,sqllite_database)
        peptide = Peptide(peptide_type = peptide_type, b_cluster = target_cluster, y_cluster=y_cluster)
        return peptide
    else:
        b_cluster = get_best_cluster_match(target_cluster, complimentary_clusters)
        peptide_type = get_peptide_type(b_cluster,target_cluster,sqllite_database)
        peptide = Peptide(peptide_type = peptide_type, b_cluster = b_cluster, y_cluster=target_cluster)
        return peptide

def get_peptides(clusters,sqllite_database):
    peptides = []
    b_clusters = [cluster for cluster in clusters if cluster.cluster_type == 'b']
    y_clusters = [cluster for cluster in clusters if cluster.cluster_type == 'y']
    for b_cluster in b_clusters:
        peptide = get_peptide(b_cluster, y_clusters,sqllite_database)
        peptides.append(peptide)
    for y_cluster in y_clusters:
        peptide = get_peptide(y_cluster,b_clusters,sqllite_database)
        peptides.append(peptide)
    return peptides        

#https://github.com/ryanlayerlab/hypedsearch/blob/48f4818c9e0b21283781c369311e63eb1e65c994/src/scoring/scoring.py#L61
def get_rescored_peptide(peptide,sqllite_database):

    return None

def get_rescored_peptides(peptides,sqllite_database):
    rescored_peptides = []
    for peptide in peptides:
        rescored_peptide = get_rescored_peptide(peptide,sqllite_database)
        rescored_peptides.append(rescored_peptide)
    return rescored_peptides

def get_aligned_peptide(rescored_peptide):
    if rescored_peptide is None:
        return None
    else:
        peptide_type = 'hybrid'
        b_clusters = rescored_peptide.b_clusters if rescored_peptide.b_clusters is not None else []
        y_clusters = rescored_peptide.y_clusters if rescored_peptide.y_clusters is not None else []
        b_scores = [len(cluster) for cluster in b_clusters]
        y_scores = [len(cluster) for cluster in y_clusters]
        kmer = rescored_peptide.b_clusters[0].longest_kmer
        total_score = sum(b_scores) + sum(y_scores)
        total_gaussian_score = total_score
        left_proteins = ["Protein1"]
        right_proteins = ["Protein2"]
        sequence = "PEPTIDESEQ"
        extensions = []
        precursor_mass = kmer.precursor_mass  
        precursor_charge = kmer.precursor_charge
        total_mass_error = 0.0 
        total_count = len(rescored_peptide.b_clusters) + len(rescored_peptide.y_clusters)    
        return None

def construct_aligned_peptides(rescored_peptides):
    aligned_peptides = []
    for rescored_peptide in rescored_peptides:
        aligned_spectrum = get_aligned_peptide(rescored_peptide)
        aligned_peptides.append(aligned_spectrum)
    return aligned_peptides

def create_aligned_peptides(experiment_parameters):
    precursors = experiment_parameters.precursors
    ppm_tolerance = experiment_parameters.ppm_tolerance    
    filtered_fragments = get_filtered_fragments(precursors, ppm_tolerance)
    sqllite_database = experiment_parameters.sqllite_database
    number_decimal_places = experiment_parameters.number_decimal_places
    matched_fragments = get_matched_fragments(filtered_fragments,sqllite_database,ppm_tolerance, number_decimal_places)
    matched_proteins = get_matched_proteins(matched_fragments,sqllite_database)
    clusters = get_clusters(matched_proteins)
    extended_clusters = get_extended_clusters(clusters)
    native_peptides = get_native_peptides(clusters)
    hybrid_peptides = get_hybrid_peptides(cluster)
    print(extended_clusters[0])
    return None
    # rescored_peptides = get_rescored_peptides(peptides,sqllite_database)
    # aligned_peptides = construct_aligned_peptides(rescored_peptides)
    # return aligned_peptides

def create_aligned_peptides_with_target(experiment_parameters):
    return None

def get_aligned_peptides(experiment_parameters):
    target_seq = experiment_parameters.target_seq
    if len(target_seq) > 0:
        aligned_peptides = create_aligned_peptides_with_target(experiment_parameters)
        return aligned_peptides
    else:
        aligned_peptides = create_aligned_peptides(experiment_parameters)
        return aligned_peptides