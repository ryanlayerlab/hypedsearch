from objects import ExperimentParameters, Precursor, KMer, MatchedFragment, Protein, MatchedProtein, Cluster, Peptide
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
    matched_fragment = MatchedFragment(fragment=fragment, b_kmers=b_kmers, y_kmers=y_kmers)
    return matched_fragment    

def get_native_matched_fragments(fragments,sqllite_database,ppm_tolerance,number_decimal_places):
    matched_fragments = []
    for fragment in fragments:
        matched_fragment = get_matched_fragment(fragment,sqllite_database,ppm_tolerance,number_decimal_places)
        if len(matched_fragment.b_kmers) + len(matched_fragment.y_kmers) > 0:
            matched_fragments.append(matched_fragment)
    return matched_fragments

def calculate_mass(sequence):
    return sum(AMINO_ACIDS.get(aa, 0) for aa in sequence)

def find_precursor_by_id(precursor_list, target_id):
    for precursor in precursor_list:
        if precursor.id == target_id:
            return precursor
    return None

def get_synthetic_b_kmer(matched_precursor, kmer, protein_sequence):
    bad_chars = ['B', 'X', 'U', 'Z', 'O', 'J']
    synthetic_kmers = []
    kmer_start = kmer.location_start
    index = kmer.location_end
    precursor_mass = matched_precursor.mass
    i = 0

    while index + i < len(protein_sequence):
        current_char = protein_sequence[index + i]
        
        if current_char in bad_chars:
            i += 1
            continue
        
        new_fragment = protein_sequence[kmer_start:index + 1]
        new_mass = calculate_mass(new_fragment)
        
        if new_mass <= precursor_mass:
            new_kmer = KMer(
                fragment_id = kmer.fragment_id,
                precursor_mass = kmer.precursor_mass,
                precursor_charge = kmer.precursor_charge,
                protein_id=kmer.protein_id,
                kmer_mass=new_mass,
                location_start=kmer_start,
                location_end=index,
                ion=kmer.ion,
                charge=kmer.charge,
                subsequence=new_fragment,
                kmer_type='S'
            )
            synthetic_kmers.append(new_kmer)
        
        if new_mass >= precursor_mass:
            break
            
        i += 1

    if len(synthetic_kmers) > 0:
        return synthetic_kmers[-1]
    else:
        return None

def get_synthetic_y_kmer(matched_precursor, kmer, protein_sequence):
    bad_chars = ['B', 'X', 'U', 'Z', 'O', 'J']
    synthetic_kmers = []
    index = kmer.location_end
    kmer_start = kmer.location_start
    precursor_mass = matched_precursor.mass
    i = 0

    while index - i > 0:
        current_char = protein_sequence[index - i]
        
        if current_char in bad_chars:
            i -= 1
            continue
        
        new_fragment = protein_sequence[kmer_start:index]
        new_mass = calculate_mass(new_fragment)
            
        if new_mass <= precursor_mass:
            new_kmer = KMer(
                fragment_id = kmer.fragment_id,
                precursor_mass = kmer.precursor_mass,
                precursor_charge = kmer.precursor_charge,
                protein_id=kmer.protein_id,
                kmer_mass=new_mass,
                location_start=kmer_start,
                location_end=index,
                ion=kmer.ion,
                charge=kmer.charge,
                subsequence=new_fragment,
                kmer_type='S'
            )
            synthetic_kmers.append(new_kmer)
        
        if new_mass >= precursor_mass:
            break
            
        i += 1

    if len(synthetic_kmers) > 0:
        return synthetic_kmers[-1]
    else:
        return None

def get_protein_sequence(protein_id,sqllite_database):
    protein_row = sqllite_database.get_protein(protein_id)
    protein = Protein(*protein_row)
    return protein.sequence

def add_synthetic_kmers_to_fragment(matched_fragment,precursors,sqllite_database):
    precursor_id = matched_fragment.fragment.precursor_id
    precursor = find_precursor_by_id(precursors, precursor_id)
    synthetic_b_kmers = []
    synthetic_y_kmers = []
    for b_kmer in matched_fragment.b_kmers:
        protein_id = b_kmer.protein_id
        protein_sequence = get_protein_sequence(protein_id,sqllite_database)
        synthetic_b_kmer =  get_synthetic_b_kmer(precursor, b_kmer, protein_sequence)
        synthetic_b_kmers.append(synthetic_b_kmer)
    for y_kmer in matched_fragment.y_kmers:
        protein_id = y_kmer.protein_id
        protein_sequence = get_protein_sequence(protein_id,sqllite_database)
        synthetic_y_kmer = get_synthetic_y_kmer(precursor, y_kmer, protein_sequence)
        synthetic_y_kmers.append(synthetic_y_kmer)
    matched_fragment.b_kmers.extend(synthetic_b_kmers)
    matched_fragment.y_kmers.extend(synthetic_y_kmers)
    return matched_fragment

def get_all_matched_fragments(native_matched_fragments, precursors, sqllite_database):
    matched_fragments = []
    for native_matched_fragment in native_matched_fragments:
        matched_fragment = add_synthetic_kmers_to_fragment(native_matched_fragment,precursors,sqllite_database)
        matched_fragments.append(matched_fragment)
    return matched_fragments

def get_matched_protein(protein_id, group, sqllite_database):
    protein_row = sqllite_database.get_protein(protein_id)
    protein = Protein(*protein_row)    
    for fragment_id, mf in enumerate(group):
        matched_protein = MatchedProtein(
            protein=protein,
            b_kmers=mf.b_kmers,
            y_kmers=mf.y_kmers
        )
    return matched_protein

def get_protein_id(matched_fragments):
    if matched_fragments.b_kmers and matched_fragments.b_kmers[0]:
        return matched_fragments.b_kmers[0].protein_id
    elif matched_fragments.y_kmers and matched_fragments.y_kmers[0]:
        return matched_fragments.y_kmers[0].protein_id
    return None

def get_matched_proteins(matched_fragments,sqllite_database):
    matched_fragments.sort(key=get_protein_id) 
    grouped_fragments = groupby(matched_fragments, key=get_protein_id)
    matched_proteins = []
    for protein_id, group in grouped_fragments:
        group = list(group)
        matched_protein = get_matched_protein(protein_id, group,sqllite_database)
        matched_proteins.append(matched_protein)    
    return matched_proteins

def get_clusters(matched_proteins):
    clusters = []
    
    for matched_protein in matched_proteins:
        protein = matched_protein.protein
        kmer_groups = defaultdict(list)
        for kmer in filter(None, matched_protein.b_kmers):
            kmer_groups[('b', protein, kmer.location_start)].append(kmer)
        
        for kmer in filter(None, matched_protein.y_kmers):
            kmer_groups[('y', protein, kmer.location_end)].append(kmer)

        for (cluster_type, protein, start), kmers in kmer_groups.items():
            longest_kmer = max(kmers, key=lambda k: k.location_end - k.location_start)
            ending_index = longest_kmer.location_end
            score = len(kmers)
            cluster = Cluster(
                protein=protein,
                cluster_type=cluster_type,
                starting_index=start,
                ending_index=ending_index,
                longest_kmer=longest_kmer, 
                score=score  
            )
            clusters.append(cluster)
    
    return clusters

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
        peptide = Peptide(peptide_type = peptide_type, b_clusters = [target_cluster], y_clusters=[y_cluster])
        return peptide
    else:
        b_cluster = get_best_cluster_match(target_cluster, complimentary_clusters)
        peptide_type = get_peptide_type(b_cluster,target_cluster,sqllite_database)
        peptide = Peptide(peptide_type = peptide_type, b_clusters = [b_cluster], y_clusters=[target_cluster])
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
    native_matched_fragments = get_native_matched_fragments(filtered_fragments,sqllite_database,ppm_tolerance, number_decimal_places)
    matched_fragments = get_all_matched_fragments(native_matched_fragments, precursors, sqllite_database)
    matched_proteins = get_matched_proteins(matched_fragments,sqllite_database)
    clusters = get_clusters(matched_proteins)
    peptides = get_peptides(clusters,sqllite_database)
    rescored_peptides = get_rescored_peptides(peptides,sqllite_database)
    aligned_peptides = construct_aligned_peptides(rescored_peptides)
    return aligned_peptides

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