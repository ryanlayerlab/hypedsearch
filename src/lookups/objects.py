from collections import namedtuple

Protein = namedtuple('Protein', ['description','sequence'])
Database = namedtuple('Database', ['fasta_file', 'proteins', 'kmers'])

Spectrum = namedtuple(
    'Spectrum', [
        'num',
        'mz_values',
        'abundance',
        'precursor_mass',
        'precursor_charge',
        'file_name',
        'id', 
        'retention_time',
        'precursor_abundance'
    ],
    defaults=[-1,[], [], 0, '', '', '', '', 0]
)

SequenceAlignment = namedtuple(
    'SequenceAlignment', 
    ['proteins', 'sequence', 'b_score', 'y_score', 'total_score', 'precursor_distance', 'total_mass_error'],
    defaults=[[], '', 0.0, 0.0, 0.0, 100, 100]
)

HybridSequenceAlignment = namedtuple(
    'HybridSequenceAlignment', 
    ['left_proteins', 'right_proteins', 'sequence', 'hybrid_sequence', 
        'b_score', 'y_score', 'total_score', 'precursor_distance', 'total_mass_error'],
    defaults=[[], [], '', '', 0.0, 0.0, 0.0, 100, 100]
)

Alignment_Instrumentation = namedtuple(
    'Alignment_Instrumentation', 
    ['B_and_Y_full_bipartite_alignment', 'average_dataset_size_1',
    'seconds_op_1','removing_ambiguous_hybrids_time','average_dataset_size_2',
    'seconds_op_2','matching_precursor_masses_time','average_dataset_size_3',
    'seconds_op_3','turning_matches_into_objects_time','average_dataset_size_4',
    'seconds_op_4',
    'initial_sequences_with_too_many_or_few_amino_acids_to_try_to_precursor_match',
    'avg_b_score','avg_y_score','avg_total_score'],
    defaults=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
)

Alignments = namedtuple(
    'Alignments', 
    ['spectrum', 'alignments', 'Alignment_Instrumentation'], 
    defaults=[Spectrum([], [], 0, 0, 0.0, ''), [], None]
)

MPSpectrumID = namedtuple(
    'MPSpectrumID', 
    ['b_hits', 'y_hits', 'spectrum', 'ppm_tolerance', 'precursor_tolerance'],
    defaults=[[], [], None, -1, 0]
)

DEVFallOffEntry = namedtuple(
    'DEVFallOffEntry', 
    ['hybrid', 'truth_sequence', 'fall_off_operation', 'meta_data'], 
    defaults=[False, '', '', {}]
)

Identification_Instrumentation = namedtuple(
    'Identification_Instrumentation', 
    ['average_b_scoring_time', 'average_y_scoring_time', 'time_to_filter_out_top_50_kmers',
     'average_extension_time', 'average_non_hybrid_refinement_time', 'average_non_hybrid_scoring_time',
      'average_hybrid_refinement_time', 'average_hybrid_scoring_time', 'average_alignment_time'],
    defaults=[0,0,0,0,0,0,0,0,0]
)

Aligned_Spectrum = namedtuple(
    'Aligned_Spectrum',
    ['hybrid','left_proteins','right_proteins','sequence','b_scores','y_scores','total_score',
     'total_gaussian_score','extensions','precursor_mass','precursor_charge','total_mass_error','total_count'])

AlignedSpectrumsParams = namedtuple(
    'AlignedSpectrumParams',
    ['spectrums', 'sqllite_database', 'max_peptide_length', 'ppm_tolerance', 'precursor_tolerance', 'number_hybrids', 
     'number_natives', 'target_seq'])

TargetDataParams = namedtuple('TargetDataParams', ['target_seq', 'mz_values', 'ppm_tolerance', 'precursor_mass', 'prec_tol', 'precursor_charge'])
#original_target_seq, built_database, spectrum.mz_values, ppm_tolerance, spectrum.precursor_mass, prec_tol, spectrum.precursor_charge
AlignmentParams = namedtuple('AlignmentParams', ['spectrum', 'sqllite_database', 'ppm_tolerance'])
#spectrum, sqllite_database, max_peptide_length, built_database, ppm_tolerance
CheckMatchedMassesParams = namedtuple('CheckMatchedMassesParams', ['matched_masses_b', 'matched_masses_y', 'target_left_pids', 'target_left_indices', 'target_right_pids', 'target_right_indices'])
#matched_masses_b, matched_masses_y, target_left_pids, target_left_indices, target_right_pids, target_right_indices
FindFromPrecursorParams = namedtuple('FindFromPrecursorParams', ['converted_b', 'matched_masses_b', 'spectrum', 'ppm_tolerance', 'proteins_database'])
#converted_b, matched_masses_b, spectrum, ppm_tolerance, built_database.proteins
CreateHitsParams = namedtuple('CreateHitsParams', ['spectrum', 'converted_b', 'converted_y', 'matched_masses_b', 'matched_masses_y', 'num'])
#spectrum, converted_b, converted_y, matched_masses_b, matched_masses_y, spectrum.num
SortedClustersParams = namedtuple('SortedClustersParams', ['built_database', 'converted_b', 'converted_y', 'b_hits', 'y_hits', 'precursor_charge', 'ppm_tolerance', 'num'])
#built_database, converted_b, converted_y, b_hits, y_hits, spectrum.precursor_charge, ppm_tolerance, spectrum.num
CheckSortedClustersParams = namedtuple('CheckSortedClustersParams', ['b_sorted_clusters', 'y_sorted_clusters', 'good_b_entries', 'good_y_entries'])
#b_sorted_clusters, y_sorted_clusters, good_b_entries, good_y_entries
SearchSpaceParams = namedtuple('SearchSpaceParams', ['b_sorted_clusters', 'y_sorted_clusters', 'precursor_charge'])
#b_sorted_clusters, y_sorted_clusters, spectrum.precursor_charge
CheckSearchesParams = namedtuple('CheckSearchesParams', ['b_search_space', 'y_search_space', 'target_left_pids', 'target_right_pids', 'target_left_indices', 'target_right_indices', 'target_seq', 'precursor_charge', 'ppm_tolerance'])
#b_search_space, y_search_space, target_left_pids, target_right_pids, target_left_indices, target_right_indices, target_seq, spectrum.precursor_charge, ppm_tolerance
PairNativesParams = namedtuple('PairNativesParams', ['b_search_space', 'y_search_space', 'precursor_mass', 'prec_tol'])
#b_search_space, y_search_space, spectrum.precursor_mass, prec_tol
PairIndicesParams = namedtuple('PairIndicesParams', ['b_search_space', 'y_search_space', 'precursor_mass', 'prec_tol', 'precursor_charge', 'score_filter'])
#b_search_space, y_search_space, spectrum.precursor_mass, prec_tol, spectrum.precursor_charge, score_filter
CheckMergesParams = namedtuple('CheckMergesParams', ['unique_hybrid_merged_seqs', 'unique_native_merged_seqs', 'good_b_searches', 'good_y_searches'])
#unique_hybrid_merged_seqs, unique_native_merged_seqs, good_b_searches, good_y_searches
RescoreMergesParams = namedtuple('RescoreMergesParams', ['unique_merges', 'spectrum', 'ppm_tolerance'])
#unique_merges, spectrum, ppm_tolerance
CheckRescoredMergesParams = namedtuple('CheckRescoredMergesParams', ['unique_rescored', 'good_merged_seqs'])
#unique_rescored, good_merged_seqs
PostprocessedAlignmentsParams = namedtuple('PostprocessedAlignmentsParams', ['built_database', 'unique_rescored', 'spectrum', 'num_hybrids', 'num_natives'])

TargetData = namedtuple('TargetData', ['target_seq', 'target_left_pids', 'target_right_pids', 'target_left_indices', 'target_right_indices', 'target_score'])
AlignmentData = namedtuple('AlignmentData', ['b_precursor', 'y_precursor', 'matched_masses_b', 'matched_masses_y'])


#converted_b, converted_y, matched_masses_b, matched_masses_y 
GoodEntries = namedtuple('GoodEntries', ['good_b_entries', 'good_y_entries'])
#good_b_entries, good_y_entries
PrecursorHitResult = namedtuple('PrecursorHitResult', ['best_precursor_hit', 'score_filter'])
#best_prec_hit, score_filter 
Hits = namedtuple('Hits', ['b_hits', 'y_hits'])
#b_hits, y_hits
SortedClusters = namedtuple('SortedClusters', ['b_sorted_clusters', 'y_sorted_clusters'])
#b_sorted_clusters, y_sorted_clusters 
SearchSpace = namedtuple('SearchSpace', ['b_search_space', 'y_search_space'])
#b_search_space, y_search_space
GoodSearches = namedtuple('GoodSearches', ['good_b_searches', 'good_y_searches'])
#good_b_searches, good_y_searches

