from collections import namedtuple

Database = namedtuple(
    'Database', 
    ['fasta_file', 'proteins', 'kmers'], 
    defaults=['', {}, {}]
)

DatabaseEntry = namedtuple(
    'DatabaseEntry', 
    ['sequence', 'description'],
    defaults=['', '']
)

Spectrum = namedtuple(
    'Spectrum', [
        'num',
        'mz_values',
        'abundance',
        'precursor_mass',
        'precursor_charge',
        'file_name',
        'id', 
        'other_metadata'
    ],
    defaults=[-1,[], [], 0, '', '', '', {}]
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

