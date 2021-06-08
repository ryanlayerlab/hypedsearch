from collections import namedtuple

Database = namedtuple(
    'Database', 
    ['fasta_file', 'proteins', 'kmers'], 
    defaults=['', {}, {}]
)
'''Holds proteins, fasta file, protein tree, and kmer masses 

:ivar fasta_file: The name of the input fasta file
:type fasta_file: str
:ivar proteins: A dictionary of proteins where keys are the entry
    name and the value is a list of DatabaseEntry objects
:type proteins: dict
:ivar kmers: A dictionary mapping kmers to a list of source protein names
:type kmers: dict
'''

DatabaseEntry = namedtuple(
    'DatabaseEntry', 
    ['sequence', 'description'],
    defaults=['', '']
)
'''Contains protein information

:ivar sequence: the full protein sequence
:type sequence: str
:ivar description: the name of the protein
:type description: str
'''

Spectrum = namedtuple(
    'Spectrum', [
        'spectrum',
        'abundance',
        'total_intensity',
        'ms_level',
        'scan_number',
        'precursor_mass',
        'precursor_charge',
        'file_name',
        'id', 
        'other_metadata'
    ],
    defaults=[[], [], 0, 0, -1, 0, '', '', '', {}]
)
'''Holds information regarding an MS or MS/MS spectrum

:ivar spectrum: m/z float values of an MS run
:type spectrum: list
:ivar abundance: floats describing the abundance of each peak value. Index *i*
    is the abundance of the m/z value at index *i* of the spectrum
:type abundance: list
:ivar ms_level: MS experiment level
:type ms_level: int
:ivar scan_number: scan number of spectrum in the MS run
:type scan_number: int
:ivar precursor_mass: precursor mass of the MS run (i.e. the mass of the whole sequence)
:type precursor_mass: float
:ivar precursor_charge: charge of the precursor mass
:type precursor_charge: int
:ivar file_name: name of the source file of the spectrum
:type file_name: str
:ivar other_metadata: other metadata associated with the spectrum not in the above
:type other_metadata: dict
'''

SequenceAlignment = namedtuple(
    'SequenceAlignment', 
    ['proteins', 'sequence', 'b_score', 'y_score', 'total_score', 'precursor_distance', 'total_mass_error'],
    defaults=[[], '', 0.0, 0.0, 0.0, 100, 100]
)
'''Alignment information for a non-hybrid sequence alignment

:ivar proteins: proteins where the aligned sequence is found
:type proteins: list
:ivar sequence: the string of amino acids that were found as the alignment
:type sequence: str
:ivar b_score: b ion score of the sequence
:type b_score: float
:ivar y_score: y ion score of the sequence
:type y_score: float
:ivar total_score: the score given to the sequence
:type total_score: float
:ivar precursor_distance: the absolute value of the difference between the observed precursor
    mass and the calculated precursor mass of the aligned sequence
:type precursor_distance: float
:ivar total_mass_error: the sum of the absolute values of the error between an aligned
    amino acid mass and the matched observed mass
:type total_mass_error: float
'''

HybridSequenceAlignment = namedtuple(
    'HybridSequenceAlignment', 
    ['left_proteins', 'right_proteins', 'sequence', 'hybrid_sequence', 
        'b_score', 'y_score', 'total_score', 'precursor_distance', 'total_mass_error'],
    defaults=[[], [], '', '', 0.0, 0.0, 0.0, 100, 100]
)
'''Alignment information for a non-hybrid sequence alignment

:ivar left_proteins: proteins that contain the sequence of amino acids that contribute
    to the left side of the hybrid peptide
:type left_proteins: list
:ivar right_proteins: proteins that contain the sequence of amino acids that contribute 
    to the right side of the hybrid peptide
:type right_proteins: list
:ivar sequence: the string of amino acids that were found as the alignment
:type sequence: str
:ivar hybrid_sequence: the string of amino acids that were found as the alignment
    with special characters [(), -] where - denotes a hybrid sequence with no 
    overlap (left-right) and () denotes a hybrid with an overlap (left(overlap)right)
:type hybrid_sequence: str
:ivar b_score: b ion score of the sequence
:type b_score: float
:ivar y_score: y ion score of the sequence
:type y_score: float
:ivar total_score: the score given to the sequence
:type total_score: float
:ivar precursor_distance: the absolute value of the difference between the observed precursor
    mass and the calculated precursor mass of the aligned sequence
:type precursor_distance: float
:ivar total_mass_error: the sum of the absolute values of the error between an aligned
    amino acid mass and the matched observed mass
:type total_mass_error: float
'''

Alignments = namedtuple(
    'Alignments', 
    ['spectrum', 'alignments'], 
    defaults=[Spectrum([], [], 0, 0, 0.0, ''), []]
)
'''Contains the spectrum with SequenceAlignments and HybridSequenceAlignments

:ivar spectrum: the observed spectrum
:type spectrum: Spectrum
:ivar alignments: SequenceAlignment and HybridSequenceAlignment objects
:type alignments: list
'''

MPSpectrumID = namedtuple(
    'MPSpectrumID', 
    ['b_hits', 'y_hits', 'spectrum', 'ppm_tolerance', 'precursor_tolerance', 
        'n', 'digest_type'],
    defaults=[[], [], None, -1, 0, 0, '']
)
'''Holds information to pass to processes during multiprocessing (MP)

:ivar b_hits: k-mers found from the b ion search
:type b_hits: list
:ivar y_hits: k-mers found from the y ion search
:type y_hits: list
:ivar spectrum: observed spectrum
:type spectrum: Spectrum
:ivar ppm_tolerance: parts per million error allowed when matching masses
:type ppm_tolerance: int
:ivar precursor_tolerance: parts per million error allowed when matching precursor mass
:type precursor_tolerance: int
:ivar n: the number of aligments to keep 
:type n: int
:ivar digest_type: the digest performed on the sample
:type digest_type: str
'''

DEVFallOffEntry = namedtuple(
    'DEVFallOffEntry', 
    ['hybrid', 'truth_sequence', 'fall_off_operation', 'meta_data'], 
    defaults=[False, '', '', {}]
)
'''DEVELOPMENT USE ONLY

Holds data about when the components that make up the desired overlapping sequence
falls off and can no longer make the correct alignment

:ivar hybrid: whether or not the desired alignment is a hybrid
:type hybrid: bool
:ivar truth_sequence: the desired string alignment
:type truth_sequence: str
:ivar fall_off_operation: which operation the sequence was no longer attainable
:type fall_off_operation: str
:ivar meta_data: any extra information pertaining to the operation
:type meta_data: dict
'''