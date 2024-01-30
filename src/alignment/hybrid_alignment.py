from lookups.objects import Database, Spectrum
from lookups.utils import all_perms_of_s
from scoring import scoring
from alignment import alignment_utils

import computational_pipeline.database as database

def __replace_ambiguous_hybrid(
    hybrid: tuple, 
    db: Database, 
    observed: Spectrum
    ) -> (str, str):
    '''Attempt to replace a hybrid with a sequence from the database. 

    :param hybrid: tuple of (nonhybrid sequence, hybrid sequence)
    :type hybrid: tuple
    :param db: source proteins
    :type db: Database
    :param observed: observed spectrum
    :type observed: Spectrum

    :returns: input or (nonhybrid, None)
    :rtype: tuple
    '''
    # get the sequence without the hybrid characters -()
    nonhyb = hybrid[0]

    # see if the sequence exists as a non hybrid. If so, return that
    if len(database.get_proteins_with_subsequence(db, nonhyb)):
        return ((nonhyb, None))

    # Try replacing all L with I and vice versa
    possible = all_perms_of_s(nonhyb, 'LI')

    # if we had no other permutations, return the hybrid
    if len(possible) == 0:
        return hybrid

    # try and find a sequnce that could be explained by an LI switch
    for p in possible:

        # if this permutation exists as a nonhybrid, return it
        if len(database.get_proteins_with_subsequence(db, p)):
            return ((p, None))

    # if we didn't find a sequence that could be found in the database, 
    # just return the hybrid input 
    return hybrid

#################### Public functions ####################

def replace_ambiguous_hybrids(
    hybrid_alignments: list, 
    db: Database, 
    observed: Spectrum
    ) -> list:
    '''Remove any ambiguous hybrid alignments that can be explained by non hybrid sequences.
    The returned list has the sequences or their replacements in the same order that 
    they were in on entry. 

    Amino acids L and I are swapped and tried in the search due to the ambiguity 
    in their mass

    :param hybrid_alignments: tuples of attemted hybrid alignments of 
        (non hybrid sequence, hybrid sequence)
    :type hybrid_alignments: list
    :param db: source proteins
    :type db: Database
    :param observed: observed spectrum
    :type observed: Spectrum

    :returns: If no replacements are found, the input. Otherwise a tuple of 
        (non hybrid, None) is inserted in its position
    :rtype: list
    '''

    return [
        __replace_ambiguous_hybrid(hybrid_alignment, db, observed)\
         for hybrid_alignment in hybrid_alignments
    ]

def hybrid_alignment(seq1: str, seq2: str) -> (str, str):
    '''Create a hybrid alignment from 2 sequences. If an overlap between these two sequences
    is found, () are placed around the ambiguous section. If there is no overlap, then 
    seq2 is appended to seq1 with a - at the junction

    :param seq1: left sequence
    :type seq1: str
    :param seq2: right sequence
    :type seq2: str

    :returns: hybrid sequence without special characters, hybrid with special 
        characters
    :rtype: tuple
        
    :Example:

    >>> hybrid_alignment('ABCDE', 'DEFGH')
    >>> ('ABCDEFGH', 'ABC(DE)FGH')

    :Example:

    >>> hybrid_alignment('ABCD', 'EFGH')
    >>> ('ABCDEFGH', 'ABCD-EFGH')
    '''
    alignment = ''
    hybalignment = ''
    attempted_overlap = alignment_utils.align_overlaps(seq1, seq2)

    # If an alignment was made (the shorter length means combined sequences is shorter than
    # just appending them), identify the ambiguous area
    if attempted_overlap is not None and 0 < len(attempted_overlap) < len(seq1) + len(seq2):
        
        # there is an overlap and some ambiguity
        # get the starting point of seq2
        rightstart = attempted_overlap.index(seq2)
        leftend = len(seq1) - 1

        # range between leftend and rigth start is ambiguous
        middle_sec = attempted_overlap[rightstart:leftend + 1]
        alignment = attempted_overlap
        hybalignment = attempted_overlap[:rightstart] \
                        + '(' + middle_sec + ')' \
                        + attempted_overlap[leftend+1:]
    
    # therwise just append right seq to the left with a - to seperate them 
    else:
        alignment = seq1 + seq2
        hybalignment = seq1 + '-' + seq2
        
    return (alignment, hybalignment)