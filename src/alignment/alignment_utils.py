from objects import Database, Spectrum
from scoring import scoring

import database
import gen_spectra
import utils

import re
import math

#################### Private functions ####################

def __get_surrounding_amino_acids(
    parent_sequence: str, 
    sequence: str, 
    count: int
    ) -> list:
    '''Get the amino acids that surround a sequence. Return the (left, right) 
    *count* number of amino acids

    :param parent_sequence: protein sequence to pull from 
    :type parent_sequence: str
    :param sequence: subsequence we are looking for
    :type sequence: str
    :param count: the number of surrounding amino acids to get per side
    :type count: int

    :returns: tuples of (left amino acids, right amino acids) from all occurances
        in the parent sequence
    :rtype: list
    '''
    # keep track of the pairs
    flanking_pairs = []

    # get all the positions of sequence in the parent sequence
    occurances = [m.start() for m in re.finditer(sequence, parent_sequence)]

    for o in occurances:
        flanking_pairs.append(
            (parent_sequence[max(0, o - count): o], 
            parent_sequence[o + len(sequence): min(len(parent_sequence), o + len(sequence) + count)])
        )

    return flanking_pairs

def __add_amino_acids(
    spectrum: Spectrum, 
    sequence: str, 
    db: Database, 
    gap: int = 3, 
    tolerance: float = 1.0
    ) -> list:
    '''Try and add amino acids to get the closest calculated precursor mass to 
    the observed precursor mass

    :param spectrum: observed spectrum
    :type spectrum: Spectrum
    :param sequence: the alignment to add amino acids to
    :type sequence: str
    :param db: holds to source proteins
    :type db: Database
    :param gap: number of allowed amino acids try and add to each side
        (default is 3)
    :type gap: int
    :param tolerance: the tolerance (in Daltons) to allow when comparing 
        precursor masses. 
        (default is 1)
    :type tolerance: float

    :returns: all perturbed str sequences with a calculated precursor mass that 
        is within the tolerance of the observed precursor mass
    :rtype: list
    '''

    filled_in  = []

    # get the parents of the sequence(s) and add or subtract amino acids
    parents = get_parents(sequence, db, 'b')
    parents += get_parents(sequence, db, 'y')

    # if its hybrid, we should fill it in in all possible ways
    if utils.HYBRID_ALIGNMENT_PATTERN.findall(sequence):

        # go through each set of parents
        for l_p in parents[0]:
            for r_p in parents[1]:
                
                # since its a hybrid, we only want to add to the middle. We know that
                # the most left and most right SHOULD be the best scoring bits, so adding 
                # to the left of the left will only hurt, and adding to the right of the right
                # will only hurt
                left_seq, right_seq = utils.split_hybrid(sequence)
                
                # get the left and right protein sequences
                left_seqs = database.get_entry_by_name(db, l_p)
                right_seqs = database.get_entry_by_name(db, r_p)

                for left_prot in left_seqs:
                    for right_prot in right_seqs:

                        # get the aas to the right of the left sequence
                        left_append = [x[1] for x in __get_surrounding_amino_acids(left_prot.sequence, left_seq, gap)]
                        
                        # get the aas to the left of the right sequence
                        right_prepend = [x[0] for x in __get_surrounding_amino_acids(right_prot.sequence, right_seq, gap)]
                        
                        # slowly add each set of amino acids
                        for to_append in left_append:
                            for to_prepend in right_prepend:
                                
                                # slowly add each
                                for i in range(len(to_append) + 1):
                                    for j in range(len(to_prepend) + 1):
                                        
                                        new_left = left_seq + to_append[:i]
                                        new_right = ('' if j == 0 else to_prepend[-j:]) + right_seq
                                        
                                        # create the new sequence and get the new precursor mass
                                        new_seq = align_overlaps(new_left, new_right)
                                        
                                        new_prec = gen_spectra.get_precursor(
                                            new_seq.replace('(', '').replace(')', '').replace('-', ''), 
                                            spectrum.precursor_charge
                                        )
                                        
                                        # find the precursor distance, and if its close enough, keep it 
                                        pd = scoring.precursor_distance(spectrum.precursor_mass, new_prec)
                                        
                                        if pd <= tolerance:
                                            filled_in.append(new_seq)
        
    # if its nonhybrid, try left and right side
    else:

        for p in parents[0]:

            # get the parent sequence
            entries = database.get_entry_by_name(db, p)

            for entry in entries:

                p_seq = entry.sequence

                # get the flanking amino acid pairs
                for flanking_pair in __get_surrounding_amino_acids(p_seq, sequence, gap):

                    # go through all possible combinations of the flanking pairs
                    for i in range(gap + 1):
                        for j in range(gap - i + 1):

                            # get the new sequence and its precursor mass. Add it to filled in list
                            new_seq = flanking_pair[0][gap-i:] + sequence + flanking_pair[1][:j]

                            new_prec = gen_spectra.get_precursor(new_seq, spectrum.precursor_charge)

                            # get the precursor distance, and if it is close enough, keep it
                            p_d = scoring.precursor_distance(spectrum.precursor_mass, new_prec)

                            if p_d <= tolerance:
                                filled_in.append(new_seq)
    return filled_in

def __remove_amino_acids(
    spectrum: Spectrum, 
    sequence: str, 
    gap: int = 3, 
    tolerance: float = 1
    ) -> list:
    '''Remove up to gap number of amino acids to try and match precursor mass

    :param spectrum: observed spectrum
    :type spectrum: Spectrum
    :param sequence: the alignment to add amino acids to
    :type sequence: str
    :param gap: number of allowed amino acids try and add to each side
        (default is 3)
    :type gap: int
    :param tolerance: the tolerance (in Daltons) to allow when comparing 
        precursor masses. 
        (default is 1)
    :type tolerance: float

    :returns: all perturbed str sequences with a calculated precursor mass that 
        is within the tolerance of the observed precursor mass
    :rtype: list
    '''

    attempted = []
    
    # treat hybrids different than non-hybrids
    if '-' in sequence or '(' in sequence or ')' in sequence:
        
        # get the left and right seperately
        left_seq, right_seq = utils.split_hybrid(sequence)
        
        # since this is a hybrid, we assume that the left is ~correct and the right is ~correct
        # so we only want to remove amino acids from the middle section
        for i in range(gap + 1):
            for j in range(gap - i + 1):
                
                # take left to :-i and right j:
                new_left = left_seq[:-i] if i > 0 else left_seq
                new_right = right_seq[j:]
                
                # create a new hybrid
                new_seq = align_overlaps(new_left, new_right)
                
                # get the new precursor
                new_prec = gen_spectra.get_precursor(
                    new_seq.replace('-', '').replace('(', '').replace(')', ''),
                    spectrum.precursor_charge
                )
                
                # get the new precursor distance
                pd = scoring.precursor_distance(spectrum.precursor_mass, new_prec)
                
                # if the precursor distance is within our tolerance, append it
                if pd <= tolerance:
                    attempted.append(new_seq)

    # otherwise, just take up to gap off from the left and the right
    else:
        for i in range(gap + 1):
            new_seq1 = sequence[i:]
            new_seq2 = sequence[:-i]
            
            # cacluate the new precurosrs and add to attempted if within the tolerance
            new_prec1 = gen_spectra.get_precursor(new_seq1, spectrum.precursor_charge)
            new_prec2 = gen_spectra.get_precursor(new_seq2, spectrum.precursor_charge)
            
            pd1 = scoring.precursor_distance(spectrum.precursor_mass, new_prec1)
            pd2 = scoring.precursor_distance(spectrum.precursor_mass, new_prec2)
            
            if pd1 <= tolerance:
                attempted.append(new_seq1)
            
            if pd2 <= tolerance:
                attempted.append(new_seq2)

    return list(set(attempted))
     

#################### Public functions ####################

def align_overlaps(seq1: str, seq2: str) -> str:
    '''Attempt to align two string sequences. It will look at the right side of 
    seq1 and left side of seq2 to overlap the two strings. If no overlap is 
    found, seq2 is appended to seq1
    
    :param seq1: the left sequence
    :type seq1: str
    :param seq2: the right sequence
    :type seq2: str

    :returns:
    :rtype: str

    :Example: 

    >>> align_overlaps('ABCD', 'CDEF')
    >>> 'ABCDEF'

    :Example:

    >>> align_overlaps('ABCD', 'EFGH')
    >>> 'ABCD-EFGH'
    '''

    alignment = None
    # if we have a perfect overlap, return it
    if seq1 == seq2:
        return seq1
    
    # if the first bit of seq2 == all of seq1 or if the last bit of seq1 == all of seq2, return
    # (ex: ABC, ABCDE -> ABCDE)
    # (ex: ABCDE, CDE -> ABCDE)
    idx_len = min(len(seq1), len(seq2))

    if seq1[:idx_len] == seq2[:idx_len]:
        return seq1 if len(seq1) > len(seq2) else seq2

    if seq1[-idx_len:] == seq2[-idx_len:]:
        return seq1 if len(seq1) > len(seq2) else seq2
    
    # try and find an alignment. seq2 should overlap as much of the right of seq1 as possible
    # get the starting points. 
    # Starting points means we found the first character in seq2 in seq1
    start_points = [i for i in range(len(seq1)) if seq2[0] == seq1[i]]
    for sp in start_points:

        # try and see if extending it makes it match
        # a correct overlap should mean we run out of characters in 
        # seq1 before we hit the end of seq2
        for i in range(sp, len(seq1)):

            # if the next value i in seq1 does not equal the next 
            # characeter i-sp in seq2
            if i-sp < 0 or i-sp > len(seq2)-1 or seq1[i] != seq2[i-sp]:
                i -= 1
                break

        # if i hits the the end of seq1, we have an overlap
        if i == len(seq1) - 1:
            s2_start = len(seq1) - sp
            right_seq = seq2[s2_start:] if s2_start < len(seq2) else ''
            alignment = seq1 + right_seq
            break
  
    # if no overlpa exists, just make append seq2 to seq1
    if alignment is None:
        alignment = seq1 + '-' + seq2

    return alignment

def match_precursor(
    spectrum: Spectrum, 
    sequence: str, 
    db: Database, 
    gap: int = 3, 
    tolerance: float = 1
    ) -> list:
    '''Try and fill in the gaps of an alignment. This is primarily focused on 
    filling in the gaps left by the difference in precursor mass. If we find that
    the difference is more than GAP amino acids, then an empty list is returned

    :param spectrum: observed spectrum
    :type spectrum: Spectrum
    :param sequence: the attempted alignment
    :type sequence: str
    :param db: source of proteins
    :type db: Database
    :param gap: the maximum number of allowed number of amino acids to 
        add/subtract to/from the original sequence. 
        (default is 3)
    :type gap: int
    :param tolerance: the error (in Daltons) allowed when trying to match
        a calculated precursor to the observed precursor mass. 
        (default is 1)
    :type tolerance: float

    :returns: sequences with a precursor within the tolerance of the observed
        precursor
    :rtype: list
    '''

    # remove special characters for generating sequences
    clean_seq = sequence.replace('-', '').replace('(', '').replace(')', '')

    # get the theoretical precursor mass
    theory_precrusor = gen_spectra.get_precursor(clean_seq, spectrum.precursor_charge)

    # determine the number of amino acids we could be off
    estimated_off = abs(utils.predicted_len_precursor(spectrum, clean_seq) - len(clean_seq))

    # if there are too many to add or subtract, return None
    if gap < estimated_off:
        return [None]

    # add amino acids
    if spectrum.precursor_mass > theory_precrusor:

        return __add_amino_acids(spectrum, sequence, db, gap, tolerance)

    # subtract amino acids:
    else:

        return __remove_amino_acids(spectrum, sequence, gap, tolerance)

def get_parents(
    seq: str, 
    db: Database, 
    ion: str = None
    ) -> (list, list):
    ''' Get the parents of a sequence. If the sequence is a hybrid sequence, 
    then the second entry of the tuple holds a list of proteins for the right 
    contributor, otherwise the right entry is empty.

    :param seq: the sequence to look for the parents of
    :type seq: str
    :param db: the source proteins
    :type db: Database
    :param ion: if left as None, look for the full string. Otherwise try 
        to recursivley look by taking off sides of the kmer depending on the ion 
        type. 
        (default is None)
    :type ion: str

    :returns: if hybrid, (left source proteins, right source proteins)
        else, (source proteins, None)
    :rtype: (list, list)

    :Example:

    >>> # non hybrid peptide
    >>> get_parents('ABCDE', db, None)
    >>> ([protein1, protein2], None)

    :Example:

    >>> # hybrid peptide
    >>> get_parents('ABC(DE)FGH', db, None)
    >>> ([protein1, protein2], [protein3])
    '''
    get_sources = lambda s: database.get_proteins_with_subsequence(db, s)
    get_sources_ion = lambda s, i: database.get_proteins_with_subsequence_ion(db, s, i)

    # If the sequence is hybrid, split it to find each parent
    if utils.HYBRID_ALIGNMENT_PATTERN.findall(seq):

        # get the left and right sequnces
        left_seq, right_seq = utils.split_hybrid(seq)
        
        return (get_sources_ion(left_seq, 'b'), get_sources_ion(right_seq, 'y'))

    # otherwise just look up the one for the full sequence
    if ion is not None and ion in 'by':
        return (get_sources_ion(seq, ion), None)

    return (get_sources(seq), None)

def extend_non_hybrid(kmer, spectrum: Spectrum, ion: str, db: Database) -> list:
    seq = kmer[2][2]
    extensions = []
    extension_len = utils.predicted_len_precursor(spectrum, seq) - len(seq)
    if extension_len <= 0:
        return [seq]
    parents, _ = get_parents(seq, db)
    for parent in parents:
        entries = database.get_entry_by_name(db, parent)
        for entry in entries:
            seq_idxes = [m.start() for m in re.finditer(seq, entry.sequence)]
            for seq_idx in seq_idxes:
                if 'y' in ion:
                    min_idx = max(0, seq_idx - extension_len)
                    extensions.append(entry.sequence[min_idx:len(seq) + seq_idx])
                else:
                    max_idx = min(len(entry.sequence), seq_idx + len(seq) + extension_len)
                    extensions.append(entry.sequence[seq_idx:max_idx])

    return extensions