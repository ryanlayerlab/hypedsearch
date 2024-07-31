from lookups.objects import Database, Spectrum
from scoring import scoring
import computational_pipeline.gen_spectra
import lookups.utils
import re
import math

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
    filled_in  = []
    parents = get_parents(sequence, db, 'b')
    parents += get_parents(sequence, db, 'y')
    if lookups.utils.HYBRID_ALIGNMENT_PATTERN.findall(sequence):
        for l_p in parents[0]:
            for r_p in parents[1]:
                left_seq, right_seq = lookups.utils.split_hybrid(sequence)
                left_seqs = preprocessing.database_generator.get_entry_by_name(db, l_p)
                right_seqs = preprocessing.database_generator.get_entry_by_name(db, r_p)

                for left_prot in left_seqs:
                    for right_prot in right_seqs:
                        left_append = [x[1] for x in __get_surrounding_amino_acids(left_prot.sequence, left_seq, gap)]
                        right_prepend = [x[0] for x in __get_surrounding_amino_acids(right_prot.sequence, right_seq, gap)]
                        for to_append in left_append:
                            for to_prepend in right_prepend:
                                for i in range(len(to_append) + 1):
                                    for j in range(len(to_prepend) + 1):
                                        new_left = left_seq + to_append[:i]
                                        new_right = ('' if j == 0 else to_prepend[-j:]) + right_seq
                                        new_seq = align_overlaps(new_left, new_right)
                                        
                                        new_prec = computational_pipeline.gen_spectra.get_precursor(
                                            new_seq.replace('(', '').replace(')', '').replace('-', ''), 
                                            spectrum.precursor_charge
                                        )
                                        pd = scoring.precursor_distance(spectrum.precursor_mass, new_prec)
                                        
                                        if pd <= tolerance:
                                            filled_in.append(new_seq)
        
    else:

        for p in parents[0]:
            entries = preprocessing.database_generator.get_entry_by_name(db, p)
            for entry in entries:
                p_seq = entry.sequence
                for flanking_pair in __get_surrounding_amino_acids(p_seq, sequence, gap):
                    for i in range(gap + 1):
                        for j in range(gap - i + 1):
                            new_seq = flanking_pair[0][gap-i:] + sequence + flanking_pair[1][:j]
                            new_prec = computational_pipeline.gen_spectra.get_precursor(new_seq, spectrum.precursor_charge)
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

    attempted = []
    if '-' in sequence or '(' in sequence or ')' in sequence:
        left_seq, right_seq = lookups.utils.split_hybrid(sequence)
        for i in range(gap + 1):
            for j in range(gap - i + 1):
                new_left = left_seq[:-i] if i > 0 else left_seq
                new_right = right_seq[j:]
                new_seq = align_overlaps(new_left, new_right)
                new_prec = computational_pipeline.gen_spectra.get_precursor(
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
            new_prec1 = computational_pipeline.gen_spectra.get_precursor(new_seq1, spectrum.precursor_charge)
            new_prec2 = computational_pipeline.gen_spectra.get_precursor(new_seq2, spectrum.precursor_charge)
            
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

    clean_seq = sequence.replace('-', '').replace('(', '').replace(')', '')
    theory_precrusor = computational_pipeline.gen_spectra.get_precursor(clean_seq, spectrum.precursor_charge)
    estimated_off = abs(lookups.utils.predicted_len_precursor(spectrum, clean_seq) - len(clean_seq))

    if gap < estimated_off:
        return [None]

    if spectrum.precursor_mass > theory_precrusor:

        return __add_amino_acids(spectrum, sequence, db, gap, tolerance)
    else:

        return __remove_amino_acids(spectrum, sequence, gap, tolerance)

def get_parents(
    seq: str, 
    db: Database, 
    ion: str = None
    ):
    get_sources = lambda s: preprocessing.database_generator.get_proteins_with_subsequence(db, s)
    get_sources_ion = lambda s, i: preprocessing.database_generator.get_proteins_with_subsequence_ion(db, s, i)

    if lookups.utils.HYBRID_ALIGNMENT_PATTERN.findall(seq):
        left_seq, right_seq = lookups.utils.split_hybrid(seq)
        return (get_sources_ion(left_seq, 'b'), get_sources_ion(right_seq, 'y'))

    if ion is not None and ion in 'by':
        return (get_sources_ion(seq, ion), None)

    return (get_sources(seq), None)

def extend_non_hybrid(kmer, spectrum: Spectrum, ion: str, database: Database) -> list:
    seq = kmer
    extensions = []
    extension_len = lookups.utils.predicted_len_precursor(spectrum, seq) - len(seq)
    if extension_len <= 0:
        return [seq]
    parents, _ = get_parents(seq, database)
    for parent in parents:
        entries = database.get_entry_by_name(database, parent)
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