from scoring import scoring
from objects import Alignment_Instrumentation, Spectrum, SequenceAlignment, HybridSequenceAlignment, Database, Alignments, DEVFallOffEntry
from alignment import alignment_utils, hybrid_alignment
import objects
import utils
import database
import gen_spectra

import math
import re

import time

FIRST_ALIGN_TIME = 0
AMBIGUOUS_REMOVAL_TIME = 0
PRECURSOR_MASS_TIME = 0
OBJECTIFY_TIME = 0

FIRST_ALIGN_COUNT = 0
AMBIGUOUS_REMOVAL_COUNT = 0
PRECURSOR_MASS_COUNT = 0
OBJECTIFY_COUNT = 0

OUT_OF_RANGE_SEQS = 0
TOTAL_ITERATIONS = 0

global extension_times
extension_times = []
global initial_alignment_times
initial_alignment_times = []
global Non_hybrid_refine_time
Non_hybrid_refine_time = []
global non_hybrid_scoring_times
non_hybrid_scoring_times = []
global Hybrid_refine_times
Hybrid_refine_times = []
global hybrid_scoring_times
hybrid_scoring_times = []

def same_protein_alignment(
    seq1: str, 
    seq2: str, 
    parent_sequence: str
    ) -> (str, str):
    '''Attempt to create a non-hybrid alignment from two sequences from the same 
    protein. If the two sequences do not directly overlap but are close enough 
    and from the same protein, make the alignment. If not, create a hybrid 
    alignment from the two input halves. If one compeletely overlaps the other, 
    use the larger sequence as the alignment.
    
    :param seq1: left sequence
    :type seq1: str
    :param seq2: right sequence
    :type seq2: str
    :param parent_sequence: parent sequence of seq1 and seq2
    :type parent_sequence: str

    :returns: if hybrid sequence (sequence without special charcters, sequence
        with hybrid sequence)
        else (sequence, None)
    :rtype: (str, str or None)

    :Example:
    
    >>> same_protein_alignment('ABC', 'CDE', 'ABCDEFG')
    >>> ('ABCDE', None)

    :Example:
    
    >>> same_protein_alignment('ABC', 'FGH', 'ABCDEFHI')
    >>> ('ABCFGH', 'ABC-FGH')
    ''' 

    # check to see if they are equal or one covers the entirety of the other
    if seq1 == seq2:
        return (seq1, None)
    
    if seq1 in seq2:
        return (seq2, None)
    
    if seq2 in seq1: 
        return (seq1, None)
    
    # get the number of gap amino acids allowed
    gap_aa = max(1, len(parent_sequence) // 100)
    
    # get the positions of the left sequence from the protein
    left_start = [m.start() for m in re.finditer(seq1, parent_sequence)]
    
    # get the positions of the right sequence from the protein
    right_start = [m.start() for m in re.finditer(seq2, parent_sequence)]
    
    # if EVERY position of seq1 is to the RIGHT of ALL positions of seq2, make a hybrid 
    if all([r < l for r in right_start for l in left_start]):
        return hybrid_alignment.hybrid_alignment(seq1, seq2)
    
    # check to see if any of the points are within the gap_aa limit
    nonhybrid_alignments = []

    for l in left_start:
        for r in right_start:
            
            # if the right is to the left of left, continue
            if r < l:
                continue
                
            # if the right start - left start plus the length of the subsequence 
            # is less than the gap, just take the starting position of l and ending position 
            # of seq2 to make the full alignment
            if r - (l + len(seq1)) <= gap_aa:
                overlapped = parent_sequence[l: r + len(seq2)]
                nonhybrid_alignments.append(overlapped)
        
    # if no nonhybrids could be made, return a hybrid alignment
    if len(nonhybrid_alignments) == 0:
        return hybrid_alignment.hybrid_alignment(seq1, seq2)
    
    # we have at least one good one. Return the shortest one
    nonhybrid_alignments.sort(key=lambda x: len(x))
    return (nonhybrid_alignments[0], None)

def align_b_y(
    b_kmers: list, 
    y_kmers: list, 
    spectrum: Spectrum, 
    db: Database
    ) -> list:
    '''Try and connect all b and y k-mers and try and make either hybrid 
    or non hybrid string alignments from them.

    :param b_kmers: kmers from b ion masses
    :type b_kmers: list
    :param y_kmers: kmers from y ion masses
    :type y_kmers: list
    :param spectrum: observed spectrum
    :type spectrum: Spectrum
    :param db: source proteins
    :type db: Database

    :results: tuples of alignments. If hybrid, (sequence, sequence with 
        special hybrid characters), otherwise (sequence, None)
    :rtype: list
    '''

    # try and create an alignment from each b and y ion sequence pairs
    spec_alignments = []

    for b_seq in b_kmers:

        # get all the b proteins
        b_proteins = database.get_proteins_with_subsequence(db, b_seq)

        for y_seq in y_kmers:

            # ge the y proteins
            y_proteins = database.get_proteins_with_subsequence(db, y_seq)
            
            # the sequence is from the same protein, try and overlap it
            if any([x in y_proteins for x in b_proteins]):

                # get each protein they have in common
                shared_prots = [x for x in y_proteins if x in b_proteins]
                
                # try each of them 
                for sp in shared_prots:

                    # get the sequence from the entry for alignment
                    prot_seqs = database.get_entry_by_name(db, sp)

                    for prot_entry in prot_seqs:

                        # append any alignments made from these 2 sequences
                        spec_alignments.append(
                            same_protein_alignment(b_seq, y_seq, prot_entry.sequence)
                        )
                
                # try just a dumb hybrid too to make sure
                spec_alignments.append((f'{b_seq}{y_seq}', f'{b_seq}-{y_seq}'))

            # otherwise try hybrid alignment
            else: 
                spec_alignments.append(hybrid_alignment.hybrid_alignment(b_seq, y_seq))
        
    # remove repeats
    return list(set([x for x in spec_alignments if x is not None]))

def extend_base_kmers(
    b_kmers: list, 
    y_kmers: list, 
    spectrum: Spectrum, 
    db: Database
    ) -> (list,list):
    '''Extend all the base b and y ion matched k-mers to the predicted length
    to try and find a non-hybrid alignment

    :param b_kmers: kmers from b ion masses
    :type b_kmers: list
    :param y_kmers: kmers from y ion masses
    :type y_kmers: list
    :param spectrum: observed spectrum
    :type spectrum: Spectrum
    :param db: source proteins
    :type db: Database

    :results: extended ion kmers (strings)
    :rtype: list
    '''
    # try and create an alignment from each extended b and y ion sequence
    extended_b = []
    extended_y = []
    #[item for sublist in t for item in sublist]
    for seq in b_kmers:
        extended_b += [x for x in alignment_utils.extend_non_hybrid(seq, spectrum, 'b', db)]
    for seq in y_kmers:
        extended_y += [x for x in alignment_utils.extend_non_hybrid(seq, spectrum, 'y', db)]

    # file1 = open("metadata.txt", "a")
    # file1.write("b_kmers: \n")
    # file1.writelines([(L + '\n') for L in b_kmers])
    # file1.write("extended_b_kmers: \n")
    # file1.write("extended_b_size: ")
    # file1.write(str(len(extended_b)) + '\n')
    # for seq in extended_b:
    #     file1.write(seq + '\n')
    # file1.write("y_kmers: \n")
    # file1.writelines([(L + '\n') for L in y_kmers])
    # file1.write("extended_y_kmers: \n")
    # file1.write("extended_y_size: ")
    # file1.write(str(len(extended_y)) + '\n')
    # for seq in extended_y:
    #     file1.write(seq + '\n')

    # spec_align = extended_y + extended_b
    # #Checking for duplicates
    # duplicate_list = []
    # for seq in extended_b:
    #     for y_seq in extended_y:F
    #         if seq == y_seq:
    #             duplicate_list += seq

    # file1.write("There are " + str(len(duplicate_list)) + " duplicates out of " + str(len(spec_align)) + " spectra \n")
    # file1.close()

    return extended_b, extended_y

def refine_alignments(
    spectrum: Spectrum, 
    db: Database, 
    alignments: list, 
    precursor_tolerance: int = 10,
    DEV: bool = False, 
    truth: dict = None, 
    fall_off: dict = None
) -> list:
    '''Refine the rough alignmnets made. This includes precursor matching and 
    ambiguous hybrid removals/replacements

    :param spectrum: observed spectrum in question
    :type spectrum: Spectrum
    :param db: Holds all the source sequences
    :type db: Database
    :param alignments: tuples of ('nonhybrid_sequence', None or 'hybrid_sequence') alignments
    :type alignments: list
    :param precursor_tolerance: the parts per million error allowed when trying 
        to match precursor masses. 
        (default is 10)
    :type percursor_tolerance: int
    :param DEV: set to True if truth is a valid dictionary and fall off detection is 
        desired
        (default is False)
    :type DEV: bool
    :param truth: a set of id keyed spectra with the desired spectra. A better 
        description of what this looks like can be 
        seen in the param.py file. If left None, the program will continue normally
        (default is None)
    :type truth: dict
    :param fall_off: only works if the truth param is set to a dictionary. This 
        is a dictionary (if using multiprocessing, needs to be process safe) 
        where, if a sequence loses the desired sequence, a key value pair of 
        spectrum id, DevFallOffEntry object are added to it. 
        (default is None)
    :type fall_off: dict

    :returns: tuples of refined alignments ('nonhybrid_sequence', None or 'hybrid_sequence')
    :rtype: list
    '''

    global PRECURSOR_MASS_COUNT, AMBIGUOUS_REMOVAL_COUNT, OUT_OF_RANGE_SEQS
    global PRECURSOR_MASS_TIME, AMBIGUOUS_REMOVAL_TIME

    # get the predicted length of the sequence and allow for a 25% gap to be filled in
    predicted_len = utils.predicted_len(spectrum.precursor_mass, spectrum.precursor_charge)
    allowed_gap = math.ceil(predicted_len * .25)

    # Limit our search to things that match our precursor mass
    # try and fill in the gaps that are in any alignments
    st = time.time()
    precursor_matches = []

    for sequence_pairs in alignments:
        
        # take the sequence. If hybrid, take the hybrid, otherwise the non hybrid
        sequence = sequence_pairs[0] if sequence_pairs[1] is None else sequence_pairs[1]

        # add the closer precursors to the list
        p_ms = [
            x for x in \
            alignment_utils.match_precursor(spectrum, sequence, db, gap=allowed_gap, tolerance=precursor_tolerance)
        ]

        if len(p_ms) and p_ms[0] is None:
            OUT_OF_RANGE_SEQS += 1
            continue

        precursor_matches += p_ms

    PRECURSOR_MASS_COUNT += len(alignments)
    PRECURSOR_MASS_TIME += time.time() - st

    # check to see if we no longer have the match. At this point we should
    if DEV:

        # get the id, the sequnce, and if its a hybrid
        _id = spectrum.id
        is_hybrid = truth[_id]['hybrid']
        truth_seq = truth[_id]['sequence']

        if not utils.DEV_contains_truth_exact(truth_seq, is_hybrid, precursor_matches):

            # add metadata about wwhat we had before filling in precursor and 
            # what we ended up after
            metadata = {
                'sequences_before_precursor_filling': alignments, 
                'sequences_after_precursor_filling': precursor_matches, 
                'observed_precursor_mass': spectrum.precursor_mass, 
                'observed_percursor_charge': spectrum.precursor_charge, 
                'allowed_gap': allowed_gap
            }

            fall_off[_id] = DEVFallOffEntry(
                is_hybrid, 
                truth_seq, 
                'precursor_filling', 
                metadata
            )

            # exit the alignment
            return Alignments(spectrum, [])

    # seperate the hybrids from the non hybrids for later analysis
    nonhyba, hyba = [], []
    for p_m in precursor_matches:

        if '-' in p_m or '(' in p_m or ')' in p_m:
            hyba.append((p_m.replace('-', '').replace('(', '').replace(')', ''), p_m))
        else:
            nonhyba.append((p_m, None))
    
    # replace any hybrid alignments that are seen that can be explained by non 
    # hybrid sequences
    st = time.time()
    updated_hybrids = [] if len(hyba) == 0 else hybrid_alignment.replace_ambiguous_hybrids(hyba, db, spectrum)

    # check to see if we lost the match, but only if the sequence is a hybrid
    if DEV and truth[spectrum.id]['hybrid']:

        # get the id, the sequnce, and if its a hybrid
        _id = spectrum.id
        is_hybrid = truth[_id]['hybrid']
        truth_seq = truth[_id]['sequence']

        # updated hybrids is a list of [(nonhyb, hybrid)]. Do the first value because sometimes the 
        # second value is none because its no longer a hybrid
        if not utils.DEV_contains_truth_exact(truth_seq, is_hybrid, [x[0] for x in updated_hybrids]):

            # add some metadata about what we had before and after ambiguous changing
            metadata = {
                'before_ambiguous_removal': hyba, 
                'after_ambiguous_removal': updated_hybrids
            }

            fall_off[_id] = DEVFallOffEntry(
                is_hybrid, 
                truth_seq, 
                'removing_ambiguous_hybrids', 
                metadata
            )

            # exit the alignment
            return Alignments(spectrum, [])

    AMBIGUOUS_REMOVAL_COUNT += len(hyba)
    AMBIGUOUS_REMOVAL_TIME += time.time() - st


    return nonhyba + updated_hybrids

def attempt_alignment_dev(    
    spectrum: Spectrum, 
    truth: bool = None, 
    fall_off: bool = None,
    a: list = None
)-> Alignments:
    b_seqs = [x[0] for x in a]
    y_seqs = [x[0] for x in a]
    _id = spectrum.id
    is_hybrid = truth[_id]['hybrid']
    truth_seq = truth[_id]['sequence']

    if not utils.DEV_contains_truth_parts(truth_seq, is_hybrid, b_seqs, y_seqs):
        metadata = {
            'alignments': a, 
            'before_alignments_b': b_seqs, 
            'before_alignments_y': y_seqs
        }
        fall_off[_id] = DEVFallOffEntry(
            is_hybrid, 
            truth_seq, 
            'first_alignment_round', 
            metadata
        )
        return Alignments(spectrum, [], None)

def attempt_alignment_first_pass(
    spectrum: Spectrum, 
    db: Database, 
    n: int = 3, 
    ppm_tolerance: int = 20, 
    precursor_tolerance: int = 10,
    digest_type: str = '',
    truth: bool = None, 
    fall_off: bool = None,
    DEV: bool = False,
    OBJECTIFY_COUNT: int = 0,
    OBJECTIFY_TIME: int = 0,
    a: list = None,
    is_last: bool = False
)-> Alignments:
    refine_start = time.time()
    non_hybrid_refined = refine_alignments(
        spectrum, 
        db, 
        [x for x in a if x[1] is None], 
        precursor_tolerance=precursor_tolerance, 
        DEV=DEV, 
        truth=truth, 
        fall_off=fall_off
    )
    refine_time = time.time() - refine_start
    Non_hybrid_refine_time.append(refine_time)
    non_hybrid_alignments = []
    tracker = {}
    st = time.time()
    for nhr, _ in non_hybrid_refined:
        if nhr in tracker: 
            continue
        tracker[nhr] = True 
        scoring_start = time.time()
        p_d = scoring.precursor_distance(
            spectrum.precursor_mass, 
            gen_spectra.get_precursor(nhr, spectrum.precursor_charge)
        )
        b_score = scoring.score_sequence(
            spectrum.mz_values, 
            sorted(gen_spectra.gen_spectrum(nhr, ion='b')['spectrum']), 
            ppm_tolerance
        )

        y_score = scoring.score_sequence(
            spectrum.mz_values, 
            sorted(gen_spectra.gen_spectrum(nhr, ion='y')['spectrum']), 
            ppm_tolerance
        )

        total_error = scoring.total_mass_error(spectrum, nhr, ppm_tolerance)

        t_score = b_score + y_score + scoring.digest_score(nhr, db, digest_type)

        non_hybrid_scoring_times.append(time.time() - scoring_start)
        parents = alignment_utils.get_parents(nhr, db)

        non_hybrid_alignments.append(   
            SequenceAlignment(
                parents[0], 
                nhr, 
                b_score, 
                y_score, 
                t_score, 
                p_d, 
                total_error
            )
        )

    OBJECTIFY_COUNT += len(non_hybrid_refined)
    OBJECTIFY_TIME += time.time() - st
    if any([x.total_score >= 1.5 * len(x.sequence) for x in non_hybrid_alignments]):
        sorted_alignments = sorted(
            non_hybrid_alignments, 
            key=lambda x: (
                x.total_score, 
                math.inf if x.total_mass_error <= 0 else 1/x.total_mass_error, 
                math.inf if x.precursor_distance <= 0 else 1/x.precursor_distance, 
                x.b_score, 
                x.y_score
            ), 
            reverse=True
        )
        top_n_alignments = sorted_alignments[:n]
        if is_last:
            alignment_instrumentation = objects.Alignment_Instrumentation(
            B_and_Y_full_bipartite_alignment = FIRST_ALIGN_TIME,
            average_dataset_size_1 = FIRST_ALIGN_COUNT/TOTAL_ITERATIONS,
            seconds_op_1 = FIRST_ALIGN_TIME/FIRST_ALIGN_COUNT,
            removing_ambiguous_hybrids_time = AMBIGUOUS_REMOVAL_TIME,
            average_dataset_size_2 = AMBIGUOUS_REMOVAL_COUNT/TOTAL_ITERATIONS,
            seconds_op_2 = AMBIGUOUS_REMOVAL_TIME/AMBIGUOUS_REMOVAL_COUNT,
            matching_precursor_masses_time = PRECURSOR_MASS_TIME,
            average_dataset_size_3 = PRECURSOR_MASS_COUNT/TOTAL_ITERATIONS,
            seconds_op_3 = PRECURSOR_MASS_TIME/PRECURSOR_MASS_COUNT,
            turning_matches_into_objects_time = OBJECTIFY_TIME,
            average_dataset_size_4 = OBJECTIFY_COUNT/TOTAL_ITERATIONS,
            seconds_op_4 = OBJECTIFY_TIME/OBJECTIFY_COUNT,
            initial_sequences_with_too_many_or_few_amino_acids_to_try_to_precursor_match = OUT_OF_RANGE_SEQS/PRECURSOR_MASS_COUNT
            )
            current_alignments = Alignments(spectrum, top_n_alignments,alignment_instrumentation)
        return current_alignments,None
    else:
        return None, non_hybrid_alignments        

def attempt_alignment_second_pass(    
    spectrum: Spectrum, 
    db: Database, 
    n: int = 3, 
    ppm_tolerance: int = 20, 
    precursor_tolerance: int = 10,
    digest_type: str = '',
    truth: bool = None, 
    fall_off: bool = None,
    DEV: bool = False,
    OBJECTIFY_COUNT: int = 0,
    OBJECTIFY_TIME: int = 0,
    a: list = [],
    non_hybrid_alignments: list = [],
    is_last: bool = False
) -> Alignments:
    refine_start = time.time()
    hybrid_refined = refine_alignments(
        spectrum, 
        db, 
        [x for x in a if x[1] is not None], 
        precursor_tolerance=precursor_tolerance, 
        DEV=DEV, 
        truth=truth, 
        fall_off=fall_off
    )
    refine_time = time.time() - refine_start
    Hybrid_refine_times.append(refine_time)
    hybrid_alignments = []
    tracker = {}
    st = time.time()
    for hr, special_hr in hybrid_refined:
        if hr in tracker: 
            continue
        tracker[hr] = True 
        scoring_start = time.time()
        p_d = scoring.precursor_distance(
            spectrum.precursor_mass, 
            gen_spectra.get_precursor(hr, spectrum.precursor_charge)
        )
        b_score = scoring.score_sequence(
            spectrum.mz_values, 
            sorted(gen_spectra.gen_spectrum(hr, ion='b')['spectrum']), 
            ppm_tolerance
        )
        y_score = scoring.score_sequence(
            spectrum.mz_values, 
            sorted(gen_spectra.gen_spectrum(hr, ion='y')['spectrum']), 
            ppm_tolerance
        )
        total_error = scoring.total_mass_error(spectrum, hr, ppm_tolerance)
        hybrid_scoring_times.append(time.time() - scoring_start)
        parents = alignment_utils.get_parents(hr, db)
        t_score = None
        if special_hr is None:
            t_score = b_score + y_score + scoring.digest_score(hr, db, digest_type)
            non_hybrid_alignments.append(
                 SequenceAlignment(
                    parents[0], 
                    hr, 
                    b_score, 
                    y_score, 
                    t_score, 
                    p_d, 
                    total_error
                )
            )
            #instrumentation
            i = 0
            avg_b_score = 0
            avg_y_score = 0
            avg_total_score = 0
            for score in non_hybrid_alignments:
                 avg_b_score = avg_b_score + score[2]
                 avg_y_score = avg_y_score + score[3]
                 avg_total_score = avg_total_score + score[4]
                 i = i + 1
            alignment_instrumentation = objects.Alignment_Instrumentation(
            avg_b_score = avg_b_score/i,
            avg_y_score = avg_y_score/i,
            avg_total_score = avg_total_score/i
            )
        else: 

            t_score = scoring.hybrid_score(spectrum, special_hr, ppm_tolerance)\
            + scoring.digest_score(special_hr, db, digest_type)
            hybrid_alignments.append(
                HybridSequenceAlignment(
                    parents[0], 
                    parents[1], 
                    hr, 
                    special_hr, 
                    b_score, 
                    y_score, 
                    t_score, 
                    p_d, 
                    total_error
                )
            )
            #instrumentation
            i = 0
            avg_b_score = 0
            avg_y_score = 0
            avg_total_score = 0
            for score in hybrid_alignments:
                avg_b_score = avg_b_score + score[4]
                avg_y_score = avg_y_score + score[5]
                avg_total_score = avg_total_score + score[6]
                i = i + 1
            objects.Alignment_Instrumentation(
            avg_b_score = avg_b_score/i,
            avg_y_score = avg_y_score/i,
            avg_total_score = avg_total_score/i
            )


    OBJECTIFY_COUNT += len(hybrid_refined)
    OBJECTIFY_TIME += time.time() - st
    all_alignments = non_hybrid_alignments + hybrid_alignments
    sorted_alignments = sorted(
            all_alignments, 
            key=lambda x: (
                x.total_score, 
                math.inf if x.total_mass_error <= 0 else 1/x.total_mass_error, 
                math.inf if x.precursor_distance <= 0 else 1/x.precursor_distance, 
                1/len(x),
                x.b_score, 
                x.y_score
            ), 
            reverse=True
        )

    top_n_alignments = sorted_alignments[:n]
    #instrumentation
    if is_last:
        alignment_instrumentation = objects.Alignment_Instrumentation(
        B_and_Y_full_bipartite_alignment = FIRST_ALIGN_TIME,
        average_dataset_size_1 = FIRST_ALIGN_COUNT/TOTAL_ITERATIONS,
        seconds_op_1 = FIRST_ALIGN_TIME/FIRST_ALIGN_COUNT,
        removing_ambiguous_hybrids_time = AMBIGUOUS_REMOVAL_TIME,
        average_dataset_size_2 = AMBIGUOUS_REMOVAL_COUNT/TOTAL_ITERATIONS,
        seconds_op_2 = AMBIGUOUS_REMOVAL_TIME/AMBIGUOUS_REMOVAL_COUNT,
        matching_precursor_masses_time = PRECURSOR_MASS_TIME,
        average_dataset_size_3 = PRECURSOR_MASS_COUNT/TOTAL_ITERATIONS,
        seconds_op_3 = PRECURSOR_MASS_TIME/PRECURSOR_MASS_COUNT,
        turning_matches_into_objects_time = OBJECTIFY_TIME,
        average_dataset_size_4 = OBJECTIFY_COUNT/TOTAL_ITERATIONS,
        seconds_op_4 = OBJECTIFY_TIME/OBJECTIFY_COUNT,
        initial_sequences_with_too_many_or_few_amino_acids_to_try_to_precursor_match = OUT_OF_RANGE_SEQS/PRECURSOR_MASS_COUNT
        )
    return Alignments(spectrum, top_n_alignments,alignment_instrumentation) 

def attempt_alignment(
    spectrum: Spectrum, 
    db: Database, 
    b_hits: list,
    y_hits: list, 
    n: int = 3, 
    ppm_tolerance: int = 20, 
    precursor_tolerance: int = 10,
    digest_type: str = '',
    DEBUG: bool = False, 
    is_last: bool = False, 
    truth: bool = None, 
    fall_off: bool = None
) -> Alignments:
    global FIRST_ALIGN_TIME, AMBIGUOUS_REMOVAL_TIME, PRECURSOR_MASS_TIME, OBJECTIFY_TIME
    global FIRST_ALIGN_COUNT, AMBIGUOUS_REMOVAL_COUNT, PRECURSOR_MASS_COUNT, OBJECTIFY_COUNT
    global TOTAL_ITERATIONS
    TOTAL_ITERATIONS += 1
    DEV = truth is not None and fall_off is not None
    extension_time = time.time()
    b_non_hybrids, y_non_hybrids = extend_base_kmers(b_hits, y_hits, spectrum, db)
    extension_times.append(time.time() - extension_time)
    non_hybrids = b_non_hybrids + y_non_hybrids
    st = time.time()
    initial_alignment_start = time.time()
    a = align_b_y(b_hits, y_hits, spectrum, db) + [(kmer, None) for kmer in non_hybrids]
    initial_alignment_times.append(time.time() - initial_alignment_start)
    FIRST_ALIGN_COUNT += len(b_hits) + len(y_hits)
    FIRST_ALIGN_TIME += time.time() - st
    if DEV:
        return attempt_alignment_dev(spectrum,truth,fall_off,a)
    alignments, non_hybrid_alignments = attempt_alignment_first_pass(spectrum,db,n,ppm_tolerance,precursor_tolerance,digest_type,truth,fall_off,DEV,OBJECTIFY_COUNT,OBJECTIFY_TIME,a,is_last)
    if alignments is not None:
        return alignments
    else:
        return attempt_alignment_second_pass(spectrum,db,n,ppm_tolerance,precursor_tolerance,digest_type,truth,fall_off,DEV,OBJECTIFY_COUNT,OBJECTIFY_TIME,a,non_hybrid_alignments,is_last)
   
