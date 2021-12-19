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

def same_protein_alignment(seq1: str, seq2: str, parent_sequence: str):
    if seq1 == seq2:
        return (seq1, None)    
    if seq1 in seq2:
        return (seq2, None)    
    if seq2 in seq1: 
        return (seq1, None)    
    gap_aa = max(1, len(parent_sequence) // 100)
    left_start = [m.start() for m in re.finditer(seq1, parent_sequence)]
    right_start = [m.start() for m in re.finditer(seq2, parent_sequence)]
    if all([r < l for r in right_start for l in left_start]):
        return hybrid_alignment.hybrid_alignment(seq1, seq2)
    nonhybrid_alignments = []
    for l in left_start:
        for r in right_start:
            if r < l:
                continue
            if r - (l + len(seq1)) <= gap_aa:
                overlapped = parent_sequence[l: r + len(seq2)]
                nonhybrid_alignments.append(overlapped)
    if len(nonhybrid_alignments) == 0:
        return hybrid_alignment.hybrid_alignment(seq1, seq2)
    nonhybrid_alignments.sort(key=lambda x: len(x))
    return (nonhybrid_alignments[0], None)

def align_b_y(b_kmers: list, y_kmers: list, spectrum: Spectrum, db: Database):
    spec_alignments = []
    for b_kmer in b_kmers:
        b_seq = b_kmer[2][2]
        b_proteins = database.get_proteins_with_subsequence(db, b_seq)
        for y_kmer in y_kmers:
            y_seq = y_kmer[2][2]
            y_proteins = database.get_proteins_with_subsequence(db, y_seq)
            if any([x in y_proteins for x in b_proteins]):
                shared_prots = [x for x in y_proteins if x in b_proteins]
                for sp in shared_prots:
                    prot_seqs = database.get_entry_by_name(db, sp)
                    for prot_entry in prot_seqs:
                        spec_alignments.append(
                            same_protein_alignment(b_seq, y_seq, prot_entry.sequence)
                        )
                spec_alignments.append((f'{b_seq}{y_seq}', f'{b_seq}-{y_seq}'))
            else: 
                spec_alignments.append(hybrid_alignment.hybrid_alignment(b_seq, y_seq))
    return list(set([x for x in spec_alignments if x is not None]))

def extend_base_kmers(b_kmers: list, y_kmers: list, spectrum: Spectrum, db: Database):
    extended_b = []
    extended_y = []
    for seq in b_kmers:
        extended_b += [x for x in alignment_utils.extend_non_hybrid(seq, spectrum, 'b', db)]
    for seq in y_kmers:
        extended_y += [x for x in alignment_utils.extend_non_hybrid(seq, spectrum, 'y', db)]
    return extended_b, extended_y

def refine_alignments(spectrum: Spectrum, db: Database, alignments: list, precursor_tolerance: int = 10,
    DEV: bool = False, truth: dict = None, fall_off: dict = None):
    global PRECURSOR_MASS_COUNT, AMBIGUOUS_REMOVAL_COUNT, OUT_OF_RANGE_SEQS
    global PRECURSOR_MASS_TIME, AMBIGUOUS_REMOVAL_TIME
    predicted_len = utils.predicted_len(spectrum.precursor_mass, spectrum.precursor_charge)
    allowed_gap = math.ceil(predicted_len * .25)
    st = time.time()
    precursor_matches = []
    for sequence_pairs in alignments:
        sequence = sequence_pairs[0] if sequence_pairs[1] is None else sequence_pairs[1]
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
    if DEV:
        _id = spectrum.id
        is_hybrid = truth[_id]['hybrid']
        truth_seq = truth[_id]['sequence']
        if not utils.DEV_contains_truth_exact(truth_seq, is_hybrid, precursor_matches):
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
            return Alignments(spectrum, [])
    nonhyba, hyba = [], []
    for p_m in precursor_matches:
        if '-' in p_m or '(' in p_m or ')' in p_m:
            hyba.append((p_m.replace('-', '').replace('(', '').replace(')', ''), p_m))
        else:
            nonhyba.append((p_m, None))
    st = time.time()
    updated_hybrids = [] if len(hyba) == 0 else hybrid_alignment.replace_ambiguous_hybrids(hyba, db, spectrum)
    if DEV and truth[spectrum.id]['hybrid']:
        _id = spectrum.id
        is_hybrid = truth[_id]['hybrid']
        truth_seq = truth[_id]['sequence']
        if not utils.DEV_contains_truth_exact(truth_seq, is_hybrid, [x[0] for x in updated_hybrids]):
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
            return Alignments(spectrum, [])
    AMBIGUOUS_REMOVAL_COUNT += len(hyba)
    AMBIGUOUS_REMOVAL_TIME += time.time() - st
    return nonhyba + updated_hybrids

def attempt_alignment_dev(spectrum: Spectrum, truth: bool = None, fall_off: bool = None,a: list = None):
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

def attempt_alignment_first_pass(spectrum: Spectrum, db: Database, n: int = 3, ppm_tolerance: int = 20, 
    precursor_tolerance: int = 10,digest_type: str = '',truth: bool = None, fall_off: bool = None,
    DEV: bool = False,OBJECTIFY_COUNT: int = 0,OBJECTIFY_TIME: int = 0,a: list = None,is_last: bool = False):
    refine_start = time.time()
    non_hybrid_refined = refine_alignments(spectrum, db, [x for x in a if x[1] is None], 
        precursor_tolerance=precursor_tolerance, DEV=DEV, truth=truth, fall_off=fall_off)
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

def attempt_alignment_second_pass(spectrum: Spectrum, db: Database, n: int = 3, 
    ppm_tolerance: int = 20, precursor_tolerance: int = 10,digest_type: str = '',truth: bool = None, 
    fall_off: bool = None,DEV: bool = False,OBJECTIFY_COUNT: int = 0,OBJECTIFY_TIME: int = 0,
    a: list = [],non_hybrid_alignments: list = [],is_last: bool = False):
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
            i = 0
            avg_b_score = 0
            avg_y_score = 0
            avg_total_score = 0
            for score in hybrid_alignments:
                avg_b_score = avg_b_score + score[4]
                avg_y_score = avg_y_score + score[5]
                avg_total_score = avg_total_score + score[6]
                i = i + 1
            alignment_instrumentation = objects.Alignment_Instrumentation(
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

def attempt_alignment(spectrum: Spectrum, db: Database, b_hits: list,y_hits: list, n: int = 3, 
    ppm_tolerance: int = 20, precursor_tolerance: int = 10,digest_type: str = '',DEBUG: bool = False, 
    is_last: bool = False, truth: bool = None, fall_off: bool = None):
    global FIRST_ALIGN_TIME, AMBIGUOUS_REMOVAL_TIME, PRECURSOR_MASS_TIME, OBJECTIFY_TIME
    global FIRST_ALIGN_COUNT, AMBIGUOUS_REMOVAL_COUNT, PRECURSOR_MASS_COUNT, OBJECTIFY_COUNT
    global TOTAL_ITERATIONS
    TOTAL_ITERATIONS += 1
    DEV = truth is not None and fall_off is not None
    b_non_hybrids, y_non_hybrids = extend_base_kmers(b_hits, y_hits, spectrum, db)
    non_hybrids = b_non_hybrids + y_non_hybrids
    a = align_b_y(b_hits, y_hits, spectrum, db) + [(kmer, None) for kmer in non_hybrids]
    FIRST_ALIGN_COUNT += len(b_hits) + len(y_hits)
    if DEV:
        return attempt_alignment_dev(spectrum,truth,fall_off,a)
    alignments, non_hybrid_alignments = attempt_alignment_first_pass(spectrum,db,n,ppm_tolerance,precursor_tolerance,digest_type,truth,fall_off,DEV,OBJECTIFY_COUNT,OBJECTIFY_TIME,a,is_last)
    if alignments is not None:
        return alignments
    else:
        return attempt_alignment_second_pass(spectrum,db,n,ppm_tolerance,precursor_tolerance,digest_type,truth,fall_off,DEV,OBJECTIFY_COUNT,OBJECTIFY_TIME,a,non_hybrid_alignments,is_last)
   
