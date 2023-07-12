from scoring import scoring
from objects import Spectrum, SequenceAlignment, HybridSequenceAlignment, Database, Alignments, DEVFallOffEntry
from constants import PROTON_MASS, WATER_MASS
from alignment import alignment_utils, hybrid_alignment
from gen_spectra import get_precursor
from preprocessing.clustering import calc_from_sequences
import objects
import utils
import database
import gen_spectra
from sqlite import database_file

import math
import collections
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
        b_seq = b_kmer
        b_proteins = database.get_proteins_with_subsequence(db, b_seq)
        for y_kmer in y_kmers:
            y_seq = y_kmer
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
        return top_n_alignments
    else:
        return None, non_hybrid_alignments        

def attempt_alignment_second_pass(spectrum: Spectrum, db: Database, n: int = 3, 
    ppm_tolerance: int = 20, precursor_tolerance: int = 10,digest_type: str = '',truth: bool = None, 
    fall_off: bool = None,DEV: bool = False,OBJECTIFY_COUNT: int = 0,OBJECTIFY_TIME: int = 0,
    a: list = [],non_hybrid_alignments: list = []):
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
    return Alignments(spectrum, top_n_alignments) 

def attempt_alignment(spectrum: Spectrum, db: Database, b_hits: list,y_hits: list, n: int = 3, 
    ppm_tolerance: int = 20, precursor_tolerance: int = 10,digest_type: str = '',DEBUG: bool = False, 
    truth: bool = None, fall_off: bool = None):
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
    alignments, non_hybrid_alignments = attempt_alignment_first_pass(spectrum,db,n,ppm_tolerance,precursor_tolerance,digest_type,truth,fall_off,DEV,OBJECTIFY_COUNT,OBJECTIFY_TIME,a)
    if alignments is not None:
        return alignments
    else:
        return attempt_alignment_second_pass(spectrum,db,n,ppm_tolerance,precursor_tolerance,digest_type,truth,fall_off,DEV,OBJECTIFY_COUNT,OBJECTIFY_TIME,a,non_hybrid_alignments)

def make_merge(b, y, b_seq, y_seq):
    new_b = (b[0], b[1], b[2], b[3], b_seq)
    new_y = (y[0], y[1], y[2], y[3], y_seq)
    return (b[3] + y[3], b[1] - y[2], y[2]-b[1], new_b, new_y)

def native_get_extensions(precursor_mass,prec_charge,b_side,y_side,prec_tol):
    extended_cluster = collections.namedtuple('sorted_cluster', 'score pid start end mz charge')
    tol = utils.ppm_to_da(precursor_mass, prec_tol)
    extensions = []
    #if both the clusters simply glue together
    if b_side.end -1 == y_side.start: # Does it have to be -1?
        combined_prec = gen_spectra.calc_precursor_as_disjoint(b_side.mz, y_side.mz, b_side.charge, y_side.charge, prec_charge)
        if abs(combined_prec - precursor_mass) <= tol:
            b_cluster = extended_cluster(b_side.score, b_side.pid, b_side.start, b_side.end, b_side.mz, b_side.charge)
            y_cluster = extended_cluster(y_side.score, y_side.pid, y_side.start, y_side.end, y_side.mz, y_side.charge)
            return [(b_side.score + y_side.score, b_cluster, y_cluster)]
    
    #if the clusters are not h
    for b in b_side.components:
        if b[2]-1 == y_side.start:
            total_precursor = gen_spectra.calc_precursor_as_disjoint(b[0], y_side.mz, b[4], y_side.charge, prec_charge)
            if abs(total_precursor - precursor_mass) < prec_tol:
                b_half = extended_cluster(score = b_side.score, pid=b[5], start=b[1], end = b[2], mz = b[0], charge= b[4])
                y_half = extended_cluster(score = y_side.score, pid=y_side.pid, start=y_side.start, end = y_side.end, mz = y_side.mz, charge=y_side.charge)
                extensions.append((b[0] + y_side.score, b_half, y_half))
    return extensions

def get_extensions(precursor_mass, precursor_charge, b_side, y_side, prec_tol):
    extended_cluster = collections.namedtuple('sorted_cluster', 'score pid start end mz charge')
    tol = utils.ppm_to_da(precursor_mass, prec_tol)
    extensions = []
    
    #Check if the two pieces just go together already  
    combined_prec = gen_spectra.calc_precursor_as_disjoint(b_side.mz, y_side.mz, b_side.charge, y_side.charge, precursor_charge)
    if abs(combined_prec - precursor_mass) <= tol:
        b_cluster = extended_cluster(b_side.score, b_side.pid, b_side.start, b_side.end, b_side.mz, b_side.charge)
        y_cluster = extended_cluster(y_side.score, y_side.pid, y_side.start, y_side.end, y_side.mz, y_side.charge)
        return [(b_side.score + y_side.score, b_cluster, y_cluster)]
        
    #Only extend b
    for b in b_side.components:
        b_mass = b[0]
        this_prec = gen_spectra.calc_precursor_as_disjoint(b_mass, y_side.mz, b[4], y_side.charge, precursor_charge)
        if this_prec > precursor_mass + tol:
            break
        elif abs(this_prec - precursor_mass) <= tol:
            extended_b_cluster = extended_cluster(score = b_side.score, pid=b[5], start=b[1], end = b[2], mz = b[0], charge= b[4])
            extended_y_cluster = extended_cluster(score = y_side.score, pid=y_side.pid, start=y_side.start, end = y_side.end, mz = y_side.mz, charge = y_side.charge)
            extensions.append((b_side.score + y_side.score, extended_b_cluster, extended_y_cluster))
    
    #Only extend y
    for y in y_side.components:
        y_mass = y[0]
        this_prec = gen_spectra.calc_precursor_as_disjoint(b_side.mz, y_mass, b_side.charge, y[4], precursor_charge)
        if this_prec > precursor_mass + tol:
            break
        elif abs(this_prec - precursor_mass) <= tol:
            extended_b_cluster = extended_cluster(score = b_side.score, pid=b_side.pid, start=b_side.start, end = b_side.end, mz = b_side.mz, charge = b_side.charge)
            extended_y_cluster = extended_cluster(score = y_side.score, pid=y[5], start=y[1], end = y[2], mz = y[0], charge= y[4])
            extensions.append((b_side.score + y_side.score, extended_b_cluster, extended_y_cluster))
    
    #Either side needs to be extended
    for b in b_side.components:
        for y in y_side.components:
            b_mass, y_mass = b[0], y[0]
            this_prec = gen_spectra.calc_precursor_as_disjoint(b_mass, y_mass, b[4], y[4], precursor_charge)
            if this_prec > precursor_mass + tol:
                break
            elif abs(this_prec - precursor_mass) <= tol:
                extended_b_cluster = extended_cluster(score = b_side.score, pid=b[5], start=b[1], end = b[2], mz = b[0], charge= b[4])
                extended_y_cluster = extended_cluster(score = y_side.score, pid=y[5], start=y[1], end = y[2], mz = y[0], charge= y[4])
                extensions.append((b_side.score + y_side.score, extended_b_cluster, extended_y_cluster))
            
    return extensions

def find_alignments(native_merged, hybrid_merged, obs_prec, prec_charge, tol, max_len, prec_tol):
    extended_cluster = collections.namedtuple('sorted_cluster', 'score pid start end mz charge')
    natural_alignments, hybrid_alignments = [], []
    for i, comb_seq in enumerate(native_merged): #Maybe read the top 50 of these, look for matches and go until we find something
        b_side = comb_seq[1]
        y_side = comb_seq[2]
        # print("Natural", i)
        if y_side.start >= b_side.end: #no overlap but b before y
            natural_alignments = natural_alignments + native_get_extensions(obs_prec, prec_charge, b_side, y_side, prec_tol)
        elif b_side.start <= y_side.start and b_side.end <= y_side.end and y_side.start < b_side.end: #some overlap
            combined_precursor = calc_from_sequences(b_side.start, y_side.end, b_side.pid, max_len, prec_charge)
            if abs(combined_precursor - obs_prec) < tol:
                b_cluster = extended_cluster(b_side.score, b_side.pid, b_side.start, b_side.end, b_side.mz, b_side.charge)
                y_cluster = extended_cluster(y_side.score, y_side.pid, y_side.start, y_side.end, y_side.mz, y_side.charge)
                natural_alignments.append((b_side.score + y_side.score, b_cluster, y_cluster))
        else:
            hybrid_alignments = hybrid_alignments + get_extensions(obs_prec, prec_charge, b_side, y_side, prec_tol)

    total_extension_time = 0
    for i, comb_seq in enumerate(hybrid_merged):
        b_side, y_side = comb_seq[1], comb_seq[2]
        extension_time = time.time()
        hybrid_alignments = hybrid_alignments + get_extensions(obs_prec, prec_charge, b_side, y_side, prec_tol)
        total_extension_time = total_extension_time + (time.time() - extension_time)
    # print("\n Average extension time:",total_extension_time/len(hybrid_merged))
    return natural_alignments, hybrid_alignments

def pair_indices(b_search_space, y_search_space, prec_mass, prec_tol, prec_charge, score_filter):
    pairs = []
    tol = utils.ppm_to_da(prec_mass, prec_tol)
    sorted_b_keys = sorted(b_search_space) #might want to produce this outside this function so we don't do this again
    sorted_y_keys = sorted(y_search_space, reverse = True)
    
    good_b_prec, good_y_prec = gen_spectra.get_precursor("DLQTLAL", prec_charge), gen_spectra.get_precursor("EVE", prec_charge)
        
    b_end, y_end = len(sorted_b_keys), len(sorted_y_keys)
    b_ctr, y_ctr = 0, 0
    
    for b_prec in sorted_b_keys:
        missing_mass_upper = prec_mass - b_prec + (prec_charge * PROTON_MASS)/prec_charge + WATER_MASS/prec_charge + tol
        missing_mass_low = prec_mass - b_prec + (prec_charge * PROTON_MASS)/prec_charge + WATER_MASS/prec_charge - tol
        for y_prec in sorted_y_keys:
            if y_prec > missing_mass_upper:
                continue
            elif y_prec < missing_mass_low:
                break
            else:
                for b in b_search_space[b_prec]:
                    for y in y_search_space[y_prec]:
                        if b[7] + y[7] > score_filter:
                            pairs.append((b, y))
    
    # while b_ctr < b_end and y_ctr < y_end: #(188.58901117418367, (376.1707458496094, 56, 59, 1, 1, 274, 'EVE', 2))
    #     b_prec, y_prec = sorted_b_keys[b_ctr], sorted_y_keys[y_ctr]
    #     if b_prec == 387.223811816879:
    #         print('here')
    #     missing_mass_upper = prec_mass - b_prec + (prec_charge * PROTON_MASS)/prec_charge + WATER_MASS/prec_charge + tol
    #     if y_prec > missing_mass_upper:
    #         y_ctr += 1
    #         continue
            
    #     missing_mass_low = prec_mass - b_prec + (prec_charge * PROTON_MASS)/prec_charge + WATER_MASS/prec_charge - tol
    #     if y_prec < missing_mass_low:
    #         b_ctr += 1
    #         continue
        
    #     for b in b_search_space[b_prec]:
    #         for y in y_search_space[y_prec]:
    #             if "DLQTLAL" == b[6] and y[6] == "EVE":
    #                 print('here')
    #             if b[7] + y[7] > score_filter:
    #                 pairs.append((b, y))
    #     b_ctr += 1

    return pairs
        
def find_from_prec(converted_b, matched_masses_b, input_spectrum, ppm_tolerance, protein_list):
    prec_matches = []
    prec_hits = matched_masses_b[converted_b]
    for hit in prec_hits:
        if hit[4] == 2:
            hit_score, hit_abundance = scoring.prec_score(hit, input_spectrum, ppm_tolerance, protein_list)
            prec_matches.append((hit_score, hit_abundance, hit))
        
    prec_matches = sorted(prec_matches, key = lambda x: (x[0], x[1]), reverse=True)
    if len(prec_matches) != 0:
        return prec_matches[0], prec_matches[0][0]
    else:
        return [], 0
    
def make_native_pair(b, ion):
    y_cluster = (gen_spectra.max_mass(b[6], 'b' if ion == 0 else 'y', b[4]), b[1], b[2], ion, b[4], b[5], b[6], 0)
    return y_cluster
    
def pair_natives(b_search_space, y_search_space, prec_mass, prec_tol, score_filter):
    # in the case of a native, the b must be before the y, the two seqs must come from the same protein, and the two seqs together cannot be less than the max_pep_len. Actually, 
    # Actually, since we have all the extensions we could just find the good extension
    native_pairs = []
    tol = utils.ppm_to_da(prec_mass, prec_tol)
    for b_prec in sorted(b_search_space):
        if abs(b_prec - prec_mass) < tol:
            for b in b_search_space[b_prec]:
                # if b[7] >= score_filter
                y_pair = make_native_pair(b, 0)
                native_pairs.append((b,y_pair))
    
    for y_prec in sorted(y_search_space):
        if abs(y_prec - prec_mass) < tol:
            for y in y_search_space[y_prec]:
                # if y[7] >= score_filter:
                b_pair = make_native_pair(y, 1)
                native_pairs.append((b_pair, y))
                
    return native_pairs