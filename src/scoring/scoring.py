from scoring import mass_comparisons
from lookups.objects import Spectrum, Database
from lookups.utils import ppm_to_da
from preprocessing import clustering
from lookups.constants import WATER_MASS, AMMONIUM
from math import exp
import computational_pipeline.gen_spectra
import lookups.utils
import preprocessing.database_generator
import re
import json 
import os

script_dir = os.path.dirname(__file__)
json_dir = '/'.join(script_dir.split('/')[:-1])

def calc_mass_given_other_explanations(unique_m, seq, mz):
    oEXPnum = (len(unique_m[mz]) - 1)/ len(unique_m[mz])
    if oEXPnum == 0:
        return 0

    p = 0
    for seq2 in unique_m[mz]:
        if seq != seq2: 
            p = p + 1/len(seq2)
    return p

def Bayes_given_mass(pH, seq, mz, unique_m):
    pEH = 1/len(seq)
    pnH = 1-pH
    pEnH = calc_mass_given_other_explanations(unique_m, seq, mz)
    prob = (pH * pEH)/((pH*pEH)+(pnH*pEnH))
    return prob

def calc_bayes_score(seq, mz, unique_m, indices, kmer_set):
    pH = len(seq)/len(kmer_set)
    prob = Bayes_given_mass(pH, seq, mz, unique_m)
    pH = prob
    return prob

def parse_indices(index_set):
    indices = []
    for index in index_set:
        string = str(index)
        A = string.rstrip().split(',')
        start = A[0]
        end = A[1]
        seq = A[2]
        mz = A[3]
        disallowed_characters = " ()\'"
        for character in disallowed_characters:
            start = start.replace(character, "")
            end = end.replace(character, "")
            seq = seq.replace(character, "")
            mz = mz.replace(character, "")
        
        target_tuple = (int(start), int(end), seq, float(mz))
        indices.append(target_tuple)
    
    return indices

def rescore_with_seq(sequence, ppm_tolerance, input_masses):
    spectrum = computational_pipeline.gen_spectra.gen_spectrum(sequence)
    masses = sorted(spectrum['spectrum'])
    input_masses = sorted(input_masses)
    o_ctr, t_ctr = 0, 0
    observed = input_masses[o_ctr]
    theoretical = masses[t_ctr]
    total_score = 0
    while (o_ctr < len(input_masses) and t_ctr < len(masses)):
        tol = ppm_to_da(observed, ppm_tolerance)
        if theoretical < observed - tol:
            t_ctr = t_ctr + 1
            if t_ctr < len(masses):
                theoretical = masses[t_ctr]
        elif observed + tol < theoretical: #The bug is with 810 around here
            o_ctr = o_ctr + 1
            if o_ctr < len(input_masses):
                observed = input_masses[o_ctr]
        elif abs(observed-theoretical) <= tol:
            total_score = total_score + 1
            o_ctr = o_ctr + 1
            t_ctr = t_ctr + 1
            if o_ctr < len(input_masses) and t_ctr < len(masses):
                observed = input_masses[o_ctr]
                theoretical = masses[t_ctr]
                
    return(total_score)
    
def score_by_dist(b_side, y_side, obs_prec, prec_charge, max_len, hybrid, protein_list):
    if hybrid:
        full_seq = clustering.find_sequence(b_side.pid, b_side.start, b_side.end, protein_list) + clustering.find_sequence(y_side.pid, y_side.start, y_side.end, protein_list)
    else:
        full_seq = clustering.find_sequence(b_side.pid, b_side.start, y_side.end, protein_list)
    combined_precursor = computational_pipeline.gen_spectra.get_precursor(full_seq, prec_charge)
    dist = abs(combined_precursor - obs_prec)
    return dist

def calc_overlap(masses, input_masses, ppm_tolerance):
    total_score = 0
    o_ctr, t_ctr = 0, 0
    observed = input_masses[o_ctr]
    theoretical = masses[t_ctr]
    while (o_ctr < len(input_masses) and t_ctr < len(masses)):
        tol = ppm_to_da(observed, ppm_tolerance)
        if theoretical < observed - tol:
            t_ctr = t_ctr + 1
            if t_ctr < len(masses):
                theoretical = masses[t_ctr]
        elif observed + tol < theoretical:
            o_ctr = o_ctr + 1
            if o_ctr < len(input_masses):
                observed = input_masses[o_ctr]
        elif abs(observed-theoretical) <= tol:
            total_score = total_score + 1
            o_ctr = o_ctr + 1
            t_ctr = t_ctr + 1
            if o_ctr < len(input_masses) and t_ctr < len(masses):
                observed = input_masses[o_ctr]
                theoretical = masses[t_ctr]
                
    return(total_score)

def overlap_scoring(b_side, y_side, input_masses, ppm_tolerance, hybrid, proteins):
    if hybrid:
        b_sequence, y_sequence = clustering.find_sequence(b_side.pid, b_side.start, b_side.end, proteins), clustering.find_sequence(y_side.pid, y_side.start, y_side.end, proteins)
        sequence = b_sequence + y_sequence
    else:
        sequence = clustering.find_sequence(b_side.pid, b_side.start, y_side.end, proteins)
    spectrum = computational_pipeline.gen_spectra.gen_spectrum(sequence)
    masses = sorted(spectrum['spectrum'])
    input_masses = sorted(input_masses)
    score = calc_overlap(masses, input_masses, ppm_tolerance)
    return score

def losing_water(b_side, y_side, input_masses, ppm_tolerance, hybrid, proteins):
    if hybrid:
        b_sequence, y_sequence = clustering.find_sequence(b_side.pid, b_side.start, b_side.end, proteins), clustering.find_sequence(y_side.pid, y_side.start, y_side.end, proteins)
        sequence = b_sequence + y_sequence
    else:
        sequence = clustering.find_sequence(b_side.pid, b_side.start, y_side.end, proteins)
    spectrum = computational_pipeline.gen_spectra.gen_spectrum(sequence)
    masses = sorted(spectrum['spectrum'])
    minus_water = []
    for mass in masses:
        minus_water.append(mass - WATER_MASS)
    input_masses = sorted(input_masses)
    score = calc_overlap(minus_water, input_masses, ppm_tolerance)
    return score    

def second_scoring(natural_alignments, hybrid_alignments, input_spectrum, tol, proteins, max_len):
    rescored_naturals, rescored_hybrids = [], []
    for comb_seq in natural_alignments:
        b_side, y_side = comb_seq[1], comb_seq[2]
        dist = score_by_dist(b_side, y_side, input_spectrum.precursor_mass, input_spectrum.precursor_charge, max_len, False, proteins)
        score = overlap_scoring(b_side, y_side, input_spectrum.mz_values, tol, False, proteins)
        score += losing_water(b_side, y_side, input_spectrum.mz_values, tol, False, proteins)
        rescored_naturals.append((score, 1/dist, comb_seq, 0))
    for comb_seq in hybrid_alignments:
        b_side, y_side = comb_seq[1], comb_seq[2]
        dist = score_by_dist(b_side, y_side, input_spectrum.precursor_mass, input_spectrum.precursor_charge, max_len, True, proteins)
        score = overlap_scoring(b_side, y_side, input_spectrum.mz_values, tol, True, proteins)
        score += losing_water(b_side, y_side, input_spectrum.mz_values, tol, True, proteins)
        rescored_hybrids.append((score, 1/dist, b_side, y_side, 1))
    return rescored_naturals, rescored_hybrids

def calc_overlap(masses, input_masses, ppm_tolerance):
    total_score = 0
    tiebreaker = 1
    ppm_sum = 0
    o_ctr, t_ctr = 0, 0
    observed = input_masses[o_ctr]
    theoretical = masses[t_ctr]
    while (o_ctr < len(input_masses) and t_ctr < len(masses)):
        tol = ppm_to_da(observed, ppm_tolerance)
        if theoretical < observed - tol:
            t_ctr = t_ctr + 1
            if t_ctr < len(masses):
                theoretical = masses[t_ctr]
        elif observed + tol < theoretical:
            o_ctr = o_ctr + 1
            if o_ctr < len(input_masses):
                observed = input_masses[o_ctr]
        elif abs(observed-theoretical) <= tol:
            total_score = total_score + 1
            ppm_mass_error = ((observed-theoretical)/theoretical) * 1000000 #converting to ppm
            ppm_sum += ppm_mass_error
            similarityScore_gauss = exp(-0.5 * abs(ppm_mass_error/ppm_tolerance)**2)
            tiebreaker = tiebreaker * similarityScore_gauss
            o_ctr = o_ctr + 1
            t_ctr = t_ctr + 1
            if o_ctr < len(input_masses) and t_ctr < len(masses):
                observed = input_masses[o_ctr]
                theoretical = masses[t_ctr]
                
    return total_score, tiebreaker, ppm_sum

def modified_overlap_scoring(sequence, input_masses, ppm_tolerance):
    spectrum = computational_pipeline.gen_spectra.generate_spectrum(sequence)
    masses = sorted(spectrum['spectrum'])
    input_masses = sorted(input_masses)
    score, tiebreaker, mass_error_sum = calc_overlap(masses, input_masses, ppm_tolerance)
    return score, tiebreaker, mass_error_sum

def modified_losing_water(sequence, input_masses, ppm_tolerance):
    masses,_ = computational_pipeline.gen_spectra.calculate_masses(sequence)
    input_masses = sorted(input_masses)
    score, tiebreaker, mass_error_sum = calc_overlap(masses, input_masses, ppm_tolerance)
    return score, tiebreaker, mass_error_sum

def modified_losing_ammonium(sequence, input_masses, ppm_tolerance):
    masses,_ = computational_pipeline.gen_spectra.calculate_masses(sequence)
    input_masses = sorted(input_masses)
    score, tiebreaker, mass_error_sum = calc_overlap(masses, input_masses, ppm_tolerance)
    return score, tiebreaker, mass_error_sum

def rescore_merges(rescore_merges_params):
    unique_merge_space = rescore_merges_params.unique_merge_space
    input_spectrum = rescore_merges_params.input_spectrum
    ppm_tol = rescore_merges_params.ppm_tol
    rescored_unique = dict()
    for key, hyb in unique_merge_space:
        score, tiebreaker, ppm_sum = modified_overlap_scoring(key, input_spectrum.mz_values, ppm_tol) # counts peaks
        minus_water_score, minus_water_tiebreaker, minus_water_ppm_sum = modified_losing_water(key, input_spectrum.mz_values, ppm_tol) # counts again but strip water
        minus_ammonium_score, minus_ammonium_tiebreaker, minus_ammonium_ppm_sum = modified_losing_ammonium(key, input_spectrum.mz_values, ppm_tol) # counts again but strip ammonium
        score += minus_water_score + minus_ammonium_score
        ppm_sum += minus_water_ppm_sum + minus_ammonium_ppm_sum
        tiebreaker * minus_water_tiebreaker * minus_ammonium_tiebreaker
        if (score/len(key), tiebreaker, key, hyb) not in rescored_unique.keys():
            rescored_unique[(score, tiebreaker, key, hyb)] = []
        for b, y in unique_merge_space[(key, hyb)]:
            rescored_unique[(score, tiebreaker, key, hyb)].append((score, tiebreaker, (b,y), ppm_sum))
    return rescored_unique
               
def prec_overlap_scoring(input_masses, ppm_tolerance, pid, start, end, protein_list):
    sequence = clustering.find_sequence(pid, start, end, protein_list)
    spectrum = computational_pipeline.gen_spectra.generate_spectrum(sequence)
    masses = sorted(spectrum['spectrum'])
    input_masses = sorted(input_masses)
    score, tiebreaker, _ = calc_overlap(masses, input_masses, ppm_tolerance)
    return score, tiebreaker

def prec_score(hit, input_spectrum, ppm_tolerance, protein_list):
    score, tiebreaker = prec_overlap_scoring(input_spectrum.mz_values, ppm_tolerance, hit[5], hit[1], hit[2], protein_list)
    return score, tiebreaker
