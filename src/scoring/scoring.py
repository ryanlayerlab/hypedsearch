from scoring import mass_comparisons
from objects import Spectrum, Database
from utils import ppm_to_da
from preprocessing import clustering
from constants import WATER_MASS
import gen_spectra
import utils
import database

import re

# load digests json for digest scoring
import json 
import os

script_dir = os.path.dirname(__file__)
json_dir = '/'.join(script_dir.split('/')[:-1])
digest_file = os.path.join(json_dir, 'digests.json')

digests = json.load(open(digest_file, 'r'))

def calc_mass_given_other_explanations(unique_m, seq, mz):
    oEXPnum = (len(unique_m[mz]) - 1)/ len(unique_m[mz])
    if oEXPnum == 0:
        return 0
    else:
        p = 0
        for i, seq2 in enumerate(unique_m[mz]):
            if seq == seq2:
                continue
            else:
                p = p + 1/len(seq2)
        return p

def Bayes_given_mass(pH, seq, mz, unique_m):
    pEH = 1/len(seq)
    pnH = 1-pH
    pEnH = calc_mass_given_other_explanations(unique_m, seq, mz)
    prob = (pH * pEH)/((pH*pEH)+(pnH*pEnH))
#     print(seq,pH,pEH,(pH*pEH),(pnH*pEnH),pnH,pEnH)
    return prob

def calc_bayes_score(seq, mz, unique_m, indices, kmer_set):
    pH = len(seq)/len(kmer_set)
    for index in reversed(indices):
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
    spectrum = gen_spectra.gen_spectrum(sequence)
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
    
def score_by_dist(b_side, y_side, obs_prec, prec_charge, max_len, hybrid, protein_list): #Change bullet points to be a query
    #((bmass,bstart,bend,ion,charge,pid)(ymass,ystart,yend,ion,charge,pid))
    if hybrid:
        full_seq = clustering.find_sequence(b_side.pid, b_side.start, b_side.end, protein_list) + clustering.find_sequence(y_side.pid, y_side.start, y_side.end, protein_list)
    else:
        full_seq = clustering.find_sequence(b_side.pid, b_side.start, y_side.end, protein_list)
    combined_precursor = gen_spectra.get_precursor(full_seq, prec_charge)
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
    spectrum = gen_spectra.gen_spectrum(sequence)
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
    spectrum = gen_spectra.gen_spectrum(sequence)
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
