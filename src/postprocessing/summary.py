from lookups.utils import make_dir, make_valid_dir_string
from file_io import JSON
from lookups.objects import Alignments

import pandas as pd
import json
import os

SUMMARY_NAME = 'summary'
HYBRID_PREFIX = 'hybrid_'

def json_file(results: dict, output_dir: str) -> None:
    json_file_name = os.path.join(output_dir, f'{SUMMARY_NAME}.json')
    dictified = {}
    for name, alignment in results.items():
        dictified[name] = {
            'spectrum': alignment.spectrum._asdict(), 
            'alignments': [x._asdict() for x in alignment.alignments]
        }
    JSON.save_dict(json_file_name, dictified)
    
def get_protein_strings(left_proteins, right_proteins):
    left_protein_string, right_protein_string = '',''
    for protein in left_proteins:
        left_protein_string += protein + "/"
    for protein in right_proteins:
        right_protein_string += protein + "/"
    
    return left_protein_string[:-1], right_protein_string[:-1]

def get_score_strings(b_scores, y_scores):
    b_score_string, y_score_string = '',''
    for score in b_scores:
        b_score_string += str(score) + "/"
    for score in y_scores:
        y_score_string += str(score) + "/"
        
    return b_score_string[:-1], y_score_string[:-1]
        
def get_extensions_strings(extensions):
    b_extension_string = ''
    y_extension_string = ''
    left = extensions[0]
    right = extensions[1]
    for b in left:
        b_extension_string += b + "/"
    for y in right:
        y_extension_string += y + "/"
    return b_extension_string[:-1], y_extension_string[:-1]

def get_output_row_from_aligned_spectrums(spectrum_id, aligned_spectrum):
            hybrid = aligned_spectrum.hybrid
            left_proteins,right_proteins = aligned_spectrum.left_proteins, aligned_spectrum.right_proteins
            sequence = aligned_spectrum.sequence
            b_scores = aligned_spectrum.b_scores
            y_scores = aligned_spectrum.y_scores
            total_score = aligned_spectrum.total_score
            total_gaussian_score = str(aligned_spectrum.total_gaussian_score)
            extensions = aligned_spectrum.extensions
            precursor_mass, precursor_charge = str(aligned_spectrum.precursor_mass), str(aligned_spectrum.precursor_charge)
            total_mass_error = str(aligned_spectrum.total_mass_error)
            left_protein_string, right_protein_string = get_protein_strings(left_proteins, right_proteins)
            b_score_string, y_score_string = get_score_strings(b_scores, y_scores)
            b_extension_strings, y_extension_strings = get_extensions_strings(extensions)
            total_count = int(total_score * len(sequence))
            fraction_form = str(total_count) + "/" + str(len(sequence))
            output_row = str(spectrum_id) + '\t' + hybrid + '\t' + sequence + '\t' + str(total_score) + '\t' + fraction_form + "\t" + total_gaussian_score + '\t' + total_mass_error + "\t" + precursor_mass + '\t' + precursor_charge + '\t' + left_protein_string + '\t' + right_protein_string + '\t' + b_score_string + '\t' + y_score_string + '\t' + b_extension_strings + '\t' + y_extension_strings + '\n' 
            return output_row

def tsv_file(results: dict, output_dir: str) -> None:
    mac = 0
    hybrids, nonhybrids = [], []
    alignment: Alignments
    for name, alignment in results.items():
        if len(alignment.alignments) == 0:
            mac += 1
            continue
        topalignment = alignment.alignments[0]._asdict()
        topalignment['entry name'] = name
        topalignment['id'] = alignment.spectrum.id
        if 'hybrid_sequence' in topalignment:
            hybrids.append(topalignment)
        else:
            nonhybrids.append(topalignment)
    hybridresults = pd.DataFrame(hybrids)
    with open(f'{output_dir + HYBRID_PREFIX + SUMMARY_NAME}.tsv', 'w') as ho:
        ho.write(hybridresults.to_csv(sep='\t'))
    del hybridresults
    del hybrids
    nonhybridresults = pd.DataFrame(nonhybrids)
    output_file = os.path.join(output_dir, f'{SUMMARY_NAME}.tsv')
    with open(output_file, 'w') as nho:
        nho.write(nonhybridresults.to_csv(sep='\t'))

def generate(alignments: dict, output_dir='./') -> None:
    output_dir = make_valid_dir_string(output_dir)
    make_dir(output_dir)
    json_file(alignments, output_dir)
    tsv_file(alignments, output_dir)

def write_aligned_spectrums_to_disk(aligned_spectrums, output_folder_path, output_file_name ):
        file_name = os.path.basename(output_file_name)
        A = file_name.split(".")
        base_file_name = "HS_"+ A[0] + ".txt"
        output_file_path = os.path.join(output_folder_path,base_file_name)   
        with open(output_file_path, 'w') as t:
            t.write("spectrum_id" + '\t' + "hybrid" + '\t' + "sequence" + '\t' + "total score" + '\t' + "#peaks/length" + '\t' + "total gaussian score" + '\t' + "total mass error" + '\t' + "precursor mass" + '\t' + "precursor charge" + '\t' + "left kmer" + '\t' + "right kmer" + '\t' + "b score" + '\t' + "y score" + '\t' + "prev aa" + '\t' + "next aa" + '\n')
            spectrum_id = 0
            for aligned_spectrum in aligned_spectrums:
                output_row = get_output_row_from_aligned_spectrums(spectrum_id,aligned_spectrum)
                t.write(output_row)
                spectrum_id+=1
