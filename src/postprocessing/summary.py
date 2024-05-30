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

def create_text_file(aligned_spectrums: dict, txt_file_name: str) -> None:
    with open(txt_file_name, 'w') as t:
        t.write("spectrum_id" + '\t' + "hybrid" + '\t' + "sequence" + '\t' + "total score" + '\t' + "#peaks/length" + '\t' + "total gaussian score" + '\t' + "total mass error" + '\t' + "precursor mass" + '\t' + "precursor charge" + '\t' + "left kmer" + '\t' + "right kmer" + '\t' + "b score" + '\t' + "y score" + '\t' + "prev aa" + '\t' + "next aa" + '\n')
        if aligned_spectrums is None:
            pass
        else:
            for index, target_alignments in enumerate(aligned_spectrums):
                for alignments in target_alignments:
                    alignment = alignments[0]
                    spec_num = str(index)
                    t.write(alignment)

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
        filename = os.path.basename(output_file_name)
        A = filename.split(".")
        base_file = "HS_"+ A[0] + ".txt"
        output_file = os.path.join(output_folder_path,base_file)        
        create_text_file(aligned_spectrums, output_file)