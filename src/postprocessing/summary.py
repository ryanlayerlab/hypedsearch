from utils import make_dir, make_valid_dir_string
from file_io import JSON
from objects import Alignments

import pandas as pd
import json
import os

SUMMARY_NAME = 'summary'
HYBRID_PREFIX = 'hybrid_'

def json_file(results: dict, output_dir: str) -> None:
    '''
    Generate a summary json file for the results made

    Inputs:
        results: (dict) containing the results made. The key of each entry is the name
                        and each entry should be an Alignments namedtuple
        output_dir: (str) path to the output directory
    Outputs:
        None
    '''
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


def text_file(results: dict, txt_file_name: str) -> None:
    with open(txt_file_name, 'w') as t:
        t.write("spectrum_id" + '\t' + "hybrid" + '\t' + "sequence" + '\t' + "total score" + '\t' + "total abundance" + '\t' + "left kmer" + '\t' + "right kmer" + '\t' + "b score" + '\t' + "y score" + '\t' + "prev aa" + '\t' + "next aa" + '\n')
        for i, x in enumerate(results):
            target_alignments = x
            for alignment in target_alignments:
                spec_num = str(i)
                hybrid = alignment[0]
                left_proteins,right_proteins = alignment[1], alignment[2]
                left_protein_string, right_protein_string = get_protein_strings(left_proteins, right_proteins)
                sequence = alignment[3]
                b_scores = alignment[4]
                y_scores = alignment[5]
                b_score_string, y_score_string = get_score_strings(b_scores, y_scores)
                total_score = str(alignment[6])
                total_abundance = str(alignment[7])
                extensions = alignment[8]
                b_extension_strings, y_extension_strings = get_extensions_strings(extensions)
                
                t.write(spec_num + '\t' + hybrid + '\t' + sequence + '\t' + total_score + '\t' + total_abundance + '\t' + left_protein_string + '\t' + right_protein_string + '\t' + b_score_string + '\t' + y_score_string + '\t' + b_extension_strings + '\t' + y_extension_strings + '\n')

def tsv_file(results: dict, output_dir: str) -> None:
    '''
    Write the results of the experiment to 2 tsv files. One tsv file is for 
    non hybrid identified sequences, the other is for hybrid identified sequences.

    Inputs: 
        results:    (dict) results of the search. The key of each entry is the name
                            of each entry and the value is an Alignments namedtuple 
        output_dir: (str) path to the directory to save the tsvs
    Outputs:
        None
    '''
    mac = 0

    # seperate the hybrids from the nonhybrids
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

    # move to pandas dataframe for easy writing
    hybridresults = pd.DataFrame(hybrids)
    with open(f'{output_dir + HYBRID_PREFIX + SUMMARY_NAME}.tsv', 'w') as ho:
        ho.write(hybridresults.to_csv(sep='\t'))

    del hybridresults
    del hybrids

    nonhybridresults = pd.DataFrame(nonhybrids)
    output_file = os.path.join(output_dir, f'{SUMMARY_NAME}.tsv')
    
    with open(output_file, 'w') as nho:
        nho.write(nonhybridresults.to_csv(sep='\t'))

    print(f'Could not make an alignment for {mac}/{len(results)} spectra ({int(100 * mac / len(results))}%)')

def generate(alignments: dict, output_dir='./') -> None:
    '''
    Generate a summary text and json file for the alignments made

    Inputs:
        alignments: dict containing the alignments made. The key of each entry is the name
                    of the file appended with scan number, and the values should be Alignments
    kwargs:
        output_dir: str path to the output directory. Default=./
    Outputs:
        None
    '''
    output_dir = make_valid_dir_string(output_dir)
    make_dir(output_dir)

    json_file(alignments, output_dir)
    tsv_file(alignments, output_dir)

def generate_to_txt(matched_spectra, spectra_files, output_dir) -> None:

    for i, file in enumerate(spectra_files):
        alignments = matched_spectra[i]
        filename = os.path.basename(file)
        A = filename.split(".")
        base_file = "HS_"+ A[0] + ".txt"
        output_file = os.path.join(output_dir, "Hypedsearch_outputs" ,base_file)
        text_file(alignments, output_file)