from lookups.utils import make_dir, make_valid_dir_string
from file_io import JSON
from lookups.objects import Alignments

import pandas as pd
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
    left_protein_string = '/'.join(left_proteins)
    right_protein_string = '/'.join(right_proteins)
    
    return left_protein_string, right_protein_string

def get_score_strings(b_scores, y_scores):
    b_score_string = '/'.join(b_scores)
    y_score_string = '/'.join(y_scores)
        
    return b_score_string, y_score_string
        
def get_extensions_strings(extensions):
    left, right = extensions[:2]
    b_extension_string = '/'.join(left)
    y_extension_string = '/'.join(right)

    return b_extension_string, y_extension_string

def create_text_file(results: dict, txt_file_name: str) -> None:
    with open(txt_file_name, 'w') as t:
        t.write("spectrum_id" + '\t' + "hybrid" + '\t' + "sequence" + '\t' + "total score" + '\t' + "#peaks/length" + '\t' + "total gaussian score" + '\t' + "total mass error" + '\t' + "precursor mass" + '\t' + "precursor charge" + '\t' + "left kmer" + '\t' + "right kmer" + '\t' + "b score" + '\t' + "y score" + '\t' + "prev aa" + '\t' + "next aa" + '\n')
        if results is not None:
            for i, x in enumerate(results):
                target_alignments = x
                for alignment in target_alignments:
                    spec_num = str(i)
                    hybrid = alignment[0]

                    left_proteins,right_proteins = alignment[1:3]
                    left_protein_string, right_protein_string = get_protein_strings(left_proteins, right_proteins)

                    sequence, b_scores, y_scores = alignment[3:6]
                    b_score_string, y_score_string = get_score_strings(b_scores, y_scores)

                    total_score = alignment[6]
                    total_gaussian_score = str(alignment[7])
                    extensions = alignment[8]
                    b_extension_strings, y_extension_strings = get_extensions_strings(extensions)
                    precursor_mass, precursor_charge = str(alignment[9]), str(alignment[10])
                    total_mass_error = str(alignment[11])
                    total_count = int(total_score * len(sequence))
                    fraction_form = str(total_count) + "/" + str(len(sequence))
                    t.write(spec_num + '\t' + hybrid + '\t' + sequence + '\t' + str(total_score) + '\t' + fraction_form + "\t" + total_gaussian_score + '\t' + total_mass_error + "\t" + precursor_mass + '\t' + precursor_charge + '\t' + left_protein_string + '\t' + right_protein_string + '\t' + b_score_string + '\t' + y_score_string + '\t' + b_extension_strings + '\t' + y_extension_strings + '\n')


def tsv_file(results: dict, output_dir: str) -> None:
    mac = 0
    hybrids, nonhybrids = [], []
    alignment: Alignments
    for name, alignment in results.items():
        if len(alignment.alignments) == 0:
            mac += 1
            continue
        top_alignment = alignment.alignments[0]._asdict()
        top_alignment['entry name'] = name
        top_alignment['id'] = alignment.spectrum.id
        if 'hybrid_sequence' in top_alignment:
            hybrids.append(top_alignment)
        else:
            nonhybrids.append(top_alignment)
    hybrid_results = pd.DataFrame(hybrids)
    with open(f'{output_dir}{HYBRID_PREFIX}{SUMMARY_NAME}.tsv', 'w') as hybrid_output:
        hybrid_output.write(hybrid_results.to_csv(sep='\t'))
    non_hybrid_results = pd.DataFrame(nonhybrids)
    output_file = os.path.join(output_dir, f'{SUMMARY_NAME}.tsv')
    with open(output_file, 'w') as non_hybrid_output:
        non_hybrid_output.write(non_hybrid_results.to_csv(sep='\t'))

def generate(alignments: dict, output_dir='./') -> None:
    output_dir = make_valid_dir_string(output_dir)
    make_dir(output_dir)
    json_file(alignments, output_dir)
    tsv_file(alignments, output_dir)

def write_matched_spectrum_to_disk(alignments, output_folder_path, output_file_name ) -> None:
        filename = os.path.basename(output_file_name)
        A = filename.split(".")
        base_file = f"HS_{A[0]}.txt"
        output_file = os.path.join(output_folder_path,base_file)
        create_text_file(alignments, output_file)