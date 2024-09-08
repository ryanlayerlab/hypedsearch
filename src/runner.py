
import os
from pyteomics import fasta
from objects import FastaDatabase, ExperimentParameters, Fragment, Precursor
from lookups.sqlite_database import Sqllite_Database
import identification
import lookups.utils
import multiprocessing as mp
from datetime import datetime
from lookups.utils import file_exists
from pyteomics import mzml

def get_number_of_cores(number_of_cores):
    revised_number_of_cores = number_of_cores = max(1, number_of_cores)
    min_number_of_cores = min(revised_number_of_cores, mp.cpu_count() - 1)
    return min_number_of_cores

def create_fasta_database(database_file_path):
    fasta_database = FastaDatabase(database_file_path,{},{})
    prots = []
    for entry in fasta.read(database_file_path):
        prots.append(entry)
    fasta_database = fasta_database._replace(proteins=prots)
    return fasta_database  

def create_sqllite_database(database_file_path,max_peptide_length,digest_left,digest_right):
    fasta_database = create_fasta_database(database_file_path)
    sqllite_database = Sqllite_Database(max_peptide_length, True)
    kv_proteins = [(k, v) for k, v in fasta_database.proteins]  
    sqllite_database.populate_database(kv_proteins, max_peptide_length, digest_left, digest_right)
    return sqllite_database

def get_existing_sqllite_database(max_peptide_length):
    sqllite_database = Sqllite_Database(max_peptide_length, False)
    return sqllite_database

def get_output_file_name(spectra_file_paths):
    current_datetime = datetime.now()
    formatted_datetime = current_datetime.strftime('%Y%m%d%H%M%S')
    first_spectra_file_path = spectra_file_paths[0]
    last_token_with_extension = os.path.basename(first_spectra_file_path)
    filename, extension = os.path.splitext(last_token_with_extension)
    return_value = filename + "_" + formatted_datetime
    return return_value

def relative_abundance_filtering(masses, abundances, percentage):
    ti = sum(abundances)
    min_value = ti * percentage
    filtered_mass_abundances = [x for x in zip(masses, abundances) if x[1] >= min_value]
    masses = [float(x) for x, _ in filtered_mass_abundances]
    abundances = [float(x) for _, x in filtered_mass_abundances]
    return (masses, abundances)

def peak_filtering(masses, abundances, num_peaks):
    mass_abundances = zip(masses, abundances)
    mass_abundances = sorted(mass_abundances, key=lambda x: x[1], reverse=True)[:num_peaks]
    mass_abundances.sort(key=lambda x: x[0])
    masses = [float(x) for x, _ in mass_abundances]
    abundances = [float(x) for _, x in mass_abundances]
    return (masses, abundances)

def get_precursors(spectra_file_path,number_peaks):
    if not file_exists(spectra_file_path):
        print('File {} not found. Please make sure that this file exists'.format(spectra_file_path))
        return

    precursors = []
    filecontents = mzml.read(spectra_file_path)
    content: dict
    for precursor_id, content in enumerate(filecontents):
        precursor_mass = None
        precursor_charge = 0
        if not len(content['precursorList']['precursor']) or not len(content['precursorList']['precursor'][0]['selectedIonList']['selectedIon']):
            precursor_mass = max(masses)
            precursor_charge = 1
        else:
            precursor_mass = float(content['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['selected ion m/z'])
            if 'peak intensity' in content['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]:
                precursor_abundance = float(content['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['peak intensity'])
            else:
                precursor_abundance = None
            precursor_charge = int(content['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['charge state'])
        precursor_description = content.get('id', '')
        retention_time = content['scanList']['scan'][0]['scan start time']
        masses = list(content['m/z array'])
        abundances = list(content['intensity array'])
        if number_peaks > 0:
            adjusted_masses, adjusted_abundances = peak_filtering(masses, abundances, number_peaks)
        elif relative_abundance_filter > 0:
            while relative_abundance_filter > 1:
                relative_abundance_filter /= 100
            adjusted_masses, adjusted_abundances = relative_abundance_filtering(masses, abundances, relative_abundance_filter)
        fragments = []
        for fragment_id, adjusted_mass in enumerate(adjusted_masses):
            adjusted_abundance = adjusted_abundances[fragment_id]
            fragment = Fragment(id=fragment_id,precursor_id=precursor_id, precursor_mass=precursor_mass, mz_value=adjusted_mass,abundance=adjusted_abundance)       
            fragments.append(fragment)
        precursor = Precursor(id=precursor_id,description=precursor_description,mass=precursor_mass,charge=precursor_charge,retention_time=retention_time,abundance=precursor_abundance,fragments=fragments)
        precursors.append(precursor)
    return precursors

def get_all_precursors(spectra_file_paths,number_peaks):
    all_precursors = []
    for spectra_file_path in spectra_file_paths:
        precursors = get_precursors(spectra_file_path,number_peaks)
        all_precursors.append(precursors)
    return precursors

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

def get_output_row_from_aligned_peptide(spectrum_id, aligned_peptide):
            hybrid = aligned_peptide.hybrid
            left_proteins,right_proteins = aligned_peptide.left_proteins, aligned_peptide.right_proteins
            sequence = aligned_peptide.sequence
            b_scores = aligned_peptide.b_scores
            y_scores = aligned_peptide.y_scores
            total_score = aligned_peptide.total_score
            total_gaussian_score = str(aligned_peptide.total_gaussian_score)
            extensions = aligned_peptide.extensions
            precursor_mass, precursor_charge = str(aligned_peptide.precursor_mass), str(aligned_peptide.precursor_charge)
            total_mass_error = str(aligned_peptide.total_mass_error)
            left_protein_string, right_protein_string = get_protein_strings(left_proteins, right_proteins)
            b_score_string, y_score_string = get_score_strings(b_scores, y_scores)
            b_extension_strings, y_extension_strings = get_extensions_strings(extensions)
            total_count = int(total_score * len(sequence))
            fraction_form = str(total_count) + "/" + str(len(sequence))
            output_row = str(spectrum_id) + '\t' + hybrid + '\t' + sequence + '\t' + str(total_score) + '\t' + fraction_form + "\t" + total_gaussian_score + '\t' + total_mass_error + "\t" + precursor_mass + '\t' + precursor_charge + '\t' + left_protein_string + '\t' + right_protein_string + '\t' + b_score_string + '\t' + y_score_string + '\t' + b_extension_strings + '\t' + y_extension_strings + '\n' 
            return output_row

def write_aligned_peptides_to_disk(aligned_peptides, output_folder_path, output_file_name ):
    file_name = os.path.basename(output_file_name)
    A = file_name.split(".")
    base_file_name = "HS_"+ A[0] + ".txt"
    output_file_path = os.path.join(output_folder_path,base_file_name)   
    with open(output_file_path, 'w') as t:
        t.write("spectrum_id" + '\t' + "hybrid" + '\t' + "sequence" + '\t' + "total score" + '\t' + "#peaks/length" + '\t' + "total gaussian score" + '\t' + "total mass error" + '\t' + "precursor mass" + '\t' + "precursor charge" + '\t' + "left kmer" + '\t' + "right kmer" + '\t' + "b score" + '\t' + "y score" + '\t' + "prev aa" + '\t' + "next aa" + '\n')
        spectrum_id = 0
        for aligned_peptide in aligned_peptides:
            output_row = get_output_row_from_aligned_peptide(spectrum_id,aligned_peptide)
            t.write(output_row)
            spectrum_id+=1
 
def get_experiment_parameters(args: dict):
    all_precursors = get_all_precursors(args['spectra_file_paths'],args['number_peaks'])
    sqllite_database = None
    if args['create_sqllite_database']:
        sqllite_database = create_sqllite_database(args['database_file_path'],args['max_peptide_length'], args['digest_left'], args['digest_right'])
    else:
        sqllite_database = get_existing_sqllite_database(args['max_peptide_length'])
    max_peptide_length=args['max_peptide_length']
    ppm_tolerance=args['ppm_tolerance']
    precursor_tolerance=args['precursor_tolerance']
    number_hybrids=args['number_hybrids']
    number_natives=args['number_natives']
    target_seq = args['target_seq']
    experiment_parameters = ExperimentParameters(id=0,precursors=all_precursors,sqllite_database=sqllite_database,
                                                max_peptide_length=max_peptide_length, ppm_tolerance=ppm_tolerance,
                                                precursor_tolerance=precursor_tolerance,number_hybrids=number_hybrids,
                                                number_natives=number_natives,target_seq=target_seq)
    return experiment_parameters

def handle_results(args, aligned_peptides):
    output_file_name = get_output_file_name(args['spectra_file_paths']) 
    lookups.utils.make_dir(args['output_folder_path'])
    output_folder_path=args['output_folder_path']
    write_aligned_peptides_to_disk(aligned_peptides, output_folder_path, output_file_name)

def run(args: dict):
    experiment_parameters = get_experiment_parameters(args)
    aligned_peptides = identification.get_aligned_peptides(experiment_parameters)  
    # handle_results(args,aligned_peptides)
    # return aligned_peptides