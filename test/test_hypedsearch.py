import sys, os.path  
import unittest
import shutil

src_path = (os.path.abspath(os.path.join(os.path.dirname(__file__), '..')) + '/src/')
sys.path.append(src_path)
import runner, utils, database, gen_spectra, database_preprocessing

from alignment import alignment_utils
from objects import Spectrum
from scoring import mass_comparisons, scoring
from constants import AMINO_ACIDS, SINGLY_CHARGED_B_BASE, DOUBLY_CHARGED_B_BASE, SINGLY_CHARGED_Y_BASE, DOUBLY_CHARGED_Y_BASE, PROTON_MASS, WATER_MASS

class Test_Main(unittest.TestCase):
    
    def get_arguments(self):
        dirname = os.path.dirname(__file__)
        spectra_file_paths = [os.path.join(dirname, '../data/spectra/hybrid_nod2e3.mzML')]  
        database_file_path = os.path.join(dirname, '../data/database/sample_database.fasta')
        database_file = database.build(database_file_path)
        output_dir = '../output'
        min_peptide_len = 3
        max_peptide_len = 30
        ppm_tolerance = 20
        precursor_tolerance = 10
        verbose = False
        peak_filter = 25
        relative_abundance_filter = .01
        digest = ''
        debug = False
        cores = 1
        n = 5
        truth_set = ''

        return {'spectra_files': spectra_file_paths, 'database_file': database_file,
            'output_dir': output_dir, 'min_peptide_len': min_peptide_len,
            'max_peptide_len': max_peptide_len,'tolerance': ppm_tolerance,
            'precursor_tolerance': precursor_tolerance,'verbose': verbose, 
            'peak_filter': peak_filter, 'relative_abundance_filter': relative_abundance_filter,
            'digest': digest, 'DEBUG': debug, 'cores': cores,'n': n,
            'truth_set': truth_set}        

    def setUp(self): 
        self.sequence = 'MALWAR'

    def test_ppm_to_da(self):
        #Tests the convert ppm to da function with a mass of 100 and a tolerance of 20. 20 ppm of 100 should be .002
        mass = 100
        tol = 20
        self.assertEqual(utils.ppm_to_da(mass, tol), .002, '20 ppm of 100 should be .002')

    def test_file_exists(self):
        #Run the utils.file_exists function with two sample files.
        filename = 'requirements.txt'
        filename2 = 'thisfiledoesnotexist.txt'
        self.assertTrue(utils.file_exists(filename)) # This file exists in the directory
        self.assertFalse(utils.file_exists(filename2)) # This file does not exist in the directory

    def test_make_valid_dir_string(self):
        #Run the utils.make_valid_dir_string function with path which includes os seperator character and
        # one path which does not include an os seperator character
        path1 = '/mnt/c/Users'
        path2 = "/mnt/c/Users/"
        self.assertEqual(utils.make_valid_dir_string(path1), '/mnt/c/Users/')
        self.assertEqual(utils.make_valid_dir_string(path2), '/mnt/c/Users/')
    
    def test_make_dir(self):
        #Run the utils.make_dir function with a directory which exists and a directory which does not exist
        #Case1: Directory does not exist. Should make a new directory
        dir = os.getcwd() + '/not_a_dir'
        self.assertFalse(os.path.isdir(dir))
        utils.make_dir(dir)

        #Case2: Directory already exists. Nothing should happen
        self.assertTrue(os.path.isdir(dir))
        utils.make_dir(dir)
        self.assertTrue(os.path.isdir(dir))
        shutil.rmtree(dir)
    
    def test_make_valid_text_file(self):
        #Run the utils.make_valid_text_file function with a file which is a text file and a file which is not a text file
        filename = 'requirements.txt'
        filename2 = 'testfile.yaml'
        self.assertEqual(utils.make_valid_text_file(filename), 'requirements.txt')
        self.assertEqual(utils.make_valid_text_file(filename2), 'testfile.yaml.txt')
    
    def test_make_valid_json_file(self):
        #Run the utils.make_valid_json_file function with a file which is a json file and a file which is not a json file
        filename = 'test.json'
        filename2 = 'testfile.yaml'
        self.assertEqual(utils.make_valid_json_file(filename), 'test.json')
        self.assertEqual(utils.make_valid_json_file(filename2), 'testfile.yaml.json')

    def test_make_valid_csv_file(self):
        #Run the utils.make_valid_csv_file function with a file which is a json file and a file which is not a json file
        filename = 'test.csv'
        filename2 = 'testfile.yaml'
        self.assertEqual(utils.make_valid_csv_file(filename), 'test.csv')
        self.assertEqual(utils.make_valid_csv_file(filename2), 'testfile.yaml.csv')

    def test_make_valid_fasta_file(self):
        #Run the utils.make_valid_fasta_file function with a file which is a json file and a file which is not a json file
        filename = 'test.fasta'
        filename2 = 'testfile.yaml'
        self.assertEqual(utils.make_valid_fasta_file(filename), 'test.fasta')
        self.assertEqual(utils.make_valid_fasta_file(filename2), 'testfile.yaml.fasta')
    
    def test_is_json(self):
        #Run the utils.is_json function with a file which is a json and a file which is not a json file
        filename = 'test.json'
        filename2 = 'test.csv'
        self.assertEqual(utils.is_json(filename), True)
        self.assertEqual(utils.is_json(filename2), False)
    
    def test_is_fasta(self):
        #Run the utils.is_fasta function with a file which is a fasta and a file which is not a fasta file
        filename = 'test.fasta'
        filename2 = 'test.csv'
        self.assertEqual(utils.is_fasta(filename), True)
        self.assertEqual(utils.is_fasta(filename2), False)
    
    def test_is_dir(self): 
        #Run the utils.is_dir function with a path which is a valid path to a directory and a path which isn't
        dir = os.getcwd() + '/not_a_dir'
        self.assertFalse(utils.is_dir(dir))
        utils.make_dir(dir)
        self.assertTrue(utils.is_dir(dir))
        shutil.rmtree(dir)

    def test_is_file(self):
        #Run the utils.is_file function with a file and a file which does not exist
        filename = 'requirements.txt'
        filename2 = 'thisfiledoesnotexist.txt'
        self.assertEqual(utils.is_file(filename), True)
        self.assertEqual(utils.is_file(filename2), False)

    def test_all_perms_of_s(self): 
        #Run the utils.all_perms_of_s function with two strings and varying keywords

        string1 = 'LMWHOMP'
        keyletter1 = 'LJI' #expected result: ['LMWHOMP', 'JMWHOMP', 'IMWHOMP']
        string2 = 'MALWAR MZHL'
        keyletter2 = 'LH' #expected result: ['MALWAR MZHL', 'MAHWAR MZHL', 'MALWAR MZLL', 'MALWAR MZHH', 'MAHWAR MZLL', 
        #'MAHWAR MZHH', 'MAHWAR MZLH', 'MALWAR MZLH']
        self.assertEqual(sorted(utils.all_perms_of_s(string1, keyletter1)), sorted(['LMWHOMP', 'JMWHOMP', 'IMWHOMP']))
        self.assertEqual(sorted(utils.all_perms_of_s(string2, keyletter2)), sorted(['MALWAR MZHL', 'MAHWAR MZHL', 
            'MALWAR MZLL', 'MALWAR MZHH', 'MAHWAR MZLL', 'MAHWAR MZHH', 'MAHWAR MZLH', 'MALWAR MZLH']))
    
    def test_to_percent(self):
        #Run the utils.to_percent function with two different values to convert to a percent.
        value1 = 53
        total1 = 100 #Expected: 53%
        value2 = 234
        total2 = 3456 #Expected: 7%
        self.assertEqual(utils.to_percent(value1, total1), 53)
        self.assertEqual(utils.to_percent(value2, total2), 6)
    
    def test_predicted_len(self):
        #Run the utils.predicted_len function with two different precursor masses.
        precursor = 240 #Expected len would be 4
        charge = 1
        self.assertEqual(utils.predicted_len(precursor, charge), 4)
        precursor = precursor * 5
        charge = charge + 1 #Expected len would be 32
        self.assertEqual(utils.predicted_len(precursor, charge), 32)
        
    def test_hashable_boundaries(self):
        #run the hashable_boundaries function with two lists of boundaries
        boundaries1 = [2,10] #Expected to return 2-10
        boundaries2 = [3,5,7] #Expected to return nothing
        self.assertEqual(utils.hashable_boundaries(boundaries1), '2-10')
        self.assertEqual(utils.hashable_boundaries(boundaries2), None)

    def test_split_hybrid(self):  
        #run the split_hybrid function with two samples sequence
        Sequence1 = 'MAL-WAR' #Expected left is "MAL" and expected right is "WAR"
        Sequence2 = 'DLTQTL-B' #Expected left is "DLTQTL" and expected right is "B"
        self.assertEqual(utils.split_hybrid(Sequence1), ('MAL', 'WAR'))
        self.assertEqual(utils.split_hybrid(Sequence2), ('DLTQTL', 'B'))
    
    def test_b_ions(self):
        #test the b_ions function with the sequence, "MALWAR".
        target_sequence = [AMINO_ACIDS['M'], AMINO_ACIDS['M'] + AMINO_ACIDS['A'], 
        AMINO_ACIDS['M'] + AMINO_ACIDS['A'] + AMINO_ACIDS['L'],
        AMINO_ACIDS['M'] + AMINO_ACIDS['A'] + AMINO_ACIDS['L'] + AMINO_ACIDS['W'],
        AMINO_ACIDS['M'] + AMINO_ACIDS['A'] + AMINO_ACIDS['L'] + AMINO_ACIDS['W'] + AMINO_ACIDS['A'],
        AMINO_ACIDS['M'] + AMINO_ACIDS['A'] + AMINO_ACIDS['L'] + AMINO_ACIDS['W'] + AMINO_ACIDS['A'] + AMINO_ACIDS['R']]
        #adding singly charged base
        singly_b_seq = [n + SINGLY_CHARGED_B_BASE for n in target_sequence]
        self.assertEqual(gen_spectra.b_ions(self.sequence, 1), singly_b_seq)
        #adding doubly charged base
        doubly_b_seq = [(n + DOUBLY_CHARGED_B_BASE) / 2 for n in target_sequence]
        #testing doubly charged
        self.assertEqual(gen_spectra.b_ions(self.sequence, 2), doubly_b_seq)
    
    def test_y_ions(self):
        #Test the y_ions function with the sequence, 'MALWAR".
        target_sequence = [AMINO_ACIDS['R'], AMINO_ACIDS['A'] + AMINO_ACIDS['R'],
        AMINO_ACIDS['R'] + AMINO_ACIDS['A'] + AMINO_ACIDS['W'],
        AMINO_ACIDS['R'] + AMINO_ACIDS['A'] + AMINO_ACIDS['W'] + AMINO_ACIDS['L'],
        AMINO_ACIDS['R'] + AMINO_ACIDS['A'] + AMINO_ACIDS['W'] + AMINO_ACIDS['L'] + AMINO_ACIDS['A'],
        AMINO_ACIDS['R'] + AMINO_ACIDS['A'] + AMINO_ACIDS['W'] + AMINO_ACIDS['L'] + AMINO_ACIDS['A'] + AMINO_ACIDS['M']]
        #add in single/doubly charged
        singly_charged_seq = [n + SINGLY_CHARGED_Y_BASE for n in target_sequence]
        self.assertEqual(gen_spectra.y_ions(self.sequence, 1), singly_charged_seq)
        #testing doubly charged
        doubly_charged_seq = [(n + DOUBLY_CHARGED_Y_BASE) / 2 for n in target_sequence]
        self.assertEqual(gen_spectra.y_ions(self.sequence, 2), doubly_charged_seq)

    def test_calc_masses(self):
        #Test the calc masses function with the sequence, "MALWAR".
        #The masses should be all of the b and y ion masses in one list.
        expected = sorted((gen_spectra.b_ions(self.sequence)) + (gen_spectra.y_ions(self.sequence)))
        actual, _ = gen_spectra.calc_masses(self.sequence)
        self.assertEqual(sorted(actual), expected)
    
    def test_max_mass(self):
        #Test the max mass function with the sequence, "MALWAR"
        Expected_max = 728.379201
        charge = 1
        charge2 = 2
        self.assertEqual(gen_spectra.max_mass(self.sequence, 'b', charge), Expected_max + SINGLY_CHARGED_B_BASE)
        self.assertEqual(gen_spectra.max_mass(self.sequence, 'y', charge2), (Expected_max + DOUBLY_CHARGED_Y_BASE) / 2)
    
    def test_get_precursor(self):
        #Test the get_precursory function with the sequence, "MALWAR"
        Expected_prec = AMINO_ACIDS['M'] + AMINO_ACIDS['A'] + AMINO_ACIDS['L'] + AMINO_ACIDS['W'] + AMINO_ACIDS['A'] + AMINO_ACIDS['R'] + WATER_MASS
        self.assertEqual(gen_spectra.get_precursor(self.sequence), Expected_prec + PROTON_MASS)
    
    def test_gen_spectrum(self):
        #Test the gen_spectrum function with the sequence, "MALWAR"
        Expected_spectrum = gen_spectra.b_ions(self.sequence) + gen_spectra.y_ions(self.sequence)
        precursor_mass = gen_spectra.get_precursor(self.sequence, 2) #gen_spectrum assumes the charge is 2 if the charge is not stated
        actual = gen_spectra.gen_spectrum(self.sequence)
        self.assertEqual(actual['spectrum'], Expected_spectrum)
        self.assertEqual(actual['precursor_mass'], precursor_mass)
    
    def test_optimized_compare_masses(self):
        #Test the scoring.optimized_compare_masses function
        #easy test case:
        self.assertEqual(mass_comparisons.optimized_compare_masses([1, 2, 4], [1, 3, 4], 1, False), 2)
        #testing tolerance:
        self.assertEqual(mass_comparisons.optimized_compare_masses([1, 2, 4], [1, 3, 4], 1000000, False), 3)
        #test the optimized compare masses function with some input spectra and a database of hybrids

    def test_score_sequence(self):
        #Test the scoring.score_sequence function
        #easy test case:
        self.assertEqual(scoring.score_sequence([1, 2, 4], [1, 3, 4], 1, False), 2)
        #testing tolerance:
        self.assertEqual(scoring.score_sequence([1, 2, 4], [1, 3, 4], 1000000, False), 3)

    def test_hybrid_score(self):
        #Test the scoring.hybrid_score function
        # say our b ions found are A, C, E
        # and y ions found are D, A
        # our scoring then works like
        # .5(bA) + .5(bC) + 1(bE) + .5 (yD) + 1(yA) 
        hybrid_seq = 'ABC-DEF'
        lesser_point = .5
        greater_point = 1.0
        #Building spectrum
        hits = []
        hits.append(max(gen_spectra.b_ions('A', 1))) #bA
        hits.append(max(gen_spectra.b_ions('ABC', 1))) #bC
        hits.append(max(gen_spectra.b_ions('ABCDE', 1))) #bE
        hits.append(max(gen_spectra.y_ions('DEF', 1))) #yD
        hits.append(max(gen_spectra.y_ions('ABCDEF', 1))) #yA
        sample_spectrum = Spectrum(hits)
        expected = 3.5
        actual = scoring.hybrid_score(sample_spectrum, hybrid_seq, 20, lesser_point, greater_point)
        self.assertEqual(actual, expected)

    def test_precursor_distance(self):
        #Test the scoring.precursor_distance function
        observed_precursor = 10
        reference_precursor = 15
        expected_precursor_dist = 5
        self.assertEqual(scoring.precursor_distance(observed_precursor, reference_precursor), expected_precursor_dist)
        #test the other way
        self.assertEqual(scoring.precursor_distance(reference_precursor, observed_precursor), expected_precursor_dist)

    def test_database_dependant_functions(self):
        # Build database
        dirname = os.path.dirname(__file__)
        spectra_file_paths = [os.path.join(dirname, '../data/spectra/hybrid_nod2e3.mzML')]  
        database_file_path = os.path.join(dirname, '../data/database/sample_database.fasta')
        db = database_preprocessing.database_and_spectra_preprocessing(spectra_file_paths, database_file_path)

        def test_digest_score(self):
            #Test the scoring.digest_score function
            #Testing with a non-hybrid
            self.assertEqual(scoring.digest_score('DAP', db, 'trypsin'), 2)
            #Testing with a cis-spliced hybrid
            self.assertEqual(scoring.digest_score('DAP-GSTMYPGIADR', db, 'trypsin'), 2)

        def test_gen_extensions(self):
            # test b extensions
            seq = 'DLQTLAL'
            ion = 'b'
            dirname = os.path.dirname(__file__)
            spectra_file_paths = [os.path.join(dirname, '../data/spectra/hybrid_nod2e3.mzML')]  
            database_file_path = os.path.join(dirname, '../data/database/sample_database.fasta') 
            # Generate spectra for some peptide. Note that we are only testing this function's ability to generate extensions
            # from the database so the input spectrum is not relevant. 
            spectrum = Spectrum(gen_spectra.gen_spectrum('DLQTLALWSRM'), [], 0, 0, -1, gen_spectra.get_precursor('DLQTLALWSRM'), 1, '', '', {})
            # Import database (mouse_filtered.fasta) in NOD2_E3 and do all usual preprocessing
            # Generate kmer set. This is going to be every kmer of size len(seq)

            # Generate all extensions
            b_ext = []
            
            b_ext += [x for x in alignment_utils.extend_non_hybrid(seq, spectrum, 'b', db)]
            
            # Check there are no missing extensions
            # These extensions were found manually
            b_test_ext = ['DLQTLALEVAQQK', 'DLQTLALEVARQK']
            self.assertEqual(sorted(b_ext), sorted(b_test_ext))
            #calculate extension length
            extension_len = utils.predicted_len_precursor(spectrum, seq) - len(seq)

            # test y extensions
            ion = 'y'
            # Generate spectra for some peptide. Note that we are only testing this function's ability to generate extensions
            # from the database so the input spectrum is not relevant. 
            spectrum = Spectrum(gen_spectra.gen_spectrum('DLQTLALWSRM'), [], 0, 0, -1, gen_spectra.get_precursor('DLQTLALWSRM'), 1, '', '', {})
            # Import database (mouse_filtered.fasta) in NOD2_E3 and do all usual preprocessing
            # db = database.build(database_file_path)
            # Generate kmer set. This is going to be every kmer of size len(seq)

            # Generate all extensions
            y_ext = []
            
            y_ext += [x for x in alignment_utils.extend_non_hybrid(seq, spectrum, 'y', db)]
            
            # Check there are no missing extensions
            # These extensions were found manually
            y_test_ext = ['GGPGAGDLQTLAL', 'LGGSPGDLQTLAL']
            self.assertEqual(sorted(y_ext), sorted(y_test_ext))
            #calculate extension length
            extension_len = utils.predicted_len_precursor(spectrum, seq) - len(seq)

    def test_CICD(self):
        self.assertEqual(utils.CICD_test(), 1)

    
    def test_runner(self):  
        arguments = self.get_arguments()
        matched_spectra = runner.run(arguments)
        actual = len(matched_spectra)
        expected = 3
        self.assertEqual(expected, actual, 'matched_spectra length is three')


if __name__ == "__main__":
    unittest.main()
