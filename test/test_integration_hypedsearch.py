import sys, os.path  
import unittest
import shutil

src_path = (os.path.abspath(os.path.join(os.path.dirname(__file__), '..')) + '/src/')
sys.path.append(src_path)
import runner, utils, database, gen_spectra, database_preprocessing

from alignment import alignment_utils
from objects import Spectrum
from constants import AMINO_ACIDS, SINGLY_CHARGED_B_BASE, DOUBLY_CHARGED_B_BASE, SINGLY_CHARGED_Y_BASE, DOUBLY_CHARGED_Y_BASE, PROTON_MASS, WATER_MASS

class Test_Main(unittest.TestCase):

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
        db = database_preprocessing.database_and_spectra_preprocessing(spectra_file_paths, database_file_path)
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



if __name__ == "__main__":
    unittest.main()