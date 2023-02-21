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
        def get_arguments(self):
            dirname = os.path.dirname(__file__)
            spectra_file_paths = [os.path.join(dirname, '../data/spectra/NOD2_E3/hybrid_nod2e3.mzML')]  
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
    
        def test_runner(self):  
            arguments = self.get_arguments()
            matched_spectra = runner.run(arguments)
            actual = len(matched_spectra)
            expected = 3
            self.assertEqual(expected, actual, 'matched_spectra length is three')


if __name__ == "__main__":
    unittest.main()
