import sys, os.path  
import unittest
src_dir = (os.path.abspath(os.path.join(os.path.dirname(__file__), '..')) + '/src/')
sys.path.append(src_dir)
from src import runner, utils, database, runner, objects


class Test_Main(unittest.TestCase):
    def setUp(self):
        sys.path.insert(0, "../src")

    def get_arguments(self):
        dirname = os.path.dirname(__file__)
        spectra_file_paths = [os.path.join(dirname, '../data/spectra/hybrid_nod2e3.mzML')]  
        database_file_path = os.path.join(dirname, '../data/database/sample_database.fasta')
        database_file = src.database.build(database_file_path)
        output_dir = '../output'
        min_peptide_len = 3
        max_peptide_len = 30
        ppm_tolerance = 20
        precursor_tolerance = 10
        verbose = False
        peak_filter = 25
        relative_abundance_filter = .01
        digest = ''
        debug = True
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
        matched_spectra = src.runner.run(arguments)
        actual = len(matched_spectra)
        expected = 3
        self.assertEqual(expected, actual, 'matched_spectra length is three')

if __name__ == "__main__":
    unittest.main()
