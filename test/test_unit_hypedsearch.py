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
    def setUp(self): 
        self.sequence = 'MALWAR'

    def test_ppm_to_da(self):
        #Tests the convert ppm to da function with a mass of 100 and a tolerance of 20. 20 ppm of 100 should be .002
        mass = 100
        tol = 20
        self.assertEqual(utils.ppm_to_da(mass, tol), .002, '20 ppm of 100 should be .002')

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

    def test_overlap_intervals(self):
        #Run the utils.overlap_intervals function with two different intervals. One set will overlap and one won't
        intervals1 = [[0,3], [2,5], [3,7]] #Expected to return [0,7]
        intervals2 = [[-1,4], [6,15]] #Expected to return itself
        self.assertEqual(utils.overlap_intervals(intervals1), [[0,7]])
        self.assertEqual(utils.overlap_intervals(intervals2), intervals2)
    
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
    
    def test_predicted_len_precursor(self): 
        #Run the predicted_len_precursor function with the sequence 'MAL' and the spectrum for 'MALWAR'. 
        sequence = 'MAL'
        spectrum = Spectrum(gen_spectra.gen_spectra('MALWAR'), [], 0, 0, -1, gen_spectra.get_precursor('MALWAR'), 1)
        expected_length = 7 #Note that while "MALWAR" has a length of 6, the calculation is intentially rounded up because
        #it is better to overshoot than undershoot
        self.assertEqual(utils.predicted_len_precursor(spectrum, sequence), expected_length)
    
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



if __name__ == "__main__":
    unittest.main()
