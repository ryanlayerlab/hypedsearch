from gen_spectra import calc_masses
from utils import ppm_to_da

from bisect import bisect

# def cmp_spectra_spectra__JAN_2020(spec: list, reference: list) -> float:
#     '''
#     Score two spectra against eachother
#     Simple additive scoring algorithm with a length divider of the smaller sequence
    
#     Inputs:
#         spec: list of floats mass spectra of first sequence
#         reference: list of floats mass spectra of second sequence
#     Outputs:
#         float score from comparison
#     '''
#     if len(spec) == 0 or len(reference) == 0:
#         return
#     streak = 0
#     last = False
#     score = 0
#     max_streak = 0
#     for mass in spec:
#         if last == True:
#             streak += 1
#             max_streak = max([streak, max_streak])

#         if mass in reference:
#             score += 1
#             last = True

#         else:
#             streak = 0
#             last = False
    
#     score += max_streak
#     divider = min([len(spec), len(reference)]) / 2
#     return score / divider


# def compare_masses__FEB_2020(spectrum: list, reference: list) -> float:
#     '''
#     CREATED FEB 26 2020
#     Score two spectra against eachother. Simple additive scoring with bonuses for streaks
#     Divides by the length of the reference to make it length biased for the reference

#     Inputs:
#         spectrum:   list of floats (from mass spectra)
#         reference:  list of floats (calculated from protein sequence)
#     Outputs:
#         score:      float score 
#     '''
#     if len(spectrum) == 0 or len(reference) == 0:
#         return
#     streak = 0
#     last = False
#     score = 0
#     max_streak = 0
#     for mass in spectrum:
#         if last == True:
#             streak += 1
#             max_streak = max([streak, max_streak])

#         if mass in reference:
#             score += 1
#             last = True

#         else:
#             streak = 0
#             last = False
    
#     score += max_streak
#     score /= (len(reference) / 2)
#     return score 

# def compare_masses__APRIL_2020(spectrum: list, reference: list) -> float:
#     '''
#     CREATED APRIL 6 2020
#     Score two spectra against eachother. Simple additive scoring with bonuses for streaks
#     Divides by the length of the reference to make it length biased for the reference

#     Note:   the difference between this one and the February one is which spectrum
#             is being iterated through. This one iterates through the reference first

#     Inputs:
#         spectrum:   list of floats (from mass spectra)
#         reference:  list of floats (calculated from protein sequence)
#     Outputs:
#         score:      float score 
#     '''
#     if len(spectrum) == 0 or len(reference) == 0:
#         return 0.0
#     streak = 0
#     last = True
#     score = 0
#     max_streak = 0
#     for refmass in reference:
    
#         if refmass in spectrum:
#             if last == True:
#                 streak += 1
#                 max_streak = max([streak, max_streak])
#             score += 1
#             last = True 

#         else:
#             streak = 0
#             last = False
    
#     score += max_streak
#     score /= (float(len(reference)) / 2)
#     return score 

# def compare_masses(spectrum: list, reference: list) -> float:
#     '''
#     CREATED APRIL 27 2020
#     Score two spectra against eachother. Simple additive scoring with bonuses for streaks
#     Divides by the length of the reference to make it length biased for the reference

#     Note:   the difference between this one and the other April one is this one divides 
#             by the length of the spectrum. This is because all extended kmers will 
#             be getting longer, so we want the score to increase with the length, not stay the same

#     Inputs:
#         spectrum:   list of floats (from mass spectra)
#         reference:  list of floats (calculated from protein sequence)
#     Outputs:
#         score:      float score 
#     '''
#     if len(spectrum) == 0 or len(reference) == 0:
#         return 0.0
#     streak = 0
#     last = True
#     score = 0
#     max_streak = 0
#     for refmass in reference:
    
#         if refmass in spectrum:
#             if last == True:
#                 streak += 1
#                 max_streak = max([streak, max_streak])
#             score += 1
#             last = True 

#         else:
#             streak = 0
#             last = False
    
#     score += max_streak
#     score /= float(len(spectrum))
#     return score 

# def compare_spectra_sequence_ion_type(spectra: list, reference: str, ion: str) -> float:
#     '''
#     MARCH 11 2020
#     Generate a score by the comparison of a list of masses against a reference sequences
#     Additive scoring divided by the reference length

#     Inputs:
#         spectra:    list of masses from an mzml file
#         reference:  string of amino acids to compare spectra to
#         ion:        string ion type to compare. Possilbe are {'b', 'y'}
#     Outputs: 
#         float score of the comparision
#     '''
#     reference_ions, _ = calc_masses(reference, ion=ion)
#     return compare_masses(spectra, reference_ions)

# def compare_sequence_sequence_ion_type(spectra: str, reference: str, ion: str) -> float: 
#     '''
#     MARCH 11 2020
#     Generate a score by the comparison of two sequences
#     Additive scoring divided by the refernce length
    
#     Inputs:
#         spectra:   string amino acid sequence in question
#         reference: string reference amino acid sequence to compare to 
#         ion:       string ion type to compare. Possible types are {'b', 'y'}
#     Ouputs:
#         float score of the comparison
#     '''
#     spectra_ions, _ = calc_masses(spectra, ion=ion)
#     reference_ions , _= calc_masses(reference, ion=ion)
#     return compare_masses(spectra_ions, reference_ions)

# def optimized_compare_masses__MAY(observed: list, reference: list, ppm_tolerance=20, needs_sorted=False) -> float:
#     '''
#     CREATED MAY 19 2020
#     Score two spectra against eachother. Simple additive scoring with bonuses for streaks
#     Divides by the length of the reference to make it length biased for the reference

#     Note:   the difference between this one and the April 27 one is this one attempts
#             to be more optimized in terms of search complexity. Uses binary search
#             for the faster search. NOTE: observed should be sorted before passed in.
#             If not, pass True in the needs_sorted flag

#     Inputs:
#         observed:       (list of floats) spectrum being scored
#         reference:      (list of floats) reference spectrum to score observed againsts
#     kwargs:
#         ppm_tolerance:  (float) the mass error tolerance allowed. Default=20
#         needs_sorted:   (bool) the observed mass needs to be sorted before binary search. Default=False
#     Outputs:
#         score:      float score 
#     '''
#     if len(observed) == 0 or len(reference) == 0:
#         return 0.0

#     if needs_sorted:
#         observed.sort()
        
#     def boundaries(mass):
#         tol = ppm_to_da(mass, ppm_tolerance)
#         return [mass - tol, mass + tol]
                
#     # calculate the boundaries for each of the reference masses for binary search
#     observed_boundaries = []
#     for obs in observed:
#         observed_boundaries += boundaries(obs)
        
#     # local variables for score
#     streak = 0
#     last = True
#     score = 0
#     max_streak = 0
    
#     for ref in reference:
#         # see if an observed is in the spectrum +/- the ppm by binary searchy
#         found = bisect(observed_boundaries, ref) % 2
    
#         # increment score if found
#         if found:
#             if last == True:
#                 streak += 1
#                 max_streak = max([streak, max_streak])
#             score += 1
#             last = True 

#         else:
#             streak = 0
#             last = False
    
#     score += max_streak
#     return score

#   CREATED JULY 1 2020
def optimized_compare_masses(
    observed: list, 
    reference: list, 
    ppm_tolerance: int = 20, 
    needs_sorted: bool = False
    ) -> float:
    '''Score two spectra against eachother. Simple additive scoring of ions found

    :param observed: observed set of m/z values
    :type observed: list
    :param reference: reference set of m/z values
    :type reference: list
    :param ppm_tolerance: parts per million mass error allowed when matching masses. 
        (default is 20)
    :type ppm_tolerance: int
    :param needs_sorted: Set to true if either the observed or reference need to 
        be sorted. 
        (default is False)
    :type needs_sorted: bool

    :returns: the number of matched ions
    :rtype: int

    :Example:

    >>> optimized_compare_masses([1, 2, 4], [1, 3, 4], 1, False)
    >>> 2
    '''
    if len(observed) == 0 or len(reference) == 0:
        return 0.0

    if needs_sorted:
        observed.sort()
        reference.sort()
        
    def boundaries(mass):
        tol = ppm_to_da(mass, ppm_tolerance)
        return [mass - tol, mass + tol]
                
    # calculate the boundaries for each of the reference masses for binary search
    observed_boundaries = []
    for obs in observed:
        observed_boundaries += boundaries(obs)
        
    # local variables for score
    return sum([1 for ref in reference if bisect(observed_boundaries, ref) % 2])
    