from scoring import mass_comparisons
from objects import Spectrum, Database
from utils import ppm_to_da

import gen_spectra
import utils
import database

import re

# load digests json for digest scoring
import json 
import os

script_dir = os.path.dirname(__file__)
json_dir = '/'.join(script_dir.split('/')[:-1])
digest_file = os.path.join(json_dir, 'digests.json')

digests = json.load(open(digest_file, 'r'))

# def score_subsequence(pepspec: list, subseq: str, ppm_tolerance=20) -> (float, float):
#     '''
#     Score a mass spectrum to a substring of tagged amino acids

#     Inputs:
#         pepspec:        (list of floats) the mass spectrum to score
#         subseq:         (str) amino acids to score the spectrum against
#     kwargs:
#         ppm_tolerance:  (int) the tolerance to accepted while scoring. Default=20
#     Outputs:
#         (b_score, y_score): (float, float) the b and y ion scores generated from this comparison
#     '''
#     kmerspec_b = gen_spectra.gen_spectrum(subseq, ion='b')['spectrum']
#     kmerspec_y = gen_spectra.gen_spectrum(subseq, ion='y')['spectrum']
#     b_score = mass_comparisons.optimized_compare_masses(pepspec, kmerspec_b, ppm_tolerance=ppm_tolerance)
#     y_score = mass_comparisons.optimized_compare_masses(pepspec, kmerspec_y, ppm_tolerance=ppm_tolerance)
#     return (b_score, y_score)

def score_sequence(
    observed: list, 
    theoretical: list, 
    ppm_tolerance: int = 20, 
    needs_sorted: bool = False
    ) -> float:
    '''Score a mass spectrum to a substring of tagged amino acids

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

    >>> score_sequence([1, 2, 4], [1, 3, 4], 1, False)
    >>> 2
    '''

    return mass_comparisons.optimized_compare_masses(
        observed, 
        theoretical, 
        ppm_tolerance=ppm_tolerance, 
        needs_sorted=needs_sorted
    )


# def backbone_score(observed: Spectrum, reference: str, ppm_tolerance: int) -> int:
#     '''
#     Scoring algorithm based on backbone coverage of the reference. The scoring algorithm 
#     returns a number between 0 and 100/len(observed) + 3*(len(reference)-1). The calculation of the score is as follows:
    
#         1. A percentage is given for the number of bond sites successfully identified
#         2. For each bond site that has > 1 ion that describes it, an extra point is awarded. 
        
#     Example:
#         reference:   ABCDE, 4 junctions to describe
#         observed:    ions: b1+, y1++, y2+, b4+
        
#         ions for A, E, DE, ABCD found. Coverage = A**DE with DE described by both E and ABCD
        
#         Score = %(3/4) + 1 = 75/len(observed) + 1 = 4.75 
        
#     Inputs:
#         observed:       (Spectrum) spectrum being scored against
#         reference:      (str) reference amino acid sequence being scored against the spectrum
#         ppm_tolerance:  (int) tolerance to allow in ppm for each peak
#     Outputs:
#         (int) score according the the function described above
#     '''
#     if len(reference) < 2:
#         return 0

#     jcount = [0 for _ in range(len(reference)-1)]
    
#     for ion in ['b', 'y']:
#         for charge in [1, 2]:
#             singled_seq = reference[:-1] if ion == 'b' else reference[1:]
#             peaks = gen_spectra.gen_spectrum(singled_seq, charge=charge, ion=ion)['spectrum']
#             peaks = peaks if ion == 'b' else peaks[::-1]
#             for i in range(len(peaks)):
#                 da_tol = ppm_to_da(peaks[i], ppm_tolerance)
#                 if any([peaks[i] - da_tol <= obs_peak <= peaks[i] + da_tol for obs_peak in observed.spectrum]):
#                     jcount[i] += 1
    
#     jcoverage = int(100 * sum([1 if jc > 0 else 0 for jc in jcount]) / len(jcount))
#     # extrapoints = sum([jc - 1 if jc > 1 else 0 for jc in jcount])
    
#     # make it a function of the number of observed peaks matched
#     return jcoverage# + extrapoints

# def intensity_backbone_score(observed: Spectrum, reference: str, ppm_tolerance: int) -> int:
#     '''
#     Scoring algorithm based on backbone coverage of the reference. The scoring algorithm 
#     returns a number between 0 and 100/len(observed) + 3*(len(reference)-1). The calculation of the score is as follows:
    
#         1. A percentage is given for the number of bond sites successfully identified
#         2. For each bond site that has > 1 ion that describes it, an extra point is awarded. 
        
#     Example:
#         reference:   ABCDE, 4 junctions to describe
#         observed:    ions: b1+, y1++, y2+, b4+
        
#         ions for A, E, DE, ABCD found. Coverage = A**DE with DE described by both E and ABCD
        
#         Score = %(3/4) + 1 = 75/len(observed) + 1 = 4.75 
        
#     Inputs:
#         observed:       (Spectrum) spectrum being scored against
#         reference:      (str) reference amino acid sequence being scored against the spectrum
#         ppm_tolerance:  (int) tolerance to allow in ppm for each peak
#     Outputs:
#         (int) score according the the function described above
#     '''
#     if len(reference) < 2:
#         return 0

#     jcount = [0 for _ in range(len(reference)-1)]

#     # keep track of the abundances that contribute to our score
#     ided_abundances = 0
    
#     for ion in ['b', 'y']:
#         for charge in [1, 2]:
#             # take off the trailing or leading amino acid from the reference according to ion type
#             singled_seq = reference[:-1] if ion == 'b' else reference[1:]

#             # get the m/z peaks
#             peaks = gen_spectra.gen_spectrum(singled_seq, charge=charge, ion=ion)['spectrum']
            
#             # go through each peak and try and match it to an observed one
#             for i in range(len(peaks)):

#                 # take tolerance into account 
#                 da_tol = ppm_to_da(peaks[i], ppm_tolerance)

#                 # get hits
#                 peak_hits = list(map(lambda idx_x: idx_x[0] if peaks[i] - da_tol <= idx_x[1] <= peaks[i] + da_tol else None, enumerate(observed.spectrum)))

#                 # remove None from peak hits
#                 peak_hits = [idx for idx in peak_hits if idx is not None]

#                 # if the len > 1, increment
#                 if len(peak_hits) > 0:
#                     jcount[i] += 1
#                     ided_abundances += sum([observed.abundance[idx] for idx in peak_hits])
    
#     jcoverage = sum([1 if jc > 0 else 0 for jc in jcount])
#     #extrapoints = sum([jc - 1 if jc > 1 else 0 for jc in jcount]) // 3
    
#     # make it a function of the number of observed peaks matched
#     #jcoverage /= len(observed.spectrum)
#     return jcoverage #+ extrapoints) * (ided_abundances / observed.total_intensity)


# def ion_backbone_score(observed: Spectrum, reference: str, ion: str, ppm_tolerance: int) -> float:
#     '''
#     Scoring algorithm based on backbone coverage of the reference. The scoring algorithm 
#     returns a number between 0 and (100/+ 3*(len(reference)-1))/len(observed) . The calculation of the score is as follows:
    
#         1. A percentage is given for the number of bond sites successfully identified
#         2. For each bond site that has > 1 ion that describes it, an extra point is awarded. 
        
#     Example:
#         reference:   ABCDE, 4 junctions to describe
#         observed:    ions: b1+, y1++, y2+, b4+
        
#         ions for A, E, DE, ABCD found. Coverage = A**DE with DE described by both E and ABCD
        
#         Score = %(3/4) + 1 = 75/len(observed) + 1 = 4.75 
        
#     Inputs:
#         observed:       (Spectrum) spectrum being scored against
#         reference:      (str) reference amino acid sequence being scored against the spectrum
#         ion:            (str) the ion type to focus on. Options are 'b' or 'y'
#         ppm_tolerance:  (int) tolerance to allow in ppm for each peak
#     Outputs:
#         (float) score according to the formula 
#     '''
#     jcount = [0 for _ in range(len(reference)-1)]
    
#     for charge in [1, 2]:
#         singled_seq = reference[:-1] if ion == 'b' else reference[1:]
#         peaks = gen_spectra.gen_spectrum(singled_seq, charge=charge, ion=ion)['spectrum']
#         peaks = peaks if ion == 'b' else peaks[::-1]
#         for i in range(len(peaks)):
#             da_tol = ppm_to_da(peaks[i], ppm_tolerance)
#             if any([peaks[i] - da_tol <= obs_peak <= peaks[i] + da_tol for obs_peak in observed.spectrum]):
#                 jcount[i] += 1
#     divider = 1 if len(jcount) <= 1 else len(jcount)
#     jcoverage = int(100 * sum([1 if jc > 0 else 0 for jc in jcount]) / divider)
#     #extrapoints = sum([jc - 1 if jc > 1 else 0 for jc in jcount])
    
#     return jcoverage #+ extrapoints

# def intensity_ion_backbone_score(observed: Spectrum, reference: str, ion: str, ppm_tolerance: int) -> float:
#     '''
#     Scoring algorithm that factors in how much of the backbone is identified along with the abundance each peak 
#     contributes to the total score. The formula is:

#         count(# backbone cleavages found) + 1 for every additional hit of an already identifed cleavage * total percentage
#         of intensity covered by identfied peaks

#     Example:
#         reference:   ABCDE, 4 junctions to describe
#         observed:    ions: b1+, y1++, y2+, b4+, with relative intensities of (100, 200, 100, 50) of a total 600
        
#         ions for A, E, DE, ABCD found. Coverage = A**DE with DE described by both E and ABCD
        
#         Score = %(3/4) + 1 = 75 + 1 = 76 *(100 + 200 + 100 + 50)/600 = 57 

#     Inputs:
#         observed:       (Spectrum) spectrum being scored against
#         reference:      (str) reference amino acid sequence being scored against the spectrum
#         ion:            (str) the ion type to focus on. Options are 'b' or 'y'
#         ppm_tolerance:  (int) tolerance to allow in ppm for each peak
#     Outputs:
#         (float) score according to the formula 
#     '''
#     # check to see if observed is nothing
#     if len(observed.spectrum) == 0:
#         return 0

#     # keep track of the junction (bond) sites found
#     jcount = [0 for _ in range(len(reference)-1)]

#     # keep track of the abundances that contribute to our score
#     ided_abundances = 0
    
#     for charge in [1, 2]:

#         # take off the trailing or leading amino acid from the reference according to ion type
#         singled_seq = reference[:-1] if ion == 'b' else reference[1:]

#         # get the m/z peaks
#         peaks = gen_spectra.gen_spectrum(singled_seq, charge=charge, ion=ion)['spectrum']
        
#         # go through each peak and try and match it to an observed one
#         for i in range(len(peaks)):

#             # take tolerance into account 
#             da_tol = ppm_to_da(peaks[i], ppm_tolerance)

#             # get hits
#             peak_hits = list(map(lambda idx_x: idx_x[0] if peaks[i] - da_tol <= idx_x[1] <= peaks[i] + da_tol else None, enumerate(observed.spectrum)))

#             # remove None from peak hits
#             peak_hits = [idx for idx in peak_hits if idx is not None]

#             # if the len > 1, increment
#             if len(peak_hits) > 0:
#                 jcount[i] += 1
#                 ided_abundances += sum([observed.abundance[idx] for idx in peak_hits])

#     # if any entry has at least 1, that bond has been identified at least once
#     jcoverage = sum([1 if jc > 0 else 0 for jc in jcount])

#     # if an entry has more than 1, we give it extra points
#     #extrapoints = sum([jc - 1 if jc > 1 else 0 for jc in jcount]) // 3
    
#     return jcoverage# + extrapoints) * (ided_abundances / observed.total_intensity)

# def ion_intensity_percentage(observed: Spectrum, reference: str, ppm_tolerance: int, ion: str) -> int:
#     '''
#     The percentage the matched m/z values make up of the overall intensity. For example,
#     if 4 peaks are matched between the observed and reference and those peaks make up 10% of the
#     total intensity, 10 is returned. Spectrum m/z values should be sorted low to high with the intensity
#     values indexed in the same order as the m/z values

#     Inputs:
#         observed:       (Spectrum) observed spectrum
#         reference:      (str) the reference sequence
#         ppm_tolerance:  (int) the tolerance (in ppm) to allow when matching peaks
#     Output:
#         (int) rounded percentage value
#     '''
#     if len(reference) == 0:
#         return 0

#     # generate the spectrum for the reference sequence
#     refspec = gen_spectra.gen_spectrum(reference, ion=ion)['spectrum']
        
#     def boundaries(mass):
#         tol = ppm_to_da(mass, ppm_tolerance)
#         return [mass - tol, mass + tol]
                
#     # calculate the boundaries for each of the reference masses for binary search
#     reference_boundaries = []
#     for mz in refspec:
#         reference_boundaries += boundaries(mz)
        
#     # get all indices where the observed gets a hit
#     idxes = [i for i, mz in enumerate(observed.spectrum) if bisect(reference_boundaries, mz) % 2]

#     # take the intensities at idxes and sum them
#     score = int( 100 * sum([observed.abundance[i] for i in idxes]) / observed.total_intensity)
    
#     return score


def hybrid_score(
    observed: Spectrum, 
    hybrid_seq: str, 
    ppm_tolerance: int, 
    lesser_point: float = .5, 
    greater_point: float = 1.0
    ) -> float:
    '''A score for hybrid sequences. b ions found to the left of the hybrid 
    junction and y ions found to the right of the hybrid junctions will be 
    rewarded a point of value *lesser_point*. b ions found to the right of the 
    hybrid junction and y ions found to the left of hybrid junction will be 
    awarded a point of value *greater_point*.

    :param observed: observed spectrum
    :type observed: Spectrum
    :param hybrid_seq: hybrid string sequence
    :type hybrid_seq: str
    :param ppm_tolerance: mass error allowed in parts per million when matching masses 
    :type ppm_tolerance: int
    :param lesser_point: point awarded to ions found on their respective side of 
        the hybrid junction. 
        (default is .5)
    :type lesser_point: float
    :param greater_point: point awarded to ions found on their non respective side 
        of the hybrid junction. 
        (default is 1.0)
    :type greater_point: float

    :returns: the score 
    :rtype: float 

    :Example: 

    >>> hybrid_seq = 'ABC-DEF'
    >>> lesser_point = .5
    >>> greater_point = 1.0
    >>> # say our b ions found are A, C, E
    >>> # and y ions found are D, A
    >>> # our scoring then works like
    >>> # .5(bA) + .5(bC) + 1(bE) + .5 (yD) + 1(yA) 
    >>> hybrid_score(spectrum, hybrid_seq, 20, lesser_point, greater_point)
    >>> 3.5
    '''

    if '-' not in hybrid_seq and '(' not in hybrid_seq and ')' not in hybrid_seq:
        return 0

    score = 0

    # get a non hybrid for scoring purposes
    non_hyb = hybrid_seq.replace('-', '').replace('(', '').replace(')', '')

    # get the index to the left of which b ions will only get .5 points
    b_split = hybrid_seq.index('-') if '-' in hybrid_seq else hybrid_seq.index('(')

    # get the index to the right of which y ions will only get .5 points
    y_split = len(hybrid_seq) - hybrid_seq.index('-') if '-' in hybrid_seq else len(hybrid_seq) - hybrid_seq.index(')') 

    # generate b and y separately to be sure
    b_spec = sorted(gen_spectra.gen_spectrum(non_hyb, ion='b')['spectrum'])
    y_spec = sorted(gen_spectra.gen_spectrum(non_hyb, ion='y')['spectrum'])

    # convert the spectra into lists of tuples
    gen_range = lambda x: (x - ppm_to_da(x, ppm_tolerance), x + ppm_to_da(x, ppm_tolerance))
    b_ranges = [gen_range(x) for x in b_spec]
    y_ranges = [gen_range(x) for x in y_spec]

    # do a merge search where we linearly search each mass in the observed twice
    b_range_i = 0
    observed_i = 0
    while b_range_i < len(b_ranges) and observed_i < len(observed.spectrum):

        # if observed is larger than the range, increment range
        if observed.spectrum[observed_i] > b_ranges[b_range_i][1]:
            b_range_i += 1

        # if observed is smaller than the range, increment observed
        elif observed.spectrum[observed_i] < b_ranges[b_range_i][0]:
            observed_i += 1

        # otherwise its in the range, see what to increment score by, and increment observed
        else:
            score += 1 if b_range_i >= b_split else .5
            observed_i += 1

    y_range_i = 0
    observed_i = 0
    while y_range_i < len(y_ranges) and observed_i < len(observed.spectrum):

        # if observed is larger than the range, increment range
        if observed.spectrum[observed_i] > y_ranges[y_range_i][1]:
            y_range_i += 1

        # if observed is smaller than the range, increment observed
        elif observed.spectrum[observed_i] < y_ranges[y_range_i][0]:
            observed_i += 1

        # otherwise its in the range, see what to increment score by, and increment observed
        else:
            score += 1 if y_range_i >= y_split else .5
            observed_i += 1

    return score

def precursor_distance(observed_precursor: float, reference_precursor: float) -> float:
    '''The absolute distance between the observed precursor and reference precursor

    :param observed_precursor: the observed precursor mass
    :type observed_precursor: float
    :param reference_precursor: the precursor mass of the reference sequence
    :type reference_precursor: float

    :returns: the absolute value of the difference between the two
    :rtype: float
    '''

    return abs(observed_precursor - reference_precursor)


# def xcorr(observed: np.ndarray, reference: np.ndarray, value=50) -> float:
#     '''
#     An xcorr scoring algorithm 
    
#     Formula: dot(reference, y'), y' = observed - (sum(-75, 75)observed[i]/150)
    
#     Inputs should be sparsely populated np arrays. Indexing should observe the formula
    
#     idx = int(m/w), m is mass, w is bin width (aka tolerance)
    
#     Inputs:
#         observed:    (np.ndarray) list of peaks normalized to value
#         reference:   (np.ndarray) list of peaks normalized to value
#     kwargs: 
#         value:       (number) value given to the normalized peaks. Default=50
#     Outputs:
#         (float) the score of the formula shown above
#     '''
#     def sum_term(i):
#         min_idx = max(0, i-75)
#         max_idx = min(len(observed) - 1, i + 75)
#         return np.sum(observed[min_idx:max_idx])/(max_idx - min_idx)

#     # calculate the y prime term 
#     y_prime = np.asarray([observed[i] - sum_term(i) for i in range(len(observed))])
        
#     # fill y prime or reference to be the same size. fill with zeros
#     if len(y_prime) < len(reference):
#         y_prime = np.concatenate((y_prime, [0 for _ in range(len(reference) - len(y_prime))]))
#     else:
#         reference = np.concatenate((reference, [0 for _ in range(len(y_prime) - len(reference))]))
        
#     return np.dot(reference, y_prime)/(sum(reference)*value)


# def overlap_score(l: list, ion: str) -> float:
#     '''
#     A score based on how well the sequences that make up the list overlap

#     Example:   
#         l: [ABC, ABCDE, ABCDLMNO, ABCDEF]

#         + 3 for ABC overlapping all of them
#         + 1 for ABCDE overlapping ABCDEF
#         - .5 for ABCDEF not overlapping ABCDLMNO
#         Final score: 3.5

#     Inputs:
#         l:      (list) sequences checking for overlap
#         ion:    (str) ion type. If y, strings compared right to left
#     Outputs:
#         (float) score of the overlap
#     '''
#     l.sort(key=lambda x: len(x))
    
#     score = 0
    
#     for i, e in enumerate(l[:-1]):
        
#         # if its y, we want to go from right to left
#         if ion == 'y':
#             e = e[::-1]
        
#         for j, e2 in enumerate(l[i+1:]):
            
#             # if its y, we want to go from right to left
#             if ion == 'y':
#                 e2 = e2[::-1]
            
#             score += 1 if e2[:len(e)] == e else -.5
            
#     return score

def total_mass_error(observed: Spectrum, alignment: str, tolerance: int) -> float:
    '''The sum of all of the mass errors for every matched mass between the 
    observed and the alignment.
    
    :param observed: observed spectrum
    :type observed: Spectrum
    :param alignment: the string alignment
    :type alignment: str
    :param tolerance: parts per million tolerance allowed when matching masses
    :type tolerance: int

    :returns: sum of the absolute values of all mass errors
    :rtype: float
    '''

    # clean the input alignment
    sequence = alignment.replace('-', '').replace(')', '').replace('(', '')

    # generate the spectrum
    alignment_spectrum = gen_spectra.gen_spectrum(sequence)['spectrum']

    # sort them both
    sorted_observed = sorted(observed.spectrum)
    sorted_alignment = sorted(alignment_spectrum)

    # i is for observed, j for the str alignment
    i, j = 0, 0

    # keep track of total error
    total_error = 0

    while i < len(sorted_observed) and j < len(sorted_alignment):

        # see if the mass at j is +- the mass at i
        da_tol = ppm_to_da(sorted_alignment[j], tolerance)

        # if alignment < observed - tolerance, increment alignment
        if sorted_alignment[j] < sorted_observed[i] - da_tol:
            j += 1

        # if alignment > observed + tolerance, increment observed
        elif sorted_alignment[j] > sorted_observed[i] + da_tol:
            i += 1

        # finally ad tot total error and increment both
        else:
            total_error += abs(sorted_alignment[j] - sorted_observed[i])
            i+= 1
            j += 1

    return total_error


def digest_score(sequence: str, db: Database, digest_type: str) -> int:
    '''The additional points *sequence* gets if it follows the digest rules of 
    the specified digest type

    :param sequence: hybrid or non hybrid sequence to analyze
    :type sequence: str
    :param db: source proteins
    :type db: Database
    :param digest_type: what kind of digest was performed
    :type digest_type: str

    :returns: additional points the sequence gets by following the digest
    :rtype: int
    '''

    # if the digest type is not in digests, return 0
    if digest_type not in digests:
        return 0

    digest = digests[digest_type]

    # if left digest cuts left and sequence [0] matches, give it a point
    left_point = 1 \
        if any(
            [x['amino_acid'] == sequence[0] for x in digest['start'] \
            and x['cut_position'] == 'left']
        ) else 0

    # if right digest cuts right and sequence[-1] matches, give it a point
    right_point = 1 \
        if any(
            [x['amino_acid'] == sequence[-1] for x in digest['end'] \
            and x['cut_position'] == 'right']
        ) else 0

    if left_point + right_point == 2:
        return 2

    # check to see if its a hybrid sequence
    if utils.HYBRID_ALIGNMENT_PATTERN.findall(sequence):

        # get the left and right halves
        left, right = utils.__split_hybrid(sequence)

        # well first check to see if we can assign a point to left
        # before even looking at the source proteins. If we just look 
        # at the first amino acid and it follows the digest rule, bam we 
        # golden
        left_point = 1 \
            if left_point == 1 or any(
                [x['amino_acid'] == left[0] for x in digest['start'] \
                and x['cut_position'] == 'left']
            ) else 0

        # do the same for the right
        right_point = 1 \
            if right_point == 1 or any(
                [x['amino_acid'] == sequence[-1] for x in digest['end'] \
                and x['cut_position'] == 'right']  
            ) else 0

        if left_point == 0:

            # find source proteins
            left_source_proteins = database.get_proteins_with_subsequence_ion(
                db, left, 'b'
            )

            # look to see if there is any source protein where the amino acid
            # to the left of our amino acid follows the digest rule
            for lsp in left_source_proteins:

                for entry in database.get_entry_by_name(db, lsp):

                    # if entry.sequence at index of left has amino acid 
                    # to left, add point
                    indices = [m.start() for m in re.finditer(left, entry.sequence)]

                    for i in indices:

                        # see if any of the right cut start match sequence[i-1]
                        if i == 0:
                            continue

                        left_point = 1 if any(
                            [x['amino_acid'] == entry.sequence[i-1] for x in \
                            digest['start'] and x['cut_position'] == 'right']
                        ) else 0 

                        if left_point == 1:
                            break 

                    if left_point == 1:
                        break 

                if left_point == 1: 
                    break 

        if right_point == 0:

            # get source proteins
            right_source_proteins = database.get_proteins_with_subsequence_ion(
                db, right, 'y'
            )

            for rsp in right_source_proteins:

                for entry in database.get_entry_by_name(db, rsp):

                    indices = [m.start() + len(right) for \
                        m in re.finditer(right, entry.sequence)]

                    for i in indices:
                        
                        if i == len(entry.sequence) - 1:
                            continue

                        right_point = 1 if any(
                            [x['amino_acid'] == entry.sequence[i] for x in \
                            digest['end'] and x['cut_position'] == 'left']
                        ) else 0

                        if right_point == 1:
                            break 

                    if right_point == 1:
                        break 

                if right_point == 1:
                    break 

        return right_point + left_point

    else:
        # find source proteins
        source_proteins = database.get_proteins_with_subsequence(
            db, sequence
        )

        for sp in source_proteins:

            for entry in database.get_entry_by_name(db, sp):

                # if entry.sequence at index of left has amino acid 
                # to left, add point
                indices = [m.start() for m in re.finditer(sequence, entry.sequence)]

                for i in indices:

                    if i == 0:
                        continue
                    
                    # make left point 1 if already 1 or if the amino acid 
                    # to the left of our sequence follows a right cut in the 
                    # digest
                    left_point = 1 if left_point == 1 or any(
                        [x['amino_acid'] == entry.sequence[i-1] for x in \
                        digest['start'] and x['cut_position'] == 'right']
                    ) else 0 

                    # make right point 1 if already 1 or if the amino acid
                    # to the right of our sequence follows a left cut in 
                    # the digest
                    right_point = 1 if right_point == 1 or any(
                        [x['amino_acid'] == entry.sequence[i+len(sequence)] \
                        for x in digest['end'] and x['cut_position'] == 'left']
                    ) else 0

        return left_point + right_point
