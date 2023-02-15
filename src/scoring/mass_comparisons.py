from gen_spectra import calc_masses
from utils import ppm_to_da, hashable_boundaries

from bisect import bisect


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
        
    #hack
    #the_type = type(reference) #python = 'dict' #cpp = 'list'
    updated_reference = reference      
    if isinstance(updated_reference,dict):
        reference = updated_reference.get('spectrum')
    
    # local variables for score
    return_value = 0
    for ref in reference:
        if bisect(observed_boundaries, ref) % 2:
            return_value = return_value + 1
    # return_value = sum([1 for ref in reference if bisect(observed_boundaries, ref) % 2])
    return return_value
    


    
