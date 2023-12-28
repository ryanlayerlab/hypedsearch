from computational_pipeline.gen_spectra import calc_masses
from lookups.utils import ppm_to_da, hashable_boundaries

from bisect import bisect

def optimized_compare_masses(
    observed: list, 
    reference: list, 
    ppm_tolerance: int = 20, 
    needs_sorted: bool = False
    ) -> float:
    if len(observed) == 0 or len(reference) == 0:
        return 0.0

    if needs_sorted:
        observed.sort()
        reference.sort()
        
    def boundaries(mass):
        tol = ppm_to_da(mass, ppm_tolerance)
        return [mass - tol, mass + tol]
                
    observed_boundaries = []
    for obs in observed:
        observed_boundaries += boundaries(obs)
        
    updated_reference = reference      
    if isinstance(updated_reference,dict):
        reference = updated_reference.get('spectrum')
    
    return_value = 0
    for ref in reference:
        if bisect(observed_boundaries, ref) % 2:
            return_value = return_value + 1
    return return_value
    


    
