from file_io import spectra

def load_spectra(
    spectra_files: list, 
    ppm_tol: int, 
    peak_filter: int = 0, 
    relative_abundance_filter: float = 0.0
    ) -> (list, list, dict):
    '''Load all the spectra files into memory and merge all spectra into one 
    massive list for reduction of the search space

    :param spectra_files: full string paths the the spectra files
    :type spectra_files: list
    :param ppm_tol: parts per million mass error allowed for making boundaries 
    :type ppm_tol: int
    :param peak_filter: the top X most abundant spectra to keep. If left as 0, 
        *relative_abundance_filter* is used instead. 
        (default is 0)
    :type peak_filter: int
    :param relative_abundance_filter: the percentage of the total abundance a 
        peak must make up in order to pass the filter. Value should be between 
        [0, 1). A realistic value is .005 (.5%). If *peak_filter* is non-zero, 
        that value is used instead. 
        (default is 0.0)
    :type relative_abundance_filter: float
  

    :returns: Spectra objects from file, overlapped boundaries of [lower_bound, upper_bound], 
        mapping from a m/z value to the index of the boundaries that the m/z fits in
    :rtype: (list, list, dict)
    '''
    
    # the single dimension of all the rounded spectra to use for reducing the search space
    linear_spectra = []

    # a list of all boundaries namedtuples
    all_spectra = []

    # go through each spectra file and load them into memory
    for spectra_file in spectra_files:
        these_spectra = spectra.load(
            spectra_file, 
            peak_filter=peak_filter, 
            relative_abundance_filter=relative_abundance_filter
        )

        all_spectra += these_spectra

        # go through each mass of each s, load it into memory, 
        # round the numbers to 3 decimal places for easier, and append to linear_spectra
        linear_spectra += list(set([
            x for spectrum in these_spectra for x in spectrum.spectrum
        ]))

    # sort the linear spectra
    linear_spectra.sort()

    # turn the all spectra list into a list of boundaries
    def make_boundaries(mz):
        da_tol = ppm_to_da(mz, ppm_tol)
        return [mz - da_tol, mz + da_tol]

    boundaries = [make_boundaries(mz) for mz in linear_spectra]

    # make overlapped boundaries larger boundaries
    boundaries = overlap_intervals(boundaries)

    # make a mapping for mz -> boundaries
    b_i, s_i = 0, 0
    mz_mapping = {}
    while s_i < len(linear_spectra):
        
        # if the current boundary encapsulates s_i, add to list
        if boundaries[b_i][0] <= linear_spectra[s_i] <= boundaries[b_i][1]:
            mz_mapping[linear_spectra[s_i]] = b_i 
            s_i += 1

        # else if the s_i < boundary, increment s_i
        elif linear_spectra[s_i] < boundaries[b_i][0]:
            s_i += 1

        # else if s_i > boundary, incrment b_i
        elif linear_spectra[s_i] > boundaries[b_i][1]:
            b_i += 1

    return (all_spectra, boundaries, mz_mapping)