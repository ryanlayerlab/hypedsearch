def relative_abundance_filtering(
    masses: list, 
    abundances: list, 
    percentage: float
    ) -> (list, list):
    '''Take all peaks from the spectrum who's abundance is at least *percentage* 
    of the total abundances. It is assumed that the masses and abundances lists 
    share ordering

    :param masses: m/z values 
    :type masses: list
    :param abundances: abundance value for the m/z values. Abundance at entry 
        *i* corresponds to m/z valuat entry *i*
    :type abundances: list
    :param percentage: the minimum percentage of the total abundance a peak must
        have to pass the filter. Values are in the range [0, 1). A relatively 
        realistic value is .005 (.5%)
    :type percentage: float

    :returns: filtered masses, filtered abundaces
    :rtype: (list, list)
    '''

    # total intensity
    ti = sum(abundances)

    # find the filter value
    min_value = ti * percentage

    # zip them up and take values that pass the filter
    filtered_mass_abundances = [x for x in zip(masses, abundances) if x[1] >= min_value]

    # split them off and return
    masses = [float(x) for x, _ in filtered_mass_abundances]
    abundances = [float(x) for _, x in filtered_mass_abundances]

    return (masses, abundances)


def peak_filtering(masses: list, abundances: list, num_peaks: int) -> (list, list):
    '''Take the most abundant peaks and return the sorted masses with the abundances.
    It is assumed that the masses and abundances lists share ordering


    :param masses: m/z values 
    :type masses: list
    :param abundances: abundance value for the m/z values. Abundance at entry 
        *i* corresponds to m/z valuat entry *i*
    :type abundances: list
    :param num_peaks: the top X most abundant peaks 
    :type num_peask: int

    :returns: filtered masses, filtered abundaces
    :rtype: (list, list)
    '''

    # zip the abundance and the m/z values together
    mass_abundances = zip(masses, abundances)
    
    # sort by key 1, the abundance, and take the top peak filter results
    mass_abundances = sorted(mass_abundances, key=lambda x: x[1], reverse=True)[:num_peaks]

    # sort them now by the value of m/z
    mass_abundances.sort(key=lambda x: x[0])

    # seperate them
    masses = [float(x) for x, _ in mass_abundances]
    abundances = [float(x) for _, x in mass_abundances]

    return (masses, abundances)