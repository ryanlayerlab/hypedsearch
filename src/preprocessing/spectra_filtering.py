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

    total_intensity = sum(abundances)

    # find the filter value
    min_value = total_intensity * percentage

    # zip them up and take values that pass the filter
    filtered_mass_abundances = [x for x in zip(masses, abundances) if x[1] >= min_value]

    # split them off and return
    masses, abundances = zip(*filtered_mass_abundances)
    masses, abundances = map(float, masses), map(float, abundances)

    return list(masses), list(abundances)


def peak_filtering(masses: list, abundances: list, num_peaks: int) -> (list, list):
    mass_abundances = zip(masses, abundances)
    mass_abundances = sorted(mass_abundances, key=lambda x: x[1], reverse=True)[:num_peaks]
    mass_abundances.sort(key=lambda x: x[0])

    # seperate them
    masses, abundances = zip(*mass_abundances)
    masses, abundances = map(float, masses), map(float, abundances)

    return list(masses), list(abundances)
