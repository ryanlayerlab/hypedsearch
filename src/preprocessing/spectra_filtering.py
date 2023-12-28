def relative_abundance_filtering(
    masses: list, 
    abundances: list, 
    percentage: float
    ) -> (list, list):
    ti = sum(abundances)
    min_value = ti * percentage
    filtered_mass_abundances = [x for x in zip(masses, abundances) if x[1] >= min_value]
    masses = [float(x) for x, _ in filtered_mass_abundances]
    abundances = [float(x) for _, x in filtered_mass_abundances]
    return (masses, abundances)


def peak_filtering(masses: list, abundances: list, num_peaks: int) -> (list, list):
    mass_abundances = zip(masses, abundances)
    mass_abundances = sorted(mass_abundances, key=lambda x: x[1], reverse=True)[:num_peaks]
    mass_abundances.sort(key=lambda x: x[0])
    masses = [float(x) for x, _ in mass_abundances]
    abundances = [float(x) for _, x in mass_abundances]
    return (masses, abundances)