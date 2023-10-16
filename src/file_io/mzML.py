from utils import file_exists
from objects import Spectrum
from preprocessing import spectra_filtering
from pyteomics import mzml

# def prettyPrint(d, indent=0):
#    for key, value in d.items():
#         print('\t' * indent + str(key))
#         if isinstance(value, dict):
#             pretty(value, indent+1)
#         else:
#             print('\t' * (indent+1) + str(value))

def read(filename: str, peak_filter=0, relative_abundance_filter=0) -> list:
    '''
    read an .mzML file into memory. Filter out the peaks by the type specified.
    If both filter types are set to 0, all peaks are returned, otherwise filtering 
    will happen. If both filters are given a value, peak_filter is given preference.

    Inputs:
        filename:       (str) path to the file to import
    kwargs:
        peak_filter:                (int) the top number of peaks to keep. Default=0
        relative_abundance_filter:  (float) thee percentage of abundance a peak must have to be
                                    considered. Must be in range(0, 1) or the integer is converted
                                    to a decimal in that range. Default=0
    Outputs:
        (list) Spectrum namedtuple instances
    '''
    if not file_exists(filename):
        print('File {} not found. Please make sure that this file exists'.format(filename))
        return

    spectra = [] 
    filecontents = mzml.read(filename)

    for spec_num, content in enumerate(filecontents):
        # prettyPrint(content) # testing
        # print(yaml.dump(content))

        masses = list(content['m/z array'])
        abundances = list(content['intensity array'])

        # peak filter if the number is > 0
        if peak_filter > 0:
            masses, abundances = spectra_filtering.peak_filtering(masses, abundances, peak_filter)

        # if peak filter is not set and relative abundances is, filter by that
        elif relative_abundance_filter > 0:

            # if an integer is given, make it a float in the range (0, 1)
            while relative_abundance_filter > 1:
                relative_abundance_filter /= 100

            masses, abundances = spectra_filtering.relative_abundance_filtering(masses, abundances, relative_abundance_filter)

        # get the total intensity
        ti = sum(abundances)

        # get the precursor and its charge
        # we will assume its the first entry in the list
        precursor = None
        precursor_charge = 0
        precursor_abundance = 0

        if not len(content['precursorList']['precursor']) or not len(content['precursorList']['precursor'][0]['selectedIonList']['selectedIon']):
            precursor = max(masses)
            precursor_charge = 1

        else:
            selected_ion = content['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]
            precursor = float(selected_ion.get('selected ion m/z', 0))
            precursor_charge = int(selected_ion.get('charge state', 0))
            precursor_abundance = float(selected_ion.get('peak intensity', 0))

        # get the id
        id = content.get('id', '')

        retention_time = content['scanList']['scan'][0]['scan start time']
        
        spectra.append(
            Spectrum(
            spec_num,
            masses,
            abundances,
            precursor,
            precursor_charge,
            filename, 
            id,
            retention_time,
            precursor_abundance
        ))

    return spectra