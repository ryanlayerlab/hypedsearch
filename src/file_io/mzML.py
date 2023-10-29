from lookups.utils import file_exists
from lookups.objects import Spectrum
from preprocessing import spectra_filtering
from pyteomics import mzml
import yaml

def read(filename: str, peak_filter=0, relative_abundance_filter=0) -> list:
    if not file_exists(filename):
        print('File {} not found. Please make sure that this file exists'.format(filename))
        return

    spectra = []
    filecontents = mzml.read(filename)
    content: dict
    for spec_num, content in enumerate(filecontents):
        masses = list(content['m/z array'])
        abundances = list(content['intensity array'])
        if peak_filter > 0:
            masses, abundances = spectra_filtering.peak_filtering(masses, abundances, peak_filter)
        elif relative_abundance_filter > 0:
            while relative_abundance_filter > 1:
                relative_abundance_filter /= 100
            masses, abundances = spectra_filtering.relative_abundance_filtering(masses, abundances, relative_abundance_filter)
        ti = sum(abundances)
        precursor = None
        precursor_charge = 0

        if not len(content['precursorList']['precursor']) or not len(content['precursorList']['precursor'][0]['selectedIonList']['selectedIon']):
            precursor = max(masses)
            precursor_charge = 1

        else:
            precursor = float(content['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['selected ion m/z'])
            precursor_abundance = float(content['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['peak intensity'])
            precursor_charge = int(content['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['charge state'])
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