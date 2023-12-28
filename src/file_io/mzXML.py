from pyrsistent import s
from lookups.utils import file_exists
from lookups.objects import Spectrum
from preprocessing import spectra_filtering
from pyteomics import mzxml

def read(filename: str, peak_filter=0, relative_abundance_filter=0) -> list:
    if not file_exists(filename):
        print('File {} not found. Please make sure that this file exists'.format(filename))
        return

    spectra = []
    filecontents = mzxml.read(filename)
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
        precursor = content['precursorMz'][0]['precursorMz']
        precursor_charge = content['precursorMz'][0]['precursorCharge']
        id = content.get('id', '')
        ms_level = content['msLevel']
        scan_number = content['num']
        other_metadata = content['scanOrigin']
        retention_time = content['scan'][0]['scan start time']
        spectra.append(Spectrum(
            spec_num,
            masses,
            abundances,
            precursor,
            precursor_charge,
            filename, 
            id,
            other_metadata,
            retention_time
        ))

    return spectra