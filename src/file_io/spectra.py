from file_io import mzML, mzXML
from utils import file_exists

def load(filename: str, peak_filter: int = 0, relative_abundance_filter: float = 0.0) -> list:
    '''
    Reads one of the supported file types for spectra and returns a list of spectrum
    objects

    Inputs:
        filename:   (str) the name of the spectra file
    kwargs:
        peak_filter:                (int) the top number of peaks to keep. Default=0
        relative_abundance_filter:  (float) thee percentage of abundance a peak must have to be
                                    considered. Must be in range(0, 1) or the integer is converted
                                    to a decimal in that range. Default=0
    Outputs:
        (list) Spectrum namedtuple objects
    '''

    if not file_exists(filename):
        print(f'File {filename} not found. Please make sure that this file exists')

    # return based on file type
    ext = filename.split('.')[-1]

    if ext.lower() == 'mzxml':
        return mzXML.read(filename, peak_filter, relative_abundance_filter)

    elif ext.lower() == 'mzml':
        return mzML.read(filename, peak_filter, relative_abundance_filter)

    else:
        print(f'File {filename} is not of supported types (mzML, mzXML)')
        return []