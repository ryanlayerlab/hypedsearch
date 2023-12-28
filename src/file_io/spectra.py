from file_io import mzML, mzXML
from lookups.utils import file_exists

def load(filename: str, number_peaks: int = 0, relative_abundance: float = 0.0) -> list:
    if not file_exists(filename):
        print(f'File {filename} not found. Please make sure that this file exists')
    ext = filename.split('.')[-1]

    if ext.lower() == 'mzxml':
        mzml_file = mzXML.read(filename, number_peaks, relative_abundance)
        return mzml_file

    elif ext.lower() == 'mzml':
        return mzML.read(filename, number_peaks, relative_abundance)

    else:
        print(f'File {filename} is not of supported types (mzML, mzXML)')
        return []