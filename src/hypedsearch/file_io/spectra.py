from file_io import mzML, mzXML
from utils import file_exists

def load(filename: str, peak_filter: int = 0, relative_abundance_filter: float = 0.0) -> list:
    if not file_exists(filename):
        print(f'File {filename} not found. Please make sure that this file exists')
    ext = filename.split('.')[-1]
    if ext.lower() == 'mzxml':
        return mzXML.read(filename, peak_filter, relative_abundance_filter)
    elif ext.lower() == 'mzml':
        return mzML.read(filename, peak_filter, relative_abundance_filter)
    else:
        print(f'File {filename} is not of supported types (mzML, mzXML)')
        return []