from pyopenms import MSExperiment, MSSpectrum, MzMLFile, Peak1D, Precursor

'''write_mzml

DESC:
    create a mass spectrum file in mzml of sequences
Inputs:
    file_name: str name to save the file in
    spectra: list of dictionaries of the form [{spectrum: list[floats], precursor_mass: float}]
             This data is written to file (the spectrum and the precursor)
kwargs:
    title_prefix: str name to give as prefix to the name of each spectrum. Default=Spectrum
    output_dir: str name of the directory to save files to. Default=./
Outputs:
    list of strings of file paths
'''
def write_mzml(file_name, spectra, title_prefix='Spectrum ', output_dir='./'):
    if '.mzml' not in file_name.lower():
        file_name += '.mzML'

    output_file = output_dir + file_name

    exp = MSExperiment()
    sp_count = 0
    for spectrum in spectra:
        spec = MSSpectrum()
        spec.setMSLevel(2)
        name = str.encode(title_prefix + str(sp_count))
        spec.setName(name)
        sp_count += 1
        
        i = [500 for _ in spectrum['spectrum']]
        spec.set_peaks([spectrum['spectrum'], i])
        spec.setMSLevel(2)
        prec = Precursor()
        prec.setCharge(2)
        prec.setMZ(spectrum['precursor_mass'])
        spec.setPrecursors([prec])
        spec.sortByPosition()
        exp.addSpectrum(spec)

    MzMLFile().store(output_file, exp)

    return output_file