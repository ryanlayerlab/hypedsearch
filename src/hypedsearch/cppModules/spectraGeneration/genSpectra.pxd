from libcpp.string cimport string 
from libcpp.vector cimport vector
from libcpp cimport bool

cdef extern from "genSpectra.hpp":
    vector[float] genSpectrum(string sequence, string ion, int charge, bool sort)