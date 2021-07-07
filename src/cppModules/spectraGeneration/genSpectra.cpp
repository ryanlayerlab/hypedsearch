#include "genSpectra.hpp"

/**
 * Return the spectrum (floats) that represents the spectrum for b ions with a certain charge
 * 
 * @param sequence      string  the sequnce of amino acids 
 * @param charge        int     the charge to use for generating the seuqnce. If -1, both singly and doubly are returned
 * 
 * @return vector<float>        the spectrum
*/
std::vector<float> bIons(std::string sequence, int charge){
    std::vector<float> masses;

    if (charge == 1 || charge == -1){
        float total = SINGLY_CHARGED_B_BASE;
        for (char aa : sequence){
            total += AMINO_ACIDS[aa];
            masses.push_back(total);
        }   
    }

    if (charge == 2 || charge == -1){
        float total = DOUBLY_CHARGED_B_BASE;
        for (char aa : sequence){
            total += AMINO_ACIDS[aa];
            masses.push_back(total / 2.0);
        }
    }

    return masses;
} 

/**
 * Return the spectrum (floats) that represents the spectrum for y ions with a certain charge
 * 
 * @param sequence      string  the sequnce of amino acids 
 * @param charge        int     the charge to use for generating the seuqnce. If -1, both singly and doubly are returned
 * 
 * @return vector<float>        the spectrum
*/
std::vector<float> yIons(std::string sequence, int charge){
    std::vector<float> masses;

    if (charge == 1 || charge == -1){
        float total = SINGLY_CHARGED_Y_BASE;
        for (char aa : sequence){
            total += AMINO_ACIDS[aa];
            masses.push_back(total);
        }   
    }

    if (charge == 2 || charge == -1){
        float total = DOUBLY_CHARGED_Y_BASE;
        for (char aa : sequence){
            total += AMINO_ACIDS[aa];
            masses.push_back(total / 2.0);
        }
    }

    return masses;
} 


/**
 * Generate the spectrum of an amino acids 
 * 
 * @param sequence      string  the sequnce of an amino acids
 * @param ion           string  the ion to generate the sequnce for. Either "b", "y", or "" for both
 * @param charge        int     the charge to generate the spectrum for. Either 1, 2 or -1 for both
 * @param sort          bool    whether or not to sort the spectrum
 * 
 * @return vector<float>        the spectrum
*/
std::vector<float> genSpectrum(std::string sequence, std::string ion, int charge, bool sort){
    std::vector<float> spectrum; 
    if (ion == "b" || ion == ""){
        auto bSpec = bIons(sequence, charge);
        for (float mass: bSpec) spectrum.push_back(mass);
    }

    if (ion == "y" || ion == ""){
        auto ySpec = yIons(sequence, charge);
        for (float mass: ySpec) spectrum.push_back(mass);
    }

    if (sort){
        std::sort(spectrum.begin(), spectrum.end());
    }

    return spectrum;
}