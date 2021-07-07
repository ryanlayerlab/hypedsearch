#include <unordered_map>
#include <vector>
#include <string>
#include <algorithm>

#ifndef GENSPECTRA_H
#define GENSPECTRA_H

static std::unordered_map<char, float> AMINO_ACIDS ({
    {'A', 71.037114}, 
    {'R', 156.101111},
    {'N', 114.042927},
    {'D', 115.026943},
    {'C', 103.009185},
    {'E', 129.042593},
    {'Q', 128.058578},
    {'G', 57.021464},
    {'H', 137.058912},
    {'I', 113.084064},
    {'L', 113.084064},
    {'K', 128.094963},
    {'M', 131.040485},
    {'F', 147.068414},
    {'P', 97.052764},
    {'S', 87.032028},
    {'T', 101.047679},
    {'U', 150.95363},
    {'W', 186.079313},
    {'Y', 163.06332},
    {'V', 99.068414},
    {'X', 0.0},           
    {'B', 113.084064},  
    {'Z', 0.0}           
});

//This is the mass of water. Adding the mass of water to the sum of all the residue masses gives the mass of the peptide.
const float WATER_MASS = 2 * 1.007825035 + 15.99491463;

const float SINGLY_CHARGED_Y_BASE = 3 * 1.007825035 + 15.99491463 - 0.0005486; //for the OH to turn the residue CO on the C-terminus into COOH + 1 proton to make NH into NH2 and 1 proton make positively charged
const float DOUBLY_CHARGED_Y_BASE = 4 * 1.007825035 + 15.99491463 - 2 * 0.0005486; //another proton to make doubly charged

const float SINGLY_CHARGED_B_BASE = 1.007825035 - 0.0005486; //for the H to turn the residue NH on the N-terminus into NH2
const float DOUBLY_CHARGED_B_BASE = 2 * 1.007825035 - 2 * 0.0005486; //adding one more proton this time to make it doubly charged

/**
 * Return the spectrum (floats) that represents the spectrum for b ions with a certain charge
 * 
 * @param sequence      string  the sequnce of amino acids 
 * @param charge        int     the charge to use for generating the seuqnce. If -1, both singly and doubly are returned
 * 
 * @return vector<float>        the spectrum
*/
std::vector<float> bIons(std::string sequence, int charge=-1);

/**
 * Return the spectrum (floats) that represents the spectrum for y ions with a certain charge
 * 
 * @param sequence      string  the sequnce of amino acids 
 * @param charge        int     the charge to use for generating the seuqnce. If -1, both singly and doubly are returned
 * 
 * @return vector<float>        the spectrum
*/
std::vector<float> yIons(std::string sequence, int charge=-1);

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
std::vector<float> genSpectrum(std::string sequence, std::string ion="", int charge=-1, bool sort=true);

#endif