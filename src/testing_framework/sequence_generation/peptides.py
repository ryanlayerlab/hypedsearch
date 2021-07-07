from sequence_generation.digest import tryptic, random_digest
from math import ceil
from random import choice, randint
import numpy as np

digest_functions = {
    'random': random_digest,
    'trypsin': tryptic
}

def length_dist(num_peptides: int, dist='beta', min_length=4, max_length=35, a=2.3, b=7) -> list:
    '''
    Return a list of values that are the lengths peptides should be

    Inputs:
        num_peptides: int number of lengths to generate
    kwargs:
        dist: str the type of distribution to use. Default=beta (based on experimental data)
                other values: random
        min_length: int minium length peptide to generate. Default=4
        max_length: int maximum length peptide to generate. Default=35
        a: float alpha value. Default=2.3 (based on experimental data)
        b: float beta value. Default=7 (based on experimental data)
    Outputs:
        list of int lengths drawn from distribution
    '''
    if dist.lower() == 'beta':
        y = max_length
        x = min_length
        return x + (np.random.beta(a, b, size=(num_peptides)) * (y - x))
    else: # assume random
        return np.random.random_integers(min_length, high=max_length, size=num_peptides)

def __make_hybrid_pep(hybrid_prot: dict, min_length=4, max_length=20, dist='beta') -> dict:
    '''
    generate a hybrid peptide
    Inputs:
        hybrid_prot: a dictionary with entries 
            {
                protein: str,
                left_parent_end: int,
                name: str,
            }
    kwargs:
        min_length: int minimum length peptide to generate. Default=4
        max_length: int maximum length peptide to generate. Default=20
        dist: str name of distribution to pull from. Default='beta'
    Outputs:
        dictionary of form
        {
            left_parent_starting_position: int
            left_parent_ending_position: int (inclusive)
            right_parent_starting_position: int
            right_parent_ending_position: int (inclusive)
            left_parent_name: str
            right_parent_name: str
            sequence: str
        }
    '''
    
    l = length_dist(1, dist=dist, min_length=min_length, max_length=max_length)
    left_contr = randint(1, l-1)
    right_contr = l - left_contr
    j_site = hybrid_prot['left_parent_end_position']
    pep = hybrid_prot['protein'][j_site-left_contr + 1:j_site+right_contr + 1]
    return {
        'sequence': pep,
        'left_parent_starting_position': hybrid_prot['protein'].index(pep),
        'left_parent_ending_position': j_site,
        'right_parent_starting_position': hybrid_prot['right_parent_start_position'],
        'right_parent_ending_position': hybrid_prot['right_parent_start_position'] + right_contr - 1,
        'left_parent_name': hybrid_prot['left_parent_name'],
        'right_parent_name': hybrid_prot['right_parent_name']
    }

def __make_hybrid_peps_brute_force(hybrid_prot: dict, max_contribution=10, min_contribution=3, min_length=2) -> list:
    '''
    Generate every hybrid peptide from min length to 2*max_contribution including the hybrid junction site
    Inputs:
        hybrid_prot:    a dictionary with entries 
            {
                protein: str,
                left_parent_end: int,
                name: str,
            }
    kwargs:
        max_contribution:   int max contribution from each side. Default=10
        min_contribution:   int minimum contribution from each side. Default=3
        min_length:         mininum length peptide to create. Default=2
    Outputs:
        list of dictionarie of form
        {
            left_parent_starting_position: int
            left_parent_ending_position: int (inclusive)
            right_parent_starting_position: int
            right_parent_ending_position: int (inclusive)
            left_parent_name: str
            right_parent_name: str
            sequence: str
        }
    '''
    hyb_peps = []
    j_site = hybrid_prot['left_parent_end_position']

    for i in range(min_contribution, max_contribution):
        for j in range(min_contribution, max_contribution):
            pep = hybrid_prot['protein'][j_site - i + 1:j_site + j + 1]
            if len(pep) < min_length:
                continue
            hyb_peps.append({
                'sequence': pep,
                'hybrid_sequence': pep[:i] + '-' + pep[i:],
                'left_parent_starting_position': hybrid_prot['protein'].index(pep),
                'left_parent_ending_position': j_site,
                'right_parent_starting_position': hybrid_prot['right_parent_start_position'],
                'right_parent_ending_position': hybrid_prot['right_parent_start_position'] + j - 1,
                'left_parent_name': hybrid_prot['left_parent_name'],
                'right_parent_name': hybrid_prot['right_parent_name']
            })
    return hyb_peps

def __generate_hybrids(hybrid_prots: list, num_gen=10, peptide_name_prefix='HYBRID_PEPTIDE_', min_length=4, max_length=20) -> list:
    '''
    Generate hybrid peptides. Guaranteed to capture the junction

    Inputs:
        hybrid_prots: list of dictionaries. Should be of the form 
            [{
                left_parent_name: str,
                right_parent_name: str,
                left_parent_sequence: str,
                right_parent_sequence: str,
                left_parent_start: int,
                right_parent_start: int, 
                left_parent_contribution: int,
                right_parent_contribution: int,
                protein: str, 
                name: str
            }]
    kwargs:
        num_gen: int number of hybrid peptides to genereate. if num_gen < len(hybrid_prots), then 
            the remaining will be ignored. After one is generated from. Default=10
        peptide_name_prefix: str prefix to add as the name to each protein. Default=HYBRID_PEPTIDE
        min_length: int minimum length peptide to generate. Default=4
        max_length: int maximum length peptide to generate. Default=20
    Outputs:
        list of dictionaries of the form 
            {
                left_parent_starting_position: int
                left_parent_ending_position: int (inclusive)
                right_parent_starting_position: int
                right_parent_ending_position: int (inclusive)
                left_parent_name: str
                right_parent_name: str
                sequence: str
            }
    '''
    first_round = num_gen if num_gen <= len(hybrid_prots) else len(hybrid_prots)
    second_round = num_gen if num_gen > len(hybrid_prots) else 0
    hybrid_peps = []
    name_c = 0

    fill_zeros = len(str(num_gen))

    for i in range(first_round):
        print('Generating hybrid peptide {}/{}[{}%]\r'.format(name_c, num_gen, int(float(name_c)/float(num_gen) * 100)), end="")
        hyb_pep = __make_hybrid_pep(hybrid_prots[i], min_length=min_length, max_length=max_length)
        hyb_pep['peptide_name'] = peptide_name_prefix + str(name_c).zfill(fill_zeros)
        hybrid_peps.append(hyb_pep)
        name_c += 1

    if second_round == 0: 
        print('\nFinshed generating hybrid peptides')
        return hybrid_peps

    for _ in range(second_round):
        print('Generating hybrid peptide {}/{}[{}%]\r'.format(name_c, num_gen, int(float(name_c)/float(num_gen) * 100)), end="")
        hyb_pep = __make_hybrid_pep(choice(hybrid_prots), min_length=min_length, max_length=max_length)
        hyb_pep['peptide_name'] = peptide_name_prefix + str(name_c).zfill(fill_zeros)
        hybrid_peps.append(hyb_pep)
        name_c += 1

    print('\nFinished generating hybrid peptides')
    return hybrid_peps

def __generate_hybrids_brute_force(hybrid_prots: list, peptide_name_prefix='HYBRID_PEPTIDE', min_length=2, max_contribution=10, min_contribution=3) -> list:
    '''
    Generate all possible hybrid peptides from a window of a hybrid 

    Inputs:
        hybrid_prots: list of dictionaries. Should be of the form 
            [{
                left_parent_name: str,
                right_parent_name: str,
                left_parent_sequence: str,
                right_parent_sequence: str,
                left_parent_start: int,
                right_parent_start: int, 
                left_parent_contribution: int,
                right_parent_contribution: int,
                protein: str, 
                name: str
            }]
    kwargs:
        peptide_name_prefix: str prefix to add as the name to each protein. Default=HYBRID_PEPTIDE
        min_length: int minimum length peptide to generate. Default=2
        max_contribution: int maximum to allow each side of the hybrid to allow. Default=10
    Outputs:
        list of dictionaries of the form 
            {
                left_parent_starting_position: int
                left_parent_ending_position: int (inclusive)
                right_parent_starting_position: int
                right_parent_ending_position: int (inclusive)
                left_parent_name: str
                right_parent_name: str
                sequence: str
            }
    '''
    hybrid_peps = []
    name_c = 0

    for prot in hybrid_prots:
        hybrid_peps += __make_hybrid_peps_brute_force(prot, min_contribution=min_contribution, max_contribution=max_contribution, min_length=min_length)
    fill_zeros = len(str(len(hybrid_peps)))
    for pep in hybrid_peps:
        pep['peptide_name'] = peptide_name_prefix + str(name_c).zfill(fill_zeros)
        name_c += 1
    
    return hybrid_peps

############################################################################################################
#          END PRIVATE FUNCTIONS
############################################################################################################

def gen_peptides(proteins: list, n: int, min_length=3, max_length=20, digest='random', min_contribution=3, hybrid_list=False, dist='beta') -> list:
    '''
    Generates peptides from proteins given
    
    Inputs:
        proteins: list of dictionaries of the form {'name': str, 'sequence': str}
                NOTE: if hybrid_list set to true, form is 
                {
                    left_parent_name: str,
                    right_parent_name: str,
                    left_parent_sequence: str,
                    right_parent_sequence: str,
                    left_parent_start: int,
                    right_parent_start: int, 
                    left_parent_contribution: int,
                    right_parent_contribution: int,
                    protein: str, 
                    name: str
                }
        n: number of peptides to generate
    kwargs:
        min_length: int minimum length peptide to generate. Default=3
        max_length: int maximum length peptide to generate. Default=20
        digest: str type of digest to perform. Default=random
        hybrid_list: bool if the proteins passed in are hybrids and you wish to capture a junction point, set to True. Default=False
        dist: str name of distribution to use for peptide length. Default=beta (based on experimental data)
    Outputs:
        list of dictionaries of form 
        {
            left_parent_starting_position: int
            left_parent_ending_position: int (inclusive)
            right_parent_starting_position: int
            right_parent_ending_position: int (inclusive)
            left_parent_name: str
            right_parent_name: str
            sequence: str
        }
    '''
    digest = 'random' if digest.lower() not in digest_functions.keys() else digest.lower()
    if not hybrid_list:
        return digest_functions[digest](proteins, n, min_length=min_length, max_length=max_length, dist=dist)

    else:
        # return __generate_hybrids(proteins, num_gen=n, min_length=min_length, max_length=max_length)
        return __generate_hybrids_brute_force(proteins, min_length=min_length, min_contribution=min_contribution)


    