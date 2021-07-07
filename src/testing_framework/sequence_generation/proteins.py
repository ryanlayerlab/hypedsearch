from random import randint

def generate_hybrids(prots, num_gen, min_contribution=10, name_prefix='HYBRID_'):
    '''
    Generate a number of parent proteins with certain constraints

    Inputs:
        num_gen: int number of hybrid_proteins to generate
        prots: list of dictionaries of the form
            [{
                'name': str,
                'sequence': str
            }]
    kwargs:
        min_contribution: int minimum number of AAs to use from each parent. Default=10
        name_prefix: str prefix to add before every name of protein. Default=HYBRID_
    Outputs:
        list of dictionaries. Form is
        [{
            left_parent_name: str,
            right_parent_name: str,
            left_parent_sequence: str,
            right_parent_sequence: str,
            left_parent_end: int,
            right_parent_start: int, 
            left_parent_contribution: int,
            right_parent_contribution: int,
            protein: str, 
            name: str, 
        }]
    '''
    hybrids = []
    fill_zeros = len(str(num_gen))

    for i in range(num_gen):
        print('Generating hybrid protein {}/{}[{}%]\r'.format(i, num_gen, int(float(i)/float(num_gen) * 100)), end="")
        left_parent = prots[randint(0, len(prots)-1)]
        right_parent = prots[randint(0, len(prots)-1)]
        # make sure both parents are long enough
        while len(left_parent['sequence']) < min_contribution:
            left_parent = prots[randint(0, len(prots)-1)]
        while len(right_parent['sequence']) < min_contribution:
            right_parent = prots[randint(0, len(prots)-1)]

        left_end = randint(0, len(left_parent['sequence'])-1)
        right_start = randint(0, len(right_parent['sequence'])-1)
        while left_end < min_contribution:
            left_end = randint(0, len(left_parent['sequence'])-1)
        while len(right_parent['sequence']) - right_start < min_contribution:
            right_start = randint(0, len(right_parent['sequence'])-1)
        hybrid = left_parent['sequence'][:left_end] + right_parent['sequence'][right_start:]
        name = name_prefix + str(i).zfill(fill_zeros)
        hybrid_d = {
            'left_parent_name': left_parent['name'],
            'right_parent_name': right_parent['name'],
            'left_parent_sequence': left_parent['sequence'],
            'right_parent_sequence': right_parent['sequence'],
            'left_parent_end_position': left_end - 1,        # inclusive index of the end
            'right_parent_start_position': right_start,
            'left_parent_contribution': left_end,
            'right_parent_contribution': len(right_parent['sequence']) - right_start,
            'protein': hybrid,
            'name': name
        }
        hybrids.append(hybrid_d)
    print('\nFinished generating hybrid proteins')
    return hybrids

