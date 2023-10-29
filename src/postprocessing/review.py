from src.constants.objects import SequenceAlignment, Alignments

def __is_swap_up_to_dist(a: str, b: str, i: int, j: int, dist: int, d: list) -> list:
    '''
    Helper function to edit distance long swaps. Helps identify potential swaps
    and if two chars could be swapped then the entry to d is returned
    '''
    
    # first index not worth looking at
    if i == 0 or j == 0:
        return []
    
    # iterate up to our current index or the max swap distance
    iter_len = min([i, j, dist])
    
    # keep track of swaps we find
    swaps = []
    
    for k in range(1, iter_len + 1):
        
        # if it is a swap then keep track of it
        if a[i] == b[j - k] and a[i - k] == b[j]:
            swaps.append((i-k, j-k))
        
    return swaps

def __edit_distance_long_swaps(a: str, b: str, dist: int = 0) -> int:
    '''
    Find the edit distance between two strings allowing for swaps up to a distance dist. 
    
    Example:
        a: 'ABCDE', b: 'ADCBE'
        
        edit distance with swap dist of 0: 2 (substitution of B for D and D for B)
        edit distance with swap dist of 1: 1 (swap of B and D)
        
    Limitations: if a position has 2 swaps, it will not be picked up. Example:
        a: 'CBEDA', b: 'ABCDE'
        
        A has been swapped with C then with E, but the final output will be edit distance of 3
        for all of the swaps
        
    Inputs:
        a:   (str) the first string
        b:   (str) the second string
        dist:(int) the swapping distance allowed. Default=0
        
    Outputs:
        (int) the minimum edit distance
    '''
    
    d = [[0 for _ in range(len(b))] for _ in range(len(a))]
    
    for i in range(len(a)):
        d[i][0] = i
        
    for j in range(len(b)):
        d[0][1] = j
        
    for i in range(len(a)):
        for j in range(len(b)):
            
            # look for swaps
            swaps = __is_swap_up_to_dist(a, b, i, j, dist, d)
            
            if a[i] == b[j]:
                d[i][j] = d[i-1][j-1]
                
            elif len(swaps):
                
                # get all swaps possible 
                swap_values = [d[x][y] for x, y in swaps]
                                
                d[i][j] = min(swap_values + [
                    d[i-1][j] + 1,   # deletion
                    d[i][j-1] + 1,   # substitution
                ])
                
            else:
                d[i][j] = min([
                    d[i-1][j] + 1,  # deletion
                    d[i][j-1] + 1,  # insertion
                    d[i-1][j-1] + 1,# substitution
                ])
    
    return d[len(a)-1][len(b)-1]

# def prioritize_non_hybrids(results: dict) -> dict:
#     '''
#     Look through all of the results and try to prioritize non hybrid results

#     Inputs:
#         results:   (dict) a value is a list of SequenceAlignments
#     Outputs:
#         (dict) updated sorted 
#     '''
#     for _id, sas in results.itesm():

#         # sas is a list of sequence alignments

def __digest_score(sequence: str, digest_type: str) -> int:
    '''
    The additional points a sequence gets if it follows the digest rule. 


    '''
    pass

def tie_breaker(results: dict, digest_type: str, n: int) -> dict:
    '''
    Look through all the results and try to break as many ties as possible

    Inputs:
        results:        (dict) values are Alignments
        digest_type:    (str) the digest type used
        n:              (int) return only the top n results
    Outputs:
        (dict) updated results
    '''

    for _id, alignments in results.items():

        # go through every sequence and see if they follow the digest rule
        # hard code for test rn 

        new_sas = []
        for sa in alignments.alignments:
            
            if sa.sequence[0] == 'D':
                new_sas.append(sa._replace(total_score=sa.total_score + 1))

            else:
                new_sas.append(sa)

        a = Alignments(
            alignments.spectrum, 
            sorted(new_sas, key=lambda x: (x.total_score, 1/x.precursor_distance), reverse=True)[:n]
        )

        results[_id] = a

    return results