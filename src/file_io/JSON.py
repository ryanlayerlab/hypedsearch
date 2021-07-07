import json
from utils import make_valid_json_file

def save_dict(file: str, d: dict) -> None:
    '''
    Save a dictionary in a file

    Inputs:
        file:       string a file path to save all the data in 
        d:     a dictionary of scores to save
    Outputs:
        None
    '''
    file = make_valid_json_file(file)
    if type(d) != dict:
        raise Exception('d should be type dict. Type of d: {}'.format(type(d)))
    with open(file, 'w') as o:
        json.dump(d, o)