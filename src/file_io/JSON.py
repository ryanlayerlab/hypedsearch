import json
from lookups.utils import make_valid_json_file

def save_dict(file: str, d: dict) -> None:
    file = make_valid_json_file(file)
    if type(d) != dict:
        raise Exception('d should be type dict. Type of d: {}'.format(type(d)))
    with open(file, 'w') as o:
        json.dump(d, o)