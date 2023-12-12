import json
from lookups.utils import make_valid_json_file

def save_dict(file: str, d: dict) -> None:
    file = make_valid_json_file(file)
    if not isinstance(d, dict):
        raise Exception(f'{d} should be type dict. Type of {d}: {type(d)}')
    with open(file, 'w') as json_file:
        json.dump(d, json_file)