# from lookups.utils import make_dir, make_valid_dir_string
# from file_io import JSON
# from lookups.objects import Alignments

# import pandas as pd
# import json
# import os

# SUMMARY_NAME = 'summary'
# HYBRID_PREFIX = 'hybrid_'

# def json_file(results: dict, output_dir: str) -> None:
#     json_file_name = os.path.join(output_dir, f'{SUMMARY_NAME}.json')
#     dictified = {}
#     for name, alignment in results.items():
#         dictified[name] = {
#             'spectrum': alignment.spectrum._asdict(), 
#             'alignments': [x._asdict() for x in alignment.alignments]
#         }
#     JSON.save_dict(json_file_name, dictified)
    



# def tsv_file(results: dict, output_dir: str) -> None:
#     mac = 0
#     hybrids, nonhybrids = [], []
#     alignment: Alignments
#     for name, alignment in results.items():
#         if len(alignment.alignments) == 0:
#             mac += 1
#             continue
#         topalignment = alignment.alignments[0]._asdict()
#         topalignment['entry name'] = name
#         topalignment['id'] = alignment.spectrum.id
#         if 'hybrid_sequence' in topalignment:
#             hybrids.append(topalignment)
#         else:
#             nonhybrids.append(topalignment)
#     hybridresults = pd.DataFrame(hybrids)
#     with open(f'{output_dir + HYBRID_PREFIX + SUMMARY_NAME}.tsv', 'w') as ho:
#         ho.write(hybridresults.to_csv(sep='\t'))
#     del hybridresults
#     del hybrids
#     nonhybridresults = pd.DataFrame(nonhybrids)
#     output_file = os.path.join(output_dir, f'{SUMMARY_NAME}.tsv')
#     with open(output_file, 'w') as nho:
#         nho.write(nonhybridresults.to_csv(sep='\t'))

# def generate(alignments: dict, output_dir='./') -> None:
#     output_dir = make_valid_dir_string(output_dir)
#     make_dir(output_dir)
#     json_file(alignments, output_dir)
#     tsv_file(alignments, output_dir)


