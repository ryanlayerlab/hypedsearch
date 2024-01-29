from lookups.objects import Database, DatabaseEntry
from collections import defaultdict
import json, os

script_dir = os.path.dirname(__file__)
json_dir = '/'.join(script_dir.split('/')[:-1])
digest_file = os.path.join(json_dir, 'digests.json')
digests = json.load(open(digest_file, 'r'))

def digest(db: Database, digest_type: str, missed_cleavages: int) -> Database:
    if digest_type not in digests:
        return db
    digest_rules = digests[digest_type]
    starts = {s['amino_acid']: s['cut_position'] for s in digest_rules['start']}
    ends = {s['amino_acid']: s['cut_position'] for s in digest_rules['end']}
    new_prots = defaultdict(list)
    for p_name, entries in db.proteins.items():
        digested = []
        for entry in entries:
            for pos, aa in enumerate(entry.sequence):
                if aa in starts:
                    s = pos if starts[aa] == 'left' else pos + 1
                    allowed_misses = missed_cleavages
                    for j in range(pos, len(entry.sequence)):
                        if allowed_misses < 0:
                            e = j-1 if ends[entry.sequence[j-1]] == 'left' else j 
                            
                            digested.append((entry.sequence[s:e], s, e))
                            break
                        if j == len(entry.sequence) - 1:
                            digested.append((entry.sequence[s:], s, len(entry.sequence)))
                            break
                        if entry.sequence[j] in ends:
                            allowed_misses -= 1
        for d in digested:
            new_prots[p_name].append(DatabaseEntry(d[0], entry.description))
    db = db._replace(proteins=new_prots)
    return db

