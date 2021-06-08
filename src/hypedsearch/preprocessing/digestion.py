from objects import Database, DatabaseEntry
from collections import defaultdict

import json 

import os

script_dir = os.path.dirname(__file__)
json_dir = '/'.join(script_dir.split('/')[:-1])
digest_file = os.path.join(json_dir, 'digests.json')

digests = json.load(open(digest_file, 'r'))

def digest(db: Database, digest_type: str, missed_cleavages: int) -> Database:
    '''
    Digest each protein in the database. If no digest is done, then 
    the original database is returned. 
    NOTE: 
    The entires in the database after digestion are the names of the form
    <protein_name>_<start_position>_<end_position>

    Inputs:
        db:                 (Database) the input source
        digest_type:        (str) the digestion to perform
        missed_cleavages:   (int) the number of missed cleavages allowed
    Outputs:
        (Database) updated protein entries
    '''
    if digest_type not in digests:
        return db
    
    digest_rules = digests[digest_type]
    starts = {s['amino_acid']: s['cut_position'] for s in digest_rules['start']}
    ends = {s['amino_acid']: s['cut_position'] for s in digest_rules['end']}
    
    new_prots = defaultdict(list)
    
    for p_name, entries in db.proteins.items():
        
        # keep track of what gets digested
        digested = []

        for entry in entries:
        
            for pos, aa in enumerate(entry.sequence):
                
                if aa in starts:
                    
                    # get the starting position for this cut based on rule
                    s = pos if starts[aa] == 'left' else pos + 1
                    
                    allowed_misses = missed_cleavages
                    
                    # find all of the next ends. we will keep track of them for up to missed_cleavages
                    for j in range(pos, len(entry.sequence)):
                        
                        # if we're out of missed cleavages, break
                        if allowed_misses < 0:
                            e = j-1 if ends[entry.sequence[j-1]] == 'left' else j 
                            
                            digested.append((entry.sequence[s:e], s, e))
                            break
                        
                        # check if we're at the end
                        if j == len(entry.sequence) - 1:
                            
                            # get the cut sequence
                            digested.append((entry.sequence[s:], s, len(entry.sequence)))
                            break
                        
                        # check of this aa is an end
                        if entry.sequence[j] in ends:
                            
                            # first reduce allowed
                            allowed_misses -= 1
                            
                            # determine if we do j or j+1 based on the rule
                            # e = j if ends[entry.sequence[j]] == 'left' else j + 1
                            
                            # digested.append((entry.sequence[s:e], s, e))
                        
        for d in digested:
            new_prots[p_name].append(DatabaseEntry(d[0], entry.description))
            
    db = db._replace(proteins=new_prots)
    return db

def digestion_filtering(sequences: list, digest_type: str, missed_cleavages: int) -> list:
    '''
    Take in a list of sequences, remove any of those that do not follow 
    one of the rules of the digest, and return. 

    Example:
        sequences = [DLAT, MALW, SSWQTK, RRRTRK]
        digest_type = trypsin
        missed_cleavages = 2

        outputs = [DLAT, SSWQTK]

    If a sequence is passed in that has more cleavages in it than we can 
    allow, it is also filtered out

    Inputs:
        sequences:          (list) peptide sequences
        digest_type:        (str) the digest type we are looking at
        missed_cleavages:   (int) the number of missed cleavages to allow
    Outputs:    
        (list) filtered down peptide sequences
    '''

    # first check if the digest type is a valid one
    if digest_type not in digests:
        return sequences

    filtered = []

    # an easier way for us to check if the rules are followed
    digest_rules = digests[digest_type]

    # keep track of the digest rules 
    starts = {s['amino_acid']: s['cut_position'] for s in digest_rules['start']}
    ends = {s['amino_acid']: s['cut_position'] for s in digest_rules['end']}

    # we only care about the starts if cut_position is left 
    # similarly, we only care about the ends of cut position is right
    leading = {
        s['amino_acid']: s['cut_position'] \
        for s in digest_rules['start'] if s['cut_position'] == 'left'
    }
    ending = {
        s['amino_acid']: s['cut_position'] \
        for s in digest_rules['end'] if s['cut_position'] == 'right'
    }

    for sequence in sequences:

        # first check if the start or end follow the rules
        # set to true if we have nothing in the respective dict
        follows_first = len(leading) == 0 
        follows_last = len(ending) == 0

        # lets first check left position
        if len(leading) and sequence[0] in leading:
            follows_first = True

        # check the end
        if len(ending) and sequence[-1] in ending:
            follows_last = True 

        # if neither tule is followed, continue
        if not (follows_first or follows_last):
            continue

        # make sure that we only have at most missed missed_cleavages number
        # of both starts and ends
        num_missed = len([aa for aa in sequence if aa in starts])
        num_missed += len([aa for aa in sequence if aa in ends])

        # final filtering: if we have missed <= allowed, then add it to filtered
        if num_missed <= missed_cleavages:
            filtered.append(sequence)

    return filtered