from pyteomics import fasta
from collections import namedtuple, defaultdict
from lookups.objects import Database

def extract_protein_name(prot_entry: namedtuple) -> str:
    if '|' in prot_entry.description:
        return prot_entry.description.split('|')[-1].split(' ')[0]
    return prot_entry.description.split(' ')[0]

def build(fasta_file: str) -> Database:
    db = Database(fasta_file)
    prots = [entry for entry in fasta.read(fasta_file)]
    db = db._replace(proteins=prots)
    return db

def get_proteins_with_subsequence(db: Database, sequence: str) -> list:
    return list(set(db.kmers[sequence]))

def get_proteins_with_subsequence_ion(db: Database, sequence: str, ion: str) -> list:
    hits = []
    subseq = sequence
    while len(hits) == 0 and len(subseq) > 0:
        hs = get_proteins_with_subsequence(db, subseq)
        for h in hs:
            for entry in get_entry_by_name(db, h):
                if sequence in entry.sequence:
                    hits.append(h)
        if len(hits) == 0:
            subseq = subseq[1:] if ion == 'y' else subseq[:-1]
    return hits

def get_entry_by_name(db: Database, name: str) -> namedtuple:
    return db.proteins[name]