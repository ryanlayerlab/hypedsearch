from pyteomics import fasta
from collections import namedtuple, defaultdict

from objects import Database

def extract_protein_name(prot_entry: namedtuple) -> str:
    '''Extract the protein name from a protein entry namedtuple from pyteomics 
    fasta read. 
    
    :param prot_entry: a namedtuple with a value of 'description'
    :type prot_entry: namedtuple

    :returns: the name of the protein
    :rtype: str
    '''

    if '|' in prot_entry.description:
        return prot_entry.description.split('|')[-1].split(' ')[0]

    return prot_entry.description.split(' ')[0]

def build(fasta_file: str) -> Database:
    '''Create a Database namedtuple from a fasta file

    :param fasta_file: the full path to a fasta database file 
    :type fasta_file: str

    :returns: a Database object with the fasta file and protein fields filled in
    :rtype: Database
    '''

    db = Database(fasta_file)
    prots = defaultdict(list)
    get_name = lambda x: x.split('|')[-1].split()[0]
    for entry in fasta.read(fasta_file):
        p_name = get_name(entry.description)
        prots[p_name].append(entry)
    db = db._replace(proteins=prots)
    return db


def get_proteins_with_subsequence(db: Database, sequence: str) -> list:
    '''Find the name of all proteins that have the subsequence provided. A 
    list of these names are returned

    :param db: source of the proteins
    :type db: Database
    :param sequence: the subsequence to look for
    :type sequence: str

    :returns: all protein names of source proteins
    :rtype: list
    '''

    return list(set(db.kmers[sequence]))

def get_proteins_with_subsequence_ion(db: Database, sequence: str, ion: str) -> list:
    '''Find all protein names that have the subsequence. Recursivley search
    if the full sequence is not found immediately

    :param db: source of the proteins
    :type db: Database
    :param sequence: subsequence to look for 
    :type sequence: str
    :param ion: the ion type. Either 'b' or 'y'
    :type ion: str

    :returns: names of the source protein(s)
    :rtype: list
    '''

    hits = []
    subseq = sequence

    while len(hits) == 0 and len(subseq) > 0:

        # get some hits
        hs = get_proteins_with_subsequence(db, subseq)

        # make sure that these hits HAVE the full sequence
        for h in hs:
            for entry in get_entry_by_name(db, h):
                if sequence in entry.sequence:
                    hits.append(h)

        # if no hits, reduce the subseq according to the ion
        if len(hits) == 0:
            subseq = subseq[1:] if ion == 'y' else subseq[:-1]

    return hits

def get_entry_by_name(db: Database, name: str) -> namedtuple:
    '''Get a namedtuple of the protein entry from the database. 

    :param db: source of proteins
    :type db: Database
    :param name: the name of the protein to look for 
    :type name: str

    :returns: namedtuple with fields 'description' and 'sequence'
    :rtype: namedtuple
    '''
    return db.proteins[name]