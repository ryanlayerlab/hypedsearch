from objects import FastaDatabase, Spectrum
from lookups.utils import all_perms_of_s
from scoring import scoring
from alignment import alignment_utils

def __replace_ambiguous_hybrid(hybrid: tuple, db: FastaDatabase, observed: Spectrum):
    nonhybrid = hybrid[0]
    if len(preprocessing.database_generator.get_proteins_with_subsequence(db, nonhybrid)):
        return ((nonhybrid, None))
    possible = all_perms_of_s(nonhybrid, 'LI')
    if len(possible) == 0:
        return hybrid
    for p in possible:
        if len(computational_pipeline.database.get_proteins_with_subsequence(db, p)):
            return ((p, None))
    return hybrid

def replace_ambiguous_hybrids(
    hybrid_alignments: list, 
    db: FastaDatabase, 
    observed: Spectrum
    ) -> list:

    return [
        __replace_ambiguous_hybrid(hybrid_alignment, db, observed)\
         for hybrid_alignment in hybrid_alignments
    ]

def hybrid_alignment(seq1: str, seq2: str):
    alignment = ''
    hybalignment = ''
    attempted_overlap = alignment_utils.align_overlaps(seq1, seq2)

    if attempted_overlap is not None and 0 < len(attempted_overlap) < len(seq1) + len(seq2):
        rightstart = attempted_overlap.index(seq2)
        leftend = len(seq1) - 1
        middle_sec = attempted_overlap[rightstart:leftend + 1]
        alignment = attempted_overlap
        hybalignment = attempted_overlap[:rightstart] \
                        + '(' + middle_sec + ')' \
                        + attempted_overlap[leftend+1:]
    else:
        alignment = seq1 + seq2
        hybalignment = seq1 + '-' + seq2
    return (alignment, hybalignment)