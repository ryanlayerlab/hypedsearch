import logging
import shutil
import tempfile
from pathlib import Path
from time import time
from typing import List, Optional, Union

import click

from src.constants import (
    DEFAULT_MAX_KMER_LEN,
    DEFAULT_MIN_KMER_LEN,
    DEFAULT_NUM_PSMS,
    DEFAULT_PPM_TOLERANCE,
    PROTON_MASS,
    WATER_MASS,
    HybridFormingMethods,
)
from src.hybrids_via_clusters import HybridPeptide, SeqWithMass
from src.hypedsearch_utils import (
    HybridPeptide,
    create_hybrids_fasta,
    hypedsearch_output_stem,
    remove_native_hybrids,
)
from src.mass_spectra import Spectrum
from src.peptides_and_ions import (
    Peptide,
    compute_peptide_mz,
    get_proteins_by_name,
    get_uniq_kmer_to_protein_map,
)
from src.sql_database import Sqlite3Database
from src.utils import (
    ClickOptions,
    get_time_in_diff_units,
    log_params,
    log_time,
    relative_ppm_tolerance_in_daltons,
    setup_logger,
)

logger = logging.getLogger(__name__)


@log_time(level=logging.DEBUG)
def form_all_hybrids(
    proteins: List[Peptide],
    spectrum: Optional[Spectrum] = None,
    precursor_charge: Optional[int] = None,
    precursor_mz: Optional[float] = None,
    precursor_mz_ppm_tol: float = DEFAULT_PPM_TOLERANCE,
    min_k: int = DEFAULT_MIN_KMER_LEN,
    max_k: int = DEFAULT_MAX_KMER_LEN,
) -> List[HybridPeptide]:

    if spectrum is not None:
        precursor_charge = spectrum.precursor_charge
        precursor_mz = spectrum.precursor_mz

    logger.info("For given proteins, forming the unique kmer to protein map")
    uniq_kmer_to_protein_map = get_uniq_kmer_to_protein_map(
        proteins=proteins,
        min_k=min_k,
        max_k=max_k,
        protein_attr="name",
    )
    db_rows = [
        SeqWithMass.from_seq(seq=seq, charge=precursor_charge)
        for seq in uniq_kmer_to_protein_map
    ]
    db = Sqlite3Database()
    table_name = "kmers"
    db.create_table_from_dataclass(table_name=table_name, obj=SeqWithMass)
    db.insert_dataclasses(table_name=table_name, data_classes=db_rows)
    db.add_index(table_name=table_name, index_name="mass", colms_to_index=["mz"])

    # For each b-sequence search the y-sequence database to find the hybrids
    # that would produce a peptide within the given PPM tolerance of the precursor m/z
    mz_tolerance = relative_ppm_tolerance_in_daltons(
        ppm=precursor_mz_ppm_tol, ref_mass=precursor_mz
    )
    adjusted_precursor_mz = precursor_mz + (WATER_MASS / precursor_charge) + PROTON_MASS
    potench_hybrids = []
    for b_seq in uniq_kmer_to_protein_map:
        if b_seq == "EVEDPQVEQLE":
            pass
        b_seq_mz = compute_peptide_mz(aa_seq=b_seq, charge=precursor_charge)
        lower_bdd = adjusted_precursor_mz - mz_tolerance - b_seq_mz
        upper_bdd = adjusted_precursor_mz + mz_tolerance - b_seq_mz
        query = f"""
            SELECT
                *
            FROM {table_name} as ion
            WHERE ion.mz BETWEEN {lower_bdd} AND {upper_bdd}
        """
        matches = [match["seq"] for match in db.read_query(query=query)]
        for y_seq in matches:
            potench_hybrids.append(
                HybridPeptide(
                    b_seq=b_seq,
                    y_seq=y_seq,
                    b_prot_names=uniq_kmer_to_protein_map[b_seq],
                    y_prot_names=uniq_kmer_to_protein_map[y_seq],
                )
            )
    # Remove hybrids that are native
    potench_hybrids = remove_native_hybrids(
        proteins=proteins, hybrids=potench_hybrids, min_k=min_k, max_k=max_k
    )

    return potench_hybrids
