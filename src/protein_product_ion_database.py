import logging
from dataclasses import dataclass, field, fields
from pathlib import Path
from time import time
from typing import Any, Dict, List, Literal, Optional

from src.constants import (
    AMINO_ACID_MASSES,
    B_ION_TYPE,
    DEFAULT_MAX_K,
    DEFAULT_MIN_K,
    DEFAULT_PPM_TOLERANCE,
    EXCLUSIVE_END,
    INCLUSIVE_START,
    ION_INT_TO_TYPE,
    ION_TYPE_TO_INT,
    MEMORY,
    NEUTRAL_MASS,
    PRODUCT_ION_TABLE,
    PROTEIN_ID,
    PROTEIN_TABLE,
    PROTON_MASS,
    WATER_MASS,
    Y_ION_TYPE,
    IonTypes,
)
from src.mass_spectra import Peak, Spectrum
from src.peptides_and_ions import (
    Peptide,
    ProductIon,
    compute_peptide_mz,
    get_uniq_kmer_to_protein_map,
)
from src.sql_database import (
    PrimaryKey,
    SqlColumn,
    Sqlite3Database,
    SqlTableRow,
    add_data_classes_to_table,
)
from src.utils import (
    Kmer,
    Position,
    flatten_list_of_lists,
    get_positions_of_subseq_in_seq,
    get_time_in_diff_units,
    log_params,
    relative_ppm_tolerance_in_daltons,
)

logger = logging.getLogger(__name__)


@dataclass
class DbKmer(SqlTableRow):
    seq: str
    aa_mass: float
    protein_ids: str

    @classmethod
    def seq_to_ion(
        cls,
        seq: str,
        protein_ids: List[int],
        amino_acid_mass_lookup: Dict[str, float] = AMINO_ACID_MASSES,
    ) -> "DbKmer":
        aa_mass = sum([amino_acid_mass_lookup[aa] for aa in seq])
        return cls(
            seq=seq, aa_mass=aa_mass, protein_ids=",".join(map(str, protein_ids))
        )

    def set_charge(self, charge: int):
        self._charge = charge

    @property
    def charge(self):
        return self._charge

    @property
    def z(self):
        return self._charge

    def set_ion_type(self, ion_type: IonTypes):
        self._ion_type = ion_type.value

    @property
    def ion_type(self):
        return self._ion_type

    @property
    def ion(self):
        return self._ion_type

    def ionize(self, charge: int, ion_type: IonTypes):
        self.set_ion_type(ion_type=ion_type)
        self.set_charge(charge=charge)

    def __str__(self):
        return f"{self.ion}-{self.z}-{self.seq}"


@dataclass
class PositionedIon:
    seq: str
    charge: int
    ion_type: str
    protein_id: int
    inclusive_start: int
    exclusive_end: int


# @dataclass
# class DbProtein(SqlTableRow):
#     id: PrimaryKey[int]
#     seq: str

#     @classmethod
#     def from_peptide(cls, peptide: Peptide) -> "DbProtein":
#         assert peptide.id is not None
#         return cls(id=peptide.id, seq=peptide.seq)


@dataclass
class DbProtein(SqlTableRow):
    id: PrimaryKey[int]
    seq: str
    name: str

    @classmethod
    def from_peptide(cls, peptide: Peptide) -> "DbProtein":
        assert peptide.id is not None
        return cls(id=peptide.id, seq=peptide.seq, name=peptide.name)


def get_aa_seq_corresponding_to_protein_region(
    peptide_id: int, inclusive_start: int, exclusive_end: int, peptides: List[Peptide]
):
    peptide = list(filter(lambda peptide: peptide.id == peptide_id, peptides))
    assert len(peptide) == 1, "There should only be one matching peptide"
    peptide = peptide[0].seq
    return peptide[inclusive_start:exclusive_end]


# @dataclass
# class DbProductIon(SqlTableRow):
#     protein_id: int
#     inclusive_start: int
#     exclusive_end: int
#     charge: int
#     neutral_mass: float
#     ion_type: int

#     @classmethod
#     def kmer_to_ion(
#         cls, kmer: Kmer, charge: int, ion_type: IonTypes, protein_id: int
#     ) -> "DbProductIon":
#         product_ion = ProductIon(seq=kmer.seq, charge=charge, ion_type=ion_type)
#         return cls(
#             protein_id=protein_id,
#             # position=kmer.position,
#             inclusive_start=kmer.position.inclusive_start,
#             exclusive_end=kmer.position.exclusive_end,
#             charge=product_ion.charge,
#             neutral_mass=product_ion.neutral_mass,
#             ion_type=ION_TYPE_TO_INT[product_ion.ion_type.value],
#         )

#     @classmethod
#     def from_peptide(
#         cls,
#         peptide: Peptide,
#         charges: List[int],
#         ion_types: List[IonTypes],
#         max_k: int = DEFAULT_MAX_K,
#         min_k: int = DEFAULT_MIN_K,
#     ) -> List["DbProductIon"]:
#         """ """
#         kmers = peptide.kmers(min_k=min_k, max_k=max_k)
#         ions = []
#         for kmer in kmers:
#             for charge in charges:
#                 for ion_type in ion_types:
#                     ions.append(
#                         cls.kmer_to_ion(
#                             kmer=kmer,
#                             charge=charge,
#                             ion_type=ion_type,
#                             protein_id=peptide.id,
#                         )
#                     )
#         return ions

#     def set_aa_seq(self, db):
#         aa_seq = get_aa_seq_from_db(
#             protein_id=self.protein_id,
#             inclusive_start=self.inclusive_start,
#             exclusive_end=self.exclusive_end,
#             db=db,
#         )
#         self._aa_seq = aa_seq
#         return aa_seq

#     @property
#     def aa_seq(self):
#         return self._aa_seq

#     def get_aa_seq_from_db_peptides(self, db_peptides: List[Peptide]):
#         return get_aa_seq_corresponding_to_protein_region(
#             peptide_id=self.protein_id,
#             inclusive_start=self.inclusive_start,
#             exclusive_end=self.exclusive_end,
#             peptides=db_peptides,
#         )


def b_bounds(ppm_tolerance, query_mass, charge):
    mz_tolerance = relative_ppm_tolerance_in_daltons(
        ppm=ppm_tolerance, ref_mass=query_mass
    )

    lower_mz_bound = query_mass - mz_tolerance
    upper_mz_bound = query_mass + mz_tolerance

    lower_bound = charge * (lower_mz_bound - PROTON_MASS)
    upper_bound = charge * (upper_mz_bound - PROTON_MASS)

    return lower_bound, upper_bound


def y_bounds(ppm_tolerance, query_mass, charge):
    b_lower, b_upper = b_bounds(
        ppm_tolerance=ppm_tolerance, query_mass=query_mass, charge=charge
    )
    y_lower, y_upper = b_lower - WATER_MASS, b_upper - WATER_MASS
    return y_lower, y_upper


class ProteinProductIonDb(Sqlite3Database):
    def __init__(
        self,
        db_path: str = MEMORY,
        protein_table_name: str = PROTEIN_TABLE,
        product_ion_table_name: str = PRODUCT_ION_TABLE,
        protein_obj: type[DbProtein] = DbProtein,
        product_ion_obj: type[DbKmer] = DbKmer,
        overwrite: bool = False,
        # charges: List[int] = [1, 2, 3],
        # ion_types: List[IonTypes] = [IonTypes.B_ION_TYPE, IonTypes.Y_ION_TYPE],
        # min_k: int = DEFAULT_MIN_K,
        # max_k: int = DEFAULT_MAX_K,
    ):
        super().__init__(db_path=db_path, overwrite=overwrite)
        self.protein_table_name = protein_table_name
        self.product_ion_table_name = product_ion_table_name
        self.protein_obj = protein_obj
        self.product_ion_obj = product_ion_obj
        # self.min_k = min_k
        # self.max_k = max_k
        # self.charges = charges
        # self.ion_types = ion_types

    # def create_protein_table(self):
    #     self.create_table_from_dataclass(
    #         table_name=self.protein_table_name, obj=self.protein_obj
    #     )

    # def create_product_ion_table(self):
    #     self.create_table_from_dataclass(
    #         table_name=self.product_ion_table_name, obj=self.product_ion_obj
    # )

    # def add_peptide_and_product_ions(
    #     self,
    #     peptide: Peptide,
    # ) -> None:
    #     # Add peptide to DB
    #     db_protein = self.protein_obj.from_peptide(peptide=peptide)
    #     self.insert_dataclasses(
    #         table_name=self.protein_table_name, data_classes=[db_protein]
    #     )

    #     # Add product ions to DB
    #     db_ions = self.product_ion_obj.from_peptide(
    #         peptide=peptide,
    #         charges=self.charges,
    #         ion_types=self.ion_types,
    #         min_k=self.min_k,
    #         max_k=self.max_k,
    #     )
    #     self.insert_dataclasses(
    #         table_name=self.product_ion_table_name, data_classes=db_ions
    #     )

    # def get_ions_within_mass_tolerance(
    #     self,
    #     query_mass: float,
    #     mz_tolerance: Optional[float] = None,
    #     ppm_tolerance: Optional[float] = None,
    # ) -> List[DbProductIon]:

    #     if ppm_tolerance is not None:
    #         mz_tolerance = relative_ppm_tolerance_in_daltons(
    #             ppm=ppm_tolerance, ref_mass=query_mass
    #         )

    #     upper_bound = query_mass + mz_tolerance
    #     lower_bound = query_mass - mz_tolerance

    #     query = f"""
    #         SELECT
    #             *
    #         FROM {self.product_ion_table_name} as product_ion
    #         WHERE product_ion.{NEUTRAL_MASS} BETWEEN {lower_bound} AND {upper_bound}
    #         ORDER BY {PROTEIN_ID}, {INCLUSIVE_START}, {EXCLUSIVE_END};
    #     """
    #     matching_ions = self.read_query(query=query)
    #     return [self.product_ion_obj(**dict(ion)) for ion in matching_ions]

    def get_proteins(self):
        rows = self.all_table_rows(table_name=self.protein_table_name)
        return [self.protein_obj(**row) for row in rows]

    def get_protein_by_id(self, protein_id):
        query = f"""
            SELECT * FROM {self.protein_table_name} WHERE id = {protein_id}
        """
        matching_proteins = self.read_query(query=query)
        assert (
            len(matching_proteins) == 1
        ), f"There should only be one protein with id = {protein_id}. There are {len(matching_proteins)}"
        return self.protein_obj(**matching_proteins[0])

    def get_product_ions(self):
        rows = self.all_table_rows(table_name=self.product_ion_table_name)
        return [self.product_ion_obj(**row) for row in rows]

    def get_matching_product_ions(
        self,
        peak_mz: float,
        precursor_charge: int,
        precursor_mz: float,
        ppm_tolerance: float,
        ion_types: List[IonTypes],
    ):
        # Constant only used in this function
        bounds_fcns = {IonTypes.B_ION_TYPE: b_bounds, IonTypes.Y_ION_TYPE: y_bounds}

        # Get the product ions that match the given peak
        results = []
        for charge in range(1, precursor_charge + 1):
            for ion_type in ion_types:
                # logger.info(f"Getting charge={charge}, {ion_type} ions...")
                bounds_fcn = bounds_fcns[ion_type]
                lower_bound, upper_bound = bounds_fcn(
                    ppm_tolerance=ppm_tolerance,
                    query_mass=peak_mz,
                    charge=charge,
                )
                query = f"""
                    SELECT
                        *
                    FROM {self.product_ion_table_name} as product_ion
                    WHERE product_ion.aa_mass BETWEEN {lower_bound} AND {upper_bound}
                """
                matching_ions = self.read_query(query=query)
                matching_ions = [
                    self.product_ion_obj(**dict(ion)) for ion in matching_ions
                ]
                # logger.info(f"Matching ions PRE-filtering: {matching_ions}")
                matching_ions = list(
                    filter(
                        lambda ion: compute_peptide_mz(
                            aa_seq=ion.seq, charge=precursor_charge
                        )
                        <= precursor_mz,
                        matching_ions,
                    ),
                )

                # Add charge and ion type to ions-to-be-returned
                for ion in matching_ions:
                    ion.ionize(charge=charge, ion_type=ion_type)
                # logger.info(f"Matching ions POST-filtering: {matching_ions}")
                results.extend(matching_ions)
        return results


# def create_and_populate_protein_and_product_ion_database(
#     charges: List[int] = [1, 2, 3],
#     ion_types: List[IonTypes] = [IonTypes.B_ION_TYPE, IonTypes.Y_ION_TYPE],
#     protein_obj: DbProtein = DbProtein,
#     product_ion_obj: DbProductIon = DbProductIon,
#     protein_table_name: str = PROTEIN_TABLE,
#     product_ion_table_name: str = PRODUCT_ION_TABLE,
#     db_path: str = MEMORY,
#     min_k: int = DEFAULT_MIN_K,
#     max_k: int = DEFAULT_MAX_K,
#     peptides: Optional[List[Peptide]] = None,
#     fasta_path: Optional[str] = None,
# ) -> ProteinProductIonDb:
#     logger.info("Creating protein-product ion database...")
#     if fasta_path is not None:
#         peptides = Peptide.from_fasta(fasta_path=fasta_path)

#     db = ProteinProductIonDb(
#         charges=charges,
#         ion_types=ion_types,
#         db_path=db_path,
#         protein_table_name=protein_table_name,
#         product_ion_table_name=product_ion_table_name,
#         min_k=min_k,
#         max_k=max_k,
#         protein_obj=protein_obj,
#         product_ion_obj=product_ion_obj,
#     )

#     db.create_protein_table()
#     db.create_product_ion_table()
#     num_peptides = len(peptides)
#     for p_idx, peptide in enumerate(peptides):
#         logger.info(f"Adding protein {p_idx + 1} / {num_peptides}...")
#         db.add_peptide_and_product_ions(peptide=peptide)

#     return db


# def protein_product_ion_database_file_name(
#     charges: List[int],
#     ion_types: List[IonTypes],
#     file_name_prefix: str,
#     min_k: int,
#     max_k: int,
# ):
#     ion_types = ",".join([ion_type.value for ion_type in ion_types])
#     charges = ",".join([str(charge) for charge in charges])
#     return f"{file_name_prefix}_ionTypes={ion_types}_charges={charges}_minK={min_k}_maxK={max_k}.db"


@dataclass
class IonWithSeq:
    ion: DbKmer
    seq: str


@dataclass
class PeakWithMatchingProductIons:
    peak: Peak
    ions: List[DbKmer]


def get_aa_seq_from_db(
    protein_id: int, inclusive_start: int, exclusive_end: int, db: ProteinProductIonDb
) -> str:
    # Get protein from database
    protein = db.get_protein_by_id(protein_id=protein_id)

    # Get AA sequence from the protein
    aa_seq = protein.seq[inclusive_start:exclusive_end]

    return aa_seq


def get_product_ions_matching_spectrum(
    spectrum: Spectrum,
    db: ProteinProductIonDb,
    peak_product_ion_ppm_tolerance: float = DEFAULT_PPM_TOLERANCE,
    ion_types: List[IonTypes] = [IonTypes.B_ION_TYPE, IonTypes.Y_ION_TYPE],
) -> List[PeakWithMatchingProductIons]:
    peaks_and_matches = []

    start_time = time()
    for peak_idx, peak in enumerate(spectrum.peaks):
        matching_product_ions = db.get_matching_product_ions(
            peak_mz=peak.mz,
            precursor_charge=spectrum.precursor_charge,
            precursor_mz=spectrum.precursor_mz,
            ion_types=ion_types,
            ppm_tolerance=peak_product_ion_ppm_tolerance,
        )
        peaks_and_matches.append(
            PeakWithMatchingProductIons(peak=peak, ions=matching_product_ions)
        )
    logger.info(
        f"Peak-to-product-ion matching took {get_time_in_diff_units(time()-start_time)}"
    )
    return peaks_and_matches


def create_db(
    db_proteins: List[Peptide],
    uniq_kmer_to_protein_map: Dict[str, List[int]],
    db_path: Path = MEMORY,
    overwrite: bool = True,
) -> ProteinProductIonDb:
    # Determine whether or not the database already exists and needs to be created
    if db_path != MEMORY:
        create_db = (not db_path.exists()) or overwrite
    else:
        create_db = True
    db = ProteinProductIonDb(
        db_path=db_path,
        overwrite=overwrite,
        product_ion_obj=DbKmer,
        protein_obj=DbProtein,
    )
    if create_db:
        logger.info("Creating protein-product ion database...")
        start_time = time()
        db.create_table_from_dataclass(
            table_name=db.protein_table_name, obj=db.protein_obj
        )
        db.create_table_from_dataclass(
            table_name=db.product_ion_table_name, obj=db.product_ion_obj
        )

        # Add proteins
        logger.info("Adding proteins to DB...")
        prots = [
            db.protein_obj.from_peptide(peptide=peptide) for peptide in db_proteins
        ]
        db.insert_dataclasses(table_name=db.protein_table_name, data_classes=prots)

        # Add product ions
        logger.info("Adding kmers to DB...")
        ions = [
            DbKmer.seq_to_ion(seq=seq, protein_ids=protein_ids)
            for seq, protein_ids in uniq_kmer_to_protein_map.items()
        ]
        db.insert_dataclasses(table_name=db.product_ion_table_name, data_classes=ions)

        # Add index
        db.add_index(
            table_name=db.product_ion_table_name,
            index_name="mass",
            colms_to_index=["aa_mass"],
        )
        logger.info(
            f"Creating and indexing DB took {round(time()-start_time, 2)} seconds"
        )

    else:
        logger.info(
            f"Protein-product ion database already exists and overwrite is {overwrite}. Not creating it."
        )
    return db


def process_peak_matching_ions(
    peaks_with_matches: List[PeakWithMatchingProductIons],
    kmer_to_protein_map: Dict[str, int],
    db: ProteinProductIonDb,
):

    # Filter out the peaks with no product ion matches
    peaks_with_matches = list(
        filter(
            lambda peak_with_matches: len(peak_with_matches.ions) > 0,
            peaks_with_matches,
        )
    )

    # Just look at the ions, ignoring the peaks
    ions = flatten_list_of_lists(
        [peak_and_matches.ions for peak_and_matches in peaks_with_matches]
    )

    # Get the position of each matched product ion
    positioned_ions = []
    for ion in ions:
        for p_id in kmer_to_protein_map[ion.seq]:
            prot = db.get_protein_by_id(protein_id=p_id)
            new_ions = [
                PositionedIon(
                    seq=ion.seq,
                    charge=ion.charge,
                    ion_type=ion.ion_type,
                    protein_id=p_id,
                    inclusive_start=pos.inclusive_start,
                    exclusive_end=pos.exclusive_end,
                )
                for pos in get_positions_of_subseq_in_seq(subseq=ion.seq, seq=prot.seq)
            ]
            positioned_ions.extend(new_ions)

    return positioned_ions
