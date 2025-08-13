import logging
from collections import defaultdict
from dataclasses import dataclass, field
from functools import cached_property
from pathlib import Path
from typing import Dict, List, Literal, Optional, Tuple, Union

import yaml

from src.comet_utils import CometPSM
from src.hybrids_via_clusters import HybridPeptide
from src.mass_spectra import Spectrum
from src.protein_abundance import (
    get_and_save_prefix_counts_by_length,
    get_protein_counts_from_comet_results,
    load_comet_psms,
)
from src.utils import load_json, number_greater_than, to_json

logger = logging.getLogger(__name__)


@dataclass
class HybridPSM:
    psm: CometPSM
    hybrids: List[HybridPeptide]

    def __post_init__(self):
        for hybrid in self.hybrids:
            assert hybrid.scan == self.psm.scan
            assert hybrid.sample == self.psm.sample

    @property
    def sample(self) -> str:
        return self.psm.sample

    @property
    def scan(self) -> int:
        return self.psm.scan

    @property
    def seq(self) -> str:
        return self.psm.seq


@dataclass
class HypedsearchResults:
    results_dir: Union[str, Path]
    config_file: Union[str, Path]
    q_value_threshold: Optional[float] = None
    config: Dict = field(init=False)

    def __post_init__(self):
        # Convert strings to Paths
        self.results_dir = Path(self.results_dir)
        self.config_file = Path(self.config_file)
        # Load config
        with open(self.config_file, "r") as file:
            self.config = yaml.safe_load(file)

        # Set q_value_threshold if not provided
        if self.q_value_threshold is None:
            self.q_value_threshold = self.config.get("q_value_threshold")

    @property
    def native_run_dir(self):
        return self.results_dir / "native_run"

    @property
    def db_dir(self):
        return self.results_dir / "db"

    @property
    def hybrids_dir(self):
        return self.results_dir / "hybrids"

    @property
    def hybrid_run_dir(self):
        return self.results_dir / "hybrid_run"

    @property
    def protein_abundance_json(self):
        return self.db_dir / "protein_abundances.json"

    @property
    def prefix_abundance_json(self):
        return self.native_run_dir / "prefix_abundances.json"

    @cached_property
    def hybrids(self) -> List[HybridPeptide]:
        # Loop over sample/MZML directories
        hybrids = []
        for sample_dir in (self.results_dir / "hybrids").iterdir():
            for json_file in sample_dir.glob("*.json"):
                hybrids.extend(HybridPeptide.from_json(json=json_file))
        return hybrids

    @cached_property
    def seq_to_hybrids(self) -> Dict[str, List[HybridPeptide]]:
        seq_to_hybrids = defaultdict(list)
        for hybrid in self.hybrids:
            seq_to_hybrids[hybrid.seq].append(hybrid)
        return seq_to_hybrids

    @cached_property
    def high_confidence_native_psms(self) -> List[CometPSM]:
        return self._get_high_confidence_psms(
            q_value_threshold=self.q_value_threshold,
            psm_type="native",
        )

    @cached_property
    def high_confidence_hybrid_psms(self) -> List[HybridPSM]:
        psms = self._get_high_confidence_psms(
            q_value_threshold=self.q_value_threshold,
            psm_type="hybrid",
        )
        return self._get_hybrid_psms(psms=psms)

    @property
    def protein_abundances(self) -> Dict[str, int]:
        return load_json(self.protein_abundance_json)

    @property
    def prefix_abundances(self) -> Dict[str, int]:
        return load_json(self.prefix_abundance_json)

    def get_abundance_and_rank(
        self,
        query_value: str,
        abundance_type: Literal["protein", "prefix"],
    ) -> Tuple[int, int, int]:
        if abundance_type == "prefix":
            abundances = self.prefix_abundances[str(len(query_value))]
        elif abundance_type == "protein":
            abundances = self.protein_abundances
        else:
            raise ValueError(
                f"Invalid abundance_type: {abundance_type}. "
                "Must be either 'protein' or 'prefix'."
            )
        abundance = abundances.get(query_value, 0)
        rank = number_greater_than(
            values=list(abundances.values()),
            query_value=abundance,
        )
        total_num = len(abundances.keys())
        return (abundance, rank, total_num)

    def _get_high_confidence_psms(
        self, q_value_threshold: float, psm_type: Literal["native", "hybrid"]
    ) -> List[CometPSM]:
        psm_file = self.results_dir / f"{psm_type}_run/assign-confidence.target.txt"
        psms = CometPSM.from_txt(txt_path=psm_file)
        high_confidence_psms = [psm for psm in psms if psm.q_value <= q_value_threshold]
        return high_confidence_psms

    def _get_hybrid_psms(self, psms: List[CometPSM]) -> List[HybridPSM]:
        # Remove PSMs that don't correspond to a hybrid
        hybrid_psms = [psm for psm in psms if psm.is_hybrid]

        # Remove PSMs that correspond to a hybrid sequence that was not formed for the given scan
        true_hybrid_psms = []
        for psm in hybrid_psms:
            matching_hybrids_for_scan = list(
                filter(
                    lambda hybrid: (hybrid.scan == psm.scan)
                    and (hybrid.sample == psm.sample),
                    self.seq_to_hybrids[psm.seq],
                )
            )
            if len(matching_hybrids_for_scan) > 0:
                true_hybrid_psms.append(
                    HybridPSM(
                        psm=psm,
                        hybrids=matching_hybrids_for_scan,
                    )
                )
        return true_hybrid_psms

    def get_spectrum(
        self,
        psm: Optional[Union[HybridPSM, CometPSM]] = None,
        sample: Optional[str] = None,
        scan: Optional[int] = None,
    ) -> Optional[Path]:
        if psm is not None:
            sample = psm.sample
            scan = psm.scan

        return Spectrum.get_spectrum(
            scan=scan, mzml=Path(self.config["mzml_dir"]) / f"{sample}.mzML"
        )

    def create_protein_and_prefix_abundance_jsons(
        self,
    ):
        if not self.prefix_abundance_json.exists():
            get_and_save_prefix_counts_by_length(
                comet_results_dir=self.native_run_dir,
                q_value_threshold=self.config.get("q_value_threshold"),
                out_path=self.protein_abundance_json,
            )

        if not self.protein_abundance_json.exists():
            psms = load_comet_psms(
                comet_results_dir=self.native_run_dir,
                q_value_threshold=self.config.get("q_value_threshold"),
            )
            prot_counts = get_protein_counts_from_comet_results(psms=psms)
            to_json(data=prot_counts, out_path=self.protein_abundance_json)

    def get_scan_psm(
        self, sample: str, scan: int, psm_type: Literal["native", "hybrid"]
    ) -> Dict[str, Union[CometPSM, HybridPSM]]:
        if psm_type == "native":
            txt_path = self.native_run_dir / "assign-confidence.target.txt"
        elif psm_type == "hybrid":
            txt_path = self.hybrid_run_dir / "assign-confidence.target.txt"
        else:
            raise ValueError(
                f"Invalid psm_type: {psm_type}. Must be 'native' or 'hybrid'."
            )
        psms = CometPSM.from_txt(txt_path=txt_path)
        scan_psm = list(
            filter(lambda psm: psm.sample == sample and psm.scan == scan, psms)
        )
        assert len(scan_psm) <= 1, (
            "Should be 0 or 1 PSM for each scan in the 'assign-confidence' output. "
            f"Found {len(scan_psm)} PSMs for sample {sample} and scan {scan},  psm_type = {psm_type}."
        )
        return scan_psm

    def get_scan_data(self, sample: str, scan: int):
        return {
            "native_psm": self.get_scan_psm(
                sample=sample, scan=scan, psm_type="native"
            ),
            "hybrid_psm": self.get_scan_psm(
                sample=sample, scan=scan, psm_type="hybrid"
            ),
            "hybrids": load_json(self.hybrids_dir / f"{sample}/{scan}.json"),
            "spectrum": self.get_spectrum(sample=sample, scan=scan),
        }

    def get_data_for_scan(self, scan: int, sample: str):
        spectrum = self.get_spectrum(sample=sample, scan=scan)
        native_psm = get_psm(
            psms=self.high_confidence_native_psms,
            sample=sample,
            scan=scan,
            should_only_be_one=False,
        )
        hybrid_psm = get_psm(
            psms=self.high_confidence_hybrid_psms,
            sample=sample,
            scan=scan,
            should_only_be_one=False,
        )
        return spectrum, native_psm, hybrid_psm


def get_psm(
    psms: Union[List[HybridPSM], List[CometPSM]],
    sample: str,
    scan: int,
    should_only_be_one: bool = True,
) -> Union[HybridPSM, CometPSM]:
    psm = list(filter(lambda psm: psm.sample == sample and psm.scan == scan, psms))
    if should_only_be_one:
        assert len(psm) == 1, (
            "Expected exactly one PSM for each (sample, scan) tuple. "
            f"For ({sample}, {scan}) found {len(psm)} PSMs."
        )
    return psm[0]
