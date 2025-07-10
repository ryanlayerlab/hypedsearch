from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

from src.comet_utils import CometPSM
from src.hypedsearch_utils import HybridPeptide
from src.mass_spectra import Spectrum
from src.peptide_spectrum_comparison import PeptideSpectrumComparison
from src.run_hypedsearch import HSResult
from src.utils import dataclass_list_to_df


@dataclass
@dataclass
class HybridResults:
    hybrid_psms: List[CometPSM]
    hybrid_seqs: Dict[str, List[HybridPeptide]]

    @property
    def best_xcorr(self):
        best_hybrid_psms, xcorr = get_max_objs_and_val(
            obj_list=self.hybrid_psms, obj_attr="xcorr"
        )


@dataclass
class NativeVsHybridComparison:
    spectrum: Spectrum
    native_psms: List[CometPSM]
    hybrid_psms: List[CometPSM]
    seq_to_hybrid_map: Optional[Dict[str, List[HybridPeptide]]] = None
    native_comp: Optional[PeptideSpectrumComparison] = field(init=False, default=None)
    hybrid_comp: Optional[PeptideSpectrumComparison] = field(init=False, default=None)
    native_xcorr: Optional[float] = field(init=False, default=None)
    hybrid_xcorr: Optional[float] = field(init=False, default=None)
    native_seq: Optional[str] = field(init=False, default=None)
    hybrid_seq: Optional[str] = field(init=False, default=None)

    def __post_init__(self):
        if len(self.native_psms) > 0:
            best_native_psms, native_xcorr = get_max_objs_and_val(
                obj_list=self.native_psms, obj_attr="xcorr"
            )
            best_native_psm = best_native_psms[0]
            self.native_comp = PeptideSpectrumComparison.from_str(
                spectrum=self.spectrum, peptide=best_native_psm.seq
            )
            self.native_xcorr = native_xcorr
            self.native_seq = best_native_psm.seq

        if len(self.hybrid_psms) > 0:
            # Get best hybrid PSM
            best_hybrid_psms, hybrid_xcorr = get_max_objs_and_val(
                obj_list=self.hybrid_psms, obj_attr="xcorr"
            )
            hybrid_seqs = []
            for psm in best_hybrid_psms:
                hybrid_seqs.extend(
                    [hy.seq_with_hyphen for hy in self.seq_to_hybrid_map[psm.seq]]
                )
            best_hybrid_seq = self.get_best_hybrid_seq(
                spectrum=self.spectrum, hybrid_seqs=hybrid_seqs
            )
            self.hybrid_seq = best_hybrid_seq
            self.hybrid_comp = PeptideSpectrumComparison.from_str(
                spectrum=self.spectrum, peptide=best_hybrid_seq
            )
            self.hybrid_xcorr = hybrid_xcorr

    @staticmethod
    def get_best_hybrid_seq(
        spectrum: Spectrum,
        hybrid_seqs: List[str],
        best_attr: str = "prop_supported_hybrid_seqs",
    ) -> str:
        comps = [
            PeptideSpectrumComparison.from_str(spectrum=spectrum, peptide=hy_seq)
            for hy_seq in hybrid_seqs
        ]
        best_comps, _ = get_max_objs_and_val(obj_list=comps, obj_attr=best_attr)
        best_comp = best_comps[0]
        return best_comp.hybrid_seq

    @property
    def mzml(self):
        return self.spectrum.mzml.stem

    @property
    def scan(self):
        return self.spectrum.scan

    @property
    def native_b_ion_seqs_supported(self):
        if self.native_comp is not None:
            return self.native_comp.b_ion_seqs_supported
        else:
            return 0

    @property
    def hybrid_b_ion_seqs_supported(self):
        if self.hybrid_comp is not None:
            return self.hybrid_comp.b_ion_seqs_supported
        else:
            return 0

    @property
    def native_y_ion_seqs_supported(self):
        if self.native_comp is not None:
            return self.native_comp.y_ion_seqs_supported
        else:
            return 0

    @property
    def hybrid_y_ion_seqs_supported(self):
        if self.hybrid_comp is not None:
            return self.hybrid_comp.y_ion_seqs_supported
        else:
            return 0

    @property
    def native_prop_intensity_supported(self):
        if self.native_comp is not None:
            return self.native_comp.proportion_of_intensity_supported_by_peptide
        else:
            return 0

    @property
    def hybrid_prop_intensity_supported(self):
        if self.hybrid_comp is not None:
            return self.hybrid_comp.proportion_of_intensity_supported_by_peptide
        else:
            return 0

    @property
    def native_ion_seqs_supported(self):
        if self.native_comp is not None:
            return self.native_comp.num_ion_seqs_supported
        else:
            return 0

    @property
    def hybrid_ion_seqs_supported(self):
        if self.hybrid_comp is not None:
            return self.hybrid_comp.num_ion_seqs_supported
        else:
            return 0

    @property
    def native_prop_ion_seqs_supported(self):
        if self.native_comp is not None:
            return self.native_comp.proportion_of_seqs_supported
        else:
            return 0

    @property
    def hybrid_prop_ion_seqs_supported(self):
        if self.hybrid_comp is not None:
            return self.hybrid_comp.proportion_of_seqs_supported
        else:
            return 0

    @property
    def hybrid_supported_seqs(self):
        if self.hybrid_comp is not None:
            return self.hybrid_comp.num_hybrid_seqs_supported
        else:
            return 0

    @property
    def prop_hybrid_seqs_supported(self):
        if self.hybrid_comp is not None:
            return self.hybrid_comp.prop_supported_hybrid_seqs
        else:
            return 0

    @classmethod
    def from_hs_results(
        cls, hs_results: List[HSResult]
    ) -> List["NativeVsHybridComparison"]:
        return [
            cls(
                spectrum=hs_result.spectrum,
                native_psms=hs_result.native_psms,
                hybrid_psms=hs_result.hybrid_psms,
                seq_to_hybrid_map=hs_result.seqs_to_hybrid_peptides,
            )
            for hs_result in hs_results
        ]


HS_RESULTS_FIELDS = [
    "mzml",
    "scan",
    "native_xcorr",
    "hybrid_xcorr",
    "native_seq",
    "hybrid_seq",
    "native_ion_seqs_supported",
    "hybrid_ion_seqs_supported",
    "native_prop_ion_seqs_supported",
    "hybrid_prop_ion_seqs_supported",
    "native_prop_intensity_supported",
    "hybrid_prop_intensity_supported",
    "hybrid_supported_seqs",
    "prop_hybrid_seqs_supported",
    "native_b_ion_seqs_supported",
    "hybrid_b_ion_seqs_supported",
    "native_y_ion_seqs_supported",
    "hybrid_y_ion_seqs_supported",
]


def create_hs_results_df(
    hs_results: List[HSResult],
    fields: List[str] = HS_RESULTS_FIELDS,
):
    comps = NativeVsHybridComparison.from_hs_results(hs_results=hs_results)

    return dataclass_list_to_df(dataclass_list=comps, fields=fields)


def get_max_objs_and_val(obj_list: List[Any], obj_attr: str):
    max_val = max(getattr(obj, obj_attr) for obj in obj_list)
    max_objs = [obj for obj in obj_list if getattr(obj, obj_attr) == max_val]
    return (max_objs, max_val)
