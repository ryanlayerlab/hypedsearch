from src.comet_utils import CometPSM
from src.hypedsearch import HypedsearchRunConfig
from src.mass_spectra import Mzml
from src.utils import load_json


def test_smoke(tmp_path):
    mzml = Mzml(
        mzml="data/251028_RP_Islet_Spikes_Crashout/3_6_Meoh_RAT_Islet_B35spike.mzML"
    )
    scan = 2306594
    data = load_json(
        path="results/251028_RP_Islet_Spikes_Crashout/no_native_hybrid_comp/3_6_Meoh_RAT_Islet_B35spike/hs.config.json"
    )
    # data["parent_out_dir"] = tmp_path
    data["parent_out_dir"] = "tmp"
    hs_config = HypedsearchRunConfig(**data)
    native_outs, hybrid_outs, seq_to_hybrids = hs_config.run_hypedsearch_on_spectrum(
        spectrum=mzml.get_spectrum(scan=scan)
    )

    # hy_psms = CometPSM.from_txt(txt=hybrid_outs.target)
    # CometPSM.
    # CometPSM.save()
    pass
