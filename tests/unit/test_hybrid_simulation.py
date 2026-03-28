from src.constants import MAC_CRUX_EXECUTABLE
from src.hybrid_simulation import (
    HALF,
    RANDOM,
    HybridSimulationExperiment,
    HybridSimulator,
    run_hybrid_simulation_on_spectrum,
)
from src.kmer_database import KmerDatabase
from src.mass_spectra import Mzml
from src.psm import CometPSM
from src.utils import load_json, to_json


class Test_HybridSimulator:
    @staticmethod
    def test_smoke(tmp_path, mouse_mzml_path, test_data_dir):
        aa_seq = "SAAPAAGSAPAAAEEKK"
        run_hybrid_simulation_on_spectrum(
            cut_method=HALF,
            crux_path=MAC_CRUX_EXECUTABLE,
            out_dir=tmp_path,
            aa_seq=aa_seq,
            scan=7,
            mzml=mouse_mzml_path,
            hybrid_run_params=test_data_dir / "hybrid_run_params.json",
        )
        # Assert
        native_psms = CometPSM.from_txt(
            txt=tmp_path / "native.BMEM_AspN_Fxn4_scans1-20.comet.7-7.target.txt"
        )
        assert len(native_psms) > 0
        assert native_psms[0].seq == aa_seq

        denative_psms = CometPSM.from_txt(
            txt=tmp_path / "denativized.BMEM_AspN_Fxn4_scans1-20.comet.7-7.target.txt"
        )
        assert len(denative_psms) > 0
        assert denative_psms[0].seq != aa_seq

        hybrid_psms = CometPSM.from_txt(
            txt=tmp_path / "hybrid.BMEM_AspN_Fxn4_scans1-20.comet.7-7.target.txt"
        )
        assert len(hybrid_psms) > 0
        assert hybrid_psms[0].seq == aa_seq


class Test_HybridSimulationExperiment:
    @staticmethod
    def test_run_experiment(tmp_path, mouse_mzml_path, test_data_dir):
        sim_exp = HybridSimulationExperiment(
            hybrid_run_params=load_json(path=test_data_dir / "hybrid_run_params.json"),
            hybrid_simulator=HybridSimulator(cut_method=HALF),
            scan_to_aa_seq={7: "SAAPAAGSAPAAAEEKK"},
            mzml=mouse_mzml_path,
            out_dir=tmp_path,
        )
        sim_exp.run_experiment(n_cores=4, crux_path=MAC_CRUX_EXECUTABLE)

        # Assert
        txts = list(tmp_path.glob("*.txt"))
        assert len(txts) > 0
        for txt in txts:
            assert len(CometPSM.from_txt(txt)) > 0

    @staticmethod
    def test_load_from_json_and_run(tmp_path, test_data_dir, caplog):
        name = "hybrid_simulation.config.json"
        config_path = tmp_path / name
        data = load_json(path=test_data_dir / name)
        data["out_dir"] = str(tmp_path)
        to_json(data=data, path=config_path)
        exp = HybridSimulationExperiment.load(path=config_path)
        exp.run_experiment(n_cores=4, crux_path=MAC_CRUX_EXECUTABLE)
        # Assert
        txts = list(tmp_path.glob("*.txt"))
        assert len(txts) > 0
        for txt in txts:
            assert len(CometPSM.from_txt(txt)) > 0
