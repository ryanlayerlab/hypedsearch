import logging
from pathlib import Path

from src.constants import MAC_CRUX_EXECUTABLE, MOUSE_PROTEOME
from src.hypedsearch import (
    HybridRunParams,
    HypedsearchRunConfig,
    hybrid_run_on_spectrum,
    run_hypedsearch,
)
from src.mass_spectra import Mzml, Spectrum
from src.psm import CometPSM
from tests.conftest import update_hs_config_out_dir_and_save_json


class Test_HybridRunParams:
    @staticmethod
    def test_load(test_data_dir):
        path = test_data_dir / "hybrid_run_params.json"
        params = HybridRunParams.load(path=path)


class Test_HypedsearchRunConfig:
    class Test_load:
        @staticmethod
        def test_load_everything_in_one_json(test_data_dir):
            HypedsearchRunConfig.from_json(
                path=test_data_dir / "ex1.hypedsearch.config.json"
            )

        @staticmethod
        def test_load_where_hybrid_run_params_are_in_another_file(test_data_dir):
            HypedsearchRunConfig.from_json(
                path=test_data_dir / "ex2.hypedsearch.config.json"
            )

    @staticmethod
    def test_native_comet_run(tmp_path, test_data_dir):
        config = HypedsearchRunConfig.from_json(
            path=update_hs_config_out_dir_and_save_json(
                old_config_path=test_data_dir / "ex1.hypedsearch.config.json",
                out_dir=tmp_path,
            )
        )
        outputs = config.native_comet_run_on_all_spectra(
            crux_path=MAC_CRUX_EXECUTABLE,
            # dry_run=True
        )
        assert len(outputs) == len(config.mzml_names)
        assert outputs[0].target.exists()
        assert outputs[0].decoy.exists()
        assert len(CometPSM.from_txt(txt=outputs[0].target)) > 0

    @staticmethod
    def test_hybrid_run_on_spectrum(tmp_path, test_data_dir):
        config = HypedsearchRunConfig.from_json(
            path=update_hs_config_out_dir_and_save_json(
                old_config_path=test_data_dir / "ex1.hypedsearch.config.json",
                out_dir=tmp_path,
            )
        )
        result = config.hybrid_run_on_spectrum(
            spectrum=Spectrum.get_spectrum(
                mzml=test_data_dir / "BMEM_AspN_Fxn4_scans1-20.mzML",
                scan=7,
            ),
            crux_path=MAC_CRUX_EXECUTABLE,
            fasta_dir=tmp_path,
            delete_hybrids_fasta=False,
        )
        assert len(CometPSM.from_txt(result[1].standardized_comet_outputs.target)) > 0


class Test_run_hypedsearch:
    @staticmethod
    def test_in_parallel(tmp_path, test_data_dir, caplog):
        config_path = update_hs_config_out_dir_and_save_json(
            old_config_path=test_data_dir / "ex1.hypedsearch.config.json",
            out_dir=tmp_path,
        )
        config = HypedsearchRunConfig.from_json(path=config_path)
        config.native_comet_run_on_all_spectra(crux_path=MAC_CRUX_EXECUTABLE)
        with caplog.at_level(logging.INFO):
            run_hypedsearch(
                config=config_path,
                n_cores=4,
                crux_path=MAC_CRUX_EXECUTABLE,
                run_in_parallel=True,
            )
        assert len(list(config.hybrid_run_scan_results_dir.glob("*"))) == 12
        assert "Running HypedSearch in parallel" in caplog.text

    @staticmethod
    def test_in_serial(tmp_path, test_data_dir, caplog):
        config_path = update_hs_config_out_dir_and_save_json(
            old_config_path=test_data_dir / "ex1.hypedsearch.config.json",
            out_dir=tmp_path,
        )
        config = HypedsearchRunConfig.from_json(path=config_path)
        config.native_comet_run_on_all_spectra(crux_path=MAC_CRUX_EXECUTABLE)
        with caplog.at_level(logging.INFO):
            run_hypedsearch(
                config=config_path,
                n_cores=4,
                crux_path=MAC_CRUX_EXECUTABLE,
                run_in_parallel=False,
            )
        assert len(list(config.hybrid_run_scan_results_dir.glob("*"))) == 12
        assert "Running HypedSearch in serial" in caplog.text


class Test_april_15_2026:
    @staticmethod
    def test_smoke():
        # Arrange
        data_dir = Path("data/260119_HuIslet_TimeCourse_Procal_Spiked")
        mzml = Mzml(path=data_dir / "HuIslet_AspN_06_IL1B_2hr.mzML")
        hybrid_run_params = HybridRunParams(
            kmer_db_path="results/021726_260119_HuIslet_TimeCourse_Procal_Spiked/kmer_dbs/HuIslet_AspN_06_IL1B_2hr.kmers.db",
            fasta="fastas/uniprotkb_proteome_UP000005640_AND_revi_2025_04_29.fasta",
            fasta_fm_index="fastas/uniprotkb_proteome_UP000005640_AND_revi_2025_04_29.mfmindex",
            min_hybrid_side_len=3,
            crux_comet_params="results/021726_260119_HuIslet_TimeCourse_Procal_Spiked/inputs/crux.comet.params",
        )
        hs_config = HypedsearchRunConfig.from_json(
            path=f"results/04-15-26-260119_HuIslet_TimeCourse_Procal_Spiked/configs/{mzml.name}.json"
        )
        hs_config.get_spectra_with_no_hybrid_results()

        pass
