from src.constants import MAC_CRUX_EXECUTABLE, MOUSE_PROTEOME
from src.hybrids_via_clusters import HybridPeptide
from src.hypedsearch import (
    HybridRunParams,
    HypedsearchRunConfig,
    create_hybrids_fasta,
    hybrid_run_on_spectrum,
    run_hypedsearch_in_parallel,
)
from src.mass_spectra import Spectrum
from src.peptides_and_ions import Fasta
from src.psm import CometPSM
from src.utils import from_pickle, load_json, setup_logger


class Test_HypedsearchRunConfig:
    @staticmethod
    def test_native_run(tmp_path, test_hs_config_path):
        data = load_json(path=test_hs_config_path)
        data["parent_output_dir"] = str(tmp_path)
        config = HypedsearchRunConfig(**data)
        outputs = config.run_native_comet(
            crux_path=MAC_CRUX_EXECUTABLE,
            # dry_run=True
        )
        assert len(outputs) == len(config.mzml_names)
        assert outputs[0].target.exists()
        assert outputs[0].decoy.exists()
        assert len(CometPSM.from_txt(txt=outputs[0].target)) > 0


class Test_hybrid_run_on_spectrum:
    @staticmethod
    def test_smoke(tmp_path, test_hs_config_path, mouse_mzml_path):
        data = load_json(path=test_hs_config_path)
        data["parent_output_dir"] = str(tmp_path)
        config = HypedsearchRunConfig(**data)
        params = HybridRunParams(
            kmer_db=config.kmer_db_path,
            fasta=config.fasta,
            fasta_fm_index=from_pickle(path=config.fasta_fm_index),
            crux_comet_params=config.crux_comet_params,
            out_dir=tmp_path,
        )
        cmd_result, run = hybrid_run_on_spectrum(
            spectrum=Spectrum.get_spectrum(scan=7, mzml=mouse_mzml_path),
            params=params,
            crux_path=MAC_CRUX_EXECUTABLE,
            fasta_dir=tmp_path,
        )
        assert len(CometPSM.from_txt(txt=run.standardized_comet_outputs.target)) > 0
        assert run.standardized_comet_outputs.decoy is None

    @staticmethod
    def test_no_output_file(tmp_path, test_hs_config_path, mouse_mzml_path):
        data = load_json(path=test_hs_config_path)
        data["parent_output_dir"] = str(tmp_path)
        config = HypedsearchRunConfig(**data)
        params = HybridRunParams(
            kmer_db=config.kmer_db_path,
            fasta=config.fasta,
            fasta_fm_index=from_pickle(path=config.fasta_fm_index),
            crux_comet_params=config.crux_comet_params,
            out_dir=tmp_path,
        )
        cmd_result, run = hybrid_run_on_spectrum(
            spectrum=Spectrum.get_spectrum(scan=2, mzml=mouse_mzml_path),
            params=params,
            crux_path=MAC_CRUX_EXECUTABLE,
            fasta_dir=tmp_path,
        )


class Test_run_in_parallel:
    @staticmethod
    def test_smoke(tmp_path, test_hs_config_path):
        setup_logger()
        data = load_json(path=test_hs_config_path)
        data["parent_output_dir"] = str(tmp_path)
        config = HypedsearchRunConfig(**data)
        config_path = tmp_path / "config.json"
        config.save(path=config_path)
        run_hypedsearch_in_parallel(
            config=config_path,
            n_cores=4,
            crux_path=MAC_CRUX_EXECUTABLE,
        )
        pass
