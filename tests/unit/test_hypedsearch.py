from src.constants import MAC_CRUX_EXECUTABLE, MOUSE_PROTEOME
from src.hybrids_via_clusters import HybridPeptide
from src.hypedsearch import (
    HybridRunParams,
    HypedsearchRunConfig,
    create_hybrids_fasta,
    hybrid_run_on_spectrum,
)
from src.mass_spectra import Spectrum
from src.peptides_and_ions import Fasta
from src.psm import CometPSM
from src.utils import from_pickle, load_json


class Test_HypedsearchRunConfig:
    @staticmethod
    def test_native_run(tmp_path, test_hs_config_path):
        data = load_json(path=test_hs_config_path)
        data["parent_output_dir"] = str(tmp_path)
        config = HypedsearchRunConfig(**data)
        outputs = config.native_comet_run(
            # dry_run=True
        )
        assert len(outputs) == len(config.mzml_names)
        assert outputs[0].target.exists()
        assert outputs[0].decoy.exists()
        assert len(CometPSM.from_txt(txt=outputs[0].target)) > 0

    @staticmethod
    def test_hybrid_run(tmp_path, test_hs_config_path, mouse_spectrum):
        data = load_json(path=test_hs_config_path)
        data["parent_output_dir"] = str(tmp_path)
        config = HypedsearchRunConfig(**data)
        seq_to_hybrids, comet_outputs = config.hybrid_run_on_spectrum(
            spectrum=mouse_spectrum
        )
        assert comet_outputs.target.exists()
        hybrid_psms = CometPSM.from_txt(txt=comet_outputs.target)
        assert len(hybrid_psms) == 5
        assert hybrid_psms[0].seq == "SAAPGSAAAPAAAEEKK"

    @staticmethod
    def test_hybrid_run_on_spectrum_with_multiple_right_parental_proteins(
        tmp_path, test_data_dir
    ):
        hs_config = HypedsearchRunConfig(
            mzml_to_scans={test_data_dir / "mouse_BMEM_AspN_Fxn4.mzML": [2508]},
            parent_output_dir=tmp_path,
            crux_comet_params=test_data_dir / "crux.comet.params",
            name="test",
            min_hybrid_side_len=3,
            kmer_db=test_data_dir / "mouse_samples.kmers.db",
            fasta=MOUSE_PROTEOME,
            fasta_fm_index="fastas/SwissProt.TAW_mouse_w_NOD_IAPP.mfmindex",
        )
        _, comet_outputs = hs_config.hybrid_run_on_spectrum(
            spectrum=Spectrum.get_spectrum(
                scan=2508, mzml=list(hs_config.mzml_to_scans.keys())[0]
            )
        )
        psms = CometPSM.from_txt(txt=comet_outputs.target)
        hybrid = HybridPeptide.parse_hybrid_peptide_str(hybrid_str=psms[1].proteins[0])
        assert len(hybrid.right_proteins) == 2


class Test_hybrid_run_on_spectrum:
    @staticmethod
    def test_smoke(tmp_path, test_hs_config_path, mouse_mzml_path):
        data = load_json(path=test_hs_config_path)
        data["parent_output_dir"] = str(tmp_path)
        config = HypedsearchRunConfig(**data)
        params = HybridRunParams(
            kmer_db=config.kmer_db,
            fasta=config.fasta,
            fasta_fm_index=from_pickle(path=config.fasta_fm_index),
            crux_comet_params=config.crux_comet_params,
            out_dir=tmp_path,
        )
        cmd_result, run = hybrid_run_on_spectrum(
            spectrum=Spectrum.get_spectrum(scan=7, mzml=mouse_mzml_path),
            params=params,
            crux_path=MAC_CRUX_EXECUTABLE,
        )
        assert len(CometPSM.from_txt(txt=run.nonstandardized_outputs.target)) > 0
        assert run.nonstandardized_outputs.decoy is None

    @staticmethod
    def test_no_output_file(tmp_path, test_hs_config_path, mouse_mzml_path):
        data = load_json(path=test_hs_config_path)
        data["parent_output_dir"] = str(tmp_path)
        config = HypedsearchRunConfig(**data)
        params = HybridRunParams(
            kmer_db=config.kmer_db,
            fasta=config.fasta,
            fasta_fm_index=from_pickle(path=config.fasta_fm_index),
            crux_comet_params=config.crux_comet_params,
            out_dir=tmp_path,
        )
        cmd_result, run = hybrid_run_on_spectrum(
            spectrum=Spectrum.get_spectrum(scan=2, mzml=mouse_mzml_path),
            params=params,
            crux_path=MAC_CRUX_EXECUTABLE,
        )
