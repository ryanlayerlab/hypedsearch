
from fiji.fiji_utils import create_config_for_fiji_run


class Test_create_config_for_fiji_run:
    @staticmethod
    def test_smoke(tmp_path, test_data_dir, test_hs_config_path):
        fiji_config = create_config_for_fiji_run(
            node_data_dir=tmp_path,
            config=test_hs_config_path,
        )
        assert fiji_config.parent_output_dir.exists()
        assert fiji_config.fasta.exists()
        assert fiji_config.kmer_db.exists()
        assert fiji_config.crux_comet_params.exists()
        for mzml in fiji_config.mzml_to_scans.keys():
            assert mzml.exists()
