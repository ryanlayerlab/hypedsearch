import pytest

from fiji.fiji_utils import HypedsearchOnFijiConfig
from src.constants import MAC_CRUX_EXECUTABLE
from tests.unit.test_hypedsearch import default_test_config


class Test_HypedsearchOnFijiConfig:
    @staticmethod
    def test_init(tmp_path, test_data_dir):
        hs_config = default_test_config(test_data_dir=test_data_dir, out_dir=tmp_path)
        fiji_config = HypedsearchOnFijiConfig(hs_config=hs_config)
        assert fiji_config.hs_config
        config_path = tmp_path / "config.json"
        hs_config.to_json(path=config_path)
        HypedsearchOnFijiConfig.from_json(path=config_path)

    class Test_prepare_files_on_fiji:
        @staticmethod
        def test_smoke(tmp_path, test_data_dir):
            hs_config = default_test_config(
                test_data_dir=test_data_dir, out_dir=tmp_path
            )
            fiji_prep = HypedsearchOnFijiConfig(hs_config=hs_config)
            fiji_config = fiji_prep.prepare_files_for_hybrid_run(node_data_dir=tmp_path)
            assert fiji_config.hybrid_former.kmer_db.exists()
