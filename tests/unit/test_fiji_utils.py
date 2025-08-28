from src.fiji_utils import HypedsearchOnFijiConfig


class Test_RunHypedsearchOnFiji:
    class Test_prepare_files_on_fiji_node:
        @staticmethod
        def test_smoke(test_data_dir, tmp_path):
            config = HypedsearchOnFijiConfig.prepare_files_on_fiji_node(
                run_config=test_data_dir / "hypedsearch_config.yaml",
                node_data_dir=tmp_path,
            )

            assert config.fiji_config.out_dir
