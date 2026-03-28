from slurm.fiji_utils import SbatchConfig, create_config_for_fiji_run


class Test_create_config_for_fiji_run:
    @staticmethod
    def test_smoke(tmp_path, hs_config_path):
        fiji_config = create_config_for_fiji_run(
            node_data_dir=tmp_path,
            config=hs_config_path,
        )
        assert fiji_config.parent_output_dir.exists()
        assert fiji_config.fasta.exists()
        assert fiji_config.kmer_db_path.exists()
        assert fiji_config.crux_comet_params.exists()
        for mzml in fiji_config.mzml_to_scans.keys():
            assert mzml.exists()


class Test_SbatchConfig:
    @staticmethod
    def test_smoke():
        config = SbatchConfig(name="test", time="24:00:00")
        config.get_sbatch_directives_lines()
