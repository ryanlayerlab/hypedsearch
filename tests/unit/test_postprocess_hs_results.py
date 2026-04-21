from src.postprocess_hs_results import NativeVsHybridComparison
from src.psm import CometPSM


def default_native_vs_hybrid_comparsion(test_data_dir):
    scans = list(range(1, 21))
    native_targets = list(
        filter(
            lambda psm: psm.scan in scans,
            CometPSM.from_txt(txt=test_data_dir / "BMEM_AspN_Fxn4.target.txt"),
        )
    )
    native_decoys = list(
        filter(
            lambda psm: psm.scan in scans,
            CometPSM.from_txt(txt=test_data_dir / "BMEM_AspN_Fxn4.decoy.txt"),
        )
    )
    native_assign_conf = list(
        filter(
            lambda psm: psm.scan in scans,
            CometPSM.from_txt(
                txt=test_data_dir / "BMEM_AspN_Fxn4.assign-confidence.txt"
            ),
        )
    )
    hybrid_targets = list(
        filter(
            lambda psm: psm.scan in scans,
            CometPSM.from_txt(
                txt=test_data_dir / "BMEM_AspN_Fxn4.comet.hybrid.target.txt"
            ),
        )
    )
    # return NativeVsHybridComparison.from_txts_and_mzml(
    #     native_targets=test_data_dir / "BMEM_AspN_Fxn4.target.txt",
    #     native_decoys=test_data_dir / "BMEM_AspN_Fxn4.decoy.txt",
    #     native_assign_conf=test_data_dir / "BMEM_AspN_Fxn4.assign-confidence.txt",
    #     hybrid_targets=test_data_dir / "BMEM_AspN_Fxn4.comet.hybrid.target.txt",
    #     mzml=test_data_dir / "BMEM_AspN_Fxn4.mzML",
    # )
    return NativeVsHybridComparison.from_comet_psms(
        native_targets=native_targets,
        native_decoys=native_decoys,
        native_assign_conf=native_assign_conf,
        hybrid_targets=hybrid_targets,
        mzml=test_data_dir / "BMEM_AspN_Fxn4.mzML",
    )


class Test_NativeVsHybridComparison:
    @staticmethod
    def test_from_comet_psms(test_data_dir):
        default_native_vs_hybrid_comparsion(test_data_dir)

        # @staticmethod
        # def test_from_txts_and_mzml(test_data_dir):
        # default_native_vs_hybrid_comparsion(test_data_dir)

    @staticmethod
    def test_create_psm_dataframes(test_data_dir):
        # Arrange
        comp = default_native_vs_hybrid_comparsion(test_data_dir)
        expected_keys = ["native_targets", "native_decoys", "hybrid_targets"]
        # Act
        psm_type_to_df = comp.create_psm_dataframes()
        # Assert
        for key in expected_keys:
            assert key in psm_type_to_df


# def test_smoke(test_data_dir):
#     comp = default_native_vs_hybrid_comparsion(test_data_dir)
#     spectrum = comp.uid_to_spectrum["mzml=BMEM_AspN_Fxn4;scan=2799"]
#     spectrum.get_total_intensity()
#     comp.uid_to_spectrum_results
#     pass
