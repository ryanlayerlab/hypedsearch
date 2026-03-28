from src.postprocess_hs_results import process_native_and_hybrid_runs_via_config


def test_smoke():
    process_native_and_hybrid_runs_via_config(
        hs_config="results/022426_it_ot/configs/1_OT_CID_0h_stressFree_AspN.config.json",
        perform_native_analysis=False,
    )
