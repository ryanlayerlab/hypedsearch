from src.analyze_hypedsearch_results import HypedsearchResults


class Test_HypedsearchResults:
    class Test_initialization:
        @staticmethod
        def test_smoke():
            results_dir = "results/mouse_samples"
            config_file = "snakefiles/configs/mouse_samples.yaml"
            hsr = HypedsearchResults(results_dir=results_dir, config_file=config_file)

    class Test_hybrids:
        @staticmethod
        def test_smoke():
            # results_dir = "results/mouse_samples"
            # config_file = "snakefiles/configs/mouse_samples.yaml"
            results_dir = "results/samples_3a_and_3b"
            config_file = "snakefiles/configs/samples_3a_and_3b.yaml"
            hsr = HypedsearchResults(results_dir=results_dir, config_file=config_file)
            assert len(hsr.hybrids) > 0
            assert len(hsr.seq_to_hybrids) > 0

    class Test_get_high_confidence_psms:
        @staticmethod
        def test_native():
            results_dir = "results/samples_3a_and_3b"
            config_file = "snakefiles/configs/samples_3a_and_3b.yaml"
            hsr = HypedsearchResults(results_dir=results_dir, config_file=config_file)
            q_thresh, psm_type = 0.05, "native"
            psms = hsr.get_high_confidence_psms(
                q_value_threshold=q_thresh, psm_type=psm_type
            )
            assert all(psm.q_value <= q_thresh for psm in psms)

        @staticmethod
        def test_hybrid():
            results_dir = "results/samples_3a_and_3b"
            config_file = "snakefiles/configs/samples_3a_and_3b.yaml"
            hsr = HypedsearchResults(results_dir=results_dir, config_file=config_file)
            q_thresh, psm_type = 0.05, "hybrid"
            psms = hsr.get_high_confidence_psms(
                q_value_threshold=q_thresh, psm_type=psm_type
            )
            assert all(psm.q_value <= q_thresh for psm in psms)

    class Test_get_hybrid_psms:
        @staticmethod
        def test_smoke():
            results_dir = "results/samples_3a_and_3b"
            config_file = "snakefiles/configs/samples_3a_and_3b.yaml"
            hsr = HypedsearchResults(results_dir=results_dir, config_file=config_file)
            q_thresh, psm_type = 0.05, "hybrid"
            psms = hsr.get_high_confidence_psms(
                q_value_threshold=q_thresh, psm_type=psm_type
            )
            hybrid_psms = hsr.get_hybrid_psms(psms=psms)
            for hybrid_psm in hybrid_psms:
                assert hybrid_psm.psm.is_hybrid
                for hybrid in hybrid_psm.hybrids:
                    assert hybrid_psm.psm.scan == hybrid.scan
                    assert hybrid_psm.psm.sample == hybrid.sample


# class Test_get_hybrid_psms_from_psms:
#     @staticmethod
#     def test_smoke():
#         results_dir = "results/samples_3a_and_3b"
#         config_file = "snakefiles/configs/samples_3a_and_3b.yaml"
#         hsr = HypedsearchResults(results_dir=results_dir, config_file=config_file)
#         q_thresh, psm_type = 0.05, "hybrid"
#         psms = hsr.get_high_confidence_psms(
#             q_value_threshold=q_thresh, psm_type=psm_type
#         )
#         hybrid_psms = get_hybrid_psms_from_psms(psms=psms)
#         assert len(hybrid_psms) < len(psms)
#         assert all(psm.is_hybrid for psm in hybrid_psms)
