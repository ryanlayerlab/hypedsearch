from collections import Counter
from pathlib import Path
from unittest.mock import MagicMock, patch

from src.comet_utils import CometPSM
from src.protein_abundance import (
    get_and_plot_most_common_proteins,
    get_prefix_counts_by_length,
    get_protein_counts_from_comet_results,
    load_comet_psms,
)


class Test_get_protein_counts_from_comet_results:
    @staticmethod
    def test_from_assign_confidence_txt(test_data_dir):
        comet_results_dir = test_data_dir / "comet_results"
        psms = load_comet_psms(
            comet_results_dir=comet_results_dir,
            q_value_threshold=0.5,
        )
        protein_counts = get_protein_counts_from_comet_results(psms=psms)
        assert len(protein_counts) > 0

    @staticmethod
    def test_from_comet_target_txts(test_data_dir):
        comet_results_dir = test_data_dir / "comet_results"
        psms = load_comet_psms(
            comet_results_dir=comet_results_dir,
            top_n_psms=1,
        )
        protein_counts = get_protein_counts_from_comet_results(psms=psms)
        assert len(protein_counts) > 0


class Test_get_prefix_counts_by_length:
    @staticmethod
    def test_smoke():
        seqs = ["ACD", "ACDE", "CF"]
        prefix_counts_by_length = get_prefix_counts_by_length(seqs=seqs)
        expected = {
            1: {"A": 2, "C": 1},
            2: {"AC": 2, "CF": 1},
            3: {"ACD": 2},
            4: {"ACDE": 1},
        }
        assert prefix_counts_by_length == expected


class Test_get_and_plot_most_common_proteins:
    @staticmethod
    def test_smoke(tmp_path, test_data_dir):
        comet_results_dir = test_data_dir / "comet_results"
        q_value_threshold = 0.5
        get_and_plot_most_common_proteins(
            comet_results_dir=comet_results_dir,
            q_value_threshold=q_value_threshold,
            out_path=tmp_path / "top_10_proteins.txt",
        )
        assert (tmp_path / "top_10_proteins.txt").exists()
        assert (tmp_path / "protein_abundances.png").exists()

    @staticmethod
    def test_too_few_proteins(tmp_path, test_data_dir):
        get_and_plot_most_common_proteins(
            # comet_results_dir=test_data_dir / "too_few_proteins",
            comet_results_dir=Path("results/samples_3a_and_3b/native_run"),
            q_value_threshold=0.01,
            top_n_proteins=10,
            out_path=tmp_path,
        )
