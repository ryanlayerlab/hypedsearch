from collections import Counter
from pathlib import Path
from unittest.mock import MagicMock, patch

from src.protein_abundance import get_protein_comet_counts


class Test_get_protein_comet_counts:
    @staticmethod
    def test_no_filtering():
        with patch("pathlib.Path.glob") as mock_glob, patch(
            "src.protein_abundance.CometPSM.from_txt"
        ) as mock_from_txt:
            mock_glob.return_value = ["test", "test"]
            mock_psm_1 = MagicMock()
            mock_psm_1.proteins = ["P1", "P2"]
            mock_psm_2 = MagicMock()
            mock_psm_2.proteins = ["P2", "P3"]
            mock_psm_3 = MagicMock()
            mock_psm_3.proteins = ["P3", "P4"]
            mock_from_txt.side_effect = [[mock_psm_1, mock_psm_2], [mock_psm_3]]
            result = get_protein_comet_counts(comet_results_dir=Path("dummy"))
            assert result == Counter(["P1", "P2", "P2", "P3", "P3", "P4"])

    @staticmethod
    def test_filtering():
        with patch("pathlib.Path.glob") as mock_glob, patch(
            "src.protein_abundance.CometPSM.from_txt"
        ) as mock_from_txt:
            mock_glob.return_value = ["test", "test"]
            mock_psm_1 = MagicMock()
            mock_psm_1.proteins = ["P1", "P2"]
            mock_psm_1.num = 1
            mock_psm_2 = MagicMock()
            mock_psm_2.proteins = ["P2", "P3"]
            mock_psm_2.num = 3
            mock_psm_3 = MagicMock()
            mock_psm_3.proteins = ["P3", "P4"]
            mock_psm_3.num = 2
            mock_from_txt.side_effect = [[mock_psm_1, mock_psm_2], [mock_psm_3]]
            result = get_protein_comet_counts(
                comet_results_dir=Path("dummy"), top_n_psms=2
            )
            assert result == Counter(["P1", "P2", "P3", "P4"])
