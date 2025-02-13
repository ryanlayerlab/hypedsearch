from src.fasta_utils import get_proteins_from_fasta
from src.lookups.data_classes import Protein
from tests.fixtures_and_helpers import create_fasta


class TestGetProteinsFromFasta:
    @staticmethod
    def test_basic(tmp_path):
        fasta_path = tmp_path / "test.fasta"
        seqs = ["ATG", "CGT"]
        create_fasta(folder=tmp_path, seqs=seqs)
        expected_out = [
            Protein(protein_id=0, sequence="ATG", desc="seq1"),
            Protein(protein_id=1, sequence="CGT", desc="seq2"),
        ]

        actual = list(get_proteins_from_fasta(fasta_path=fasta_path))
        assert actual == expected_out
