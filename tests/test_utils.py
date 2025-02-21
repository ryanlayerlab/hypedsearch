from src.utils import Kmer, Position, generate_aa_kmers, get_positions_of_subseq_in_seq


class TestGenerateKmers:
    @staticmethod
    def test_default_min_k():
        aa_seq = "ACD"  # B isn't an allowed character
        max_k = 3
        expected_out = [
            Kmer(seq="A", position=Position(inclusive_start=0, exclusive_end=1)),
            Kmer(seq="C", position=Position(inclusive_start=1, exclusive_end=2)),
            Kmer(seq="D", position=Position(inclusive_start=2, exclusive_end=3)),
            Kmer(seq="AC", position=Position(inclusive_start=0, exclusive_end=2)),
            Kmer(seq="CD", position=Position(inclusive_start=1, exclusive_end=3)),
            Kmer(seq="ACD", position=Position(inclusive_start=0, exclusive_end=3)),
        ]
        actual = generate_aa_kmers(aa_seq=aa_seq, max_k=max_k)
        assert actual == expected_out

    @staticmethod
    def test_set_min_k():
        aa_seq = "ACD"  # B isn't an allowed character
        min_k, max_k = 2, 3
        expected_out = [
            Kmer(seq="AC", position=Position(inclusive_start=0, exclusive_end=2)),
            Kmer(seq="CD", position=Position(inclusive_start=1, exclusive_end=3)),
            Kmer(seq="ACD", position=Position(inclusive_start=0, exclusive_end=3)),
        ]
        actual = generate_aa_kmers(aa_seq=aa_seq, min_k=min_k, max_k=max_k)
        assert actual == expected_out


class TestGetPositionsOfSubseqInSeq:
    @staticmethod
    def test_multiple_instances():
        subseq = "AB"
        seq = "CABDDABDB"
        expected = [
            Position(inclusive_start=1, exclusive_end=3),
            Position(inclusive_start=5, exclusive_end=7),
        ]
        actual = get_positions_of_subseq_in_seq(subseq=subseq, seq=seq)
        assert actual == expected

    @staticmethod
    def test_no_instances():
        subseq = "XYZ"
        seq = "CABDDABDB"
        expected = []
        actual = get_positions_of_subseq_in_seq(subseq=subseq, seq=seq)
        assert actual == expected
