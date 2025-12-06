from pathlib import Path

from pydantic import BaseModel

from src.comet_utils import CometPSM
from src.peptides_and_ions import Fasta


class HybridFindingSimulation(BaseModel):
    

def remove_native_seq_and_run_hypedsearch(
    fasta: Path,
    psm: CometPSM,
    db_path: Path,
    min_side_len: int,
):
    fasta = Fasta(path=fasta)
    spectrum = psm.get_spectrum()
    seq = psm.seq

    # Make sure the PSM sequence appears in only one protein and in one location
    prot_name = list(fasta.proteins_that_contain_seqs([seq])[seq])
    assert len(prot_name) == 1
    prot_name = prot_name[0]
    prot_seq = fasta.protein_name_to_seq_map[prot_name]
    assert prot_seq.count(seq) == 1

    # Create the hybrid sequence by splitting the protein sequence at the PSM sequence location
    left_hy_seq, right_hy_seq = seq[: int(len(seq) / 2)], seq[int(len(seq) / 2) :]
    if len(left_hy_seq) < min_side_len or len(right_hy_seq) < min_side_len:
        raise RuntimeError("PSM sequence too short on one side")
    split_seq = prot_seq.split(seq)
    assert len(split_seq) == 2
    new_prot_seq = (
        split_seq[0][: int(len(split_seq[0]) / 2)]
        + left_hy_seq
        + split_seq[0][int(len(split_seq[0]) / 2) :]
        + split_seq[1][: int(len(split_seq[1]) / 2)]
        + right_hy_seq
        + split_seq[1][int(len(split_seq[1]) / 2) :]
    )
    new_prot = Peptide(seq=new_prot_seq, name=prot_name)
    new_proteins_for_fasta = list(
        filter(lambda prot: prot.name != prot_name, fasta.proteins)
    ) + [new_prot]
    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_dir = Path(tmp_dir)
        new_fasta = tmp_dir / "modified_proteome.fasta"
        Fasta.write_fasta(peptides=new_proteins_for_fasta, path=new_fasta)
        # Create new k-mer database
        new_db_path = tmp_dir / "kmers.db"
        proteins = list(KmerDatabase(db_path=db_path).proteins)
        kmer_db_prots = Fasta(path=new_fasta).get_proteins_by_name(names=proteins) + [
            new_prot
        ]
        kmer_db = KmerDatabase.create_db(
            db_path=new_db_path,
            proteins=kmer_db_prots,
            overwrite=True,
        )
        # Create Hypedsearch config
        hs_config = {
            "mzml_to_scans": {f"{spectrum.mzml}": [spectrum.scan]},
            "parent_out_dir": tmp_dir,
            "hybrid_former": {
                "kmer_db": new_db_path,
                "fasta": new_fasta,
            },
            "psm_scorer": {
                "comet_params": "results/251028_RP_Islet_Spikes_Crashout/inputs/crux.comet.params",
                "fasta": new_fasta,
            },
        }
        hs_config = HypedsearchRunConfig(**hs_config)
        n_out, h_out = hs_config.run_hypedsearch_on_spectrum(
            spectrum=spectrum,
        )
        n_psms = CometPSM.from_txt(txt=n_out.target)
        assert n_psms[0].seq != seq, "Native PSM should not match the original sequence"
        h_psms = CometPSM.from_txt(txt=h_out.target)
    return n_psms, h_psms, hs_config
