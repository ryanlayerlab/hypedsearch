from pathlib import Path
from typing import List

B_NEUTRAL_MASS_CALCULATOR = "src.peptides_and_ions.compute_b_ion_neutral_mass"
Y_NEUTRAL_MASS_CALCULATOR = "src.peptides_and_ions.compute_y_ion_neutral_mass"


def create_fasta(
    folder: Path,
    file_name: str = "test.fasta",
    seqs: List[str] = ["ATGCGTA", "CGTACGT"],
):
    fasta_path = folder / file_name
    lines = []
    for seq_num, seq in enumerate(seqs):
        lines.append(f">seq{seq_num + 1}\n")
        lines.append(f"{seq}\n")
    with fasta_path.open("w") as f:
        for line in lines:
            f.write(line)


def create_spectrum(
    scan_num: int,
    mz_array: List[float] = [1, 2, 3],
    intensity_array: List[float] = [4, 5, 6],
):
    return {
        "id": f"scan={scan_num}",
        "scanList": {"scan": [{"scan start time": 600}]},
        "m/z array": mz_array,
        "intensity array": intensity_array,
        "precursorList": {
            "precursor": [
                {
                    "selectedIonList": {
                        "selectedIon": [
                            {
                                "selected ion m/z": 100,
                                "charge state": 2,
                                "peak intensity": 200,
                            }
                        ]
                    }
                }
            ]
        },
    }
