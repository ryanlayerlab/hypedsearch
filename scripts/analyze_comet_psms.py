import sys
from pathlib import Path
from time import time

repo_dir = Path(__file__).parents[1]
sys.path.append(str(repo_dir))

from dataclasses import dataclass
from typing import Union

import pandas as pd
from tqdm import tqdm

from src.comet_utils import CometPSM, load_comet_data
from src.constants import PLAIN_PEPTIDE, RESULTS_DIR, SAMPLE, SCAN, IonTypes
from src.mass_spectra import get_specific_spectrum_by_sample_and_scan_num
from src.peptide_spectrum_comparison import (
    PeptideSpectrumComparison,
    compare_peptide_to_spectrum,
)
from src.peptides_and_ions import Peptide
from src.utils import (
    decompress_and_unpickle,
    make_directory,
    pickle_and_compress,
    setup_logger,
)

logger = setup_logger()

# Constants
testing = False
peak_filtering = False
top_n_peaks = 50
ppm_tolerance = 10
ion_types = [IonTypes.B_ION_TYPE, IonTypes.Y_ION_TYPE]


@dataclass
class PSM:
    comet_row: Union[CometPSM, pd.Series]
    psm: PeptideSpectrumComparison


def main():

    # comet_df = load_comet_data()
    comet_df = decompress_and_unpickle(
        "results/152_high_scoring_psms_eValQuantile<0.005_dCnQuantile>0.995.pkl"
    )
    # comet_rows = CometRow.from_dataframe(comet_df)
    # dir_name = f"comet_psms_ppmTol={ppm_tolerance}_peakFilter={int(peak_filtering)}"
    dir_name = f"high_scoring_comet_psms_ppmTol={ppm_tolerance}_peakFilter={int(peak_filtering)}"

    if peak_filtering:
        dir_name = f"{dir_name}_topNPeaks={top_n_peaks}"

    output_dir = RESULTS_DIR / dir_name
    make_directory(output_dir)

    psms = []
    # for row_num, row in enumerate(comet_rows):
    #     proposed_peptide = row.proposed_peptide
    # spectrum = row.get_corresponding_spectrum()
    num_comet_rows = comet_df.shape[0]
    for row_idx, row in comet_df.iterrows():
        start_time = time()
        logger.info(f"Processing comet row {row_idx + 1} of {num_comet_rows}")
        output_file = output_dir / f"row_idx={row_idx}.pkl"
        if output_file.exists():
            logger.info(f"{output_file} already exists! Skipping!")
            continue
        proposed_peptide = row[PLAIN_PEPTIDE]
        spectrum = get_specific_spectrum_by_sample_and_scan_num(
            sample=row[SAMPLE], scan_num=row[SCAN]
        )
        peptide = Peptide(seq=proposed_peptide)

        # Peak filtering
        if peak_filtering:
            spectrum.filter_to_top_n_peaks(n=top_n_peaks)
            assert spectrum.peaks_preprocessed

        comparison = compare_peptide_to_spectrum(
            peptide=peptide,
            spectrum=spectrum,
            ppm_tolerance=ppm_tolerance,
            ion_types=ion_types,
        )
        psm = PSM(comet_row=row, psm=comparison)
        pickle_and_compress(obj=psm, file_path=output_file)
        # psms.append(
        #     PSM(
        #         comet_row=row,
        #         psm=comparison
        #     )
        # )
        logger.info(f"Took {round(time()-start_time, 2)} seconds")
        if testing:
            if row_idx >= 2:
                break


if __name__ == "__main__":
    main()
