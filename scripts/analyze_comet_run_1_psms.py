import sys
from dataclasses import dataclass
from pathlib import Path

repo_dir = Path(__file__).parents[1]
sys.path.append(str(repo_dir))

from src.comet_utils import CometRow, load_comet_data
from src.constants import RESULTS_DIR, IonTypes
from src.peptide_spectrum_comparison import (
    PeptideSpectrumComparison,
    compare_peptide_to_spectrum,
)
from src.peptides_and_ions import Peptide
from src.utils import pickle_and_compress, setup_logger

logger = setup_logger()

# Constants
testing = True
peak_filtering = True
top_n_peaks = 50
ppm_tolerance = 10
ion_types = [IonTypes.B_ION_TYPE, IonTypes.Y_ION_TYPE]

# Load comet data
comet_df = load_comet_data()
# comet_df.head(2)
comet_rows = CometRow.from_dataframe(comet_df)


@dataclass
class PSM:
    comet_row: CometRow
    psm: PeptideSpectrumComparison


psms = []
num_comet_rows = len(comet_rows)
for row_num, row in enumerate(comet_rows):
    logger.info(f"Processing comet row {row_num + 1} of {num_comet_rows}")
    # for row_num, row in df.iterrows():
    peptide = Peptide(seq=row.proposed_peptide)
    spectrum = row.get_corresponding_spectrum()

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
    psms.append(PSM(comet_row=row, psm=comparison))
    if testing:
        if row_num >= 2:
            break

if testing:
    output_path = (
        RESULTS_DIR
        / f"TESTING_comet_run_1_psms_ppmTol={ppm_tolerance}_peakFiltering={peak_filtering}_topN={top_n_peaks}.pkl"
    )
else:
    output_path = (
        RESULTS_DIR
        / f"comet_run_1_psms_ppmTol={ppm_tolerance}_peakFiltering={peak_filtering}_topN={top_n_peaks}.pkl"
    )

logger.info(f"Pickling/Saving results to {output_path}")
pickle_and_compress(obj=psms, file_path=output_path)
