import ast
from pathlib import Path

import click
import pandas as pd

from src.hypedsearch_run_analysis import SpectrumPSMs
from src.utils import PathType, setup_logger


@click.command(
    name="process",
    context_settings={"help_option_names": ["-h", "--help"], "max_content_width": 200},
    help="""
    """,
)
@click.option(
    "--jct_csv",
    "-j",
    type=PathType(),
    required=True,
    help="",
)
@click.option(
    "--psms_json",
    "-p",
    type=PathType(),
    required=True,
    help="",
)
@click.option(
    "--results_dir",
    "-r",
    type=PathType(),
    required=True,
    help="",
)
def cli_process(jct_csv: Path, psms_json: Path, results_dir: Path):
    results_dir.mkdir(parents=True, exist_ok=True)
    psms = {psms.spectrum.uid: psms for psms in SpectrumPSMs.load(path=psms_json)}
    jct_df = pd.read_csv(jct_csv)
    # jct_df["mzml"] = jct_csv.name.split(".")[0]
    data = []
    for row_idx, row in jct_df[jct_df.num_uniq_seqs >= 3].iterrows():
        # row
        spectra_supporting = ast.literal_eval(row.spectra_supporting)
        uniq_seqs = list(ast.literal_eval(row.psm_seq_cnter).keys())
        spectra_supporting
        max_intensity = -1
        rts = []
        for spectrum_uid in spectra_supporting:
            rts.append(round(psms[spectrum_uid].spectrum.retention_time, 5))
            if psms[spectrum_uid].spectrum.precursor_intensity > max_intensity:
                max_intensity = psms[spectrum_uid].spectrum.precursor_intensity
                max_intensity_psm = psms[spectrum_uid]
        rt_range = f"{min(rts)}-{max(rts)}"
        data.append(
            (
                row.jct,
                row.num_psm_supporting,
                row.num_uniq_seqs,
                uniq_seqs,
                max_intensity_psm.hybrid_target.seq,
                max_intensity_psm.spectrum.uid,
                max_intensity_psm.spectrum.precursor_intensity,
                max_intensity_psm.spectrum.retention_time,
                rts,
                rt_range,
            )
        )

    data = pd.DataFrame(
        data,
        columns=[
            "jct",
            "num_spectra_supporting",
            "num_uniq_seqs",
            "uniq_seqs",
            "top_hybrid_seq_by_precursor_intensity",
            "top_spectrum_by_precursor_intensity",
            "top_spectrums_precursor_intensity",
            "top_spectrums_retention_time",
            "retention_times",
            "retention_time_range",
        ],
    )
    data.to_csv(
        results_dir / f"{jct_csv.stem}.processed.csv",
        index=False,
    )


@click.group(
    context_settings={"help_option_names": ["-h", "--help"], "max_content_width": 200}
)
def cli():
    pass


if __name__ == "__main__":
    setup_logger()
    cli.add_command(cli_process)
    cli()
