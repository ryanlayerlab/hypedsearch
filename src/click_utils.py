from pathlib import Path

import click

from src.constants import (
    COMET_PARAMS,
    DEFAULT_MAX_KMER_LEN,
    DEFAULT_MIN_KMER_LEN,
    DEFAULT_PPM_TOLERANCE,
)
from src.utils import get_default_comet_executable_path


class PathType(click.ParamType):
    name = "path"

    def convert(self, value, param, ctx):
        try:
            return Path(value).resolve()  # Convert to absolute Path object
        except Exception as e:
            self.fail(f"{value} is not a valid path: {e}", param, ctx)


class ClickOptions:
    def protein_names(self):
        return click.option(
            "--protein_names",
            "-pn",
            type=PathType(),
            required=True,
            help="Path to file containing the names of proteins to form hybrids from",
        )

    def fasta_path(self):
        return click.option(
            "--fasta_path",
            "-f",
            type=PathType(),
            required=True,
            help="Path to the FASTA file.",
        )

    def mzml(self):
        return click.option(
            "--mzml",
            "-m",
            type=PathType(),
            required=True,
            help="Path to the MZML file.",
        )

    def output_dir(self):
        return click.option(
            "--output_dir",
            "-o",
            type=PathType(),
            required=True,
            help="The outputs will be saved here",
        )

    def comet_exe(self):
        return click.option(
            "--comet_exe",
            "-ce",
            type=PathType(),
            default=get_default_comet_executable_path(),
            show_default=True,
            help="Path to the comet executable.",
        )

    def db_path(self):
        return click.option(
            "--db_path",
            "-d",
            type=PathType(),
            required=True,
            help="Path to the product-ion database file.",
        )

    def comet_params(self):
        return click.option(
            "--comet_params",
            "-cp",
            type=PathType(),
            default=COMET_PARAMS,
            show_default=True,
            help="Path to the comet.params file.",
        )

    def num_peaks(self):
        return click.option(
            "--num_peaks",
            type=int,
            default=0,
            show_default=True,
            help=(
                "Only the num_peaks most intense peaks in the spectra will be considered. "
                "If num_peaks=0, all peaks will be considered"
            ),
        )

    def precursor_mz_ppm_tol(self):
        return click.option(
            "--precursor_mz_ppm_tol",
            "-pmpt",
            default=DEFAULT_PPM_TOLERANCE,
            show_default=True,
            help=(
                "Precursor m/z PPM tolerance. Hybrids will be within this PPM of the precursor m/z"
            ),
        )

    def peak_to_ion_ppm_tol(self):
        return click.option(
            "--peak_to_ion_ppm_tol",
            "-pipt",
            default=DEFAULT_PPM_TOLERANCE,
            show_default=True,
            help=(
                "Precursor m/z PPM tolerance. Hybrids will be within this PPM of the precursor m/z"
            ),
        )

    def min_k(self):
        return click.option(
            "--min_k",
            "-mk",
            type=int,
            default=DEFAULT_MIN_KMER_LEN,
            show_default=True,
            help="Minimum kmer length to consider.",
        )

    def max_k(self):
        return click.option(
            "--max_k",
            "-Mk",
            type=int,
            default=DEFAULT_MAX_KMER_LEN,
            show_default=True,
            help="Maximum kmer length to consider.",
        )

    def num_psms(self):
        return click.option(
            "--num_psms",
            "-np",
            type=int,
            default=20,
            show_default=True,
            help="Number of Comet PSMs to spit out per spectrum",
        )

    def scan(self):
        return click.option(
            "--scan",
            "-s",
            type=int,
            default=0,
            show_default=True,
            help=(
                "If this is set, then the workflow will only run on the spectrum with this scan number "
                "rather than all the spectra in the MZML which is what will happen if scan=0"
            ),
        )

    def overwrite(self):
        return click.option(
            "--overwrite",
            "-ow",
            is_flag=True,
            help="If outputs already exist, this controls whether or not to overwrite them.",
        )

    def cli_options(self):
        """Return a decorator that applies all options"""

        def decorator(func):
            for opt in reversed(
                [
                    self.scan(),
                    self.num_psms(),
                    self.max_k(),
                    self.min_k(),
                    self.precursor_mz_ppm_tol(),
                    self.num_peaks(),
                    self.output_dir(),
                    self.mzml(),
                    self.fasta_path(),
                    self.protein_names(),
                ]
            ):
                func = opt(func)
            return func

        return decorator
