import os
import sys
from pathlib import Path

curr_dir = Path(os.getcwd())
assert curr_dir.stem == "hypedsearch", "Run this script from the hypedsearch folder"
sys.path.append(os.getcwd())

import pandas as pd

from src.utils import get_time_in_diff_units

folder = Path(sys.argv[1])
durations_in_seconds = []
for txt in folder.glob("*.txt"):
    df = pd.read_csv(txt, sep="\t")
    assert df.shape[0] == 1, f"Expected one row in {txt}, found {df.shape[0]}"
    durations_in_seconds.append(df["s"].iloc[0])
print(
    f"Sum of durations in folder: {get_time_in_diff_units(sum(durations_in_seconds))}"
)
print(
    f"Mean duration: {get_time_in_diff_units(sum(durations_in_seconds) / len(durations_in_seconds))}"
)
