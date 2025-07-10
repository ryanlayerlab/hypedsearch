# %%
import sys
from pathlib import Path

import pandas as pd

repo_dir = Path("/Users/erjo3868/repos/hypedsearch/hypedsearch")
sys.path.append(str(repo_dir))

# %%
# Crux Comet output
out_txt = repo_dir / "tmp/comet.1-10.txt"
crux_df = pd.read_csv(out_txt, sep="\t")


# Direct Comet output
out_txt = repo_dir / "tmp/BMEM_AspN_Fxn4.7-7.7-7.txt"
comet_df = pd.read_csv(out_txt, sep="\t", header=1)

crux_colms = set(crux_df.columns)
comet_colms = set(comet_df.columns)

print(f"Stuff in crux that's not in comet:\n{crux_colms - comet_colms}")
print(f"Stuff in crux that's not in comet:\n{comet_colms - crux_colms}")

crux_colms.union(comet_colms)
