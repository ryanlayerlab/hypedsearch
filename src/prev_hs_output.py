import pandas as pd

from src.constants import GIT_REPO_DIR, SAMPLE


def load_hs_output():
    dir_path = GIT_REPO_DIR / "data/Hypedsearch_outputs"

    dfs = []
    for f in list(dir_path.glob("*")):
        sample = f.stem[3:]
        df = pd.read_csv(f, sep="\t")
        df[SAMPLE] = sample
        dfs.append(df)
    return pd.concat(dfs, ignore_index=True)
