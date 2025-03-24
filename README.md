# Hypedsearch

## Installation

### Using Conda (recommended)

If you are using Conda for environment management, you can create and activate the environment using the provided environment.yaml file:

```bash
conda env create -f environment.yaml
conda activate hypedsearch
```

### Without using Conda (using `pip` and `venv`)

If you are not using Conda, you can:

1. ensure you have Python installed (preferably a version matching the one in `enivornment.yaml`),
2. create a virtual environment, 
3. manually extract the dependencies from `environment.yaml` into a `requirements.txt` file, and
4. install dependencies using pip:

```bash
python -m venv hypedsearch
source hypedsearch/bin/activate # Activate it (use `myenv\Scripts\activate` on Windows)
# manually extract the dependencies from environment.yaml to requirements.txt
pip install -r requirements.txt  # Install dependencies
```

## Usage

Make sure `hypedsearch` or the environment in which you installed the Hypedsearch dependencies is activated. 

```bash
python hypedsearch.py \
--mzml_dir tmp/test/spectra \
--mzml_path tmp/test/spectra/BMEM_AspN_Fxn5.mzML \
--output_dir tmp/test/output \
--db_path tmp/test/output/db.db \
--scan_num 9 \
--top_n_proteins 50 \
```

