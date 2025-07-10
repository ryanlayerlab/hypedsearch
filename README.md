# Hypedsearch <!-- omit in toc -->

- [Installation](#installation)
  - [Using Conda (recommended)](#using-conda-recommended)
  - [Without using Conda (using `pip` and `venv`)](#without-using-conda-using-pip-and-venv)
  - [Comet](#comet)
- [Usage](#usage)
- [(Potentially) Helpful commands](#potentially-helpful-commands)


## Installation

### Using Conda (recommended)

If you are using Conda for environment management, you can create and activate the environment using the provided `environment.yaml` file:

```bash
conda env create -f environment.yamlpr
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

### Comet

Hypedsearch depends on [Comet](https://comet-ms.sourceforge.net/). 
So make sure you have a Comet executable locally and can run Comet. 
If you're running Hypedsearch on a Mac, the executable `comet/comet.macos.exe` in this repo may work for you.

<!-- ## Usage

Make sure `hypedsearch` or the environment in which you installed the Hypedsearch dependencies is activated. 
Once the environment is activated, the following command should work and show the help page:

```bash 
python hypedsearch.py -h
```

Here's an example usage of Hypedsearch that should work after cloning the repo, following the installation requirements above, and, if needed, updating the paths to the Comet executable and `comet.params` file:

```bash
python hypedsearch.py \
--mzml_dir data/spectra \
--mzml_path data/spectra/BMEM_AspN_Fxn4.mzML \
--output_dir results/test \
--db_path results/test/test.db \
--scan_num 7 \
--top_n_proteins 50 \
--num_peaks 100 \
--comet_exe_path comet/comet.macos.exe \
--comet_params_path comet/comet.params \
--fasta_path fastas/Uniprot_mouse.fasta \
--cleanup False
``` -->

# Hypedsearch algorithm

1. ***Find mass spectra files***: 
find all `*.mzML` files in the `mzml_dir` directory and its subdirectories. 
As an example, say three files are found `Fxn1.mzML`, `Fxn2.mzML`, `Fxn3.mzML`.

1. ***Create file-specific subdirectories***: 
for each `<name>.mzML` found in step 1, create a subfolder called `<name>` under `output_dir`.
In our example, three subfolders under `output_dir` will be created called `Fxn1`, `Fxn2`, and `Fxn3`.

1. ***Comet run #1***: 
runs Comet on each `.mzML` file found in step 1 using the provided Comet executable `comet_exe_path` and Comet params file `comet_params_path` and using the provided FASTA file at `fasta_path` for the database of proteins that Comet searches.
The Comet result files will be named `run_1.txt` and `run_1.pep.xml` and will be saved in their `mzML` file's subfolder. 
E.g., `Fxn1.mzML`'s Comet results will be saved as `output_dir/Fxn1/run_1.txt` and `output_dir/Fxn1/run_1.pep.xml`.

1. ***Protein-Comet counts***:
next Hypedseach combines all the the Comet-found peptide-spectrum matches (PSMs) from all the `mzML` files and counts the number of times each protein appears in the "protein" column which will product data that looks like this:

  | protein | Comet count  |
  | - | - |
  | sp\|Q7TQP3\|GP119_MOUSE | 10 | 
  | sp\|P01326\|INS2_MOUSE | 25 | 
  | ... | ... | 

  We will use these Comet protein counts as a proxy for protein abundance. 

5. ***Create protein-product ion database***: 
an SQL database is created at the given `db_path`.
This database is what will be searched to find which theoretical product/fragment ions may correspond to a mass spectrum peak. 
Because we assume that hybrids only form between highly abundant proteins, only the most abundant proteins (as approximated by the Comet counts from step 4) and the product ions from those proteins will be included in the database.
Only the top N=`top_n_proteins` most abundant proteins are included in database.

6. ***Get the desired mass spectrum and perform peak filtering***: 
users specify which mass spectrum they'd like to search for hybrids for via specifying the spectrum's `mzML` file (via `mzml_path`) and scan number (`scan_num`).
Users can filter the peaks considered to the top N=`num_peaks` most intense peaks.

7. ***Find peak-matching product ions and form hybrids***: 
by searching the database created in step 5, find the product ions that may explain each peak in the spectrum. From the returned product ions, construct the possible hybrids that may correspond to the spectrum. [TODO: add more details on how this works elsewhere]

1. ***Create new FASTA file with the hybrids and most abundant proteins***:
create a FASTA file containing the proteins in the database from step 5 (i.e., the top N=`top_n_proteins` most abundant proteins) and the hybrids from step 7.
If `mzml_path=/path/to/BMEM_AspN_Fxn4.mzML` and `scan_num=7`, then the new FASTA file will be saved as `output_dir/BMEM_AspN_Fxn4/scan_num=7_hybrids.fasta`.

1. ***Comet run #2***:
run Comet again on the `.mzML` located at `mzml_path` using the new FASTA file from step 8.
The Comet run #2 result files will be named `scan_num=<scan_num>_run_2.txt` and `scan_num=<scan_num>_run_2.pep.xml` and will be saved in their `mzML` file's subfolder.
If `mzml_path=/path/to/BMEM_AspN_Fxn4.mzML` and `scan_num=7`, this step will produce these two files: `output_dir/BMEM_AspN_Fxn4/scan_num=7_run_2.txt` `output_dir/BMEM_AspN_Fxn4/scan_num=7_run_2.pep.xml`.

1.  ***Gather the spectrum-specific results***:
get the PSMs that correspond to the user-specified spectrum from the first Comet run and the second run, combine them into a new file called `scan=<scan_num>_hs_results.csv` which will be saved in the `mzML` file's subfolder.
The `.csv` file will have a `run` column.
The PSMs from Comet run #1 will have `run=1` and the PSMs from Comet run #2 will have `run=2`. 
If `mzml_path=/path/to/BMEM_AspN_Fxn4.mzML` and `scan_num=7`, this step generate the file: `output_dir/BMEM_AspN_Fxn4/scan_num=7_hs_results.csv`.

## Pipelines

- Run Comet and `assign-confidence`:
  ```
  snakemake -s snakefiles/run_comet.smk --configfile snakefiles/configs/human_samples.yaml --cores 2
  ```

## (Potentially) Helpful commands

Here are some commands that may be helpful for you. 

- Run Comet on all `.mzML` files in a directory:
  ```
  python -m src.comet_utils -h
  python -m src.comet_utils \
    -m data/spectra \
    -o results/comet_run_for_protein_abundances \
    -np 5
  ```
- Get the top N most abundant proteins:
  ```
  python -m src.protein_abundance -d results/human_samples/native_run -q 0.01 -t 10 -o results/human_samples/db/top_10_proteins.txt
  ```
- Get the PSM prefix abundances:
  ```
  python -m src.protein_abundance prefix-abundances --comet_results_dir results/human_samples/native_run --q_value_threshold 0.01 --out_path results/human_samples/native_run/prefix_abundances.json
  ```
- Create database of proteins and product ions
  ```
  python -m src.create_db -h
  python -m src.create_db \
    -d results/new_fasta/top_10_prots.db \
    -p results/new_fasta/top_10_proteins.txt \
    -f fastas/SwissProt.TAW_mouse_w_NOD_IAPP.fasta
  ```

- Run hypedsearch
  ```
  python -m src.run_hypedsearch -h 
  python -m src.run_hypedsearch \
    -m data/spectra/BMEM_AspN_Fxn4.mzML \
    -d results/new_fasta/top_10_prots.db \
    -o results/new_fasta/hs_results \
    -n 0 -p 10 -P 20
  ```

- Form all possible hybrids and then run Comet
  ```
  python -m src.form_all_hybrids -h
  python -m src.form_all_hybrids \
    -pn results/new_fasta/top_10_proteins.txt \
    -f fastas/SwissProt.TAW_mouse_w_NOD_IAPP.fasta \
    -m data/spectra/BMEM_AspN_Fxn5.mzML \
    -pt 20 \
    -s 9 \
    -o tmp
  ```

  Running on a Slurm-managed supercomputer:
  ```
  sbatch --array=9,10,15 slurm/form_all_hybrids.sbatch
  ```

- Running on Slurm-managed supercomputer:
