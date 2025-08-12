# Hypedsearch <!-- omit in toc -->
- [Installation via Conda](#installation-via-conda)
- [Overview](#overview)
- [Code organization](#code-organization)
  - [Create the kmer database (and kmer-to-protein map)](#create-the-kmer-database-and-kmer-to-protein-map)
  - [Form hybrids](#form-hybrids)
    - [Running via `snakemake`](#running-via-snakemake)
- [Running `msconvert` to convert files to MZML, subset spectra in MZML files, etc.](#running-msconvert-to-convert-files-to-mzml-subset-spectra-in-mzml-files-etc)


## Installation via Conda

If you are using Conda for environment management, you can create and activate the environment using the provided `environment.yaml` file:

```bash
conda env create -f environment.yaml
conda activate hypedsearch
```

## Overview

Blah.

## Code organization

### Create the kmer database (and kmer-to-protein map)

This code lives in `kmer_database.py`. 
To see the help message: `python -m src.kmer_database -h`.
Here's an example usage:

```
python -m src.kmer_database \
  --min_k 1 \
  --max_k 25 \
  --fasta fastas/SwissProt.TAW_mouse_w_NOD_IAPP.fasta \
  --kmer_to_protein tmp/mouse_top_10_proteins.json \
  --db_path tmp/mouse_top_10_proteins.db \
  --protein_names results/mouse_samples/db/top_10_proteins.txt
```

### Form hybrids

This code lives in `hybrids_via_clusters.py`.
To see the help message: `python -m src.hybrids_via_clusters -h`. 
Here's an example usage in which we form hybrids for a single scan:

```
python -m src.hybrids_via_clusters \
  --mzml data/spectra/mouse_samples/BMEM_AspN_Fxn4.mzML \
  --fasta fastas/SwissProt.TAW_mouse_w_NOD_IAPP.fasta \
  --scan 1 \
  --database tmp/mouse_top_10_proteins.db \
  --precursor_mz_ppm_tol 20 \
  --peak_to_ion_ppm_tol 20 \
  --out_dir tmp
```

#### Running via `snakemake`

Here's a command that you can run to make sure that you can form hybrids via snakemake:

```
snakemake -s snakefiles/form_hybrids.smk --configfile snakefiles/tests/test_form_hybrids.yaml --cores 4
```

## Running `msconvert` to convert files to MZML, subset spectra in MZML files, etc.

Run [the following Docker container](https://hub.docker.com/r/proteowizard/pwiz-skyline-i-agree-to-the-vendor-licenses), mounting any directories that you'd like available in the Docker container. 


```
docker pull proteowizard/pwiz-skyline-i-agree-to-the-vendor-licenses:skyline_daily_25.1.1.174-b787f12

docker run -it --rm \
  -v /Users/erjo3868/repos/hypedsearch/hypedsearch/data/spectra/mouse_samples:/data \
  proteowizard/pwiz-skyline-i-agree-to-the-vendor-licenses:skyline_daily_25.1.1.174-b787f12 \
  /bin/bash
```

Note that the tag=`skyline_daily_25.1.1.174-b787f12` image worked on my arm64 processor Mac.
Other tags did not work.
You may need to play around to find an image that works on your computer.

Then run [`msconvert`](https://proteowizard.sourceforge.io/tools/msconvert.html) via [`wine`](https://gitlab.winehq.org/wine/wine) (which is a program that allows one to run Windows programs on Unix).
Here are some examples:

```
# msconvert command-line help
wine msconvert --help

# generate an mzML containing a subset of another mzML's spectra
wine msconvert BMEM_AspN_Fxn4.mzML --filter "index [0,19]" --outfile BMEM_AspN_Fxn4_scans1-20.mzML -o subset_spectra_for_testing

# convert raw mass spectrometry data to mzML format 
msconvert <input_file(s)> -o <output_directory> --mzML
```
<!-- 

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

- Running on Slurm-managed supercomputer: -->