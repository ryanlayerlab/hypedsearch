# Hypedsearch <!-- omit in toc -->
- [Installation via Conda](#installation-via-conda)
- [Overview](#overview)
- [Code organization](#code-organization)
  - [Run Crux](#run-crux)
  - [Create the kmer database (and kmer-to-protein map)](#create-the-kmer-database-and-kmer-to-protein-map)
  - [Form hybrids](#form-hybrids)
    - [Running via `snakemake`](#running-via-snakemake)
  - [Run Hypedsearch](#run-hypedsearch)
- [Running `msconvert` to convert files to MZML, subset spectra in MZML files, etc.](#running-msconvert-to-convert-files-to-mzml-subset-spectra-in-mzml-files-etc)
- [Hypedsearch Docker \& Singularity container](#hypedsearch-docker--singularity-container)


## Installation via Conda

If you are using Conda for environment management, you can create and activate the environment using the provided `environment.yaml` file:

```bash
conda env create -f environment.yaml
conda activate hypedsearch
```

## Overview

Blah.

## Code organization

### Run Crux

```
python -m src.crux run-comet -m data/spectra/mouse_samples/BMEM_AspN_Fxn9.mzML -m data/spectra/mouse_samples/BMEM_AspN_Fxn8.mzML -f fastas/SwissProt.TAW_mouse_w_NOD_IAPP.fasta -c comet/crux.comet.params -n 8 -o tmp/aug18 -m data/spectra/mouse_samples
```

### Create the kmer database (and kmer-to-protein map)

This code lives in `kmer_database.py`. 
To see the help message: `python -m src.kmer_database -h`.
Here's an example usage:

```

```
python -m src.kmer_database \
  --min_k 1 \
  --max_k 25 \
  --protein_names results/mouse_samples/db/top_10_proteins.txt
  --fasta fastas/SwissProt.TAW_mouse_w_NOD_IAPP.fasta \
  --kmer_to_protein tmp/mouse_top_10_proteins.json \
  --db_path tmp/mouse_top_10_proteins.db \
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

### Run Hypedsearch

Via Python from the command-line: 

```
python -m src.hypedsearch_utils \
  --mzml tests/data/mouse_spectra.mzML \
  --scan 7 \
  --database tests/data/mouse_top_10_proteins.db \
  --fasta tests/data/mouse_proteome_SwissProt.TAW_mouse_w_NOD_IAPP.fasta \
  --crux_comet_params tests/data/crux.comet.params \
  --out_dir tmp
```

Via Snakemake:

```
snakemake -s snakefiles/run_hypedsearch.smk \
  --configfile snakefiles/configs/test_run_hypedsearch.yaml \
  --cores 4
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

## Hypedsearch Docker & Singularity container

Build the Docker container:

```bash
cd crux_docker
docker buildx build --platform linux/amd64 -t airikjohnson/hypedsearch:latest .
docker push airikjohnson/hypedsearch:latest
```

Launch a Docker container and make sure that Hypedsearch runs in the Docker container:

```bash
docker run -it --rm -v ~/repos/hypedsearch/hypedsearch/:/data airikjohnson/hypedsearch:latest

python -m src.kmer_database -d db.db -f /data/fastas/SwissProt.TAW_mouse_w_NOD_IAPP.fasta -ktp ktp.json -pn /data/results/mouse_samples/db/top_10_proteins.txt
python -m src.hypedsearch_utils -m /data/data/spectra/mouse_samples/BMEM_AspN_Fxn4.mzML -s 7 -db /tmp/db.db -f /data/fastas/SwissProt.TAW_mouse_w_NOD_IAPP.fasta -c /data/comet/crux.comet.params -o /tmp
```

If you need a singularity container: `singularity pull docker://airikjohnson/hypedsearch:latest`.

Testing the singularity container:

```bash
singularity run \
  --bind src:/hs/src \
  --bind tmp/singularity_test:/hs/results \
  --bind tests/data/mouse_spectra.mzML:/hs/spectra.mzML \
  --bind fastas/SwissProt.TAW_mouse_w_NOD_IAPP.fasta:/hs/proteome.fasta \
  --bind tests/data/mouse_top_10_proteins.db:/hs/kmers.db \
  --bind comet/crux.comet.params:/hs/crux.comet.params \
  --pwd /hs \
  hypedsearch_latest.sif /bin/bash
python -m src.hypedsearch_utils -m spectra.mzML -s 7 -db kmers.db -f proteome.fasta -c crux.comet.params -o results
```