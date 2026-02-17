# Developper notes 

So, what the hell is happening?
The function `form_spectrum_hybrids_via_clustering` forms the possible hybrids for a given spectrum.
And it returns a dictionary in which the keys are non-hyphenated hybrid sequences (e.g., "ABC") and the values are lists of `HybridPeptide` objects that correspond to the hybrid sequence (e.g., the `HybridPeptide` objects corresponding to "AB-C" and "A-BC").

The function `HypedsearchRunConfig.native_comet_run` performs the native Comet run. 

The function `HypedsearchRunConfig.hybrid_run_on_spectrum` performs the hybrid Comet run on a given spectrum. 
To run Comet on the hybrids, we need to create a FASTA file with hybrids in it which is done by the function as `create_hybrids_fasta`. 

## Python environment

Via micromamba

```bash
micromamba create -n hypedsearch python=3.11
micromamba activate hypedsearch
micromamba install -c bioconda biopython click lxml matplotlib numpy pandas pip plotly pydantic pymzml pyteomics pytest scipy seaborn snakemake -y
pip install fm-index
```

Via `uv`

```bash
uv sync
```

## Docker container

To create the image:

```bash
docker buildx build --platform linux/amd64 -t airikjohnson/hypedsearch:latest .
```

To push the image:

```bash
docker push airikjohnson/hypedsearch:latest
```

Start-up/Boot an instance of the container in interactive mode:

```bash
docker run -it --rm airikjohnson/hypedsearch:latest /bin/bash

docker run -it --rm -v ~/repos/hypedsearch/hypedsearch:/data airikjohnson/hypedsearch:latest /bin/bash
```

Testing the Docker container locally:

```bash
# Local machine
docker run -it --rm -v ~/repos/hypedsearch/hypedsearch:/hypedsearch airikjohnson/hypedsearch:latest /bin/bash

# In container
rm -rf tests/outputs
cd /hypedsearch
./src/run_hypedsearch.sh --cores 8 --config tests/data/test.hs.config.json
```

Fiji has Singularity not Docker.
To grab pull the Singularity version of a Docker image, run:

```bash
rm hypedsearch_latest.sif
singularity pull docker://airikjohnson/hypedsearch:latest
```

By default, you'll notice that inside of a container your home directory, /Users/<yourIdentikey>, is always available. In addition, the directory that you run a container from is also available. 

Then for testing run:

```bash
singularity run hypedsearch_latest.sif /bin/bash
```
