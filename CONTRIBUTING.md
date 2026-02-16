# Developper notes 

So, what the hell is happening?
The function `form_spectrum_hybrids_via_clustering` forms the possible hybrids for a given spectrum.
And it returns a dictionary in which the keys are non-hyphenated hybrid sequences (e.g., "ABC") and the values are lists of `HybridPeptide` objects that correspond to the hybrid sequence (e.g., the `HybridPeptide` objects corresponding to "AB-C" and "A-BC").

The function `HypedsearchRunConfig.native_comet_run` performs the native Comet run. 

The function `HypedsearchRunConfig.hybrid_run_on_spectrum` performs the hybrid Comet run on a given spectrum. 
To run Comet on the hybrids, we need to create a FASTA file with hybrids in it which is done by the function as `create_hybrids_fasta`. 

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
docker run -it --rm -v ~/repos/hypedsearch/hypedsearch:/data airikjohnson/hypedsearch:latest /bin/bash
```