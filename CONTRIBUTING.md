# Development notes


## Code organization

- `run_comet.py` -- stuff related to running Comet
- Create hybrids: 
  - `form_hybrids.py`
  - `create_hybrids.smk`

- Process an MZML where "process" means (1) run Comet on the MZML, (2) form hybrids, and (3) run Comet with a separate decoy run on the hybrids+FASTA.
  - `python -m src.process_mzml -m tests/data/mouse_spectra.mzML -f tests/data/mouse_proteome.fasta -o tmp/june12 -d tests/data/mouse_top_10_prots.db `

## Hybrid formation

Given a spectrum and a path to a database, the method `form_hybrids.form_hybrids_for_spectrum` (1) forms hybrids for the spectrum, (2) removes native sequences, and, (3) if `out_path` is set, saves the hybrids. 
`peptides_and_ions.Fasta.contains_seq` is used To remove native sequences.
For each hybrid sequence, the FASTA file is parsed to a `@cached_property` (so I think it's truly parsed only once) and the sequences it contains are looped over to check if the hybrid sequence is in the FASTA sequence. 
TODO: I'm not exactly sure how `@cached_property` works but, best case, for each hybrid sequence, all the sequences in the FASTA are looped over. So best, worst-case time is `O(|H|*|F|*c)` where `|H|` is the number of hybrid sequences and `|F|` is the number of FASTA sequences and `c` is the time it takes to check of a string is a substring of another string.

- Get most abundant proteins: 

## Important functions

## Running 

## Load mass spectra


- Assign confidence:

```
crux assign-confidence \
--output-dir results/mouse_samples/confident_psms \
--list-of-files T \
results/mouse_samples/hybrid_run/BMEM_AspN_Fxn4.comet.target.txt \
results/mouse_samples/hybrid_run/BMEM_AspN_Fxn5.comet.target.txt \
results/mouse_samples/hybrid_run/BMEM_AspN_Fxn6.comet.target.txt \
results/mouse_samples/hybrid_run/BMEM_AspN_Fxn7.comet.target.txt \
results/mouse_samples/hybrid_run/BMEM_AspN_Fxn8.comet.target.txt \
results/mouse_samples/hybrid_run/BMEM_AspN_Fxn9.comet.target.txt 


crux assign-confidence \
--output-dir results/mouse_samples/native_run/confident_psms \
--list-of-files T \
results/mouse_samples/native_run/BMEM_AspN_Fxn4.comet.target.txt \
results/mouse_samples/native_run/BMEM_AspN_Fxn5.comet.target.txt \
results/mouse_samples/native_run/BMEM_AspN_Fxn6.comet.target.txt \
results/mouse_samples/native_run/BMEM_AspN_Fxn7.comet.target.txt \
results/mouse_samples/native_run/BMEM_AspN_Fxn8.comet.target.txt \
results/mouse_samples/native_run/BMEM_AspN_Fxn9.comet.target.txt 
```

## To do

- [ ] Give HybridPeptide objects have a `scan` attribute, populate it and propogate it to wherever hybrids are sold
- [ ] 