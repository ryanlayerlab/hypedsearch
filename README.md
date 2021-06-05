# hypedsearch
A powerful tool for identifying both cis-spliced and trans-spliced peptides
## **hy**brid **pe**pti**d**e **search**ing tool

For more information about installation and usage instructions, please see the documentation [found here](https://hypedsearch.readthedocs.io/en/latest/). 

---

## Overview
`hypedsearch` is a tool for identifying both hybrid and non-hybrid proteins from mass spectrometry data. `hypedsearch` takes `mzML` and `fasta` files as inputs and outputs a dataframe with a quality score.  `hypedsearch` identifies sequences by using a method called "k-mer extension".  `hypedsearch` is noted for its speed and accuracy compared to its peers and is actively being used for the exploration for the causes of Type One Diabetes (T1D).  For more informaton of T1D, please see below.  `hypedsearch` does not have to be confined to just T1D research and can be applied to other analysis of mass spectrometry data.  `hypedsearch` was created in the Layer Lab at University Of Colorado Boulder and is 100% open source with a MIT license.

## Mass Spectrometry
Mass Spectrometry is used to determine the components of a protein.  Proteins are composed of a series of intermingled peptide chains and these peptide chains, in turn, are composed of a series of amino acids.  A biosample is obtained, perpared by immursing it in acid to break it down into smaller components, and then fed into a mass spectrometer.  During a mass spectrometer run, the peptide chains are bombarded with hydrogen gas to break the peptide bonds, which isolates the individual amino acids.  In addition, the hydogen gas gives the individual amino acid a charge that can be weighed - this mass to charge ratio is measured in daltons. The amino acids of the sample can then be determined by comparing the mass/charges in the result set to known mass/charges.  The output of the mass spectrometer is a `mzML` file that can be fed into `hypedsearch`.  

## K-mer Extension
`hypedsearch` identifies sequences by using a method called "k-mer extension".  A "k-mer" is a k long string, or in this case, a k-long sequence of amino acids. The process, at a high level, works like this:
1. Pre-processing of the protein database to identify all k-mers for k in the range `min_peptide_len` to `max_peptide_len`. These k-mers are generated from both the `N` and `C` termini
2. Using the k-mers generated from the `N` termninus side, attempt to identify a sequence of amino acids that describe the `b` ions in the observed spectrum.
3. Repeat step 2, but from the `C` terminus side and try to describe the `y` ions in the observed spectrum
4. Filter out the poor scoring sequences
5. For the rest of the sequences, attempt to align and overlap the two sequences to the spectrum
6. If two sequences have no overlap, or do overlap but are from different proteins, report the alignment as hybrid
7. Save all alignments

## Application Architecture
//TODO

## T1D General Background
Type One Diabetes (T1D) is a chronic autoimmune disesase where the human body no longer can produce insulin.  Without insulin, the body cannot regulate its blood sugar and glucose levels, which can lead to immediate coma-induced death as well as long-term medical conditions like heart disease, kidney failure, blindness, and limb amputation.  The mental health toll on diabetics is equally as great - depression and suicide rates are much higher in diabetics than the overall population.  There is no know cure for diabetes  - most diabetes regulate their blood sugar via finges pricks, artifical insulin injections, and monitoring their diet.  Progress has been made in the development of a "closed loop" artifical pancreas where the diabetic uses a continual glucose monitor connected to an insulin pump to better regular their insulin levels.  However, management of the disease is still difficult and the financial toll of the disease can be catastrophic for families, even those with health insurance.

## T1D Biological Background
In a Type 1 Diabetic, the body's immune system attacks a region of the pancreas called the islets of Langerhans.  The islets of Langerhans contains beta cells that produce to body's insulin - there are an estimated one to three million beta cells in an adult human.  Insulin is a peptide hormone comprised of 51 amino acids in two chains connected by two disulfide bridges.  The creation of insulin by the beta cells is a multi-step process where a single peptide chain is translated in the ribosome (called preproinsulin), cleaved into proinsulin and then folded into its correct shape in the endoplasmic reticulum (RER), and then moved via the Golgi Network to where it finally forms into mature insulin.

In 2019, it was discoved that hybrid insulin peptides (HIPs) are highly concentrated in Type One Diabetics - and seem to trigger the body's autoimmune response which ultimately attacks the beeta cells.  The fact that insulin undergos extensive posttranslational modification aligns with this theory - somewhere in the process of creating mature insulin from preproinsulin a modification is introduced which creates the HIPs.

## Installation and Usage
First clone the repository
```bash
$> git clone https://github.com/zmcgrath96/hypedsearch.git
```
Then run the setup script that will install dependencies and build the C++ code.
```bash
$hypedsearch> ./setup.sh
```
If you get a permissions error, try the following:
```bash
$> cd hypedsearch
$hypedsearch> chmod u+x setup.sh
$hypedsearch> ./setup.sh
```
## Usage
There are two ways to use hypedsearch: command line arguments or the param file

**command line arguments**
In order to see the arguments, run
```bash
$> python3 -m src.main --help
usage: main.py [-h] [--spectra-folder SPECTRA_FOLDER]
           [--database-file DATABASE_FILE] [--output-dir OUTPUT_DIR]
           [--params PARAMS] [--min-peptide-len MIN_PEPTIDE_LEN]
           [--max-peptide-len MAX_PEPTIDE_LEN] [--tolerance TOLERANCE]
           [--precursor-tolerance PRECURSOR_TOLERANCE]
           [--peak-filter PEAK_FILTER]
           [--abundance-filter REL_ABUND_FILTER] [--digest DIGEST]
           [--verbose VERBOSE] [--cores CORES] [--n N]

Tool for identifying proteins, both hybrid and non hybrid from MS/MS data

optional arguments:
-h, --help            show this help message and exit
--spectra-folder SPECTRA_FOLDER
                        Path to folder containing spectra files.
--database-file DATABASE_FILE
                        Path to .fasta file containing proteins
--output-dir OUTPUT_DIR
                        Directory to save all figures. Default=~/
--params PARAMS       Use the params.py file adjacent to main.py instead of
                        using command line arguments. Default=False
--min-peptide-len MIN_PEPTIDE_LEN
                        Minimum peptide length to consider. Default=5
--max-peptide-len MAX_PEPTIDE_LEN
                        Maximum peptide length to consider. Default=20
--tolerance TOLERANCE
                        ppm tolerance to allow in search. Deafult=20
--precursor-tolerance PRECURSOR_TOLERANCE
                        ppm tolerance to accept when matching precursor
                        masses. Default=10
--peak-filter PEAK_FILTER
                        The number of peaks to take from a spectrum. The most
                        abundant peaks will be taken. Leave blank if you want
                        no filter or to use relative abundance filter.
                        Defualt=0
--abundance-filter REL_ABUND_FILTER
                        Take only peaks from a spectrum where the abundance of
                        the peak is >= the percentage give. Leave blank if you
                        want no filter or to use peak filter. Default=0.0
--digest DIGEST       The digest performed. Default=None
--verbose VERBOSE     Extra printing to console during run. Default=True
--cores CORES         The number of cores allowed to use when searching.
                        Uses at least 1 and at most the number of available
                        cores. Default=1
--n N                 The number of alignments to keep per spectrum.
                        Default=5
```
Just run the python3 -m src.main with your spectra folder, database file, output folder and any other parameters you want to include and hit enter.

Note
If you choose to use command line arguments, do not set the params flag to true. This will read from the params file. This is discussed in the next section.

**param file**

If you open up the params.py file in the src directory, you will be presented with a many different variables and descriptions. Read the descriptions and fill in your own parameters. Once you have done this, save the params file. Finally, in order to run using this params file, run the following:
```bash
$hypedsearch> python3 -m src.main --params True
```
## References

### Diabetes
* https://en.wikipedia.org/wiki/Type_1_diabetes
* https://en.wikipedia.org/wiki/Insulin

### Hybrid Insulin Peptides
* https://diabetes.diabetesjournals.org/content/68/9/1830

### pyteomics
* Goloborodko, A.A.; Levitsky, L.I.; Ivanov, M.V.; and Gorshkov, M.V. (2013) “Pyteomics - a Python Framework for Exploratory Data Analysis and Rapid Software Prototyping in Proteomics”, Journal of The American Society for Mass Spectrometry, 24(2), 301–304. DOI: 10.1007/s13361-012-0516-6

* Levitsky, L.I.; Klein, J.; Ivanov, M.V.; and Gorshkov, M.V. (2018) “Pyteomics 4.0: five years of development of a Python proteomics framework”, Journal of Proteome Research. DOI: 10.1021/acs.jproteome.8b00717
