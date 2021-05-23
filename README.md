# hypedsearch
## **hy**brid **pe**pti**d**e **search**ing tool

For more information about hypedsearch, installation and usage instructions, please see the documentation [found here](https://hypedsearch.readthedocs.io/en/latest/). 

---

## General Background
Type One Diabetes (T1D) is a chronic autoimmune disesase where the human body no longer can produce insulin.  Without insulin, the body cannot regulate its blood sugar and glucose levels, which can lead to immediate coma-induced death as well as long-term medical conditions like heart disease, kidney failure, blindness, and limb amputation.  The mental health toll on diabetics is equally as great - depression and suicide rates are much higher in diabetics than the overall population.  There is no know cure for diabetes  - most diabetes regulate their blood sugar via finges pricks, artifical insulin injections, and monitoring their diet.  Progress has been made in the development of a "closed loop" artifical pancreas where the diabetic uses a continual glucose monitor connected to an insulin pump to better regular their insulin levels.  However, management of the disease is still difficult and the financial toll of the disease can be catastrophic for families, even those with health insurance.

## Biological Background
In a Type 1 Diabetic, the body's immune system attacks a region of the pancreas called the islets of Langerhans.  The islets of Langerhans contains beta cells that produce to body's insulin - there are an estimated one to three million beta cells in an adult human.  Insulin is a peptide hormone comprised of 51 amino acids in two chains connected by two disulfide bridges.  The creation of insulin by the beta cells is a multi-step process where a single peptide chain is translated in the ribosome (called preproinsulin), cleaved into proinsulin and then folded into its correct shape in the endoplasmic reticulum (RER), and then moved via the Golgi Network to where it finally forms into mature insulin.

In 2019, it was discoved that hybrid insulin peptides (HIPs) are highly concentrated in Type One Diabetics - and seem to trigger the body's autoimmune response which ultimately attacks the beeta cells.  The fact that insulin undergos extensive posttranslational modification aligns with this theory - somewhere in the process of creating mature insulin from preproinsulin a modification is introduced which creates the HIPs.

## Computer Science Background
The state-of-the-art hardware to analyze peptide chains is a mass-spectroitry machine that bombards samples with X to break the peptide chain into its components.  These components are small - each amino acid is about x.  Due to the imprecise nature of the machine when analyzing a sample, a computer algorithm is used to identify the amino acids and then identify hybrid peptides compared to regular peptides - hypedsearch.

## hypedsearch
`hypedsearch` is a tool for identifying both hybrid and non-hybrid proteins from mass spectrometry data. `hypedsearch` takes `mzML` and `fasta` files as inputs and outputs a dataframe with a quality score.  `hypedsearch` identifies sequences by using a method called "k-mer extension".  

A "k-mer" is a k long string, or in this case, k long sequence of amino acids. The process, at a high level, works like this:
1. Pre-processing of the protein database to identify all k-mers for k in the range `min_peptide_len` to `max_peptide_len`. These k-mers are generated from both the `N` and `C` termini
2. Using the k-mers generated from the `N` termninus side, attempt to identify a sequence of amino acids that describe the `b` ions in the observed spectrum.
3. Repeat step 2, but from the `C` terminus side and try to describe the `y` ions in the observed spectrum
4. Filter out the poor scoring sequences
5. For the rest of the sequences, attempt to align and overlap the two sequences to the spectrum
6. If two sequences have no overlap, or do overlap but are from different proteins, report the alignment as hybrid
7. Save all alignments

## Setup

## Intsallation

## References

### Diabetes
* https://en.wikipedia.org/wiki/Type_1_diabetes
* https://en.wikipedia.org/wiki/Insulin

### Hybrid Insulin Peptides
* https://diabetes.diabetesjournals.org/content/68/9/1830

### pyteomics
* Goloborodko, A.A.; Levitsky, L.I.; Ivanov, M.V.; and Gorshkov, M.V. (2013) “Pyteomics - a Python Framework for Exploratory Data Analysis and Rapid Software Prototyping in Proteomics”, Journal of The American Society for Mass Spectrometry, 24(2), 301–304. DOI: 10.1007/s13361-012-0516-6

* Levitsky, L.I.; Klein, J.; Ivanov, M.V.; and Gorshkov, M.V. (2018) “Pyteomics 4.0: five years of development of a Python proteomics framework”, Journal of Proteome Research. DOI: 10.1021/acs.jproteome.8b00717
