[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Python application](https://github.com/ASLeonard/SuBSeA/workflows/Python%20application/badge.svg)](https://github.com/ASLeonard/SuBSeA/actions?query=workflow%3A%22Python+application%22)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/3378fad4f0174fffb2170806acb68af7)](https://www.codacy.com?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=ASLeonard/SuBSeA&amp;utm_campaign=Badge_Grade)

# *Su*bunit *B*inding *Se*quence *A*lignment (SuBSeA)

Software to test qualitative hypotheses generated from the [polyomino evolution through duplication study](https://github.com/ASLeonard/duplication "Polyomino duplication repository").

This program compares the binding residues from macromolecular interfaces with optimal sequence alignment to estimate the likelihood that two protein complex interactions are related. A more detailed description of the analysis can be found in the methods [here](https://www.biorxiv.org/content/10.1101/2020.04.22.054783v1).

## Components
There are several main components to the analysis pipeline, outlined below.

  - Protein complex dataset generation
   - utility.py
  - Bioinformatic data pulling
   - domains.py
   - pisa_XML.py
  - SuBSeA analysis
   - periodic.py
   - SuBSeA.py
  - Visualisation
   - pdb_visualisation.py

## Install

This software has been tested on Python 3.6+ and several common packages listed in requirements.txt.

In addition, a working version of needle from [Emboss](http://emboss.sourceforge.net/download/) is necessary. Other implementations of needle have not been tested, but should work provided a similar output is achievable.

### Testing

Functionality can be tested by running the following command.
```shell
python -m pytest
```
Errors at this stage are likely due to a missing needle exectuable or required python packages.

## Usage examples 
A simple example can be run by providing the two subunits to compare.
```python
python binding\_alignment.py 3WJM A 3WJM B
```
Which calculates the SuBSeA confidence between the two heteromeric interfaces of the protein complex 3WJM.

If the interaction under examination is not isologous, alternate chains can be provided for comparison.
```python
python binding\_alignment.py 2IX2 A 2IX2 B --alternate_chains C A
```
Which runs the comparison of the interactions between chains A->C with the interaction between chains B->A.

Interfaces can also be compared across subunits, such as analysing homomeric precursors.
```python
python binding\_alignment.py 15C8 L 4OFD A --alternate_chains H B
```
Again which compared the interaction between 15C8 chains L->H and 4ODF chains A->B.

A larger scale analysis can be conducted with 
```python
python pipeline\_runner.py --pullINT
```
Which will compile the heteromeric comparisons needed from the dataset, and automatically download any associated files.

## Limitations
Not all protein complexes are stored in standard formats. Particularly, there are often conflicts between the PDB and PDBePISA with regards to quaternary structure and active interactions. When there are issues in compatability, it is often the case that certain interactions are calculated incorrectly, which can provide a meaningless result with no alignment.

The full analysis requires

  - FASTA sequence
  - PDBePISA macromolecular interfaces
  - CATH domains (or other homology identifier)

Any protein complex with incomplete data will struggle in this analysis, so only a subset of recorded proteins can be used correctly.
