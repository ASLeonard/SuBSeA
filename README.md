[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Python application](https://github.com/ASLeonard/SuBSeA/workflows/Python%20application/badge.svg)](https://github.com/ASLeonard/SuBSeA/actions?query=workflow%3A%22Python+application%22)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/3378fad4f0174fffb2170806acb68af7)](https://www.codacy.com?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=ASLeonard/SuBSeA&amp;utm_campaign=Badge_Grade)

# *Su*bunit *B*inding *Se*quence *A*lignment (SuBSeA)

Software to test qualitative hypotheses generated from the [polyomino evolution through duplication study](https://github.com/ASLeonard/duplication "Polyomino duplication repository").

## Purpose

## Install
This software requires Python 3.6+ and several common packages listed in requirements.txt. In addition, a working version of needle from [Emboss](http://emboss.sourceforge.net/download/) is necessary. Other implementations of needle have not been tested, but should work provided a similar output is achievable.

### Testing


## Examples
A simple example can be run through
```python
python SubSeA.py 3WJM A 3WJM B
```
which calculates the SuBSeA confidence between the two heteromeric interfaces of the protein complex 3WJM


## Scope

### Limitations
Not all protein complexes are stored in standard formats. Particularly there are often conflicts between the PDB and PDBePISA with regards to quaternary structure and active interactions. When there are issues in compatability, it is often the case that certain interactions are calculated incorrectly, which can provide a meaningless result with no alignment.
