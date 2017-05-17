# Intelligent Recombination Experiments
This repository contains scripts and results related to structure- and homology-directed protein recombination experiments.

## Current Todo List
* Add sequence identity calculation and display while selecting parents and children.

## SCHEMA-RASPP
SCHEMA-RASPP is a protocol for structure-directed recombination using a parent structure and multiple homologous sequences. We have forked the original [CalTech repository](https://github.com/mattasmith/SCHEMA-RASPP) and initialized the forked version as a submodule in this repository (origin is [hbhargava/SCHEMA-RASPP](https://github.com/hbhargava7/SCHEMA-RASPP)).

## generate_preconditions
This script is designed to generate the preconditions to run a SCHEMA-RASPP recombination experiment using any number of parent sequences.

### Usage
An input directory is required with subfolders for each potential parent sequence. At most, each subdirectory should contain one FASTA file and one pdb structure file (pdb is optional).

```
python generate_preconditions.py -i /full/path/to/input/dir -o /full/path/to/output
```

The output of this script are unaligned sequence files for all selected parents, aligned sequences for all parents, the parent-pdb unaligned fasta file, the parent-pdb aligned fasta file, and the parent pdb file. 

Subsequently, the following can be run from within the output directory to generate SCHEMA contacts:
```
python ../SCHEMA_RASPP/schemacontacts.py -pdb parent.pdb -msa allsequences_aligned.fasta -pdbal parent_aligned.fasta -o contacts.txt
```

## Dependencies
* `numpy` for miscellaneous computation
* `biopython` for FASTA parsing and more
* `clustalomega` for sequence alignment
