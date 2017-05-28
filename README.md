# Intelligent Recombination Experiments
This repository contains scripts and results related to structure- and homology-directed protein recombination experiments conducted with the [Schreiter Lab](https://www.janelia.org/lab/schreiter-lab).

SCHEMA-RASPP is a protocol for structure-directed recombination using a parent structure and multiple homologous sequences. We have forked the original [CalTech repository](https://github.com/mattasmith/SCHEMA-RASPP) and initialized the forked version as a submodule in this repository (origin is [hbhargava/SCHEMA-RASPP](https://github.com/hbhargava7/SCHEMA-RASPP)). It is important to keep the submodule up to date by periodically running `git submodule foreach git pull origin master` to fetch the latest version.

# Usage
The Python script `compute_chimeras` takes as input a single directory containing a subdirectory for each homolog to be recombined (each subdirectory contains a FASTA sequence file and, optionally, a PDB structure file). Several additional parameters are specified by the user in the course of the computation. The overall output is a list of chimeras with corresponding SCHEMA energies and mutation scores.

To begin a computation, run the following:

```
python compute_chimeras.py -i /full/path/to/input/dir -o /full/path/to/output/dir
```
This command will initiate the chimera generation process. The steps are outlined below.

# Overview of compute_chimeras
1. Search input directory for potential parent sequences and structures and ask user which sequence is to be used for structure guidance and which are to be used as normal parent structures.
2. Build FASTA files for all selected sequences (including structure parent).
3. Compute sequence alignment of all parent sequences using ClustalOmega.
4. Build FASTA file with structure parent sequence and sequence from PDB file.
5. Compute sequence alignment for both parent sequences using ClustalOmega.
6. **[In progress]:** Compute the sequence identity for all possible pairs of homologs.
7. Copy parent PDB structure file to output folder for use by SCHEMA algorithm.
8. Use `SCHEMA_RASPP.schemacontacts` via `SR_interlink` to perform radial search and determine contacts present in structure.
9. Obtain desired number of crossovers from user and use `SCHEMA_RASPP.rasppcurve` via `SR_interlink` to compute the RASPP curve.
10. Obtain desired crossover locations from user and use `SCHEMA_RASPP.schemaenergy` via `SR_interlink` to compute SCHEMA energies and mutation scores for all chimeras.

## Dependencies
* `numpy` for miscellaneous computation
* `biopython` for FASTA parsing and more
* `clustalomega` for sequence alignment
* `picker` for user multi-selection
