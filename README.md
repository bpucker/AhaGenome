# AhaGenome
Collection of scripts associated with _Arabidopsis halleri_ genome assembly.

## Genome assembly process
Trimming of short reads prior to TrioBinning with Trimmomatic v0.39 (Bolger, 2014): 2:30:10; SLIDINGWINDOW:4:15 LEADING:5 TRAILING:5 MINLEN:50 TOPHRED33



## Scripts
ONT_scaffolding.py

in_silico_genetic_map_constructor.py

jcvi_wrapper_genes.py

phasing_check.py


## Removal of spurious markers
This scripts removed single markers which are apparently misplaced on the wrong linkage group. Input and output are BED files.

remove_spurious_markers.py

## Marker colinearity check
This script performs a nummeric comparison of genetic and physical marker positions to identify outliers. Groups of X outliers with a minimal genetic size of X are flagged.

nummeric_check.py

