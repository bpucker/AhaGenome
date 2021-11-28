# AhaGenome
Collection of scripts associated with _Arabidopsis halleri_ genome assembly.

## Genome assembly process
Trimming of short reads prior to TrioBinning with Trimmomatic v0.39 (Bolger, 2014): Illumina_adapters.fa:2:30:10; SLIDINGWINDOW:4:15 LEADING:5 TRAILING:5 MINLEN:50 TOPHRED33

Canu assembly:
'genomeSize=250m' 'correctedErrorRate=0.02' 'corOutCoverage=100' 'minReadLength=2000' 'minOverlapLength=500'
'batOptions=-dg 3 -db 3 -dr 1 -ca 500 -cp 50'

## Phasing
TrioBinning was used to separate reads by parent. These reads were aligned to the assembly with minimap2. Coverage files were generated as previously described (Pucker & Brockington, 2018). Coverage per contig was calculated for the D111 and D654 read mappings. Contigs were assigned to one haplophase if the ratio between both coverage was >10.

analyze_cov_per_contig.py

cov_based_phasing.py


## Contig length distribution

contig_lenth_distr.py


## Removal of spurious markers
This scripts removed single markers which are apparently misplaced on the wrong linkage group. Input and output are BED files.

remove_spurious_markers.py

## Marker colinearity check
This script performs a nummeric comparison of genetic and physical marker positions to identify outliers. Groups of X outliers with a minimal genetic size of X are flagged.

nummeric_check.py





## Scripts
ONT_scaffolding.py

in_silico_genetic_map_constructor.py

jcvi_wrapper_genes.py

phasing_check.py
