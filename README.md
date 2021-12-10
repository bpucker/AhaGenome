# _Arabidopsis halleri_ genome sequence
Collection of scripts associated with the _Arabidopsis halleri_ genome sequence assembly.

## Genome assembly process
Trimming of short reads prior to TrioBinning with Trimmomatic v0.39 ([Bolger, 2014](https://doi.org/10.1093/bioinformatics/btu170)): Illumina_adapters.fa:2:30:10; SLIDINGWINDOW:4:15 LEADING:5 TRAILING:5 MINLEN:50 TOPHRED33

Canu assembly:
'genomeSize=250m' 'correctedErrorRate=0.02' 'corOutCoverage=100' 'minReadLength=2000' 'minOverlapLength=500'
'batOptions=-dg 3 -db 3 -dr 1 -ca 500 -cp 50'

## Phasing
TrioBinning was used to separate reads by parent. These reads were aligned to the assembly with minimap2. Coverage files were generated as previously described ([Pucker & Brockington, 2018](https://doi.org/10.1186/s12864-018-5360-z)). Coverage per contig was calculated for the D111 and D654 read mappings. Contigs were assigned to one haplophase if the ratio between both coverage was >10.



```
Usage:
  python analyze_cov_per_contig.py --cov <FILE> --out <DIR>

Mandatory:
  --cov  STR        Coverage file
  --out STR         Directory for temporary and output files.
```


`--cov` specifies the coverage input file ([Pucker & Brockington, 2018](https://doi.org/10.1186/s12864-018-5360-z)).

`--out` specifies the output folder. All result files will be placed in this folder. If this folder does not exist already, it will be created.



```
Usage:
  python cov_based_phasing.py --cov <FILE> --out <DIR>

Mandatory:
  --covinfo1  STR        Coverage file1
  --covinfo2  STR        Coverage file2
  --assembly  STR        Assembly file (FASTA)
  --haplo1    STR        Haplophase1 file (FASTA)
  --haplo2    STR        Haplophase2 file (FASTA)
  --unclass   STR        Unclassified file (FASTA)
```

`--covinfo1` specifies the coverage input file that belongs to haplophase1. This file is based on a mapping of haplophase1 reads. This file is used to calculate the average coverage per contig.

`--covinfo2` specifies the coverage input file that belongs to haplophase2. This file is based on a mapping of haplophase2 reads.. This file is used to calculate the average coverage per contig.

`--assembly` specifies the assembly input file (FASTA). Sequence names in this file need to match the sequence names in the coverage files.

`--haplo1` specifies the haplophase1 output file (FASTA).

`--haplo2` specifies the haplophase2 output file (FASTA).

`--unclass` specifies the output file of all sequences that were not assigned to haplophase1 or haplophase2 (FASTA).





## Contig length distribution

contig_lenth_distr.py


## Removal of spurious markers
This scripts removed single markers which are apparently misplaced on the wrong linkage group. Input and output are BED files. This script also removes genetic markers with an identical genetic position i.e. no recombination between two adjacent markers.

remove_spurious_markers.py

## Marker colinearity check
This script performs a nummeric comparison of genetic and physical marker positions to identify outliers. Groups of X outliers with a minimal genetic size of X are flagged.

nummeric_check.py





## Scripts
ONT_scaffolding.py

in_silico_genetic_map_constructor.py

jcvi_wrapper_genes.py

phasing_check.py




## References

Siadjeu, C.*; Pucker, B.*; Viehöver, P.; Albach, D.C.; Weisshaar, B. High Contiguity de novo Genome Sequence Assembly of Trifoliate Yam (Dioscorea dumetorum) Using Long Read Sequencing. Genes 2020, 11, 274. doi: [10.3390/genes11030274](https://doi.org/10.3390/genes11030274)

Pucker, B., Brockington, S.F. Genome-wide analyses supported by RNA-Seq reveal non-canonical splice sites in plant genomes. BMC Genomics 19, 980 (2018). doi: [10.1186/s12864-018-5360-z](https://doi.org/10.1186/s12864-018-5360-z)

