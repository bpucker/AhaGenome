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
  python cov_based_phasing.py --covinfo1 <FILE> --covinfo2 <FILE> --assembly <FILE> --haplo1 <FILE> --haplo2 <FILE> --unclass <FILE>

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

```
Usage:
  python contig_lenth_distr.py --in <FILE> --out <FILE>

Mandatory:
  --in   STR        Input file (FASTA)
  --out  STR        Output file (TXT)
```


`--in` specifies the assembly input file (FASTA).

`--out` specifies the statistics output file (TXT).


## Removal of spurious markers
This scripts removed single markers which are apparently misplaced on the wrong linkage group. Input and output are BED files. This script also removes genetic markers with an identical genetic position i.e. no recombination between two adjacent markers. This reduces the computational time required for downstream operations.


```
Usage:
  python3 remove_spurious_markers.py --in <FILE> --out <FILE>

Mandatory:
  --in   STR        Input file (TXT)
  --out  STR        Output file (TXT)
```


`--in` specifies the input file that contains all genetic markers (TXT).

`--out` specifies the output file that will contain only the clean genetic markers (TXT).



## Marker colinearity check
This script performs a nummeric comparison of genetic and physical marker positions to identify outliers. Groups of X outliers with a minimal genetic size of X are flagged.


```
Usage:
  python3 nummeric_check.py --in <FILE> --out <FILE>

Mandatory:
  --in   STR        Input file (TXT)
  --out  STR        Output file (TXT)
```


`--in` specifies the input file that contains all clean genetic markers (TXT).

`--out` specifies the output file that will contain only genetic markers that deviate from the expected position based on the assembly (TXT).




## Long read-based scaffolding
Long reads are converted into a database for BLAST. Contig ends of an assembly are searched against this database. If two contig ends hit the same read, it is possible that they should be connected.

```
Usage:
  python ONT_scaffolding.py --assembly <FILE> --reads <FILE> --out <FOLDER>

Mandatory:
  --assembly   STR        Assembly file (FASTA)
  --reads      STR        Reads file (FASTA)
  --out        STR        Output folder
  
  optional:
  --minlen     INT        Minimal alignment length
  --minsim     INT        Minimal alignment similarity
  --cpus       INT        Number of CPUs to use
```


`--assembly` specifies the assembly input file (FASTA).

`--reads` specifies the file of long reads in FASTA format.

`--out` specifies the output folder.

`--minlen` specifies the minimal alignemnt length of hits to be considered in the BLAST step.

`--minsim` specifies the minimal alignemnt similarity of hits to be considered in the BLAST step.

`--cpus` specifies the number of threads to use for BLAST.

## Construction of ALLMAPS input file
Script for the contstruction of a genetic map (BED) which serves as input for ALLMAPS.

```
Usage:
  python in_silico_genetic_map_constructor.py --assembly <FILE> --ref <FILE> --out <FOLDER>

Mandatory:
  --assembly   STR        Assembly file (FASTA)
  --ref        STR        Reads file (FASTA)
  --out        STR        Output folder
  
  optional:
  --cluster     -         Activate cluster usage
  --k           INT       Marker query length
```


`--assembly` specifies the assembly input file (FASTA).

`--ref` specifies the reference input file (FASTA).

`--out` specifies the output folder.

`--cluster` activates the compute cluster usage instead of running everything locally..

`--k` specifies the length of the query sequence that is used to identify the marker position in the new assembly.

## Analyze individual genes

jcvi_wrapper_genes.py



## Analyze phasing at individual gene level
Checks the phasing of haplophases based on sequence variants and coverage depth in a read mapping.


```
Usage:
  python phasing_check.py  --ref <FILE> --out <FOLDER> --vcf1 <VCF_FILE1> --vcf2 <VCF_FILE2>

Mandatory:
  --ref        STR        Assembly file (FASTA)
  --vcf1       STR        Variant file haplophase1 (VCF)
  --vcf2       STR        Variant file haplophase2 (VCF)
  --out        STR        Output folder
```


`--ref` specifies the assembly input file (FASTA).

`--vcf1` specifies the VCF file of haplophase1 reads mapped against the assembly.

`--vcf2` specifies the VCF file of haplophase2 reads mapped against the assembly.

`--out` specifies the output folder.


## Modify BED file for inversion of chromosomes

python3 invert_linkage_group.py




## References

Siadjeu, C.*; Pucker, B.*; Viehöver, P.; Albach, D.C.; Weisshaar, B. High Contiguity de novo Genome Sequence Assembly of Trifoliate Yam (Dioscorea dumetorum) Using Long Read Sequencing. Genes 2020, 11, 274. doi: [10.3390/genes11030274](https://doi.org/10.3390/genes11030274)

Pucker, B., Brockington, S.F. Genome-wide analyses supported by RNA-Seq reveal non-canonical splice sites in plant genomes. BMC Genomics 19, 980 (2018). doi: [10.1186/s12864-018-5360-z](https://doi.org/10.1186/s12864-018-5360-z)

