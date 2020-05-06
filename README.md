# Snipe
## Snipe:Identification of food-borne pathogens with low abundance
### Introduction
Snipe (SeNsItive Pathogen dEtection) wraps around the metagenomic strain-typing tool Pathoscope 2.0 with Strain Specific Region (SSR) based abundance rectification to improve its sensitivity and specificity in detecting pathogens at extremely low abundance.  It consists of three modules for sequencing-based metagenomic profiling. The snipeMap module the NGS reads aligns to the target and filter reference library which we build before. The snipeId module, using an Expectation-Maximization (EM) algorithm to search for the best parameter estimates for the proportion of mapped reads. The snipeRec module, by FDR control, the metagenomic reads are first aligned to the species-specific regions (SSRs) using Bowtie2 with default settings. The sum of the numbers of reads aligned to SSR with editing distance less than or equal to 3 are counted as total number of SSR reads, the identification result will be presented on both species level and strain level. Report files: 1) a summary report (*_abundance.tsv) that contains the abundance of each strain after correction and before correction. 2) corrector documented our correction factors for different species  3) result folder keep the results of Pathoscope2. 

### Pathogenic Species Currently Supported
Species name | The Number of genomes
-|-
Escherichia coli | 1028
Salmonella enterica | 836
Staphylococcus aureus | 534
Listeria monocytogenes | 201
Campylobacter jejuni | 174
Vibrio cholerae | 63
Vibrio parahaemolyticus | 55
Proteus mirabilis | 29
Yersinia enterocolitica | 18
Clostridium perfringens | 13

### Install Software Dependencies
       • bowtie2
            conda install -c bioconda bowtie2
       • pysam
            conda install -c bioconda pysam 
            

### Manual
First of all, we should:
- change directory (cd) to snipe folder
- cd into snipe directory and call snipeIndex module help for details
  ```
  cd ../snipe
  python snipe.py -h
  ```
#### map
We need the database of strains, which can be downloaded in NCBI.First you need to make sure that the index has been established otherwise the software will take a moment to build the index.

call snipeMap module help for details
```
python snipe.py map -h

python snipe.py map -1 read_1 -2 read_2 -targetRefFiles path1 -filterRefFiles path2 -o path3 -
outAlign name1 -tag name2 -indexDir path4 -t 1

Required arguments:

-1,              string                    paired-end and requires an even number of FASTQ files represented as pairs

-2,              string                    paired-end and requires an even number of FASTQ files represented as pairs

-targetRefFiles, string                    the folder which the target reference genomes index is located

-filterRefFiles, string                    the folder which the filter reference genomes index is located

-o,              string                    the folder which the output is written to

-outAlign,       string                    output alignment file name (default=outalign.sam)

-tag,            string                    experiment tag added to files generated for identification

-indexDir,       string                    index directory (default=. (current directory))

Optional arguments:

-t,              int                       number of threads to use (default: 1)
```
#### id
First you need to make sure that the map module is finished. Id module will use file .sam generated previously with map module.

call snipeId module help for details

```
python snipe.py id -h

python snipe.py id -o path1 -outAlign name1 -tag name2

Required arguments:

-o,              string                    the folder which the output is written to

-outAlign,       string                    output alignment file name (default=outalign.sam)

-tag,            string                    experiment tag added to files generated for easy identification

```
#### rec
We need  the database of the core part of strains.Attention should be paid to the naming of the core part of the strain.
.e.g.(strain named Bacillus_cereus ,core part of strain named: only_Bacillus_cereus_blastn ) 

call snipeId module help for details

```
python snipe.py rec -h

python snipe.py rec -c path1 -1 read_1 -2 read_2 -idReport file1 -dictTarget file2 -dictTemplate file3 -o path2 -t 1

Required arguments:

-c,              string                    the directory of the species specific region

-1,              string                    paired-end and requires an even number of FASTQ files represented as pairs

-2,              string                    paired-end and requires an even number of FASTQ files represented as pairs

-idReport,       string                    the directory of Report file generated by id module

-dictTarget,     string                    the dict which contain the strain to species mapping

-dictTemplate,   string                    the dict which contain the strain to name mapping

-o,              string                    the folder which the output is written to

Optional arguments:

-t,              int                       number of threads to use (default: 1)

```
