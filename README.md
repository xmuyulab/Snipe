# Snipe
## Snipe:Identification of food-borne pathogens with low abundance
### Introduction
Snipe (SeNsItive Pathogen dEtection) wraps around Pathoscope 2.0 with species-specific genomic regions (SSRs) based abundance rectification to improve its sensitivity and specificity in detecting pathogens at low abundance. 
Snipe consists of three modules:
- The snipeMap module will map unassembled metagenomic reads against a target library and remove sequences that align to the filter and host libraries. 
- The snipeId module will reassign ambiguous reads, identify microbial strains present in the sample, and estimate proportions of reads from each genome. 
- The snipeRec module will align raw reads to SSRs and generate reports containing read proportions to each genome after rectification by the a posteriori probabilities.

### Support and Contact
For any issues or concerns, please contact us at binhong@stu.xmu.edu.cn
 

### Pathogenic Species Supported
Species name | Number of SSRs
-|-
Escherichia coli | 98
Salmonella enterica | 387
Staphylococcus aureus | 169
Listeria monocytogenes | 157
Campylobacter jejuni | 132
Vibrio cholerae | 261
Vibrio parahaemolyticus | 1206
Proteus mirabilis | 1122
Yersinia enterocolitica | 141
Clostridium perfringens | 2377

### Install Software Dependencies
It is recommended to create a new conda environment:
```
conda create -n python37 python=3.7

# Activate this environment:
conda activate python37
```
       • numpy (v1.15.0)
            conda install -c conda-forge numpy
       • pandas (v0.24.2)
            conda install -c conda-forge pandas
       • bowtie2 (v2.2.5)
            conda install -c bioconda bowtie2
       • pysam (v0.15.3)
            conda install -c bioconda pysam 
            

### Manual
First of all, we should:
- change directory (cd) to snipe folder
- cd into snipe directory and call snipeIndex module help for details
  ```
  cd ../snipe
  python snipe.py -h
  ```
#### MAP
We need the database of strains, which can be downloaded from NCBI. First you need to make sure that the index has been established otherwise the software will take a moment to build the index.

call snipeMap module help for details
```
python snipe.py MAP -h

python snipe.py MAP -1 read_1 -2 read_2 -targetRefFiles path1 -filterRefFiles path2 -outDir path3 -
outAlign name1 -expTag name2 -indexDir path4 -numThreads 8

Required arguments:

-1,              string                    paired-end and requires an even number of FASTQ files represented as pairs

-2,              string                    paired-end and requires an even number of FASTQ files represented as pairs

-targetRefFiles, string                    Target Reference Genome Fasta Files Full Path (Comma Separated)

-filterRefFiles, string                    Filter Reference Genome Fasta Files Full Path (Comma Separated)

-outDir,         string                    Output Directory (Default=. (current directory))

-outAlign,       string                    Output Alignment File Name (Default=outalign.sam)

-expTag,         string                    Experiment Tag added to files generated for identification

-indexDir,       string                    index directory (default=. (current directory))

Optional arguments:

-numThreads,     int                       Number of threads to use by aligner (bowtie2) if different from default (8)
```
#### ID
First you need to make sure that the map module is finished. ID module will use file .sam generated previously with MAP module.

call snipeId module help for details

```
python snipe.py ID -h

python snipe.py ID -outDir path1 -alignFile name1 -expTag name2

Required arguments:

-outDir,         string                    Output Directory (Default=. (current directory))

-alignFile,      string                    Alignment file path

-expTag,         string                    Experiment tag added to output file for easy identification


```
#### REC
We need the database of the core part of strains. Attention should be paid to the naming of the core part of the strain.
.e.g. (strain named Bacillus_cereus, core part of strain named: only_Bacillus_cereus_blastn ) 

call snipeRec module help for details

```
python snipe.py REC -h

python snipe.py REC -ssrRef path1 -1 read_1 -2 read_2 -idReport file1 -dictTarget file2 -dictTemplate file3 -outDir path2 -numThreads 1

Required arguments:

-ssrRef,         string                    the directory of the species specific region

-1,              string                    Input Read Fastq File (Pair 1)

-2,              string                    Input Read Fastq File (Pair 2)

-idReport,       string                    Report file generated by ID module

-dictTarget,     string                    Target strain information

-dictTemplate,   string                    the dict which contain the strain to name mapping

-outDir,         string                    Output Directory (Default=.(current directory))

-expTag,         string                    Experiment tag added to output file for easy identification

Optional arguments:

-numThreads,     int                       Number of threads to use default (1)

```

### Output TSV file format
Columns in the TSV file:
#### 1.Genomes:
This is the name of the genome found in the alignment file.
#### 2.Accession ID:
Gene number used by NCBI database.
#### 3.Rectified Final Guess:
This represent the percentage of reads that are mapped to the genome in Column 1 after snipe using SSRs rectify.
#### 4.Final Guess:
This represent the percentage of reads that are mapped to the genome in Column 1 (reads aligning to multiple genomes are assigned proportionally) after snipe reassignment is performed.
#### 5.Rectified Probability:
This represent rectification coefficient.
#### 6.SSR Aligned Reads:
This represent the number of reads that are mapped to the SSRs in Column 1.
#### 7.Rectified Abundance:
This represent the abundance after using SSRs rectifiy.
#### 8.Initial Abundance:
This represent the abundance before using SSRs rectifiy.
#### 9.Final Best Hit:
This represent the percentage of reads that are mapped to the genome in Column 1 after assigning each read uniquely to the genome with the highest score and after pathoscope reassignment is performed.
#### 10.Final Best Hit Read Numbers:
This represent the number of best hit reads that are mapped to the genome in Column 1 (may include a fraction when a read is aligned to multiple top hit genomes with the same highest score) and after pathoscope reassignment is performed.


### Step-by-step example
#### 0. [Make sure you have all the ingredients](example/Ingredients.md)
#### 1. [The SnipeMap module](example/SnipeMap.md)
#### 2. [The SnipeID module](example/SnipeID.md)
#### 3. [The SnipeRec module](example/SnipeRec.md)
