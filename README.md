# LT_senterica

## Description
This pipeline uses an alignment-free similarity search algorithm (MinHash) to identify close relatives of Salmonella query sequences and subsequently construct a localized phylogenetic tree based on SNP distances. The process involves comparing the query sequences to Salmonella whole genome sequences from public sequence repositories (NCBI and BIGSdb).

## Installation

## Usage

__Required arguments__

```-i --input``` Input directory containing raw sequence .fastq files. Each sample must be organized into a separate directory under the main input directory. The name of each sample directory will be used as sample names.

Example: ```LT_main.sh -i input_sequences/```

```
input_sequences/
├── 16-4563/
│   ├── 16-4563_S1_L001_R1_001.fastq
│   └── 16-4563_S1_L001_R2_001.fastq
└── 16-4648/
    ├── 16-4648_S2_L001_R1_001.fastq
    └── 16-4648_S2_L001_R2_001.fastq

2 directories, 4 files, 16-4563 and 16-4648 are two different samples
```

```-o --output```       Output directory that will contain the tree file in .nwk format

```-s --sketch```       Mash sketch of GenomeTrakr and PubMLST Salmonella sequences

```-r --reference```    Reference Sequence for SNP tree

__Optional arguments__

```-l --list```         List of query genome paths

```-n --neighbours```   Number of neighbours to include per query [Default: 50]

```--population```      Construct the tree in the context of the entire Salmonella population structure

```-t --thread```       Number of threads used [Default: 1]

```-h --help```         Display help message