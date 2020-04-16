# LT_senterica

## Description
This pipeline uses an alignment-free similarity search algorithm (MinHash) to identify close relatives of Salmonella query sequences and subsequently construct localized phylogenetic trees based on MLST distances. The process involves comparing the query sequences to Salmonella whole genome sequences from public sequence repositories (NCBI and BIGSdb). For efficiency, the search space is significantly narrowed down by using a fixed set of reference sequences that was developed to be representative of the global Salmonella diversity.

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

```-o --output``` Output directory that will contain the tree file in .nwk format

```-s --sketch``` Path to the Salmonella reference sequences used by Mash in .msh format

```-g --gene ```  Path to the MLST scheme required for allele calling by chewBBACA

```--tmp```       Temporary file directory

__Optional arguments__

```-t --thread``` Number of threads used.
