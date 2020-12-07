# LT_senterica

## Description
This pipeline uses a threshold-free clustering approach to rapidly identify close relatives of Salmonella query sequences from Salmonella isolates deposited in public foodborne pathogen sequence repositories (GenomeTrakr and BIGSdb) by combining an alignment-free similarity search algorithm (MinHash) and density-based hiercharchical clustering algorithm (HDBSCAN). A core genome SNP tree is subsequently constructed to visualize the evolutionary relationships between the queries and the identified query neighbours. The search will require users to first download a local copy of a reduced k-mer representation of the public Salmonella whole genome sequence databases.

## Installation
Not yet available

## Usage

__Required arguments__

```-i --input``` Input directory containing raw sequence .fastq files. Each sample must be organized into a separate directory under the main input directory. The name of each sample directory will be used as sample names. [Must specify either -i or -l]

Example: ```LTS.sh molecularlinkage -i input_sequences/```

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

```-l --list```         list of query genomes, each line contains absolute path to genome assemblies. [Must specify either -i or -l]

```-o --output```       Output directory that will contain the tree file in .nwk format.

__Optional arguments__

```-t --thread```       Number of threads used [Default: 1]

```-h --help```         Display help message

## Dependencies

* mash >= 2.1
* R >= 3.6
* abricate >= 1.0.1
* shovill >= 1.1.0
* ncbi-amrfinderplus >= 3.9.3
* phame == 1.0.3