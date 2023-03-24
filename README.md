## <samp>ABACADABRA</samp> - Amplicon-Based hieRArchical Clustering Approach for for lineage Detection and ABundance estimation fRom wAstewater sequencing
![NextFlow](https://img.shields.io/badge/Nextflow-22.10.1.5828-blue.svg)![conda](https://img.shields.io/badge/Uses-conda-green.svg)
***


### DISCLAIMER
This pipeline is based on the Nextflow pipeline [ViralClust](https://github.com/klamkiew/viralclust) from  [Kevin Lamkiewicz](https://github.com/klamkiew).

Specific commit version that was reused: https://github.com/klamkiew/viralclust/commit/3203e6de334c6834877dbdffecff70df07ed80d7
***


### Overview: What is this about?
<samp>ABACADABRA</samp>
- Detect known and unknown SARS-CoV-2 lineages from amplicon-based sequencing data
- Application to wastewater sequencing data or clinical pooling samples.
***


### Installation

* Install conda 22.9.0 and mamba 1.1.0
* Next, make sure that conda is part of your <samp>$PATH</samp> variable, which is usually the case. For any issues, please refer to the respective [installation guide](https://conda.io/projects/conda/en/latest/user-guide/install/index.html).

<sub>**Warning**: Currently, <samp>ABACADABRA</samp> is running on Linux systems. Windows and MacOS have never been tested so far.</sub>
* Create a conda environment and install NextFlow within this environment:

  <details><summary>Click here to see how:</summary>

  ```bash
  conda env create -f environment.yml
  ```
   </details>


* Done!

***


### Quickstart
* `conda activate nextflow`
* `nextflow run main.nf --help`


### Parameters & Options

This paragraph is simply the help message of <samp>ABACADABRA</samp>.

<details><summary>Expand here</summary>

```
Welcome to ABRACADABRA - Amplicon-Based hieRArchical Clustering Approach for for lineage Detection and ABundance estimation fRom wAstewater sequencing

______________________________________________________________________________________________________________________________________________________

Usage example:
If you want to analyse a fastq sample from wastewater sequencing:

  nextflow run main.nf --fastq FASTQ --bed BED --primer TSV

If you already sorted your data into amplicon sets, run:
  nextflow run main.nf --sort_reads false
(Make sure, that amplicon fasta files are located at params.output/params.amplicon_data/final)

If you need to create a mixed sample from a number of fastq files: 
  nextflow run main.nf --mode simulation --prepare_mixture --mixture_fastqs CSV --mixture_lineages CSV --bed BED --primer TSV 
(Assuming that every fastq sample contains sequencing reads from only one lineage)

If you want to analyse your recently mixed sample where the reads were already sorted into amplicons, run the following command.
  nextflow run main.nf --mode simulation --sort_reads false
(Make sure that the amplicon fasta files are located in params.output//params.amplicon_data/final and that the input mixed fasta sample
and lineageDict.csv are located at params.output/input)
______________________________________________________________________________________________________________________________________________________

Mandatory Input:
--fastq FILEPATH                  Path to FASTQ file
--mixture_fastqs PATH             Path to comma-separated CSV file containing the filepaths to FASTQ files with reads to spike-in and the target abundances.
--mixture_lineages FILEPATH       Path to comma-separated CSV file mapping FASTQ file names to their lineage annotation (ONE lineage per FASTQ)
--primer                          Primer sequence file: TSV with two columns -> primer name and primer sequence
--bed                             Primer scheme (BED)

Other Input:
--reference                       SC2 index genome (FASTA)   
--usher                           Mutational barcode file (CSV)

Pipeline Options:
--mode                            In simulation mode, you can mix a sample from input fastq files at defined abundances and 
                                  annotate sampled reads with lineage labels for evaluation or load already prepared mixture samples with lineageDict.csv
                                  (STRING)
                                                      Use the "real" mode for samples where you do not have lineage annotations on the read level
                                                      [default real]
--sort_reads                      If true, input read data is sorted into amplicon fasta files. Else, the pipeline tries to load already sorted 
                                  data from params.output/params.amplicon_data/final (BOOLEAN) [default true]
--prepare_mixture                 Mix input fastq samples at predefined abundances (BOOLEAN) [default false]
--sample_size                     Target sample size when simulating a mixture. Number of reads per component is calculated 
                                  from the defiend abudnance and the sample size (INT) [default 200000]
--unknown                         If true, the pipeline uses called SNPs that don't overlap with the barcode reference to detect an "unknown" lineage (BOOLEAN) [default true]
--variant_caller                  Currently, only ivar is supported (STRING)
--threshold                       Minimal abundance cutoff. If there are lineages predicted with abundances below the threshold, lineage detection
                                  and abundance estimation is repeated without the respective lienages in the reference (FLOAT) [default 0.01]
--min_read_len                    Minimum read length allowed when filtering sorted amplicon read sets (INT) [default 100]
--max_primer_mismatch             Maximum number of nt mismatches allowed between a primer and a read durint amplicon sorting (INT) [default 3]
--umap_params                     Additional parameters for UMAP [default -k 7 --umap_metric cosine --n_components 2 --n_neighbors 15 --min_dist 0.1]
--hdbscan_params                  Additional parameters for HDBSCAN cluster analysis [default --hdbscan_metric cosine --min_samples 5 --min_cluster_size 5 --cluster_selection_epsilon 0.05]
                                  For more information and options for hdbscan and umap, please use
                                  nextflow run main.nf --ClusterHelp

Other Options:
--cores INT                       max cores per process for local use [default 1]
--max_cores INT                   max cores used on the machine for local use [default 96]
--memory INT                      max memory in GB for local use [default 16.GB]
--output PATH                     name of the result folder [default output]

Nextflow options:
____________________________________________________________________________________________

```
</details>


***
### Data
* UShER mutational barcode data frame from Freyja version with commit hash 'a62c7f7e4011eb1de0b3aa67c53b31396805b496'  
* Spike covering amplicons of the ARTIC V4 panel in example/

### Cluster Tools
<details><summary>Click here for all citations</summary>

  * UMAP:
    * `McInnes, L, Healy, J, "UMAP: Uniform Manifold Approximation and Projection for Dimension Reduction", ArXiv e-prints 1802.03426, 2018`
  * HDBscan: 
    * `Campello, Moulavi, Sander, "Density-Based Clustering Based on Hierarchical Density Estimates" In: Density-Based Clustering Based on Hierarchical Density Estimates, 2013`
</details>

***

### Graphical Workflow

![Workflow graph](/pic/workflow.png)
