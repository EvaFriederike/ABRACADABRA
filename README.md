
## <samp>ABACADABRA</samp> - Amplicon-Based Clustering Approach for for lineage Detection and ABundance estimation fRom wAstewater sequencing
![NextFlow](https://img.shields.io/badge/Nextflow-22.10.1.5828-blue.svg)![conda](https://img.shields.io/badge/Uses-conda-green.svg)
***

### DISCLAIMER
This pipeline is based on the Nextflow pipeline [ViralClust](https://github.com/klamkiew/viralclust) from  [Kevin Lamkiewicz](https://github.com/klamkiew).

Specific commit version that I reused: https://github.com/klamkiew/viralclust/commit/3203e6de334c6834877dbdffecff70df07ed80d7
***

### Overview: What is this about?
<samp>ABACADABRA</samp> - Amplicon-Based Clustering Approach for for lineage Detection and ABundance estimation fRom wAstewater sequencing
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

Let us briefly go over the most important parameters and options. There is a more detailed overview of all possible flags, parameters and additional stuff you can
do in the help of message of the pipeline - and at the [end of this file](#help-message).




***
### Data
* UShER mutational barcode data frame from Freyja version with commit hash 'a62c7f7e4011eb1de0b3aa67c53b31396805b496'  
### Cluster Tools



<details><summary>Click here for all citations</summary>

  * UMAP:
    * `McInnes, L, Healy, J, "UMAP: Uniform Manifold Approximation and Projection for Dimension Reduction", ArXiv e-prints 1802.03426, 2018`
  * HDBscan: 
    * `L. McInnes, J. Healy, S. Astels, "hdbscan: Hierarchical density based clustering" In: Journal of Open Source Software, The Open Journal, volume 2, number 11. 2017`
</details>

***

### Graphical Workflow

![Workflow graph](/pic/workflow.png)

### Help Message

This paragraph is simply the help message of <samp>ABACADABRA</samp>.

<details><summary>Expand here</summary>

```
____________________________________________________________________________________________

Welcome to ABACADABRA - Amplicon-Based Clustering Approach for for lineage Detection and ABundance estimation fRom wAstewater sequencing
____________________________________________________________________________________________

Usage example:


____________________________________________________________________________________________

Mandatory Input:

____________________________________________________________________________________________

Options:


Computing options:


Nextflow options:
____________________________________________________________________________________________

```
</details>