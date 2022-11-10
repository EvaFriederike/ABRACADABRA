#!/usr/bin/env nextflow

/*
* Source: 
*
* Clustering of viral genomes based on different algorithms and metrices
* Author: kevin.lamkiewicz@uni-jena.de
*/

nextflow.enable.dsl=2

if (params.help) {
  exit 0, helpMSG()
}

if (params.hdbscan_help) {
  exit 0, hdbscanHelp()
}

println " "
println "\u001B[32mProfile: $workflow.profile\033[0m"
println " "
println "\033[2mCurrent User: $workflow.userName"
println "Nextflow-version: $nextflow.version"
println "Starting time: $workflow.start"
println "Workdir location:"
println "  $workflow.workDir"
println "Launchdir location:"
println "  $workflow.launchDir"
println "Configuration files:"
println "  $workflow.configFiles"
println "Path for database:"
println " $params.permanentCacheDir\u001B[0m"
println " "

if (workflow.profile == 'standard' || workflow.profile.contains('local')) {
    println "\033[2mCPUs to use: $params.cores, maximal CPUs to use: $params.max_cores\u001B[0m"
    println " "
}

if ( params.profile ) {
  exit 1, "ERROR: --profile is WRONG use -profile"
}

if ( params.fasta == '') {
  exit 1, "ERROR: --fasta is a required parameter.\nMake sure it is set."
}

sortSequence = Channel.fromPath( workflow.projectDir + '/bin/sort_sequences.py', checkIfExists: true )
rc_script = Channel.fromPath( workflow.projectDir + '/bin/reverse_complement.py', checkIfExists: true )
utils = Channel.fromPath( workflow.projectDir + '/bin/utils.py', checkIfExists: true )
umap_hdbscan_script = Channel.fromPath( workflow.projectDir + '/bin/hdbscan_virus.py', checkIfExists: true )


log.info """\
    VIRALCLUST -- CLUSTER YOUR VIRUSES
    ==================================
    Input File:             $params.fasta
    Output path:            $params.output
    CPUs used:              $params.cores
    ${sw -> if (params.hdbscan_params != '') sw << "HDBSCAN parameters:     ${params.hdbscan_params}"}
    """
    .stripIndent()

if (params.fasta) {
  sequences = Channel.fromPath(params.fasta)
  reference_channel = Channel.fromPath(params.reference)
  if (params.lineageDict) {
    lineageDict = Channel.fromPath(params.lineageDict)
  }

  if (params.goi) {
    goi = Channel.fromPath(params.goi)
    include { sort_sequences as goi_sorter } from './modules/sortsequences'
    include { concat_goi } from './modules/remove_redundancy'
  }

  include { sort_sequences } from './modules/sortsequences'
  include { remove_redundancy } from './modules/remove_redundancy'
  include { hdbscan } from './modules/hdbscan'
  include { reverseComp } from './modules/reverseComp'
  include { evaluate_cluster; merge_evaluation } from './modules/evaluate'
  include { call_variants; mutation_heatmap } from './modules/cluster_profiles'
}


workflow preprocessing {
  main:
    sort_sequences(sequences)
    if (params.goi) {
      goi_sorter(goi)
      goiSorted = goi_sorter.out.sort_result
    } else {
      goiSorted = 'NO FILE'
    }
    if (params.remove_redundancy){
      remove_redundancy(sort_sequences.out.sort_result)
      non_redundant_ch = remove_redundancy.out.nr_result
    }
    else {
      non_redundant_ch = sort_sequences.out.sort_result
    }  
    if (params.goi) {
        concat_goi(remove_redundancy.out.nr_result, goiSorted)
        non_redundant_ch = concat_goi.out.nr_result
    }
  emit:
    non_redundant_ch
    goiSorted
}

workflow hdbscan_wf{
  take:
    fasta
    lineageDict
    hdbscan_params
    goiSorted
    reference_channel

  main:
    hdbscan(fasta, lineageDict, hdbscan_params, goiSorted)
    hdbscan_results = hdbscan.out.hdbscan_out
    hdbscan_fasta = hdbscan.out.cluster_fasta_files.flatten()
    call_variants(hdbscan_fasta.combine(reference_channel))
    mutation_heatmap(call_variants.out.bam_path.unique())

  emit:
    hdbscan_results
}

workflow clustering {

  take:
    non_redundant_ch
    goiSorted

  main:
    hdbscan_wf(non_redundant_ch, lineageDict, params.hdbscan_params, goiSorted, reference_channel)
    
    results_channel = Channel.value('HDBSCAN').combine(hdbscan_wf.out.hdbscan_results)

    reverseComp(results_channel)

  emit:
    results_channel
}

workflow {

  if (params.fasta) {
    preprocessing()
    clustering(preprocessing.out.non_redundant_ch, preprocessing.out.goiSorted)
  }  
  
}

def helpMSG() {
    c_reset = "\033[0m";
    c_red = "\033[1;31m"
    c_green = "\033[1;32m";
    c_yellow = "\033[0;33m";
    c_blue = "\033[0;34m";
    c_dim = "\033[2m";
    log.info """
    ____________________________________________________________________________________________

    ${c_green}Welcome to ViralClust - your pipeline to cluster viral genome sequences once and for all!${c_reset}
    ____________________________________________________________________________________________

    ${c_yellow}Usage example:${c_reset}
    nextflow run viralclust.nf --fasta "genomes.fasta"

    ____________________________________________________________________________________________

    ${c_yellow}Mandatory Input:${c_reset}
    ${c_green}--fasta PATH${c_reset}                      Path to a multiple fasta sequence file, storing all genomes that shall be clustered.
                                      Usually, this parameter has to be set, unless the parameter ${c_green}--ncbi_update${c_reset} has been set.
    ____________________________________________________________________________________________

    ${c_yellow}Options:${c_reset}
    ${c_green}--eval${c_reset}                            TODO: After clustering, calculate basic statistics of clustering results. For each
                                      tool, the minimum, maximum, average and median cluster sizes are calculated,
                                      as well as the average distance of two representative genomes.
    ${c_yellow}Cluster options:${c_reset}
    ${c_green}--hdbscan_params${c_reset}                  Additional parameters for UMAP and HDBscan cluster analysis. [default $params.hdbscan_params]
                                      For more information and options, please use
                                      ${c_green}nextflow run viralclust.nf --hdbscan_help${c_reset}.
    ${c_yellow}Computing options:${c_reset}
    ${c_green}--cores INT${c_reset}                       max cores per process for local use [default $params.cores]
    ${c_green}--max_cores INT${c_reset}                   max cores used on the machine for local use [default $params.max_cores]
    ${c_green}--memory INT${c_reset}                      max memory in GB for local use [default $params.memory]
    ${c_green}--output PATH${c_reset}                     name of the result folder [default $params.output]
    ${c_green}--permanentCacheDir PATH${c_reset}          location for auto-download data like databases [default $params.permanentCacheDir]
    ${c_green}--condaCacheDir PATH${c_reset}              location for storing the conda environments [default $params.condaCacheDir]
    ${c_green}--workdir PATH${c_reset}                    working directory for all intermediate results [default $params.workdir]

    ${c_yellow}Nextflow options:${c_reset}
    ${c_green}-with-report rep.html${c_reset}             cpu / ram usage (may cause errors)
    ${c_green}-with-dag chart.html${c_reset}              generates a flowchart for the process tree
    ${c_green}-with-timeline time.html${c_reset}          timeline (may cause errors)
    ${c_reset}____________________________________________________________________________________________

    """.stripIndent()
}

def hdbscanHelp() {
  c_reset = "\033[0m";
  c_red = "\033[1;31m"
  c_green = "\033[1;32m";
  c_yellow = "\033[0;33m";
  c_blue = "\033[0;34m";
  c_dim = "\033[2m";
  log.info """
  ____________________________________________________________________________________________
  This python program is part of ViralClust and takes several genome sequences
  from different viruses as an input.
  These will be clustered these sequences into groups (clades) based
  on their sequence similarity. For each clade, the centroid sequence is
  determined as representative genome, i.e. the sequence with the lowest
  distance to all other sequences of this clade.
  ____________________________________________________________________________________________
  Python Dependencies:
    docopt
    BioPython
    colorlog
    numpy
    scipy
    umap-learn
    hdbscan

  Contact:
    kevin.lamkiewicz@uni-jena.de

  Usage:
    hdbscan_virus.py [options] <inputSequences> <lineageDict> [<genomeOfInterest>]

  Options:
    -h, --help                              Show this help message and exits.
    -v, --verbose                           Get some extra information from viralClust during calculation. [Default: False]
    --version                               Prints the version of viralClust and exits.

    -o DIR, --output DIR                    Specifies the output directory of viralClust. [Default: pwd]
    -p PROCESSES, --process PROCESSES       Specify the number of CPU cores that are used. [Default: 1]


    -k KMER, --kmer KMER                    Length of the considered kmer. [Default: 7]
    --metric METRIC                         Distance metric applied by UMAP (if applied) and HDBSCAN.
                                            The following are supported:
                                            'euclidean', 'manhatten', 'chebyshev', 'minkwoski',
                                            'canberra', 'braycurtis',
                                            'mahalanobis', 'wminkowski', 'seuclidean',
                                            'cosine'.
                                            If an invalid metric is set, ViralClust will default back to 
                                            the cosine distance.
                                            [Default: cosine]
    --neighbors NEIGHBORS                   Number of neighbors considered by UMAP to reduce the dimension space.
                                            Low numbers here mean focus on local structures within the data, whereas 
                                            larger numbers may loose fine details. [default: 15]
    --dThreshold dTHRESHOLD                 Sets the threshold for the minimum distance of two points in the low-dimensional space.
                                            Smaller thresholds are recommended for clustering and identifying finer topological structures
                                            in the data. [Default: 0.1]
    --dimension DIMENSION                   UMAP tries to find an embedding for the input data that can be represented by a low-dimensional space.
                                            This parameter tells UMAP how many dimensions should be used for the embedding. Lower numbers may result 
                                            in loss of information, whereas larger numbers will increase the runtime. [Default: 10]

    --clusterSize CLUSTERSIZE               This parameter forces HDBSCAN to form cluster with a size larger-equal to CLUSTERSIZE.
                                            Be aware that some data points (i.e. genomes) or even whole subcluster are considered as noise, if this parameter is set too high.
                                            E.g., if a very distinct viral genus has 40 genomes and the parameter is set to anything >40, HDBSCAN will not form
                                            the genus specific cluster. [Default: 5]
    --minSample MINSAMPLE                   Intuitively, this parameter declares how conservative clustering is performed. Higher values will lead 
                                            to more points considered noise, whereas a low value causes "edge-cases" to be grouped into a cluster.
                                            The default parameter is the same as CLUSTERSIZE. [Default: CLUSTERSIZE]
    --clusterThreshold CLUSTERTHRESHOLD     A distance threshold. Clusters below this value will be merged. [Default: 0.0]
  """.stripIndent()
}