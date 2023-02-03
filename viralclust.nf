#!/usr/bin/env nextflow

/*
* Source: https://github.com/klamkiew/viralclust/commit/3203e6de334c6834877dbdffecff70df07ed80d7
* "Clustering of viral genomes based on different algorithms and metrices"
* (Author: kevin.lamkiewicz@uni-jena.de)
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
println " "

if (workflow.profile == 'standard' || workflow.profile.contains('local')) {
    println "\033[2mCPUs to use: $params.cores, maximal CPUs to use: $params.max_cores\u001B[0m"
    println " "
}

if ( params.profile ) {
  exit 1, "ERROR: --profile is WRONG use -profile"
}

if (params.mode != 'simulation' && params.mode != 'real') {
    exit 1, "ERROR: choose between simulation or real data mode!"
}
if (params.mode == 'simulation' && params.sort_reads && params.mixture_fastqs=='') {
    exit 1, "ERROR: Please input a folder of mix-in fastq samples!"
}
if (params.mode == 'simulation' && params.sort_reads && params.mixture_lineages=='') {
    exit 1, "ERROR: Please input a csv file with lineage annotation for the input mix-in samples!"
}
if (params.mode == 'real' && params.fastq=='') {
    exit 1, "ERROR: Please input fastq file!"
}

log.info """\
    ACWLA -- Amplicon-based Clustering of Wastewater Sequencing Data for SARS-CoV-2 Lineage Detection and Abundance Estimation
    ==================================
    Analysis mode:          $params.mode
    Input Files:            
    Output path:            $params.output
    CPUs used:              $params.cores
    ${sw -> if (params.hdbscan_params != '') sw << "HDBSCAN parameters:     ${params.hdbscan_params}"}
    Read length filter:     $params.min_read_len nt
    False positive abundance cutoff:  $params.threshold
    """
    .stripIndent()


if (params.sort_reads) { 
  bed = Channel.fromPath(params.bed)
  primer = Channel.fromPath(params.primer)
}
reference = Channel.fromPath(params.reference)
usher = Channel.fromPath(params.usher)


include { fastq_to_fasta } from './modules/fastq_to_fasta'
include { sort_sequences; primer_sort; read_filtering } from './modules/sortsequences'
include { primer_clipping } from './modules/primer_clipping'
include { remove_redundancy; add_redundancy } from './modules/remove_redundancy'
include { create_lineage_dict } from './modules/create_dict'
include { primerset_to_amplicon } from './modules/primerset_to_amplicon'
include { hdbscan } from './modules/hdbscan'
// include { reverseComp } from './modules/reverseComp'
include { call_variants } from './modules/cluster_profiles'
include { update_reference; lineage_prediction; amplicon_quant; aggregate_abundances } from './modules/lineage_prediction'
include { lineage_prediction as lineage_prediction_repeat; amplicon_quant as amplicon_quant_repeat; aggregate_abundances as aggregate_abundances_repeat } from './modules/lineage_prediction'


workflow prepare_mixture {

  main:
    mixed_components = Channel.fromPath("$params.mixture_fastqs/*.fastq")
    mixed_lineages = Channel.fromPath(params.mixture_lineages).splitCsv(header:false,skip:1)

    fastq_to_fasta(mixed_components)
    mixed_components_fa = fastq_to_fasta.out.map{ it -> tuple(it.baseName, it) } // [barcode, fasta]
    create_lineage_dict(mixed_components_fa.join(mixed_lineages))
    lineageDict = create_lineage_dict.out.collectFile(keepHeader: true, name: "lineageDict.csv", storeDir:"$params.output/$params.mode/input")
    fasta = fastq_to_fasta.out.collectFile(name:"query.fasta", storeDir:"$params.output/$params.mode/input")

  emit:
    lineageDict
    fasta
}


workflow get_amplicon_reads {
  take:
    fasta
  
  main:
    primer_sort(fasta.combine(primer))
    primer_reads = primer_sort.out.primer_reads.flatten().map{ it -> tuple(it.name, it) } // [ primername, fasta ]
    read_filtering(primer_reads)
    read_filtering.out.log.collectFile(storeDir: "$params.output/$params.mode/${params.runinfo}", name: "read_filtering.log")
    filtered_primer_reads = read_filtering.out.filtered.filter{ it[1].countFasta() != 0 } 

    primer_clipping(filtered_primer_reads.combine(bed).combine(reference))
    primer_clipped_ch = primer_clipping.out.primerclipped_fasta
    right_primerclipped = primer_clipped_ch.filter{ it[0] ==~ /.*RIGHT.fasta/ }
    left_primerclipped = primer_clipped_ch.filter{ it[0] ==~ /.*LEFT.fasta/ }
    
    sort_sequences(right_primerclipped)
    right_primerclipped_sorted = sort_sequences.out
    primerclipped_redundant = left_primerclipped.concat(right_primerclipped_sorted)

    remove_redundancy(primerclipped_redundant) 
    // [primername, fasta]
    nr_primersets = remove_redundancy.out.nr_fasta.map{ it -> tuple(it[0].split('_')[1], it[1]) }.groupTuple().map{ it.flatten() } // [amplicon, left.fasta, right.fasta]
    amplicon_redundancy = remove_redundancy.out.redundancy.collectFile(storeDir: "$params.output/$params.mode/$params.amplicon_data/redundant"){ it -> ["${it.name}", it] }.map{ it -> tuple(it.baseName.split('-')[0], it) }
    remove_redundancy.out.log.collectFile(storeDir: "$params.output/$params.mode/${params.runinfo}", name: "remove_redundancy.log")

    nr_primersets_handle_borders = nr_primersets.map{it -> if (it.size()==2) [it[0], it[1], file(params.dummy)] else it}
    primerset_to_amplicon(nr_primersets_handle_borders)
    amplicon_ch = primerset_to_amplicon.out    
            

  emit:
    amplicon_ch
    amplicon_redundancy
}

workflow clustering {
  take:
    amplicon_fasta
    lineageDict
    amplicon_redundancy
    amplicon_size
    //primerDict

  main:
    hdbscan(amplicon_fasta.combine(lineageDict)) // previously had primerDict as additional input
    cluster_result = hdbscan.out.amplicon_cluster
    hdbscan.out.log.collectFile(name: 'hdbscan.log', storeDir: "${params.output}/${params.mode}/${params.runinfo}")
  
    // [amplicon, ampliconsize, [cluster0.fasta,...clusterN.fasta], amplicon-duplicates.fasta]
    cluster_result_extended = amplicon_size.join(cluster_result).join(amplicon_redundancy)
    add_redundancy(cluster_result_extended.transpose())
    redundant_cluster_result = add_redundancy.out
    //results_channel = Channel.value('HDBSCAN').combine(hdbscan.out.hdbscan_out)
    // why is this still needed?
    //reverseComp(results_channel)

  emit:
    //results_channel
    redundant_cluster_result
}

// workflow evaluation {
//   take:
//     cluster_fasta_files
//     primerDict

//   main:
//     split_fasta_by_primer(cluster_fasta_files.combine(primerDict))
//     split_fasta_by_primer.out.view()
//     split_fasta_by_primer.out.combine(reference).view()
//     call_variants(split_fasta_by_primer.out.combine(reference))
//     call_variants.out.bam_path.groupTuple().view()
//     call_variants.out.bam_path.groupTuple().map{ it -> it[1] }.first().flatten().unique().view()
//     mutation_heatmap(call_variants.out.bam_path.groupTuple().map{ it -> it[1] }.first().flatten().unique())
// }

workflow prediction{
  take:
    cluster_fasta // [amplicon, ampliconsize, cluster, cluster.fasta]
    amplicon_scaler // [amplicon, ampliconsize/samplesize]

  main:
    // First, split cluster by primer?
    call_variants(cluster_fasta.combine(reference))
    // [amplicon, ampliconsize, cluster, cluster.tsv]
    allele_frequency_ch = call_variants.out.cluster_AF.filter{ it[3].baseName.contains("medaka") }
    call_variants.out.log.collectFile(storeDir: "$params.output/${params.mode}/${params.runinfo}", name: "variant_call.log")

    prediction_input = cluster_fasta.join(allele_frequency_ch, by:[0,1,2]).map{ it -> tuple(it[0], it[2], it[3].countFasta()/it[1], it[4]) }
    // [amplicon, cluster, clustersize/ampliconsize, cluster.tsv]

    update_reference(usher)
    new_usher = update_reference.out

    lineage_prediction(prediction_input.combine(new_usher))
    lineage_prediction.out.log.collectFile(name:'quantify_cluster_abundances.log', storeDir:"$params.output/${params.mode}/${params.runinfo}")
    amplicon_cluster_abundances = lineage_prediction.out.abs // [amplicon, cluster, tsv]
    amplicon_ready = amplicon_cluster_abundances.map{ it -> tuple(it[0], it[2]) }.groupTuple().map{ it ->  it[0]} // [amplicon] 

    amplicon_quant(amplicon_ready.join(amplicon_scaler)) // input [amplicon, ampliconsize/samplesize]
    amplicon_abundances = amplicon_quant.out.abs // [amplicon, tsv]
    amplicon_quant.out.log.collectFile(name:'quantify_amplicon_abundances.log', storeDir:"$params.output/${params.mode}/${params.runinfo}")

    sample_abundances = amplicon_abundances.map{ it -> it[1] }.collectFile(keepHeader: true, skip: 1, name:"scaled_amplicon_abundances.tsv") // [tsv]

    aggregate_abundances(sample_abundances, new_usher, Channel.value(params.threshold))
    fp_update = aggregate_abundances.out.csv
    if ( fp_update.splitCsv(header: false, skip: 1).map{ row -> row[-1].toInteger() }.sum() == 0 ) {
      final_output = aggregate_abundances.out.abs.collect() // [tsv]    
    }
    else{
      aggregate_abundances.out.fp_log.collectFile(name:"false_positive_detection.log", storeDir:"$params.output/${params.mode}/${params.runinfo}")
      lineage_prediction_repeat(prediction_input.combine(fp_update))
      lineage_prediction_repeat.out.log.collectFile(name:'quantify_cluster_abundances.log', storeDir:"$params.output/${params.mode}/${params.runinfo}")
      amplicon_cluster_abundances_repeat = lineage_prediction_repeat.out.abs // [amplicon, cluster, tsv]
      amplicon_ready_repeat = amplicon_cluster_abundances_repeat.map{ it -> tuple(it[0], it[2]) }.groupTuple().map{ it ->  it[0]} // [amplicon] 

      amplicon_quant_repeat(amplicon_ready_repeat.join(amplicon_scaler)) // input [amplicon, ampliconsize/samplesize]
      amplicon_abundances_repeat = amplicon_quant_repeat.out.abs // [amplicon, tsv]
      amplicon_quant_repeat.out.log.collectFile(name:'quantify_amplicon_abundances.log', storeDir:"$params.output/${params.mode}/${params.runinfo}")

      sample_abundances_repeat = amplicon_abundances_repeat.map{ it -> it[1] }.collectFile(keepHeader: true, skip: 1, name:"scaled_amplicon_abundances.tsv") // [tsv]

      aggregate_abundances_repeat(sample_abundances_repeat, fp_update, Channel.value(0))
      final_output = aggregate_abundances_repeat.out.abs.collect() // [tsv]
      // final_log_repeat = aggregate_abundances_repeat.out.log.collectFile(name:'sample_abundances.log', storeDir:"$params.output/${params.runinfo}")
    }

  emit:
    final_output

}

workflow {
  // in simulation mode either prepare a mixture or load an already prepared one
  if (params.mode == 'simulation') {
    if (params.prepare_mixture) {
      prepare_mixture()
      lineageDict = prepare_mixture.out.lineageDict
      fasta = prepare_mixture.out.fasta
    }
    else {
      lineageDict = Channel.fromPath("$params.output/$params.mode/input/lineageDict.csv")
      fasta = Channel.fromPath("$params.output/$params.mode/input/query.fasta")
    }
  }
  // in real data mode, convert fastq file to fasta
  else {
    input_ch = Channel.fromPath(params.fastq)
    fastq_to_fasta(input_ch)
    fasta = fastq_to_fasta.out
    lineageDict = Channel.fromPath(params.dummy)
  }
  lineageDict.view()
  fasta.view()

  // sort reads into amplicons or load already sorted amplicon FASTAs and redundancy data
  // from either the simulation or real data input folder
  if (params.sort_reads) {
    get_amplicon_reads(fasta)
    amplicon_fasta = get_amplicon_reads.out.amplicon_ch 
    amplicon_redundancy = get_amplicon_reads.out.amplicon_redundancy
  }
  else {
    amplicon_ch = Channel.fromPath("$params.output/$params.mode/$params.amplicon_data/final/*.fasta")
    amplicon_redundancy_ch = Channel.fromPath("$params.output/$params.mode/$params.amplicon_data/redundant/*.fasta")
    // [amplicon, amplicon.fasta]
    amplicon_fasta = amplicon_ch.map{ it -> tuple(it.baseName, it) }
    // [amplicon, amplicon-duplicates.fasta]
    amplicon_redundancy = amplicon_redundancy_ch.map{ it -> tuple(it.baseName.split('-')[0], it) }
  }
  amplicon_fasta.view()
  amplicon_size = amplicon_fasta.join(amplicon_redundancy).map{ it -> tuple(it[0], it[1].countFasta()+it[2].countFasta())} // [amplicon, size]
  amplicon_size.view()
  sample_size = amplicon_size.map{ it -> it[1] }.sum() // [size]
  sample_size.view()
  // [amplicon, ampliconsize/samplesize]
  amplicon_scaler = amplicon_size.combine(sample_size).map{ it -> tuple(it[0], it[1]/it[2]) } 

  clustering(amplicon_fasta, lineageDict, amplicon_redundancy, amplicon_size) // previously: primerdict included
  cluster_fasta = clustering.out.redundant_cluster_result // [amplicon, ampliconsize, clusterX, clusterX.fasta]
  //evaluation(clustering.out.cluster_result, get_amplicon_reads.out.primerDict)
      
  prediction(cluster_fasta, amplicon_scaler)
  abundances = prediction.out.final_output
  abundances.view()

}  

workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}  

//  TODO
def helpMSG() {
    c_reset = "\033[0m";
    c_red = "\033[1;31m"
    c_green = "\033[1;32m";
    c_yellow = "\033[0;33m";
    c_blue = "\033[0;34m";
    c_dim = "\033[2m";
    log.info """
    ____________________________________________________________________________________________

    ${c_green}Welcome to ACWLA - your pipeline to estimate SARS-CoV-2 lineage abundances from wastewater sequencing data !${c_reset}
    ____________________________________________________________________________________________

    ${c_yellow}Usage example:${c_reset}
    If you need to analyse a fastq sample from wastewater sequencing:
    nextflow run main.nf --fastq FASTQ --bed BED --primer TSV 

    If you need to create an equally mixed sample from a number of fastq files. Assuming that
    every fastq sample contains sequencing reads from only one lineage.
    nextflow run main.nf --mode simulation --prepare_mixture --mixture_fastqs PATH/TO/FASTQS --mixture_lineages CSV

    ____________________________________________________________________________________________

    ${c_yellow}Mandatory Input:${c_reset}
    ${c_green}--fastq FILEPATH${c_reset}                  Path to fastq file containing wastewater seq reads
    ${c_green}--mixture_fastqs PATH${c_reset}             Path to folder containing single-lineage fastq files
    ${c_green}--mixture_lineages FILEPATH${c_reset}       Path to csv file mapping fastq file names to their lineage annotation
                                                          (csv is comma-separated)
    ${c_green}--primer${c_reset}
    ${c_green}--bed${c_reset}

    ${c_yellow}Other Input:${c_reset}
    ${c_green}--reference${c_reset}
    ${c_green}--usher${c_reset}

    ${c_yellow}Pipeline Options:${c_reset}
    ${c_green}--mode${c_reset}                            TODO: After clustering, calculate basic statistics of clustering results. For each
                                      tool, the minimum, maximum, average and median cluster sizes are calculated,
                                      as well as the average distance of two representative genomes.
    ${c_green}--sort_reads${c_reset}
    ${c_green}--prepare_mixture${c_reset}
    ${c_green}--threshold${c_reset}
    ${c_green}--min_read_len${c_reset}
    ${c_green}--hdbscan_params${c_reset}                  Additional parameters for UMAP and HDBscan cluster analysis. [default $params.hdbscan_params]
                                      For more information and options, please use
                                      ${c_green}nextflow run viralclust.nf --hdbscan_help${c_reset}.
    
    ${c_yellow}Other Options:${c_reset}
    ${c_green}--cores INT${c_reset}                       max cores per process for local use [default $params.cores]
    ${c_green}--max_cores INT${c_reset}                   max cores used on the machine for local use [default $params.max_cores]
    ${c_green}--memory INT${c_reset}                      max memory in GB for local use [default $params.memory]
    ${c_green}--output PATH${c_reset}                     name of the result folder [default $params.output]
                                          
    ${c_yellow}Nextflow options:${c_reset}
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
  These will be UMAPped based on their kmer frequency representation and then clustered 
  into groups (clades) based on their sequence similarity. For each clade, the centroid sequence is
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
    hdbscan_virus.py [options] <inputSequences> <lineageDict> <primerDict>

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
