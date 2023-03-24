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

if (params.ClusterHelp) {
  exit 0, ClusterHelp()
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

if (params.mode == 'simulation' && params.prepare_mixture){
  if(params.mixture_fastqs=='') {
    exit 1, "ERROR: Please input a folder of mix-in fastq samples!"
  }
  if ( params.mixture_lineages=='') {
    exit 1, "ERROR: Please input a csv file with lineage annotation for the input mix-in samples!"
  }
}

if (params.mode == 'real' && params.fastq=='') {
  if (params.sort_reads){
    exit 1, "ERROR: Please input fastq file!"
  }
}
if (params.sort_reads) {
  if (params.bed == '' || params.primer == '') { exit 1, "ERROR: Please input the required bed and tsv file for the primer panel" }
}

log.info """\
    ABRACADABRA  - Amplicon-Based hieRArchical Clustering Approach for lineage Detection and ABundance estimation fRom wAstewater sequencing
    ==================================
    Analysis mode:      $params.mode
    ${sw -> if (params.mode == 'simulation') sw << "Input mixtures:     ${params.mixture_fastqs}"}
    ${sw -> if (params.mode == 'real') sw << "Input Fastq:      ${params.fastq}"}
    Primer scheme:      $params.primer
    BED file:     $params.bed
    UShER file:     $params.usher
    Read length filter:                  $params.min_read_len nt
    Primer mismatches:      $params.max_primer_mismatch
    UMAP params:      $params.umap_params
    HDBSCAN params:     $params.hdbscan_params
    Minimap2 mode:      $params.minimap_params
    False positive abundance cutoff:     $params.threshold
  
    
    CPUs used:      $params.cores
    Output path:      $params.output
    Runinfo:      $params.runinfo
    """
    .stripIndent()


if (params.sort_reads) { 
  bed = Channel.fromPath(params.bed)
  primer = Channel.fromPath(params.primer)
}
reference = Channel.fromPath(params.reference)
usher = Channel.fromPath(params.usher)
// TODO: implement check for input usher datafream: colsum != 0 for all mutations?!


include { fastq_to_fasta } from './modules/fastq_to_fasta'
include { sort_sequences; primer_sort_edit ; read_filtering; redirect_reads; log_primer_sort } from './modules/sortsequences'
include { primer_clipping } from './modules/primer_clipping'
include { create_lineage_dict; create_mixtures } from './modules/create_sample'
include { primerset_to_amplicon } from './modules/primerset_to_amplicon'
include { hdbscan } from './modules/hdbscan'
include { update_reference; call_variants; filter_variants; filter_variants_repeat; variant_metrics; variant_metrics as variant_metrics_repeat; plot_mutation_profiles; plot_mutation_profiles as plot_mutation_profiles_repeat } from './modules/cluster_profiles'
include { lineage_prediction; amplicon_quant; aggregate_abundances } from './modules/lineage_prediction'
include { lineage_prediction as lineage_prediction_repeat; amplicon_quant as amplicon_quant_repeat; aggregate_abundances as aggregate_abundances_repeat } from './modules/lineage_prediction'


workflow prepare_mixture {

  main:
    mixed_components = Channel.fromPath(params.mixture_fastqs).splitCsv(header:false).map{ it -> tuple(file(it[0]), it[1]) } // [fastq, abundance]
    mixed_lineages = Channel.fromPath(params.mixture_lineages).splitCsv(header:false,skip:1) // [samplename, lineage]

    fastq_to_fasta(mixed_components.map{ it -> it[0] })
    // [samplename, fasta]
    component_fasta = fastq_to_fasta.out.map{ it -> tuple(it.baseName, it)}
    // [samplename, abundance*sample size]
    read_number = mixed_components.map{ it -> tuple(it[0].baseName, it[1]) }.combine(Channel.value(params.sample_size)).map{ it -> tuple(it[0], (it[1].toFloat()*it[2].toFloat()).round()) }

    // input [samplename, fasta, n_reads]
    create_mixtures(component_fasta.join(read_number))
    mixed_components_fa = create_mixtures.out.component // [samplename, sub_fasta]
    create_mixtures.out.log.collectFile(name:"mix_sample.log", storeDir:"$params.output/$params.runinfo/")
    fasta = mixed_components_fa.map{ it -> it[1]}.collectFile(name:"query.fasta", storeDir:"$params.output/input")

    create_lineage_dict(mixed_components_fa.join(mixed_lineages))
    lineageDict = create_lineage_dict.out.collectFile(keepHeader: true, name: "lineageDict.csv", storeDir:"$params.output/input")
    

  emit:
    lineageDict
    fasta
}


workflow get_amplicon_reads {
  take:
    fasta
  
  main:
    redirect_reads(fasta)
    positive_fasta = redirect_reads.out
    amplicon_primers = primer.splitCsv(header: false, skip:1, sep:'\t').map{ it -> tuple(it[0].split('_')[1], it[1]) }.groupTuple().map{ it -> it.flatten()} // [amplicon, primer_L, primer_R]
    primer_sort_edit(positive_fasta.combine(amplicon_primers))
    log1 = primer_sort_edit.out.log
    primer_reads = primer_sort_edit.out.primer_reads
    log_primer_sort(primer_reads.map{ it -> it[1]}.collectFile(name:'sorted.fasta').combine(positive_fasta))
    log2 = log_primer_sort.out
    log1.collect().concat(log2).collectFile(name:"primer_sort.log", storeDir: "${params.output}/${params.runinfo}/")
    // primer_sort(fasta.combine(primer))
    // primer_reads = primer_sort.out.primer_reads.flatten().map{ it -> tuple(it.name, it) } // [ primername, fasta ]
    read_filtering(primer_reads)
    read_filtering.out.log.collectFile(storeDir: "$params.output/${params.runinfo}", name: "read_filtering.log")
    filtered_primer_reads = read_filtering.out.filtered.filter{ it[1].countFasta() != 0 } 

    primer_clipping(filtered_primer_reads.combine(bed).combine(reference))
    // primer_clipped_ch = primer_clipping.out.primerclipped_fasta
    // right_primerclipped = primer_clipped_ch.filter{ it[0] ==~ /.*RIGHT.fasta/ }
    // left_primerclipped = primer_clipped_ch.filter{ it[0] ==~ /.*LEFT.fasta/ }
    
    // sort_sequences(right_primerclipped)
    // right_primerclipped_sorted = sort_sequences.out
    // primerclipped_redundant = left_primerclipped.concat(right_primerclipped_sorted)

    // amplicon_primer_sets = primerclipped_redundant.map{ it -> tuple(it[0].split('_')[1], it[1]) }.groupTuple().map{ it.flatten() } // [amplicon, left.fasta, right.fasta]
    // primersets_handle_borders = amplicon_primer_sets.map{it -> if (it.size()==2) [it[0], it[1], file(params.dummy)] else it}
    // primerset_to_amplicon(primersets_handle_borders)
    // amplicon_ch = primerset_to_amplicon.out    
    amplicon_ch = primer_clipping.out.primerclipped_fasta    

  emit:
    amplicon_ch
}

workflow clustering {
  take:
    amplicon_fasta
    lineageDict

  main:
    hdbscan(amplicon_fasta.combine(lineageDict))
    cluster_result = hdbscan.out.amplicon_cluster
    hdbscan.out.log.collectFile(name: 'hdbscan.log', storeDir: "${params.output}/${params.runinfo}")
  
    // [amplicon, [cluster0.fasta,...clusterN.fasta], amplicon-duplicates.fasta]
    redundant_cluster_result = cluster_result.transpose().map{ it -> tuple(it[0], it[1].baseName, it[1]) }
  
  emit:
    //results_channel
    redundant_cluster_result
}

workflow prediction{
  take:
    cluster_fasta // [amplicon, cluster, cluster.fasta]
    lineageDict

  main:
    // First, split cluster by primer?
    call_variants(cluster_fasta.combine(reference))
    call_variants.out.log.collectFile(storeDir: "$params.output/${params.runinfo}", name: "variant_call.log")
    // [amplicon, cluster, cluster.fasta, cluster.tsv]
    allele_frequency_ch = call_variants.out.cluster_AF.filter{ it[2].countFasta() > 0 }
    update_reference(usher)
    new_usher = update_reference.out
    filter_variants(allele_frequency_ch.combine(new_usher).combine(reference).combine(Channel.value(params.unknown)))
    filter_variants.out.log.collectFile(name:"filter_variants.log", storeDir: "$params.output/${params.runinfo}")
    filtered_variant_ch = filter_variants.out.filtered.filter{ it[2].countFasta() > 0 }
    //sample_variants = allele_frequency_ch.map{ it -> it[3]}.collectFile(name: "sample_variants.tsv", keepHeader:true, skip: 1)
    sample_variants =  filtered_variant_ch.map{ it -> it[3]}.collectFile(name: "sample_variants.tsv", keepHeader:true, skip: 1)
    variant_metrics(sample_variants.combine(new_usher).combine(lineageDict))
    variant_metrics.out.collectFile(name:"variant_evaluation.log", storeDir: "$params.output/${params.runinfo}")
    // [amplicon, ampliconsize]
    amplicon_size = filtered_variant_ch.map{ it -> tuple(it[0], it[2].countFasta()) }.groupTuple().map{ it -> tuple(it[0], it[1].sum()) }
    plot_mutation_profiles(amplicon_size)
    amplicon_size_fwd = plot_mutation_profiles.out.fwd
    amplicon_cluster_map = filtered_variant_ch.map{ it -> tuple(it[0], it[1]) }.groupTuple() // [amplicon, [cluster1,...,clustern]]
    // [amplicon, cluster, ampliconsize]
    amplicon_size_cluster_map = amplicon_cluster_map.join(amplicon_size_fwd).transpose()
   
    sample_size = amplicon_size_fwd.map{ it -> it[1] }.sum() 
    
    //[amplicon, cluster, clustersize/ampliconsize, cluster.tsv]
    prediction_input = amplicon_size_cluster_map.join(filtered_variant_ch, by: [0,1]).map{ it -> tuple(it[0], it[1], it[3].countFasta()/it[2], it[4]) }
    // amplicon_size_cluster_map.join(filtered_variant_ch, by: [0,1]).map{ it -> tuple(it[0], it[1], it[3].countFasta(), it[2], it[3].countFasta()/it[2]) }.view()
    // [amplicon, ampliconsize/samplesize]
    amplicon_scaler = amplicon_size_fwd.combine(sample_size).map{ it -> tuple(it[0], it[1]/it[2]) } 
    // amplicon_size.combine(sample_size).join(amplicon_scaler).view()

    lineage_prediction(prediction_input.combine(new_usher))
    // new_usher = lineage_prediction.out.mod_usher
    lineage_prediction.out.log.collectFile(name:'quantify_cluster_abundances.log', storeDir:"$params.output/${params.runinfo}")
    amplicon_cluster_abundances = lineage_prediction.out.abs // [amplicon, cluster, tsv]
    amplicon_ready = amplicon_cluster_abundances.map{ it -> tuple(it[0], it[2]) }.groupTuple().map{ it ->  it[0]} // [amplicon] 

    amplicon_quant(amplicon_ready.join(amplicon_scaler)) // input [amplicon, ampliconsize/samplesize]
    amplicon_abundances = amplicon_quant.out.abs // [amplicon, tsv]
    amplicon_quant.out.log.collectFile(name:'quantify_amplicon_abundances.log', storeDir:"$params.output/${params.runinfo}")

    sample_abundances = amplicon_abundances.map{ it -> it[1] }.collectFile(keepHeader: true, skip: 1, name:"scaled_amplicon_abundances.tsv") // [tsv]

    aggregate_abundances(sample_abundances, new_usher, Channel.value(params.threshold))
    fp_update = aggregate_abundances.out.csv
    if ( fp_update.splitCsv(header: false, skip: 1).map{ row -> row[-1].toInteger() }.sum() == 0 ) {
      final_output = aggregate_abundances.out.abs.collect() // [tsv]    
    }
    else{
      aggregate_abundances.out.fp_log.collectFile(name:"false_positive_detection.log", storeDir:"$params.output/${params.runinfo}")

      // filter variants and reads using false positive annotated usher frame and the filtered variant output from before
      filter_variants_repeat(filtered_variant_ch.join(allele_frequency_ch, by: [0,1]).map{ it -> tuple(it[0], it[1], it[2], it[3], it[6]) }.combine(fp_update).combine(reference).combine(Channel.value(params.unknown)))
      // mod_new_usher = filter_variants_repeat.out.new_usher_mod
      filter_variants_repeat.out.log.collectFile(name:"filter_variants_repeat.log", storeDir: "$params.output/${params.runinfo}")
      // update downstream variant evaluation and scaling factors
      filtered_variant_repeat_ch = filter_variants_repeat.out.filtered.filter{ it[2].countFasta() > 0 }
      sample_variants_repeat =  filtered_variant_repeat_ch.map{ it -> it[3]}.collectFile(name: "sample_variants.tsv", keepHeader:true, skip: 1)
      variant_metrics_repeat(sample_variants_repeat.combine(fp_update).combine(lineageDict))
      variant_metrics_repeat.out.collectFile(name:"variant_evaluation_repeat.log", storeDir: "$params.output/${params.runinfo}")
      // [amplicon, ampliconsize]
      amplicon_size_repeat = filtered_variant_repeat_ch.map{ it -> tuple(it[0], it[2].countFasta()) }.groupTuple().map{ it -> tuple(it[0], it[1].sum()) }
      plot_mutation_profiles_repeat(amplicon_size_repeat)
      amplicon_size_fwd_repeat = plot_mutation_profiles_repeat.out.fwd
      amplicon_cluster_map_repeat = filtered_variant_repeat_ch.map{ it -> tuple(it[0], it[1]) }.groupTuple() // [amplicon, [cluster1,...,clustern]]
      // [amplicon, cluster, ampliconsize]
      amplicon_size_cluster_map_repeat = amplicon_cluster_map_repeat.join(amplicon_size_fwd_repeat).transpose()
    
      sample_size_repeat = amplicon_size_fwd_repeat.map{ it -> it[1] }.sum() 
      
      //[amplicon, cluster, clustersize/ampliconsize, cluster.tsv]
      prediction_input_repeat = amplicon_size_cluster_map_repeat.join(filtered_variant_repeat_ch, by: [0,1]).map{ it -> tuple(it[0], it[1], it[3].countFasta()/it[2], it[4]) }
      // amplicon_size_cluster_map_repeat.join(filtered_variant_repeat_ch, by: [0,1]).map{ it -> tuple(it[0], it[1], it[3].countFasta(), it[2], it[3].countFasta()/it[2]) }.view()
      // [amplicon, ampliconsize/samplesize]
      amplicon_scaler_repeat = amplicon_size_fwd_repeat.combine(sample_size_repeat).map{ it -> tuple(it[0], it[1]/it[2]) } 
      // amplicon_size_repeat.combine(sample_size_repeat).join(amplicon_scaler_repeat).view()

      // analyse clusters with updated data
      lineage_prediction_repeat(prediction_input_repeat.combine(fp_update))
      lineage_prediction_repeat.out.log.collectFile(name:'quantify_cluster_abundances_repeat.log', storeDir:"$params.output/${params.runinfo}")

      // analyse amplicons
      amplicon_cluster_abundances_repeat = lineage_prediction_repeat.out.abs // [amplicon, cluster, tsv]
      amplicon_ready_repeat = amplicon_cluster_abundances_repeat.map{ it -> tuple(it[0], it[2]) }.groupTuple().map{ it ->  it[0]} // [amplicon] 

      amplicon_quant_repeat(amplicon_ready_repeat.join(amplicon_scaler_repeat)) // input [amplicon, ampliconsize/samplesize]
      amplicon_abundances_repeat = amplicon_quant_repeat.out.abs // [amplicon, tsv]
      amplicon_quant_repeat.out.log.collectFile(name:'quantify_amplicon_abundances_repeat.log', storeDir:"$params.output/${params.runinfo}")

      // aggregate
      sample_abundances_repeat = amplicon_abundances_repeat.map{ it -> it[1] }.collectFile(keepHeader: true, skip: 1, name:"scaled_amplicon_abundances.tsv") // [tsv]
      aggregate_abundances_repeat(sample_abundances_repeat, fp_update, Channel.value(0))
      final_output = aggregate_abundances_repeat.out.abs.collect() // [tsv]
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
      lineageDict = Channel.fromPath("$params.output/input/lineageDict.csv")
      fasta = Channel.fromPath("$params.output/input/query.fasta")
    }
    lineageDict.view()
    fasta.view()
  }
  // in real data mode, convert fastq file to fasta
  else {
    if (params.sort_reads) {
      input_ch = Channel.fromPath(params.fastq)
      fastq_to_fasta(input_ch)
      fasta = fastq_to_fasta.out  
      fasta.view()
    }
    lineageDict = Channel.fromPath(params.dummy)
    lineageDict.view()
  }
  
  // sort reads into amplicons or load already sorted amplicon FASTAs and redundancy data
  // from either the simulation or real data input folder
  if (params.sort_reads) {
    get_amplicon_reads(fasta)
    amplicon_fasta = get_amplicon_reads.out.amplicon_ch 
  }
  else {
    amplicon_ch = Channel.fromPath("$params.output/$params.amplicon_data/final/*.fasta")
    // [amplicon, amplicon.fasta]
    amplicon_fasta = amplicon_ch.map{ it -> tuple(it.baseName, it) }
    // [amplicon, amplicon-duplicates.fasta]
  }
  amplicon_fasta.view()
  clustering(amplicon_fasta, lineageDict) 
  cluster_fasta = clustering.out.redundant_cluster_result // [amplicon,clusterX, clusterX.fasta]
  
  prediction(cluster_fasta, lineageDict)
  abundances = prediction.out.final_output

}  

workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}  






def helpMSG() {
    c_reset = "\033[0m";
    c_red = "\033[1;31m"
    c_green = "\033[1;32m";
    c_yellow = "\033[0;33m";
    c_blue = "\033[0;34m";
    c_dim = "\033[2m";
    log.info """
    
    ${c_green}Welcome to ABRACADABRA - Amplicon-Based hieRArchical Clustering Approach for for lineage Detection and ABundance estimation fRom wAstewater sequencing${c_reset}

    ______________________________________________________________________________________________________________________________________________________
   
    ${c_yellow}Usage example:${c_reset}
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

    ${c_yellow}Mandatory Input:${c_reset}
    ${c_green}--fastq FILEPATH${c_reset}                  Path to FASTQ file
    ${c_green}--mixture_fastqs PATH${c_reset}             Path to comma-separated CSV file containing the filepaths to FASTQ files with reads to spike-in and the target abundances.
    ${c_green}--mixture_lineages FILEPATH${c_reset}       Path to comma-separated CSV file mapping FASTQ file names to their lineage annotation (ONE lineage per FASTQ)
    ${c_green}--primer${c_reset}                          Primer sequence file: TSV with two columns -> primer name and primer sequence
    ${c_green}--bed${c_reset}                             Primer scheme (BED)

    ${c_yellow}Other Input:${c_reset}
    ${c_green}--reference${c_reset}                       SC2 index genome (FASTA)   
    ${c_green}--usher${c_reset}                           Mutational barcode file (CSV)

    ${c_yellow}Pipeline Options:${c_reset}
    ${c_green}--mode${c_reset}                            In simulation mode, you can mix a sample from input fastq files at defined abundances and 
                                      annotate sampled reads with lineage labels for evaluation or load already prepared mixture samples with lineageDict.csv
                                      (STRING)
                                                          Use the "real" mode for samples where you do not have lineage annotations on the read level
                                                          [default $params.mode]
    ${c_green}--sort_reads${c_reset}                      If true, input read data is sorted into amplicon fasta files. Else, the pipeline tries to load already sorted 
                                      data from params.output/params.amplicon_data/final (BOOLEAN) [default $params.sort_reads]
    ${c_green}--prepare_mixture${c_reset}                 Mix input fastq samples at predefined abundances (BOOLEAN) [default $params.prepare_mixture]
    ${c_green}--sample_size${c_reset}                     Target sample size when simulating a mixture. Number of reads per component is calculated 
                                      from the defiend abudnance and the sample size (INT) [default $params.sample_size]
    ${c_green}--unknown${c_reset}                         If true, the pipeline uses called SNPs that don't overlap with the barcode reference to detect an "unknown" lineage (BOOLEAN) [default $params.unknown]
    ${c_green}--variant_caller${c_reset}                  Currently, only ivar is supported (STRING)
    ${c_green}--threshold${c_reset}                       Minimal abundance cutoff. If there are lineages predicted with abundances below the threshold, lineage detection
                                      and abundance estimation is repeated without the respective lienages in the reference (FLOAT) [default $params.threshold]
    ${c_green}--min_read_len${c_reset}                    Minimum read length allowed when filtering sorted amplicon read sets (INT) [default $params.min_read_len]
    ${c_green}--max_primer_mismatch${c_reset}             Maximum number of nt mismatches allowed between a primer and a read durint amplicon sorting (INT) [default $params.max_primer_mismatch]
    ${c_green}--umap_params${c_reset}                     Additional parameters for UMAP [default $params.umap_params]
    ${c_green}--hdbscan_params${c_reset}                  Additional parameters for HDBSCAN cluster analysis [default $params.hdbscan_params]
                                      For more information and options for hdbscan and umap, please use
                                      ${c_green}nextflow run main.nf --ClusterHelp${c_reset}
    
    ${c_yellow}Other Options:${c_reset}
    ${c_green}--cores INT${c_reset}                       max cores per process for local use [default $params.cores]
    ${c_green}--max_cores INT${c_reset}                   max cores used on the machine for local use [default $params.max_cores]
    ${c_green}--memory INT${c_reset}                      max memory in GB for local use [default $params.memory]
    ${c_green}--output PATH${c_reset}                     name of the result folder [default $params.output]
                                          
    ${c_yellow}Nextflow options:${c_reset}
    ${c_reset}____________________________________________________________________________________________

    """.stripIndent()
}

def ClusterHelp() {
  c_reset = "\033[0m";
  c_red = "\033[1;31m"
  c_green = "\033[1;32m";
  c_yellow = "\033[0;33m";
  c_blue = "\033[0;34m";
  c_dim = "\033[2m";
  log.info """
  ____________________________________________________________________________________________
  The input reads are UMAPped based on their kmer frequency representation and then clustered 
  into groups (clades) based on their sequence similarity using HDBSCAN.
  ____________________________________________________________________________________________
  Python Dependencies:
    - python=3.8.15
    - numpy=1.21.0
    - biopython=1.80
    - umap-learn=0.5.3
    - hdbscan=0.8.29
    - joblib=1.1.0
    - docopt=0.6.2
    - scipy=1.4.1
    - pandas=1.5.2
    - matplotlib=3.3.2
    - datashader>=0.5.4
    - bokeh=2.4.3
    - holoviews=1.15.1
    - colorcet=3.0.1
    - scikit-image=0.19.3

  Usage:
    hdbscan_virus.py [options] <inputSequences> <lineageDict>

  Options:
    -h, --help                              Show this help message and exits.
    -v, --verbose                           Get some extra information from viralClust during calculation. [Default: False]
    -o DIR, --output DIR                    Specifies the output directory of viralClust. [Default: pwd]
    -p PROCESSES, --process PROCESSES       Specify the number of CPU cores that are used. [Default: 1]
    -k KMER, --kmer KMER                    Length of the considered kmer. [Default: 5]
    --umap_metric UMAP_METRIC               Distance metric applied by  UMAP. The following are supported:
                                            'euclidean', 'manhatten', 'chebyshev', 'minkwoski',
                                            'canberra', 'braycurtis',
                                            'mahalanobis', 'wminkowski', 'seuclidean',
                                            'cosine'.
                                            If an invalid metric is set, the default is selected. [Default: cosine] 
    --n_components N_COMPONENTS             This parameter tells UMAP how many dimensions should be used for the embedding.      [Default: 2] 
    --n_neighbors N_NEIGHBOR S              Number of neighbors considered by UMAP to reduce the dimension space. [Default: 15]
    --min_dist MIN_DIST                     Sets the threshold for the minimum distance of two points in the low-dimensional space. [Default: 0.1]
    --hdbscan_metric HDBSCAN_METRIC         Distance metric applied by  HDBSCAN.
                                            The following are supported:
                                            'euclidean', 'manhatten', 'chebyshev', 'minkwoski',
                                            'canberra', 'braycurtis',
                                            'mahalanobis', 'wminkowski', 'seuclidean',
                                            'cosine'.
                                            If an invalid metric is set, the default will be selected
                                            [Default: cosine]
    --min_samples MIN_SAMPLES               Intuitively, this parameter declares how conservative clustering is performed. Higher values will lead 
                                            to more points considered noise, whereas a low value causes "edge-cases" to be grouped into a cluster.
                                            The default parameter is the same as CLUSTERSIZE.  [Default: 5] 
    --min_cluster_size MIN_CLUSTER_SIZE     This parameter forces HDBSCAN to form cluster with a size larger-equal to CLUSTERSIZE.
                                            Be aware that some data points (i.e. genomes) or even whole subcluster are considered as noise, if this parameter is set too high.
                                            E.g., if a very distinct viral genus has 40 genomes and the parameter is set to anything >40, HDBSCAN will not form
                                            the genus specific cluster.        [Default: 5]  
    --cluster_selection_epsilon CLUSTER_SELECTION_EPSILON               A distance threshold. Clusters below this value will be merged.     [Default: 0.05] 

  """.stripIndent()
}
