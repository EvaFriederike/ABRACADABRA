/************************************************************************
Detect lineages in a read cluster based on mutation comparison. Unknown mutations
are processed as belonging to a lineage called "unknown". Cluster abundances
are scaled by #cluster_reads / #amplicon_reads
************************************************************************/
process lineage_prediction {
    label 'evaluate'
    publishDir "${params.output}/tmp/${amplicon}/", mode: 'copy', pattern: "scaled_${cluster}_abundances.tsv"
    publishDir "${params.output}/${params.abundances}/${amplicon}/", mode: 'copy', pattern: "*${cluster}_abundances.tsv"

    input:
    tuple val(amplicon), val(cluster), val(clusterscaler), path(vcf_tsv), path(usher)

    output:
    tuple val(amplicon), val(cluster), path("scaled_${cluster}_abundances.tsv"), emit: abs
    // path "mod_${usher}", emit: mod_usher, optional:true
    path "${cluster}_abundances.tsv"
    path ".command.log", emit:log

    script:
        """
        #!/bin/bash
        echo "------------------- Amplicon ${amplicon} - cluster ${cluster}-------------------"
        echo \$PWD
        python ${baseDir}/bin/detect_and_estimate.py $vcf_tsv $usher $params.unknown ${cluster}_abundances.tsv
        
        echo "Scaling for proportion of cluster reads among the whole amplicon"
        python $baseDir/bin/scale_abundances.py ${cluster}_abundances.tsv $clusterscaler scaled_${cluster}_abundances.tsv
        """
}

/************************************************************************
Sum lineages abundances across all amplicon clusters. Amplicon abundances
are scaled by #amplicon_reads / #sample_reads. Where #sample_reads are 
the number of reads that could be sorted into an amplicon set.
NOTE: remove launchdir when running on hpc
************************************************************************/
process amplicon_quant {
    label 'evaluate'
    publishDir "${params.output}/${params.abundances}/${amplicon}/", mode: 'copy', pattern: "*${amplicon}_abundances.tsv"

    input:
    tuple val(amplicon), val(amplicon_scaler)

    output:
    tuple val(amplicon), path("scaled_${amplicon}_abundances.tsv"), emit: abs
    path "${amplicon}_abundances.tsv"
    path ".command.log", emit: log

    script:
    """
    #!/bin/bash
    echo "------------------- Amplicon ${amplicon} -------------------"
    echo "\$PWD"
    echo "Sum lineage abundances across all clusters"
    
    #python ${baseDir}/bin/aggregate_cluster.py ${launchDir}/${params.output}/tmp/${amplicon}/ ${amplicon}_abundances.tsv
    python ${baseDir}/bin/aggregate_cluster.py ${params.output}/tmp/${amplicon}/ ${amplicon}_abundances.tsv

    echo "Scaling for proportion of amplicon reads among all sample reads"
    python ${baseDir}/bin/scale_abundances.py ${amplicon}_abundances.tsv $amplicon_scaler scaled_${amplicon}_abundances.tsv

    #rm -rf ${launchDir}/${params.output}/tmp/${amplicon}/
    rm -rf ${params.output}/tmp/${amplicon}/
    """

}

/************************************************************************
Sum lineages abundances across all amplicons. In case of a threshold, 
potential false positive detection are filtered and returned to initiate
a second lienage detection and quantification round with the updated 
reference data set.
NOTE_ remove launchdir when running on hpc
************************************************************************/
process aggregate_abundances {
    label 'evaluate'
    publishDir "${params.output}/final", mode: 'copy', pattern: "*_sample_abundances.tsv"
    publishDir "${params.output}/${params.runinfo}", mode: 'copy', pattern: ".command.log", saveAs: {filename -> "sample_abundances.log"}


    input:
    path tsv
    path usher
    val threshold

    output:
    path "*_sample_abundances.tsv", emit: abs
    path "barcodes_FP_annotated.csv", emit: csv, optional: true
    path ".command.log", emit: log
    path "false_positive_detection.log", emit:fp_log,  optional: true

    script:
    """
    #!/bin/bash

    echo "Aggregate amplicon abundance across the sample"

    python ${baseDir}/bin/aggregate_sample.py $tsv $threshold $usher $params.unknown
    
    #rm -rf  ${launchDir}/${params.output}/tmp/
    rm -rf  ${params.output}/tmp/
    """

}