/************************************************************************
Update reference data set. Remove all lineages that were predicted with an 
abudnance below a certain cutoff value.
************************************************************************/
process update_reference {
    label 'evaluate'
    // publishDir "${params.output}/${params.runinfo}/", mode: 'copy', pattern: ".command.log", saveAs: {filename -> "update_reference.log"}

    input:
    path usher

    output:
    path "mod_${usher}"
    //path ".command.log", emit: log

    script:
    """
    #!/bin/bash

    echo "Add column for false positive annotation"
    python ${baseDir}/bin/update_reference.py $usher mod_${usher}
    """



}

/************************************************************************
Detect lineages in a read cluster based on mutation comparison. Unknown mutations
are processed as belonging to a lineage called "unknown". Cluster abundances
are scaled by #cluster_reads / #amplicon_reads
NOTES:
- Remove $launchDir when running on HPC
************************************************************************/
process lineage_prediction {
    label 'evaluate'
    publishDir "${params.output}/${params.mode}/tmp/${amplicon}/", mode: 'copy', pattern: "scaled_${cluster}_abundances.tsv"
    publishDir "${params.output}/${params.mode}/${params.abundances}/${amplicon}/", mode: 'copy', pattern: "${cluster}_abundances.tsv"

    input:
    tuple val(amplicon), val(cluster), val(clusterscaler), path(vcf_tsv), path(usher)

    output:
    tuple val(amplicon), val(cluster), path("scaled_${cluster}_abundances.tsv"), emit: abs
    path "${cluster}_abundances.tsv"
    path ".command.log", emit:log

    script:
        """
        #!/bin/bash
        echo "------------------- Amplicon ${amplicon} - cluster ${cluster}-------------------"
        python ${baseDir}/bin/detect_and_estimate.py $vcf_tsv $usher ${cluster}_abundances.tsv
        echo "Scaling for proportion of cluster reads among the whole amplicon"
        python $baseDir/bin/scale_abundances.py ${cluster}_abundances.tsv $clusterscaler scaled_${cluster}_abundances.tsv
        """
}

/************************************************************************
Sum lineages abundances across all amplicon clusters. Amplicon abundances
are scaled by #amplicon_reads / #sample_reads. Where #sample_reads are 
the number of reads that could be sorted into an amplicon set.
NOTE:
- remove launchDir when running on HPC
************************************************************************/
process amplicon_quant {
    label 'evaluate'
    publishDir "${params.output}/${params.mode}/${params.abundances}/${amplicon}/", mode: 'copy', pattern: "${amplicon}_abundances.tsv"

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
    echo "Sum lineage abundances across all clusters"
    #python ${baseDir}/bin/aggregate_cluster.py ${launchDir}/${params.output}/tmp/${amplicon}/ ${amplicon}_abundances.tsv
    python ${baseDir}/bin/aggregate_cluster.py ${params.output}/${params.mode}/tmp/${amplicon}/ ${amplicon}_abundances.tsv
    echo "Scaling for proportion of amplicon reads among all sample reads"
    python ${baseDir}/bin/scale_abundances.py ${amplicon}_abundances.tsv $amplicon_scaler scaled_${amplicon}_abundances.tsv

    rm -rf ${params.output}/${params.mode}/tmp/${amplicon}/
    """

}

/************************************************************************
Sum lineages abundances across all amplicons. In case of a threshold, 
potential false positive detection are filtered and returned to initiate
a second lienage detection and quantification round with the updated 
reference data set.
************************************************************************/
process aggregate_abundances {
    label 'evaluate'
    publishDir "${params.output}/${params.mode}/final", mode: 'copy', pattern: "final_abundances.tsv"
    publishDir "${params.output}/${params.mode}/${params.runinfo}", mode: 'copy', pattern: ".command.log", saveAs: {filename -> "sample_abundances.log"}


    input:
    path tsv
    path usher
    val threshold

    output:
    path "final_abundances.tsv", emit: abs
    path "barcodes_FP_annotated.csv", emit: csv, optional: true
    path ".command.log", emit: log
    path "false_positive_detection.log", emit:fp_log,  optional: true

    script:
    """
    #!/bin/bash

    echo "Aggregate amplicon abundance across the sample"

    python ${baseDir}/bin/aggregate_sample.py $tsv $threshold $usher final_abundances.tsv
    
    rm -rf  ${params.output}/${params.mode}/tmp/

    """

}