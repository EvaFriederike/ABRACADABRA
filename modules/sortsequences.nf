/************************************************************************
Sort all input sequences (fasta) into positive strands (i.e. direction 
with highest ORF density).
************************************************************************/

process sort_sequences {
  label 'sortseq'
  publishDir "${params.output}/${params.amplicon_data}/positive_sorted", mode: 'copy', pattern: "*_positive.fasta"

  input:
    tuple val(primer_name), path(primer_fasta) 

  output:
    tuple val(primer_name), path("${primer_name}_positive.fasta")

  script:
  """
    #!/bin/bash

    python3 ${projectDir}/bin/sort_sequences.py $primer_fasta > ${primer_name}_positive.fasta
  """
}


/**************************************************************************************************
Filter reads by a cutoff for minimum read length.
**************************************************************************************************/
process read_filtering {
  label 'sortseq'

  input:
    tuple val(name), path(fasta)

  output:
    tuple val(name), path("filtered_${fasta}"), emit: filtered
    path ".command.log", emit: log

  script:
    """
    #!/bin/bash

    echo "--------------Filter ${name} for very short reads"
    seqkit seq -m $params.min_read_len $fasta > filtered_${fasta}
    B=\$(grep -c '>' $fasta)
    A=\$(grep -c '>' filtered_${fasta})
    echo "Number of reads before filtering: \$B"
    echo "Number of reads after filtering: \$A"
    if [ \$A -eq 0 ]; then
      echo "No remaining reads after filtering."
    fi
    """
}



/**************************************************************************************************
Sort reads into amplicons based on input primer scheme. Only reads with less than 10% N count are
considered.
If a read contains multiple primers, it is split into subreads and the subreads are sorted into the
respective amplicon. 
**************************************************************************************************/
process primer_sort {
  label 'sortseq'
  publishDir "${params.output}/${params.amplicon_data}/primer_sorted", mode: 'copy', pattern: "*.fasta"
  publishDir "${params.output}/${params.runinfo}", mode: 'copy', pattern: "primer_sort.log"
  
  input:
    tuple path(fasta), path(primer)

  output:
    path "*_{LEFT,RIGHT}.fasta", emit: primer_reads
    path "*_primer_match*.fasta"
    path "primer_sort.log"
    
  script:
  """
    #!/bin/bash
    readcount=\$(grep -c '>' ${fasta.baseName}.fasta)
    python ${baseDir}/bin/primer_sort.py $fasta $primer \$readcount $params.max_primer_mismatch
        
  """
}

/**************************************
Strategy might result in some reads being sorted into 1-3 amplicons,
depending on the range of overlap -> check the amount of concerned reads
Also, seqkit grep detected primer matches that regex finditer didn't...?
*******************************************/
process primer_sort_edit {
  label 'sortseq'
  
  input:
    tuple path(fasta), val(amplicon), val(primer1), val(primer2)

  output:
    tuple val(amplicon), path("${amplicon}.fasta"), emit: primer_reads
    path ".command.log", emit: log
    
  script:
  """
    #!/bin/bash
    echo "-------------------Amplicon $amplicon-------------------------------"
    echo $primer1 > primer.tsv
    echo $primer2 >> primer.tsv

    seqkit grep -m $params.max_primer_mismatch -f primer.tsv $fasta > matches.fasta
    echo "Primers are matching with \$(grep -c '>' matches.fasta) reads"
    python3 ${projectDir}/bin/split_reads.py matches.fasta $primer1 $primer2 $params.max_primer_mismatch

    mv split_reads.fasta ${amplicon}.fasta
       
  """
}

process log_primer_sort {
  label 'sortseq'

  input:
    tuple path(ampliconfasta), path(infasta)

  output:
    path ".command.log"

  script:
  """
    #!/bin/bash
    echo "____________________________________________________"
    echo "\$(grep -c '>' $ampliconfasta) of \$(grep -c '>' $infasta) reads were sorted into amplicon sets."
       
  """
}



process redirect_reads {
  label 'sortseq'
  
  input:
    path fasta

  output:
    path "dedup_positive_${fasta}"
    
  script:
  """
    #!/bin/bash
    echo "Sorting input reads into positive strand direction"
    python3 ${projectDir}/bin/sort_sequences.py $fasta > positive_${fasta}
    seqkit rmdup -n -o dedup_positive_${fasta} positive_${fasta}
    echo "After removing duplicates by header, \$(grep -c '>' positive_${fasta}) of\$(grep -c '>' $fasta) reads remain."
  """
}