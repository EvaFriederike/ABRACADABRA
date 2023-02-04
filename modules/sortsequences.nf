/************************************************************************
Sort all input sequences (fasta) into positive strands (i.e. direction 
with highest ORF density).
************************************************************************/

process sort_sequences {
  label 'sortseq'
  publishDir "${params.output}/${params.mode}/${params.amplicon_data}/positive_sorted", mode: 'copy', pattern: "*_positive.fasta"

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
  publishDir "${params.output}/${params.mode}/${params.amplicon_data}/primer_sorted", mode: 'copy', pattern: "*.fasta"
  publishDir "${params.output}/${params.mode}/${params.runinfo}", mode: 'copy', pattern: "primer_sort.log"
  
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