/************************************************************************
Convert fastq to fastq file
************************************************************************/

process fastq_to_fasta {
  label 'sortseq'
 
  input:
    path fastq

  output:
    path "${fastq.baseName}.fasta"

  script:
  """
    #!/bin/bash

    fastq_to_fasta -n -i $fastq -o ${fastq.baseName}.fasta
    sed -i 's/ .*//g' ${fastq.baseName}.fasta
    sed -i 's/@//g'  ${fastq.baseName}.fasta

  """
}
