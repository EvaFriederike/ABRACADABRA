/************************************************************************
Clip primers from input reads.
************************************************************************/

process primer_clipping {
  label 'medaka'
  publishDir "${params.output}/${params.amplicon_data}/primerclipped", mode: 'copy', pattern: "${primer_name}.sorted.primerclipped.fasta"
  publishDir "${params.output}/${params.amplicon_data}/final", mode: 'copy', pattern: "${primer_name}.sorted.primerclipped.fasta", saveAs: {filename -> "${primer_name}.fasta"}

  input:
    tuple val(primer_name), path(primer_reads), path(bed), path (reference)

  output:
    tuple val(primer_name), path("${primer_name}.sorted.primerclipped.fasta"), emit: primerclipped_fasta

  script:
  """
    #!/bin/bash

    minimap2 -ax map-ont $reference $primer_reads -t $task.cpus > ${primer_name}.sam
    samtools view -S -b -h ${primer_name}.sam > ${primer_name}.bam
    samtools sort ${primer_name}.bam > ${primer_name}.sorted.bam
    samtools index ${primer_name}.sorted.bam
  
    python ${baseDir}/bin/primerbed2bedpe.py $bed 
    bamclipper.sh -b ${primer_name}.sorted.bam -p primer.bedpe -n 5
    samtools fasta ${primer_name}.sorted.primerclipped.bam > ${primer_name}.sorted.primerclipped.fasta
  """
}
