/************************************************************************
* EVALUATE CLUSTER PROFILES
************************************************************************/

process call_variants {
  label 'medaka'
  publishDir "${params.output}/${params.hdbscan_output}/medaka_variants", mode: 'copy', pattern: "*.medaka.annotate.vcf"


  input:
    tuple path(fasta), path(reference)

  output:
    path "${fasta.baseName}.medaka.annotate.vcf"
    val "${params.output}/${params.hdbscan_output}/medaka_variants", emit: bam_path

  script:
  """
  #!/bin/bash

  minimap2 -ax map-ont $reference $fasta > ${fasta.baseName}.sam
  samtools view -S -b -h ${fasta.baseName}.sam > ${fasta.baseName}.bam
  samtools sort ${fasta.baseName}.bam > ${fasta.baseName}.sorted.bam
  samtools index -c ${fasta.baseName}.sorted.bam  
  
  medaka consensus --model r941_min_sup_g507 --threads 4 --chunk_len 800 --chunk_ovlp 400 ${fasta.baseName}.sorted.bam ${fasta.baseName}.medaka.consensus.hdf
  medaka variant $reference  ${fasta.baseName}.medaka.consensus.hdf   ${fasta.baseName}.medaka.vcf
  medaka tools annotate   ${fasta.baseName}.medaka.vcf $reference  ${fasta.baseName}.sorted.bam   ${fasta.baseName}.medaka.annotate.vcf
  """
}



process mutation_heatmap {
  publishDir "${params.output}/${params.hdbscan_output}", mode: 'copy', pattern: "medaka_profile.png"

  input:
    val bam_path

  output:
    path "medaka_profile.png"

  script:
  """
  Rscript /scratch/assmanne/code-dump/MutationProfile.r medaka ${projectDir}/../${bam_path}/ *.vcf medaka_profile
  """

}