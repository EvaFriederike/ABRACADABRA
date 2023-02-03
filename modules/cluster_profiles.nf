/************************************************************************
Option:
- Call variants for each primer set first before merging into amplicon variant set
- code snippets:
//  #if [ -s $left_fasta ]; then
//    # minimap2 -ax map-ont $reference $left_fasta > ${left_fasta.baseName}.sam
//    # samtools view -S -b -h ${left_fasta.baseName}.sam > ${left_fasta.baseName}.bam
//    # samtools sort ${left_fasta.baseName}.bam > ${left_fasta.baseName}.sorted.bam
//    # samtools index -c ${left_fasta.baseName}.sorted.bam  

//    # medaka consensus --model r941_min_sup_g507 --threads 4 --chunk_len 800 --chunk_ovlp 400 ${left_fasta.baseName}.sorted.bam ${left_fasta.baseName}.medaka.consensus.hdf
//    # medaka variant $reference ${left_fasta.baseName}.medaka.consensus.hdf ${left_fasta.baseName}.medaka.vcf
//    # medaka tools annotate ${left_fasta.baseName}.medaka.vcf $reference  ${left_fasta.baseName}.sorted.bam ${left_fasta.baseName}.medaka.annotate.vcf
//   #fi

//   #if [ -s $right_fasta ]; then
//    # minimap2 -ax map-ont $reference $right_fasta > ${right_fasta.baseName}.sam
//    # samtools view -S -b -h ${right_fasta.baseName}.sam > ${right_fasta.baseName}.bam
//    # samtools sort ${right_fasta.baseName}.bam > ${right_fasta.baseName}.sorted.bam
//    # samtools index -c ${right_fasta.baseName}.sorted.bam  

//    # medaka consensus --model r941_min_sup_g507 --threads 4 --chunk_len 800 --chunk_ovlp 400 ${right_fasta.baseName}.sorted.bam ${right_fasta.baseName}.medaka.consensus.hdf
//    # medaka variant $reference ${right_fasta.baseName}.medaka.consensus.hdf ${right_fasta.baseName}.medaka.vcf
//    # medaka tools annotate ${right_fasta.baseName}.medaka.vcf $reference  ${right_fasta.baseName}.sorted.bam ${right_fasta.baseName}.medaka.annotate.vcf
//   #fi

//   #if [ -s $left_fasta ] && [ -s $right_fasta ]
//   #then
//   #  bgzip -f ${left_fasta.baseName}.medaka.annotate.vcf
//   #  bgzip -f ${right_fasta.baseName}.medaka.annotate.vcf
//   #  bcftools index -f ${left_fasta.baseName}.medaka.annotate.vcf.gz
//   #  bcftools index -f ${right_fasta.baseName}.medaka.annotate.vcf.gz
//   #  bcftools merge --force-samples --missing-to-ref  -o ${right_fasta.baseName.split('_')[0]}_merged.vcf.gz -O z ${left_fasta.baseName}.medaka.annotate.vcf.gz ${right_fasta.baseName}.medaka.annotate.vcf.gz
//   #  bgzip -d ${right_fasta.baseName.split('_')[0]}_merged.vcf.gz
//   #elif [ -s $left_fasta ] && [ ! -s $right_fasta ]
//   #  then
//   #    cp ${left_fasta.baseName}.medaka.annotate.vcf ${left_fasta.baseName.split('_')[0]}_merged.vcf
//   #elif  [ ! -s $left_fasta ] && [ -s $right_fasta ]
//   #then
//   #  cp ${right_fasta.baseName}.medaka.annotate.vcf ${right_fasta.baseName.split('_')[0]}_merged.vcf
//   #else 
//   #  touch ${left_fasta.baseName.split('_')[0]}_merged.vcf
//   #fi

************************************************************************/
// process split_fasta_by_primer {
//   publishDir "${params.output}/${params.hdbscan_output}/split_clusters/${amplicon}", mode: 'copy', pattern: "${fasta.baseName}_*.fasta"


//   input:
//     tuple val(amplicon), path(fasta), path(primerDict)

//   output:
//     tuple val(amplicon), path("${fasta.baseName}_LEFT.fasta"), path("${fasta.baseName}_RIGHT.fasta")

//   script:
//   """
//   #!/bin/bash

//   python ${projectDir}/bin/split_fasta.py $fasta $primerDict 
//   mv LEFT.fasta ${fasta.baseName}_LEFT.fasta
//   mv RIGHT.fasta ${fasta.baseName}_RIGHT.fasta

//   """
// }


/************************************************************************
Call variants within an amplicon.
************************************************************************/
process call_variants {
  label 'medaka'
  publishDir "${params.output}/${params.mode}/${params.hdbscan_output}/variants/${amplicon}", mode: 'copy', pattern: "*.medaka.annotate.vcf.tsv"

  input:
    tuple val(amplicon), val(amplicon_size), val(cluster), path(fasta), path(reference)

  output:
    tuple val(amplicon), val(amplicon_size), val(cluster), path("*.tsv"), emit: cluster_AF
    path ".command.log", emit: log
   
  script:
  """
  #!/bin/bash

  echo "------------------- Amplicon ${amplicon} -------------------\n"

  minimap2 -ax map-ont $reference $fasta -t $task.cpus > ${cluster}.sam
  samtools view -S -b -h ${cluster}.sam > ${cluster}.bam
  samtools sort ${cluster}.bam > ${cluster}.sorted.bam
  samtools index -c ${cluster}.sorted.bam  

  medaka consensus --model r941_min_sup_g507 --threads $task.cpus --chunk_len 800 --chunk_ovlp 400 ${cluster}.sorted.bam ${cluster}.medaka.consensus.hdf
  medaka variant $reference ${cluster}.medaka.consensus.hdf ${cluster}.medaka.vcf

  if [ \$(grep -c -v '^#' ${cluster}.medaka.vcf) -gt 0 ]; then
    medaka tools annotate ${cluster}.medaka.vcf $reference  ${cluster}.sorted.bam ${cluster}.medaka.annotate.vcf
    medaka tools vcf2tsv ${cluster}.medaka.annotate.vcf
  else
    touch fail.tsv
    echo "--------------------------------------------------------\n"
    echo "No variants detected for ${cluster}"
    n_reads=\$(grep -c '>' $fasta)
    echo "This leads to missing information from \${n_reads} reads and reduces the overall relative abundances of detected lineages in the analysed sample!!\n"
  fi
  """
}



// process mutation_heatmap {
//   publishDir "${params.output}/${params.eval_output}", mode: 'copy', pattern: "*medaka_profile.png"

//   input:
//     val bam_path

//   output:
//     path "${bam_path.split('/')[-1]}_medaka_profile.png"

//   script:
//   """
//   Rscript ${projectDir}/bin/MutationProfile.r medaka $bam_path *.vcf medaka_profile
//   """

// }