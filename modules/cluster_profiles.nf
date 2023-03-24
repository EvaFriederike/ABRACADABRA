process update_reference {
  label 'sortseq'

  input:
   path usher
  
  output:
    path "mod_${usher}"
   
  script:
  """
  #!/bin/bash
  
  python3 ${projectDir}/bin/update_reference.py $usher

  """
}



/************************************************************************
Call variants within an amplicon.
Notes:
- maybe I need to use another alingment method or another alignment mode 
than map-ont to map Illunina reads?
- ATTENTION: medaka mode is deprecated

Also filter potential false positive variants if unknown mode is set
Notes:
- bootstrap variant filter?
- bias:
    - generate a global cutoff in case there are no known mutations in a cluster? 
        (implement filter baseline in variant calling process)
    - filter unknown variants if there are known variants i na cluster, else exclude unknown variants
      from downstream analysis
- EDIT:#  // # remove
      // # echo "filter variants for SNPs only and potential false positives"
      // # python $projectDir/bin/filter_variants.py ${cluster}.spike.variants.tsv $usher $params.unknown ${cluster}.spike.sorted.bam $reference
      // # mkdir obsolete
      // # mv ${cluster}.spike.variants.tsv obsolete
      // # VAR=\$(wc -l filtered_${cluster}.spike.variants.tsv | cut -d' ' -f1)
      // # if [[ \$VAR -le 1 ]]; then
      // #  rm filtered_${cluster}.spike.variants.tsv 
      // #  touch fail.tsv
      // #  echo "--------------------------------------------------------"
      // #  echo "No variants remaining for ${cluster} after filtering"
      // #  n_reads=\$(grep -c '>' $fasta)
      // #  echo "This leads to missing information from \${n_reads} reads and reduces the overall relative abundances of detected lineages in the analysed sample!!"
      // #  touch out_$fasta
      // #elif [ -f filtered-reads.txt ]; then
      // #  seqkit grep -v -f filtered-reads.txt $fasta -o filtered_$fasta
      // #else
      // #  cp $fasta out_$fasta
      // #fi     
************************************************************************/

process call_variants {
  label 'medaka'
  publishDir "${params.output}/${params.hdbscan_output}/variants/${amplicon}", mode: 'copy', pattern: "*.variants.tsv"
  //publishDir "${params.output}/${params.eval_output}/${amplicon}/", mode: 'copy', pattern: "AF_histogram.png", saveAs: {filename -> "${cluster}_AF_histogram.png"}

  input:
    tuple val(amplicon), val(cluster), path(fasta), path(reference)
    //, path(usher)

  output:
    //tuple val(amplicon), val(cluster), path("*_${fasta}"), path("*.tsv"), emit: cluster_AF
    tuple val(amplicon), val(cluster), path("*.spike.fasta"), path("*.tsv"), path("${cluster}.spike.sorted.bam"), emit: cluster_AF
    path ".command.log", emit: log
    //path "*.png", optional: true
   
  script:
  if ("${params.variant_caller}" == "medaka")
    """
    #!/bin/bash

    echo "------------------- Amplicon ${amplicon} - ${cluster} -------------------"
    echo \$PWD

    #minimap2 -ax map-ont $reference $fasta -t $task.cpus > ${cluster}.sam
    #samtools view -S -b -h ${cluster}.sam > ${cluster}.bam
    #samtools sort ${cluster}.bam > ${cluster}.sorted.bam
    #samtools index -c ${cluster}.sorted.bam  

    #medaka consensus --model r941_min_sup_g507 --threads $task.cpus ${cluster}.sorted.bam ${cluster}.medaka.consensus.hdf
    #medaka variant $reference ${cluster}.medaka.consensus.hdf ${cluster}.medaka.vcf

    #if [ \$(grep -c -v '^#' ${cluster}.medaka.vcf) -gt 0 ]; then
    #  medaka tools annotate ${cluster}.medaka.vcf $reference  ${cluster}.sorted.bam ${cluster}.medaka.annotate.vcf
    #  medaka tools vcf2tsv ${cluster}.medaka.annotate.vcf

    #  cp $fasta out_$fasta
    #else
    #  touch fail.tsv
    #  echo "--------------------------------------------------------"
    #  echo "No variants detected for ${cluster}"
    #  n_reads=\$(grep -c '>' $fasta)
    #  echo "This leads to missing information from \${n_reads} reads and reduces the overall relative abundances of detected lineages in the analysed sample!!"
    #  touch out_$fasta
    #fi
    """
  else if ("${params.variant_caller}"=='ivar')
  """
    #!/bin/bash

    echo "------------------- Amplicon ${amplicon} - ${cluster} -------------------"
    echo \$PWD

    minimap2 -ax $params.minimap_params $reference $fasta -t $task.cpus > ${cluster}.sam
    samtools view -S -b -h ${cluster}.sam > ${cluster}.bam
    samtools sort ${cluster}.bam > ${cluster}.sorted.bam
    samtools index -c ${cluster}.sorted.bam  

    echo "filter reads to cover spike region"
    samtools view ${cluster}.sorted.bam NC_045512.2:21563-25384 -b -h > ${cluster}.spike.sorted.bam
    echo "call variants"
    samtools mpileup -aa -A -d 600000 -B -Q 20 -q 0 -B -f $reference ${cluster}.spike.sorted.bam | ivar variants -p ${cluster}.spike.variants -q 20 -t 0 -r $reference  
    
    VAR=\$(wc -l ${cluster}.spike.variants.tsv | cut -d' ' -f1)
    if [[ \$VAR -le 1 ]]; then
      rm ${cluster}.spike.variants.tsv 
      touch fail.tsv
      echo "--------------------------------------------------------"
      echo "No variants detected for ${cluster}"
      n_reads=\$(grep -c '>' $fasta)
      echo "This leads to missing information from \${n_reads} reads and reduces the overall relative abundances of detected lineages in the analysed sample!!"
      touch empty.spike.fasta
    else
      # new
      samtools fasta ${cluster}.spike.sorted.bam > ${cluster}.spike.fasta
      # end
    fi

    """
}

process filter_variants {
  label 'medaka'
  publishDir "${params.output}/${params.eval_output}/${amplicon}/", mode: 'copy', pattern: "AF_histogram.png", saveAs: {filename -> "${cluster}_AF_histogram.png"}

  input:
     tuple val(amplicon), val(cluster), path(fasta), path(tsv), path(bam), path(usher), path(reference), val(unknown)
    

  output:
    tuple val(amplicon), val(cluster), path("*_${fasta}"), path("*.tsv"), emit: filtered
    path ".command.log", emit: log
    path "*.png", optional: true
   
  script:
    """
    #!/bin/bash
    echo "------------------- Amplicon ${amplicon} - ${cluster} -------------------"
    echo \$PWD
    
    echo "Filter variants: keep SNPs only and remove potential false positives"
    samtools index $bam
    python $projectDir/bin/filter_variants.py $tsv $usher $unknown $bam $reference

    mkdir obsolete
    mv $tsv obsolete
    
    VAR=\$(wc -l filtered_${tsv} | cut -d' ' -f1)
    if [[ \$VAR -le 1 ]]; then
      rm filtered_${tsv}
      touch fail.tsv
      echo "--------------------------------------------------------"
      echo "No variants remaining after filtering"
      n_reads=\$(grep -c '>' $fasta)
      echo "This leads to missing information from \${n_reads} reads!!"
      touch empty_$fasta
    elif [ -f filtered-reads.txt ]; then
      seqkit grep -f filtered-reads.txt $fasta -o filtered_$fasta
      echo "Keeping \$(wc -l filtered-reads.txt | cut -d' ' -f1) of \$(wc -l obsolete/$tsv | cut -d' ' -f1) reads"
    else
      echo "WARNING: Why is there no list of selected reads?"
    fi     
    """
}


process filter_variants_repeat {
  label 'medaka'
  
  input:
     tuple val(amplicon), val(cluster), path(fasta), path(tsv), path(bam), path(usher), path(reference), val(unknown)
    

  output:
    tuple val(amplicon), val(cluster), path("*_${fasta}"), path("*.tsv"), emit: filtered
    // path ("mod_${usher}"), emit: new_usher_mod
    path ".command.log", emit: log
   
  script:
    """
    #!/bin/bash
    echo "------------------- Amplicon ${amplicon} - ${cluster} -------------------"
    echo \$PWD

    echo "Filter variants: remove variants that were associated with false positive lineages only"
    samtools index $bam
    python $projectDir/bin/filter_variants_repeat.py $tsv $usher $unknown $bam $reference
    
    mkdir obsolete
    mv $tsv obsolete

    VAR=\$(wc -l filtered_${tsv} | cut -d' ' -f1)
    if [[ \$VAR -le 1 ]]; then
      rm filtered_${tsv}
      touch fail.tsv
      echo "No variants remaining for ${cluster} after filtering"
      n_reads=\$(grep -c '>' $fasta)
      echo "This leads to missing information from \${n_reads} reads!!"
      touch empty_$fasta
    elif [ -f filtered-reads.txt ]; then
      seqkit grep -f filtered-reads.txt $fasta -o filtered_$fasta
      echo "Keeping \$(wc -l filtered-reads.txt | cut -d' ' -f1) of \$(wc -l obsolete/$tsv | cut -d' ' -f1) reads"
    else
      echo "WARNING: Why is there no list of selected reads?"
    fi     
    """
}


process variant_metrics {
  label 'sortseq'
  publishDir "${params.output}/${params.eval_output}", mode: 'copy', pattern: "*.csv"

  input:
    tuple path(tsv), path(usher), path(lineageDict)

  output:
    path "variant_evaluation.log"

  script:
  """
  #!/bin/bash
  python3 ${projectDir}/bin/evaluate_variants.py $tsv $usher $lineageDict

  """

}


process plot_mutation_profiles {
  label 'sortseq'
  publishDir "${params.output}/${params.eval_output}/${amplicon}", mode: 'copy', pattern: "*.png"

  input:
    tuple val(amplicon), val(ampliconsize)

  output:
    tuple val(amplicon), val(ampliconsize), emit:fwd
    path "${amplicon}_mutation_profile.png"

  script:
  """
  #!/bin/bash
  #python3 ${projectDir}/bin/mutationprofile.py ${launchDir}/${params.output}/${params.hdbscan_output}/variants/${amplicon} $amplicon
  python3 ${projectDir}/bin/mutationprofile.py ${params.output}/${params.hdbscan_output}/variants/${amplicon} $amplicon

  """

}
