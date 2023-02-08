/************************************************************************
Represent amplicon reads by their kmer frequencies.
Reduce dimensionality of the amplicon read set for clustering. 
Then apply HDBSCAN.
Notes:
- There was the idea to provide PCA as an alternative for UMAP
- Bootstrap?
************************************************************************/

process hdbscan {
  label 'hdbscan'
  publishDir "${params.output}/${params.eval_output}/${amplicon_id}", mode: 'copy', pattern: '*.png'
  publishDir "${params.output}/${params.hdbscan_output}/${amplicon_id}", mode: 'copy', pattern: 'cluster*.fasta'


  input:
    tuple val(amplicon_id), path(amplicon_fasta), path(lineageDict)
    //path(primerDict)

  output:
    tuple val(amplicon_id), path("cluster*.fasta"), emit: amplicon_cluster, optional: true
    path "*.png", optional: true
    path ".command.log", emit: log

  script:
    """
      #!/bin/bash

      echo "------------------- Amplicon ${amplicon_id} -------------------"
      echo \$PWD

      python3 ${baseDir}/bin/hdbscan_virus.py -v -p $task.cpus $params.umap_params $params.hdbscan_params ${amplicon_id}.fasta $lineageDict 2> hdbscan.log
      cat hdbscan.log >> .command.log
      
      if [ -f cluster-1.fasta  ]; then
        touch cluster-1.fasta
        mv cluster-1.fasta  "${amplicon_id}_hdbscan_UNCLUSTERED.fasta" 
      fi 

      if [ ! -f cluster0.fasta ]; then
        echo "Only noise detected"
      fi
      

    """
}
