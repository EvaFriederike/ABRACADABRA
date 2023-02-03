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
  publishDir "${params.output}/${params.mode}/${params.eval_output}/${amplicon_id}", mode: 'copy', pattern: '*.png'
  publishDir "${params.output}/${params.mode}/${params.hdbscan_output}/${amplicon_id}", mode: 'copy', pattern: 'cluster*.fasta'


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

      touch dummy.json 
      echo "------------------- Amplicon ${amplicon_id} -------------------"

      python3 ${baseDir}/bin/hdbscan_virus.py -v -p $task.cpus $params.hdbscan_params ${amplicon_id}.fasta $lineageDict dummy.json 2> hdbscan.log
      cat hdbscan.log >> .command.log
      
      if [ -f cluster*.fasta ]; then
        touch cluster-1.fasta
        mv cluster-1.fasta  "${amplicon_id}_hdbscan_UNCLUSTERED.fasta"  
      else
        echo "No read clusters will be included into lineage detection (due to only noise or too few reads in the cluster)"
      fi
      

    """
}
