/************************************************************************
Create dictionary mapping read names to lineage annotation.
Used later on to evaluate clustering.
************************************************************************/

process create_lineage_dict {
  label 'sortseq'
 
  input:
    tuple val(name), path(fasta), val(lineage)

  output:
    path "lineageDict.csv" 

  script:
  """
    #!/bin/bash

    python ${baseDir}/bin/create_lineage_dict.py $fasta $lineage 

  """
}

/************************************************************************
Subsample for mixture sample
************************************************************************/

process create_mixtures {
  label 'sortseq'
 
  input:
    tuple val(name), path(fasta), val(n_reads)

  output:
    tuple val(name), path("sub_${fasta}"), emit: component
    path ".command.log", emit:log

  script:
  if ("${n_reads}" == params.sample_size)
  """
  #!/bin/bash
    cp $fasta sub_$fasta
  """
  else
  """
    #!/bin/bash
    echo "sampling $n_reads reads for $name"
    seqtk sample $fasta $n_reads > sub_$fasta
  """
}
