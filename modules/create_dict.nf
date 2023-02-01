/************************************************************************
Create dictoinary mapping read names to lineage annotation.
Used later on to evaluate clustering.
************************************************************************/

process create_lineage_dict {
  label 'sortseq'
 
  input:
    tuple val(barcode), path(fasta), val(lineage)

  output:
    path "lineageDict.csv" 

  script:
  """
    #!/bin/bash

    python ${baseDir}/bin/create_lineage_dict.py $fasta $lineage 

  """
}