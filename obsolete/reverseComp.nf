/************************************************************************
* Reverse Comp
*
* Generate reverse complementary sequences based on input records.
************************************************************************/

process reverseComp {
  label 'revComp'
  publishDir "${path}", mode: 'copy', pattern: '*.fasta'


  input:
    tuple val(name), val(path), path(sequences), path(cluster)


  output:
    tuple val(name), path(outpath)

  script:
  outpath = "${sequences}".reverse().replaceFirst("positive".reverse(),"negative".reverse()).reverse()
  """
  #!/bin/bash
  python3 ${baseDir}/bin/reverse_complement.py "${sequences}"

  """


}