/*******************************************************************************
Collect reads from left and right primers into a amplicon fasta file.
Note:
- there was the idea to track primer annotation for each read so that in 
downstream analysis, you can call variants for both primer pools before
merging them into amplicon variants
/******************************************************************************/

process primerset_to_amplicon {
  publishDir "${params.output}/${params.amplicon_data}/final", mode: 'copy', pattern: "${amplicon}.fasta"


  input:
    tuple val(amplicon), path(primer_fasta_x), path(primer_fasta_y)

  output:
    tuple val(amplicon), path("${amplicon}.fasta")


  script:
  if ("${primer_fasta_y}" == "dummy")
    """
    #!/bin/bash
    echo "------------------- Amplicon ${amplicon} -------------------"
    cat $primer_fasta_x > ${amplicon}.fasta
    readcount=\$(grep -c '>' $amplicon}.fasta)
    echo "Validate, the amplicon contains \${readcount} reads."

    """
  else
    """    
    #!/bin/bash
    echo "------------------- Amplicon ${amplicon} -------------------"
    cat $primer_fasta_x > ${amplicon}.fasta
    cat $primer_fasta_y >> ${amplicon}.fasta
    readcount=\$(grep -c '>' $amplicon}.fasta)
    echo "Validate, the amplicon contains \${readcount} reads."

    """
}
