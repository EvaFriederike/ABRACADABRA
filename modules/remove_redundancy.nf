/************************************************************************
Remove redundant reads and memorize them to re-introduce them after 
embedding and clustering.
************************************************************************/

process remove_redundancy {
  label 'remove'

  input:
    tuple val(primer_name), path(primer_fasta)

  output:
    tuple val(primer_name), path("${primer_name}_nr.fasta"), emit: nr_fasta
    path '*-duplicate-seqs.fasta', emit: redundancy, optional: true
    path ".command.out", emit: log

  script:
  """
    #!/bin/bash
    amplicon=\$(echo $primer_name | cut -d_ -f2)
    seqkit rmdup -s -P -d representatives-of-redundancy.fasta -D duplicate-seqs.txt -o "${primer_name}_nr.fasta" ${primer_fasta}

    if [ -f "duplicate-seqs.txt" ]; then
      python ${baseDir}/bin/write_redundancy.py $primer_fasta duplicate-seqs.txt \${amplicon}-duplicate-seqs.fasta
    else
      echo "No redundant reads in ${primer_name}"
    fi
  """
}


/************************************************************************
Re-introduce for lineage detection and abundance estimation.
************************************************************************/
process add_redundancy {
  label 'remove'

  input:
    tuple val(amplicon), val(amplicon_size), path(cluster_fasta), path(duplicates_fasta)

  output:
    tuple val(amplicon), val(amplicon_size), val("redundant_${cluster_fasta.baseName}"), path("redundant_${cluster_fasta}")
    
  script:
  """
    #!/bin/bash
    
    grep '>' $cluster_fasta > tmp.fasta
    sed -i 's/>//g' tmp.fasta
    cp $cluster_fasta redundant_${cluster_fasta}
    seqkit grep -f tmp.fasta $duplicates_fasta > redundancy.fasta

    echo \$(grep -c '>' redundancy.fasta)
    
    cat redundancy.fasta >> redundant_${cluster_fasta}
  """
}
