/************************************************************************
* HDBSCAN
*
* Utilizing UMAP and HDBSCAN to cluster sequences based on their
* k-mer vector representation and the cosine distance
************************************************************************/

process hdbscan {
  label 'hdbscan'
  publishDir "${params.output}/${params.hdbscan_output}", mode: 'copy', pattern: "*_hdbscan.fasta"
  publishDir "${params.output}/${params.hdbscan_output}", mode: 'copy', pattern: "*.clstr"
  publishDir "${params.output}/${params.hdbscan_output}", mode: 'copy', pattern: "*log"
  publishDir "${params.output}/${params.hdbscan_output}", mode: 'copy', pattern: "*_hdbscan_UNCLUSTERED.fasta"
  publishDir "${params.output}/${params.summary_output}/unclustered_sequences", mode: 'copy', pattern: '*UNCLUSTERED.fasta'
  publishDir "${params.output}/${params.summary_output}/clustered_sequences", mode: 'copy', pattern: '*_hdbscan.fasta'
  publishDir "${params.output}/${params.hdbscan_output}", mode: 'copy', pattern: '*.png'
  publishDir "${params.output}/${params.hdbscan_output}", mode: 'copy', pattern: 'cluster*.fasta'


  input:
    path(sequences)
    path(lineageDict)
    val(addParams)
    val(goi)

  output:
    tuple val("${params.output}/${params.hdbscan_output}"), path("${sequences.baseName}_hdbscan.fasta"), path("${sequences.baseName}_hdbscan.fasta.clstr"), emit: hdbscan_out
    path "${sequences.baseName}_hdbscan_UNCLUSTERED.fasta", emit: hdbscan_unclustered
    path "cluster*.fasta", emit: cluster_fasta_files
    path "*.png"
    path "hdbscan.log"

  script:
  def GOI = goi != 'NO FILE' ? "${goi}" : ''
  """
    python3 ${baseDir}/bin/hdbscan_virus.py -v -p ${task.cpus} ${addParams} ${sequences} ${lineageDict} ${GOI} 2> hdbscan.log
    ls
    touch cluster-1.fasta
    less cluster-1.fasta
    mv cluster-1.fasta  "${sequences.baseName}_hdbscan_UNCLUSTERED.fasta"

  """

}
