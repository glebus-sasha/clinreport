// The GMT list CSV names the collections; ClinReport validates all supplied
// GMT and GSEA files, regardless of how many collections the workflow emits.
// If annotations includes gene_annotations.json, ClinReport discovers it in
// --clinreport-dir automatically; no separate Nextflow input is required.
process CLINREP {
    tag "${meta.id}"
    container 'quay.io/glebussasha/clinreport:latest'

    input:
    tuple val(meta), path(samplesheet)
    tuple val(meta2), path(drug_network)
    path annotations
    path genes
    path network
    path tx2gene
    path gmt_list_file
    path gmts
    path gsea_files
    path wes
    path tf
    path dorothea
    val log2fc_cutoff
    val padj_cutoff

    output:
    tuple val(meta), path('report'), emit: app

    script:
    def wes_args = wes ? "--wes-vcf-file '${wes}' --wes-display-name 'WES · PCGR'" : ''
    def tf_args = tf ? "--tf-file '${tf}'" : ''
    def dorothea_args = dorothea ? "--dorothea-file '${dorothea}'" : ''
    def use_tf_activity = params.use_tf_activity ?: 'FALSE'
    def use_wes = params.use_wes ?: 'FALSE'
    def skip_network_processing = params.skip_network_processing ?: 'FALSE'
    def gsea_args = gsea_files.collect { "'${it}'" }.join(' ')
    def gmt_args = gmts.collect { "'${it}'" }.join(' ')
    """
    mkdir -p annotations
    for annotation in ${annotations}; do
        ln -s "\$annotation" annotations/
    done

    clinreport --prepare report \\
      --patient-id '${meta.id}' \\
      --drug-network-file '${drug_network}' \\
      --clinreport-dir annotations \\
      --gene-file '${genes}' \\
      ${wes_args} \\
      --gsea-files ${gsea_args} \\
      --gmt-files ${gmt_args} \\
      --gmt-list-file '${gmt_list_file}' \\
      ${tf_args} \\
      --interaction-network-file '${network}' \\
      --tx2gene-file '${tx2gene}' \\
      ${dorothea_args} \\
      --network-display-name 'Drug-gene interaction network' \\
      --interaction-network-name '${meta2.network_id}' \\
      --padj-cutoff ${padj_cutoff} \\
      --log2fc-cutoff ${log2fc_cutoff} \\
      --use-tf-activity '${use_tf_activity}' \\
      --use-wes '${use_wes}' \\
      --skip-network-processing '${skip_network_processing}'
    """
}
