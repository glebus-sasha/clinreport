// Full ClinReport module. Keep annotations as one directory; annotations/*
// flattens the source hierarchy during staging.
// Pass upstream paths through channels, not published output directories.
process CLINREPORT {
    container params.clinreport_container
    publishDir "${params.outdir}/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path drug_network, path annotations, path genes,
          path gsea_case, path gsea_control, path network, path tx2gene, path gmt
    path wes, optional: true
    path tf, optional: true
    path dorothea, optional: true
    val log2fc_cutoff
    val padj_cutoff

    output:
    tuple val(meta), path('report'), emit: app

    script:
    def wes_args = wes ? "--wes-vcf-file '${wes}'" : ''
    def tf_args = tf ? "--tf-file '${tf}'" : ''
    def dorothea_args = dorothea ? "--dorothea-file '${dorothea}'" : ''
    def use_tf_activity = params.use_tf_activity ?: 'FALSE'
    def use_wes = params.use_wes ?: 'FALSE'
    def skip_network_processing = params.skip_network_processing ?: 'FALSE'
    """
    clinreport --prepare report \\
      --patient-id '${meta.id}' \\
      --drug-network-file '${drug_network}' \\
      --clinreport-dir '${annotations}' \\
      --gene-file '${genes}' \\
      ${wes_args} \\
      --gsea-carcinoma-report-file '${gsea_case}' \\
      --gsea-normal-report-file '${gsea_control}' \\
      ${tf_args} \\
      --interaction-network-file '${network}' \\
      --tx2gene-file '${tx2gene}' \\
      --pathway-gmt-file '${gmt}' \\
      ${dorothea_args} \\
      --network-display-name 'Drug-gene interaction network' \\
      --interaction-network-name '${meta.network_name ?: 'STRING'}' \\
      --pathway-collection-name 'Hallmark' \\
      --pathway-id-prefix 'HALLMARK_' \\
      --padj-cutoff ${padj_cutoff} \\
      --log2fc-cutoff ${log2fc_cutoff} \\
      --use-tf-activity '${use_tf_activity}' \\
      --use-wes '${use_wes}' \\
      --skip-network-processing '${skip_network_processing}'
    """
}
