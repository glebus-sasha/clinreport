// DSL2 module: include { CLINREPORT } from './examples/clinreport'
// Pass upstream paths through channels, not published output directories.
process CLINREPORT {
    container params.clinreport_container
    publishDir "${params.outdir}/clinreport", mode: 'copy'

    input:
    path drug_network, stageAs: 'inputs/drugs.csv'
    path genes, stageAs: 'inputs/genes.tsv'
    path network, stageAs: 'inputs/network.tsv'
    path tx2gene, stageAs: 'inputs/tx2gene.tsv'
    path gsea_case, stageAs: 'inputs/gsea_case.tsv'
    path gsea_control, stageAs: 'inputs/gsea_control.tsv'
    path gmt, stageAs: 'inputs/pathways.gmt'

    output:
    path 'report', emit: app

    script:
    """
    clinreport --prepare report \\
      --drug-network-file inputs/drugs.csv \\
      --gene-file inputs/genes.tsv \\
      --interaction-network-file inputs/network.tsv \\
      --tx2gene-file inputs/tx2gene.tsv \\
      --gsea-carcinoma-report-file inputs/gsea_case.tsv \\
      --gsea-normal-report-file inputs/gsea_control.tsv \\
      --pathway-gmt-file inputs/pathways.gmt
    """
}
