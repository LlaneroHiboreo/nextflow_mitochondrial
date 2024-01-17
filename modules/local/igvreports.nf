process IGVREPORTS{
    container "/scratch_isilon/groups/dat/jgaztelu/sifffs/igvreports.sif"

    input:
    tuple val(meta), path(vcf)
    tuple val(meta), path(fasta)
    path(ideogram)

    output:
    path '*.html'      , emit: table
    
    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    create_report \\
    --ideogram $ideogram \\
    --info-columns MitoMap_Disease MitoMap_DiseaseStatus MitoMap_AF gnomAD_WG_AF_het gnomAD_WG_AF_hom \\
    --output ${prefix}.html \\
    --fasta $fasta \\
    $vcf
    """
}