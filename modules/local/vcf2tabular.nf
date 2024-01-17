process VCF2TABULAR{
    container "< container >"

    input:
    tuple val(meta), path(vcf)

    output:
    path '*'      , emit: table
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    vcf2tabular \\
    -i $vcf \\
    -o ${meta.id}.xlsx \\
    -snpEff F \\
    -format GT,AD,AF,DP,F1R2,F2R1,SB \\
    -tmpDir .

    cat <<-END_VERSIONS > versions.yml
        vcf2tabular version 1.0.1
    END_VERSIONS
    """
}
