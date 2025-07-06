process SV_visualize_Samplot {
    tag "$sv_id"
    publishDir "${params.outdir}/sv_images", mode: 'copy'

    input:
    tuple val(sv_id), path(vcf), path(bam), path(reference_genome)

    output:
    path "${sv_id}.png"

    script:
    """
    samplot plot \
        --vcf ${vcf} \
        --bam ${bam} \
        --region ${sv_id} \
        --ref ${reference_genome} \
        --out ${sv_id}.png
    """
}