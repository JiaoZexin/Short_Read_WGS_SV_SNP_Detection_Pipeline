process genotype_smoove {
    tag { sample_id }

    input:
    tuple val(sample_id), path(sorted_bam), path(ref_fa), path(sites_vcf)

    output:
    path "8.genotype/${sample_id}/*"

    script:
    """
    mkdir -p 8.genotype/${sample_id}
    start_time=\$(date +%s)
    echo "[`date`] Running smoove genotype for ${sample_id}" > 8.genotype/${sample_id}/${sample_id}_genotype.log

    smoove genotype -d -x -p 1 \\
        --name ${sample_id}-joint \\
        --outdir 8.genotype/${sample_id} \\
        --fasta ${ref_fa} \\
        --vcf ${sites_vcf} \\
        ${sorted_bam} >> 8.genotype/${sample_id}/${sample_id}_genotype.log 2>&1

    end_time=\$(date +%s)
    echo "smoove genotype finished in \$((end_time - start_time)) seconds" >> 8.genotype/${sample_id}/${sample_id}_genotype.log
    """
}
