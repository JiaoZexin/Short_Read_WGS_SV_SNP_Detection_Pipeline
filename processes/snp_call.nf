process snp_call {
    tag { sample_id }

    input:
    tuple val(sample_id), path(bam), path(ref_fa)

    output:
    path "10.snp/${sample_id}.snp.vcf"

    script:
    """
    mkdir -p 10.snp
    start_time=$(date +%s)
    echo "[`date`] SNP calling for ${sample_id}" > 10.snp/${sample_id}_snp.log

    bcftools mpileup -Ou -f ${ref_fa} ${bam} | \\
    bcftools call -mv -Ov -o 10.snp/${sample_id}.snp.vcf >> 10.snp/${sample_id}_snp.log 2>&1

    end_time=$(date +%s)
    echo "SNP calling completed in \$((end_time - start_time)) seconds" >> 10.snp/${sample_id}_snp.log
    """
}