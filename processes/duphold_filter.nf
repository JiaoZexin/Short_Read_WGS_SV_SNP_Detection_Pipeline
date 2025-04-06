process duphold_filter {
    tag { sample_id }

    input:
    tuple val(sample_id), path(vcf_file)

    output:
    path "8.duphold/${sample_id}*"

    script:
    """
    mkdir -p 8.duphold
    start_time=$(date +%s)
    echo "[`date`] Running duphold filter for ${sample_id}" > 8.duphold/${sample_id}_duphold_filter.log

    cp ${vcf_file} 8.duphold/
    cd 8.duphold
    gunzip ${sample_id}-joint-smoove.genotyped.vcf.gz

    bcftools filter -S . \
        -e '(SVTYPE = "DEL" & FMT/DHFFC[0] > 0.7 & FMT/GT[0] != "0/0") | \
            (SVTYPE = "DUP" & FMT/DHFFC[0] < 1.3 & FMT/GT[0] != "0/0") | \
            (SVTYPE = "INV" & FMT/DHFFC[0] > 1.3 & FMT/GT[0] != "0/0") | \
            (SVTYPE = "INV" & FMT/DHFFC[0] < 0.7 & FMT/GT[0] != "0/0")' \
        ${sample_id}-joint-smoove.genotyped.vcf \
        > ${sample_id}-filter-name-joint-smoove.genotyped.vcf

    # 删除原始未过滤的 VCF（保留一份用于统计）
    cp ${sample_id}-joint-smoove.genotyped.vcf original_summary_${sample_id}.vcf
    bgzip original_summary_${sample_id}.vcf
    bcftools index original_summary_${sample_id}.vcf.gz
    rm ${sample_id}-joint-smoove.genotyped.vcf

    # 对过滤后的 VCF 压缩并建立索引
    bgzip ${sample_id}-filter-name-joint-smoove.genotyped.vcf
    bcftools index ${sample_id}-filter-name-joint-smoove.genotyped.vcf.gz

    end_time=$(date +%s)
    echo "duphold filter finished in $((end_time - start_time)) seconds" >> 8.duphold/${sample_id}_duphold_filter.log
    """
}