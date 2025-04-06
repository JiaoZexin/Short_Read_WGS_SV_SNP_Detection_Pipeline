process snp_merge_filter {
    tag "snp_merge_filter"

    input:
    path vcfs from file("10.snp/*.snp.vcf")

    output:
    path "10.snp/snp.final.vcf.recode.vcf"

    script:
    """
    mkdir -p 10.snp

    for vcf in ${vcfs}; do
      bgzip -c -f \$vcf > \$vcf.gz
      bcftools index \$vcf.gz
    done

    bcftools merge --force-samples -m id -O v -o 10.snp/snp.all.vcf \\
      ${vcfs.collect{ it + '.gz' }.join(' ')}

    vcftools --vcf 10.snp/snp.all.vcf \\
        --maf 0.05 --minDP 10 --maxDP 50 --minQ 30 \\
        --remove-indels --recode --recode-INFO-all \\
        --out 10.snp/snp.final.vcf
    """
}