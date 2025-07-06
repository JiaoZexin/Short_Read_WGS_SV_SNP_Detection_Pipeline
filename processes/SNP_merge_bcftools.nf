process SNP_merge_bcftools {
    tag "snp_merge"
    label 'standard'
    conda 'envs/vcftools_bedtools.yml'
    publishDir "./Results/", mode: 'copy'

    input:
    path vcfs

    output:
    path "10.snp/snp.all.vcf", emit: merged_vcf

    script:
    """
    mkdir -p 10.snp

    # bgzip and index all input vcfs
    for vcf in ${vcfs}; do
      bgzip -c -f \$vcf > \$vcf.gz
      bcftools index \$vcf.gz
    done

    # build list of bgzipped files
    gz_vcfs=""
    for vcf in ${vcfs}; do
      gz_vcfs="\$gz_vcfs \$vcf.gz"
    done

    bcftools merge --force-samples -m id -O v -o 10.snp/snp.all.vcf \$gz_vcfs
    """
}
