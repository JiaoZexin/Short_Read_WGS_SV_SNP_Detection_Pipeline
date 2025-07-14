process SNP_filtering_bcftools {
    tag "snp_filter"
    label 'standard'
    conda 'envs/vcftools_bedtools.yml'
    publishDir "./Results/", mode: 'copy'

    input:
    path merged_vcf

    output:
    path "10.bcftools_final.snp/snp.final.vcf.recode.vcf", emit: filtered_vcf

    script:
    """
    mkdir -p 10.bcftools_final.snp

    vcftools --vcf ${merged_vcf} \\
        --maf 0.05 --minDP 10 --maxDP 50 --minQ 30 \\
        --remove-indels --recode --recode-INFO-all \\
        --out 10.bcftools_final.snp/snp.final.vcf
    """
}
