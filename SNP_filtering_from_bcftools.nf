// SV & SNP Analysis Pipeline using Nextflow
// SNP filtering
// ========== include process module 和 workflow ==========


include { SNP_merge_bcftools } from './processes/SNP_merge_bcftools.nf'
include { SNP_filtering_bcftools } from './processes/SNP_filtering_bcftools.nf'

workflow {

    snp_merge_input_ch = Channel.fromPath("${baseDir}/Results/7.bcftools/*/*.vcf").collect()
    merged_vcf_ch = SNP_merge_bcftools(snp_merge_input_ch).merged_vcf
    SNP_filtering_bcftools(merged_vcf_ch)

}