// SV & SNP Analysis Pipeline using Nextflow
// Smoove filtering
// ========== include process module 和 workflow ==========


include { SV_filtering_duphold_smoove } from './processes/SV_filtering_duphold_smoove.nf'
include { SV_filtering_paste_merge_smoove } from './processes/SV_filtering_paste_merge_smoove.nf'
include { SV_filtering_bioinformatic_smoove } from './processes/SV_filtering_bioinformatic_smoove.nf'

workflow {

    // Step 1: 逐样本输入 VCF，进入 duphold filtering
    duphold_input_ch = Channel
        .fromPath("${baseDir}/Results/8.genotype/*/*.vcf.gz")
        .map { vcf ->
            def sample = vcf.getBaseName().replaceAll("-joint-smoove\\.genotyped", "")
            tuple(sample, vcf)
        }

    // Step 1 output → Step 2 input
    sv_filtered_vcf_ch = SV_filtering_duphold_smoove(duphold_input_ch)

    // Step 2: 合并过滤后的 VCFs
    sv_merged_vcf_ch = SV_filtering_paste_merge_smoove(sv_filtered_vcf_ch)

    // Step 3: 最终过滤
    SV_filtering_bioinformatic_smoove(sv_merged_vcf_ch)

}