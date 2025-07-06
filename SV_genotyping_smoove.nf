// SV & SNP Analysis Pipeline using Nextflow
// Smoove merge
// ========== include process module 和 workflow ==========


include { SV_merge_smoove } from './processes/SV_merge_smoove.nf'
include { SV_merge_smoove_genotype } from './processes/SV_merge_smoove_genotype.nf'

workflow {
    smoove_vcf_ch = Channel.fromPath("${baseDir}/Results/6.smoove/*/*.vcf.gz")
        .collect()
    
    // SV_merge_smoove(smoove_vcf_ch)
    SV_merge_smoove_output_ch = SV_merge_smoove(smoove_vcf_ch)
    SV_merge_smoove_sites_vcf_ch = SV_merge_smoove_output_ch
        .flatten()
        .filter { it.name.endsWith('.vcf.gz') }

    // 8. Genotype输入：排序BAM、参考、sites.vcf
    SV_smoove_genoty_input_ch = Channel.fromPath(params.sample_list)
        .splitText()
        .combine(SV_merge_smoove_sites_vcf_ch)
        .map { id, sites_vcf ->
            id = id.trim()  
            def bam_path = file("${baseDir}/Results/5.sorted/${id}/")
            def ref_path = file(params.reference_path)
            tuple(id, bam_path, ref_path, sites_vcf)
        }

    SV_merge_smoove_genotype(SV_smoove_genoty_input_ch)



}