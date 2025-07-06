// SV & SNP Analysis Pipeline using Nextflow
// Optional
// ========== include process module 和 workflow ==========

include { SV_visualize_Samplot } from './processes/SV_visualize_Samplot.nf'
include { Annotation_snpeff } from './processes/Annotation_snpeff.nf'
include { Annotation_annovar } from './processes/Annotation_annovar.nf'
include { GO_KEGG_preparion } from './processes/GO_KEGG_preparion.nf'
include { GWAS_analysis } from './processes/GWAS_analysis.nf'
include { Fimpute_imputation } from './processes/Fimpute_imputation.nf'


workflow {

    // 1. SV可视化输入通道
    sv_visualize_input_ch = Channel
        .fromPath("sv_region_list.txt")
        .splitText()
        .map { line ->
            def fields = line.tokenize('\t')
            tuple(fields[0], file(fields[1]), file(fields[2]), file(fields[3]))
        }
    if (params.use_sv_visualize) {
        SV_visualize_Samplot(sv_visualize_input_ch)
    }

    // 2. snpeff注释输入通道
    snpeff_input_ch = Channel.of(
        tuple("sv", file(params.SV_vcf), file("2.index/ref.gff3")),
        tuple("snp", file(params.snp_vcf), file("2.index/ref.gff3"))
    )
    if (params.use_snpeff) {
        Annotation_snpeff(snpeff_input_ch)
    }

    // 3. annovar注释输入通道
    annovar_input_ch = Channel.of(
        tuple("sv", file(params.SV_vcf)),
        tuple("snp", file(params.snp_vcf))
    )
    if (params.use_annovar) {
        Annotation_annovar(annovar_input_ch)
    }

    // 4. GO/KEGG注释输入通道
    go_kegg_annotation_input_ch = Channel.fromPath(params.pep_fasta)
    if (params.use_go_kegg) {
        GO_KEGG_preparion(go_kegg_annotation_input_ch)
    }

    // 5. GWAS分析输入通道
    gwas_input_ch = Channel.of(
        tuple("snp", file(params.snp_vcf), file(params.phenotype)),
        tuple("sv",  file(params.SV_vcf),file(params.phenotype))
    )
    if (params.use_gwas) {
        GWAS_analysis(gwas_input_ch)
    }

    // 6. FImpute控制文件输入通道
    fimpute_input_ch = Channel.fromPath(params.control_file)
    if (params.use_imputation) {
        Fimpute_imputation(fimpute_input_ch)
    }

}