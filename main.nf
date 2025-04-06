// SV & SNP Analysis Pipeline using Nextflow

include { init_setup } from './processes/init_setup.nf'
include { build_reference_index } from './processes/build_index.nf'
include { extract_sample_ids } from './processes/extract_samples.nf'
include { align_bwa } from './processes/alignment.nf'
include { sam_to_bam } from './processes/sam_to_bam.nf'
include { sort_and_index_bam } from './processes/sort_index_bam.nf'
include { call_smoove } from './processes/sv_calling.nf'
include { merge_smoove_vcfs } from './processes/sv_merge.nf'
include { genotype_smoove } from './processes/genotyping.nf'
include { duphold_filter } from './processes/duphold_filter.nf'
include { sv_count_summary } from './processes/sv_count_summary.nf'
include { sv_paste_merge } from './processes/sv_paste_merge.nf'
include { final_filter } from './processes/final_filter.nf'
include { snp_call } from './processes/snp_call.nf'
include { snpeff_annotate } from './processes/snpeff_annotate.nf'
include { go_kegg_annotation } from './processes/go_kegg_annotation.nf'
include { snp_merge_filter } from './processes/snp_merge_filter.nf'
include { fimpute_imputation } from './processes/fimpute_imputation.nf'
include { gwas_analysis } from './processes/gwas_analysis.nf'

include { fastq_pairs; ref_fasta_ch； sam_files_ch; bam_files_ch; smoove_input_ch; smoove_vcf_ch; sites_vcf; genotype_input_ch; duphold_input_ch; sv_count_input_ch; snp_call_input_ch; snpeff_input_ch; kegg_go_input_ch; fimpute_input_ch } from './channels.nf'

workflow {
    init_setup()
    build_reference_index()
    extract_sample_ids()
    align_bwa(fastq_pairs, ref_fasta_ch)
    sam_to_bam(sam_files_ch)
    sort_and_index_bam(bam_files_ch)
    call_smoove(smoove_input_ch)
    merge_smoove_vcfs(ref_fasta_ch, smoove_vcf_ch)
    genotype_smoove(genotype_input_ch)
    duphold_filter(duphold_input_ch)
    sv_count_summary(sv_count_input_ch)
    sv_paste_merge()
    final_filter()
    snp_call(snp_call_input_ch)
    snp_merge_filter()
    snpeff_annotate(snpeff_input_ch)
    go_kegg_annotation(kegg_go_input_ch)
    fimpute_imputation(fimpute_input_ch)
    gwas_analysis()
}
