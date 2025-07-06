// SV & SNP Analysis Pipeline using Nextflow
// Detetion
// ========== include process module 和 workflow ==========


include { extract_samples_ids } from './processes/extract_samples_ids.nf'
include { fastqc } from './processes/fastqc.nf'
include { trimmomatic } from './processes/trimmomatic.nf'
include { build_reference_index } from './processes/build_reference_index.nf'
include { alignment } from './processes/alignment.nf'
include { SV_call_smoove } from './processes/SV_call_smoove.nf'
include { SV_call_delly } from './processes/SV_call_delly.nf'
include { SV_call_cnvpytor } from './processes/SV_call_cnvpytor.nf'
include { SNP_call_bcftools } from './processes/SNP_call_bcftools.nf'
include { SNP_call_deepvariant } from './processes/SNP_call_deepvariant.nf'
include { SNP_call_freebayes } from './processes/SNP_call_freebayes.nf'

workflow {
    
// ========== Workflow Execution & Channel ==========
    // Channel: extract_samples_ids
    extract_samples_ids_ch = Channel.fromPath("0.data")

    // NF: extract_samples_ids & Channel:sample.txt
    sample_txt_ch = extract_samples_ids(extract_samples_ids_ch)

    // Channel: fastq_pairs_ch
    fastq_pairs_ch = sample_txt_ch
        .splitText()
        .map { id ->
            id = id.trim()  
            def r1 = file("0.data/${id}_R1.fastq.gz")
            def r2 = file("0.data/${id}_R2.fastq.gz")
            tuple(id, r1, r2)
        }
    // fastq_pairs_ch.view()
    // NF fastqc
    // 重要 之后去掉注释
    fastq_pairs_ch.map { id, r1, r2 -> tuple(id, [r1, r2]) } | fastqc

    // NF trimmomatic
    trimmed_fastq_ch = trimmomatic(fastq_pairs_ch)

    // Channel: build_reference_index_input_ch
    build_reference_index_input_ch = Channel.of(file("0.data/ref.fa"))

    // NF: build_reference_index & Channel:ref_index_ch
    alignment_ref_index_ch = build_reference_index(build_reference_index_input_ch)
    // alignment_ref_index_ch.view()

    // NF: alignment & Channel: 3.sam/*.sam
    //alignment_input_ch = fastq_pairs_ch.combine(alignment_ref_index_ch)
    alignment_input_ch = trimmed_fastq_ch.combine(alignment_ref_index_ch)
    // alignment_input_ch.view()
    alignment_output = alignment(alignment_input_ch)


    // sorted_bam_ch = alignment_output.sorted_bam
    // bam_index_ch  = alignment_output.bam_index
    bam_path_ch = alignment_output.bam_path


    // // Step 3: 构建 smoove 输入
    // call_smoove_input_ch = sorted_bam_ch.combine(bam_index_ch)
    //     .map { bam, bai ->
    //         def sample_id = bam.getBaseName().replace('.sorted', '')
    //         tuple(sample_id, bam, bai)
    //     }

    // call_smoove_input_ch.view()
    // call_smoove_input_2_ch = call_smoove_input_ch.combine(alignment_ref_index_ch)  
    // id, bam path, ref path  
    call_input_ch = bam_path_ch.combine(alignment_ref_index_ch) 

    // Step 4: Smoove SV calling
    // 重要 之后去掉注释
    SV_call_smoove(call_input_ch)

    // Step 5: Delly SV calling
    if (params.use_delly) {
        SV_call_delly(call_input_ch)
    }
    // Step 6: CNVpytor SV calling
    if (params.use_cnvpytor) {
        SV_call_cnvpytor(call_input_ch)
    }

    // Step 7: SNP calling with BCFtools
    if (params.use_bcftools) {
        SNP_call_bcftools(call_input_ch)
    }

    // Step 8: SNP calling with DeepVariant
    if (params.use_deepvariant) {
        SNP_call_deepvariant(call_input_ch)
    }

    // Step 9: SNP calling with Freebayes
    if (params.use_freebayes) {
        SNP_call_freebayes(call_input_ch)
    }

// ========== Parameters ==========
    if (params.help) {
        log.info """
        SV/SNP Detection Pipeline 
        Please provide the necessary parameters to run the pipeline.

        Pipeline parameters:
        --adapters            Path to trim_adaptors (default: ${params.adapters})
        --sample_list         Path to sample list (default: ${params.sample_list})
        --reference_path      Path to reference index directory (default: ${params.reference_path})
        --reference_fasta     Path to reference fasta (default: ${params.reference_fasta})
        --reference_gff3      Path to reference gff3 (default: ${params.reference_gff3})
        --pep_fasta           Protein fasta (default: ${params.pep_fasta})
        --output_dir          Output directory (default: ${params.output_dir})
        --snp_vcf             SNP VCF file (default: ${params.snp_vcf})
        --SV_vcf              SV VCF file (default: ${params.SV_vcf})
        --control_file        FImpute control file (default: ${params.control_file})
        --use_delly           Use Delly for SV calling (default: false)
        --use_cnvpytor        Use CNVpytor for SV calling (default: false)
        --use_bcftools        Use BCFtools for SNP calling (default: false)
        --use_deepvariant     Use DeepVariant for SNP calling (default: false)
        --use_freebayes       Use Freebayes for SNP calling (default: false)


        Example usage:
        nextflow run main.nf -profile docker_conda --use_delly=true --use_cnvpytor=true --use_bcftools=true --use_deepvariant=true --use_freebayes=true 

        To show this help message:
        nextflow run main.nf --help
        """
        exit 0
    }

}
