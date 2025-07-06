process SNP_call_bcftools {
    tag { sample_id }
    label 'highmem'
    conda 'envs/vcftools_bedtools.yml' 
    publishDir "./Results/", mode: 'copy'

    input:
    tuple val(sample_id), path(bam_path), path(ref_index_path)

    output:
    path "7.bcftools/${sample_id}/*"

    script:
    """
    mkdir -p 7.bcftools/${sample_id}

    echo "[`date`] Running bcftools mpileup + call for ${sample_id}" > 7.bcftools/${sample_id}/${sample_id}_bcftools.log

    bcftools mpileup \\
        -f ${ref_index_path}/ref.fa \\
        ${bam_path}/${sample_id}.sorted.bam |
    bcftools call -c -v \\
        -o 7.bcftools/${sample_id}/${sample_id}.vcf 2>> 7.bcftools/${sample_id}/${sample_id}_bcftools.log

    echo "[`date`] bcftools finished for ${sample_id}" >> 7.bcftools/${sample_id}/${sample_id}_bcftools.log
    """
}
