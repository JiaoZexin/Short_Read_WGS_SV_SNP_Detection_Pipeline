process SNP_call_deepvariant {
    tag { sample_id }
    label 'bigmem'
    container 'google/deepvariant:1.6.0'
    publishDir "./Results/", mode: 'copy'

    input:
    tuple val(sample_id), path(bam_path), path(ref_index_path)

    output:
    path "7.deepvariant/${sample_id}/*"

    script:
    """
    mkdir -p 7.deepvariant/${sample_id}
    
    echo "[`date`] Running DeepVariant for ${sample_id}" > 7.deepvariant/${sample_id}/${sample_id}_deepvariant.log

    /opt/deepvariant/bin/run_deepvariant \\
        --model_type=WGS \\
        --ref=${ref_index_path}/ref.fa \\
        --reads=${bam_path}/${sample_id}.sorted.bam \\
        --output_vcf=7.deepvariant/${sample_id}/${sample_id}.vcf.gz \\
        --output_gvcf=7.deepvariant/${sample_id}/${sample_id}.g.vcf.gz \\
        --num_shards=16 >> 7.deepvariant/${sample_id}/${sample_id}_deepvariant.log 2>&1

    echo "[`date`] DeepVariant finished for ${sample_id}" >> 7.deepvariant/${sample_id}/${sample_id}_deepvariant.log
    """
}
