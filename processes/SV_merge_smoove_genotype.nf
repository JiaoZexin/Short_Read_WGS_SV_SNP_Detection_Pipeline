process SV_merge_smoove_genotype {
    tag { sample_id }
    label 'highmem'
    conda 'envs/smoove.yml'
    publishDir "./Results/", mode: 'copy'

    input:
    tuple val(sample_id), path(bam_path), path(ref_path), path(sites_vcf)

    output:
    path "8.genotype/${sample_id}/*"

    script:
    """
    mkdir -p 8.genotype/${sample_id}
    start_time=\$(date +%s)
    echo "[`date`] Running smoove genotype for ${sample_id}" > 8.genotype/${sample_id}/${sample_id}_genotype.log

    smoove genotype -d -x -p 1 \\
        --name ${sample_id}-joint \\
        --outdir 8.genotype/${sample_id} \\
        --fasta ${ref_path}/ref.fa \\
        --vcf ${sites_vcf} \\
        ${bam_path}/${sample_id}.sorted.bam 2>> 8.genotype/${sample_id}/${sample_id}_genotype.log 

    end_time=\$(date +%s)
    echo "smoove genotype finished in \$((end_time - start_time)) seconds" >> 8.genotype/${sample_id}/${sample_id}_genotype.log
    """
}

