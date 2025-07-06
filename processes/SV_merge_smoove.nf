process SV_merge_smoove {
    tag { smoove_merge }
    label 'highmem'
    conda 'envs/smoove.yml'
    publishDir "./Results/", mode: 'copy'

    input:
    path smoove_vcf_list

    output:
    path "8.smoove_merge/*"


    script:
    """
    mkdir -p 8.smoove_merge
    merge_name=\$(date +\"%Y%m%d_%H%M\")
    start_time=\$(date +%s)
    echo "[`date`] Running smoove merge" > 8.smoove_merge/merge.log

    smoove merge --name \$merge_name --fasta ${params.reference_fasta} --outdir 8.smoove_merge/ ${smoove_vcf_list.join(' ')} >> 8.smoove_merge/merge.log 2>&1

    end_time=\$(date +%s)
    echo "smoove merge finished in \$((end_time - start_time)) seconds" >> 8.smoove_merge/merge.log
    """
}
