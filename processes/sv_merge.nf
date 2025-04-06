process merge_smoove_vcfs {
    tag "smoove_merge"

    input:
    path ref_fa from file("2.index/ref.fa")
    path vcfs from smoove_vcf_ch

    output:
    path "7.merge/*"

    script:
    """
    mkdir -p 7.merge
    merge_name=\$(date +\"%Y%m%d_%H%M\")
    start_time=\$(date +%s)
    echo "[`date`] Running smoove merge" > 7.merge/merge.log

    vcf_list=\$(ls ${vcfs.join(' ')} | tr '\\n' ' ')

    smoove merge --name \$merge_name --fasta ${ref_fa} \$vcf_list >> 7.merge/merge.log 2>&1

    end_time=\$(date +%s)
    echo "smoove merge finished in \$((end_time - start_time)) seconds" >> 7.merge/merge.log
    """
}
