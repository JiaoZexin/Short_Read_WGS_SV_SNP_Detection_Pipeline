process sv_paste_merge {
    tag "sv_paste"

    input:
    path filtered_vcfs from file("8.duphold/*-filter-name-joint-smoove.genotyped.vcf.gz")

    output:
    path "8.duphold/*.smoove.square.vcf.gz"

    script:
    """
    cd 8.duphold
    merge_name=$(date +"%Y%m%d_%H%M")
    start_time=$(date +%s)
    echo "[`date`] Running smoove paste" > smoove_paste.log

    smoove paste --name $merge_name ./*-filter-name-joint-smoove.genotyped.vcf.gz >> smoove_paste.log 2>&1

    end_time=$(date +%s)
    echo "smoove paste finished in \$((end_time - start_time)) seconds" >> smoove_paste.log
    """
}
