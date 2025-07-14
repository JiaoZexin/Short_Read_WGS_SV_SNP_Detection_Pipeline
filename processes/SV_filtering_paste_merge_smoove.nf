process SV_filtering_paste_merge_smoove {
    tag "sv_paste"
    label 'standard'
    conda 'envs/smoove.yml'
    publishDir "./Results/", mode: 'copy'

    input:
    path filtered_vcfs

    output:
    path "8.duphold/*.smoove.square.vcf.gz"

    script:
    """
    mkdir -p 8.duphold
    cd 8.duphold
    merge_name=\$(date +"%Y%m%d_%H%M")
    start_time=\$(date +%s)
    echo "[`date`] Running smoove paste" > smoove_paste.log

    smoove paste --name \$merge_name ${filtered_vcfs} >> smoove_paste.log 2>&1

    end_time=\$(date +%s)
    echo "smoove paste finished in \$((end_time - start_time)) seconds" >> smoove_paste.log
    """
}
