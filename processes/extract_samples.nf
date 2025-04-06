process extract_sample_ids {
    output:
    path "sample.txt"

    script:
    """
    mkdir -p 1.fq.gz
    for fq in 0.data/*_R1.fastq.gz; do
        basename \$fq | sed 's/_R1\\.fastq\\.gz//'
    done | sort | uniq > sample.txt
    """
}
