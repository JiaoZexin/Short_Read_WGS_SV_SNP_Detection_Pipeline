process extract_samples_ids {
    label 'lowmem'
    publishDir "./Results/", mode: 'copy'
    input:
        path "0.data"
    
    output:
        path "sample.txt"

    script:
    """
    for fq in 0.data/*_R1.fastq.gz; do
        basename \$fq | sed 's/_R1\\.fastq\\.gz//'
    done | sort | uniq > sample.txt
    """
}