process build_reference_index {
    input:
    path ref_fa from file("0.data/ref.fa")

    output:
    path "2.index/*"

    script:
    """
    mkdir -p 2.index
    start_time=\$(date +%s)

    cp ${ref_fa} 2.index/ref.fa
    cd 2.index

    {
        echo "[`date`] Running: bwa index -a is ref.fa"
        bwa index -a is ref.fa
        echo "[`date`] Running: samtools faidx ref.fa"
        samtools faidx ref.fa
    } &> index.log

    end_time=\$(date +%s)
    echo "Indexing completed in \$((end_time - start_time)) seconds." >> index.log
    """
}
