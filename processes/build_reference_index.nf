process build_reference_index {
    label 'bigmem'
    conda 'envs/samtools.yml'
    publishDir "./Results/", mode: 'copy'
    input:
    path ref_fa 

    output:
    path "2.index/"


    script:
    """
    # 调试信息：检查 Conda 环境
    echo "Conda version: \$(conda --version)"
    echo "Conda path: \$(which conda)"
    echo "BWA path: \$(which bwa)"
    echo "Samtools path: \$(which samtools)"
    start_time=\$(date +%s)

    mkdir -p 2.index
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