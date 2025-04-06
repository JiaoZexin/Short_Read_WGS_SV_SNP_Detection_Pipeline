process init_setup {
    """
    mkdir -p 0.data 1.fq.gz 2.index 3.sam 4.bam 5.sorted 6.smoove 7.merge 8.genotype 8.duphold 9.final 10.snp 11.snpeff 12.GO_KEGG 13.fimpute 14.gwas

    echo "Prepare all this data and put it in 0.data："
    echo "  - FASTQ: sampleX_R1.fastq.gz 和 sampleX_R2.fastq.gz"
    echo "  - reference genome: ref.fa"
    echo "  - Annotation file: ref.gff3"
    echo "  - pep file: pep.fa"
    """
}
