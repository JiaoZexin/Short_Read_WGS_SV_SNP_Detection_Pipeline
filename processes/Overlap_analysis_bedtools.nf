process Overlap_analysis_bedtools {
    label 'standard'
    conda 'envs/vcftools_bedtools.yml'
    publishDir "./Results/", mode: 'copy'

    input:
    tuple path(sv_beds), path(bed)

    output:
    path "Results/bedtools_overlap/*"

    script:
    """
    set -e
    set -x

    mkdir -p Results/bedtools_overlap

    # 文件名列表
    SV_FILES=(22537.bed 22537.bed_freq_maf_SVPlaudit_anno.txt 8717.new.bed 8717.new.bed_freq_maf_anno.txt 22537.bed_freq_maf_SVPlaudit_anno.withoutLine22188.txt)
    ALLINFO=8717.final_allinfo.bed

    for file in \${SV_FILES[@]}; do
        prefix=\$(basename \$file .txt)
        prefix=\$(basename \$prefix .bed)

        bedtools intersect -a \$file -b ${bed} -wo > Results/bedtools_overlap/\${prefix}.GERP.bed
        bedtools intersect -a \$file -b ${bed} -wo -f 0.50 -r > Results/bedtools_overlap/\${prefix}.GERP.r.50.bed
        bedtools intersect -a \$file -b ${bed} -wo -f 0.90 -r > Results/bedtools_overlap/\${prefix}.GERP.r.90.bed
        bedtools intersect -a \$file -b ${bed} -wo -f 0.99 -r > Results/bedtools_overlap/\${prefix}.GERP.r.99.bed

        # Swap -a and -b for reverse analysis
        bedtools intersect -a ${bed} -b \$file -wo -f 0.50 > Results/bedtools_overlap/\${prefix}.GERP.50.bed
        bedtools intersect -a ${bed} -b \$file -wo -f 0.90 > Results/bedtools_overlap/\${prefix}.GERP.90.bed
        bedtools intersect -a ${bed} -b \$file -wo -f 0.99 > Results/bedtools_overlap/\${prefix}.GERP.99.bed
    done

    # 单独处理 final_allinfo.bed
    bedtools intersect -a ${bed} -b \$ALLINFO -wo >> Results/bedtools_overlap/8717.final_allinfo.GERP.bed
    """
}
