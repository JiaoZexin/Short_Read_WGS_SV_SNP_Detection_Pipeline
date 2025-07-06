process Annotation_annovar {
    tag "annovar_annotate_${label}"
    label 'standard'
    conda 'envs/annovar.yml'
    publishDir "./Results/13.annovar", mode: 'copy'

    input:
    tuple val(label), path(input_vcf)

    output:
    path "13.annovar/${label}.multianno.txt"
    path "13.annovar/${label}_annovar.log"

    script:
    """
    mkdir -p 13.annovar
    start_time=\$(date +%s)
    echo "[\$(date)] Running ANNOVAR annotation for ${label}" > 13.annovar/${label}_annovar.log

    convert2annovar.pl -format vcf4 ${input_vcf} \\
        -outfile ${label}.avinput \\
        -includeinfo >> 13.annovar/${label}_annovar.log 2>&1

    if [[ "${label}" == "snp" ]]; then
        table_annovar.pl ${label}.avinput \\
            /path/to/annovar/humandb/ \\
            -buildver hg19 \\
            -out ${label} \\
            -remove \\
            -protocol refGene,cytoBand,gnomad211_exome \\
            -operation g,r,f \\
            -nastring . \\
            -vcfinput >> 13.annovar/${label}_annovar.log 2>&1
    elif [[ "${label}" == "sv" ]]; then
        table_annovar.pl ${label}.avinput \\
            /path/to/annovar/humandb/ \\
            -buildver hg19 \\
            -out ${label} \\
            -remove \\
            -protocol refGene,cytoBand,genomicSuperDups \\
            -operation g,r,r \\
            -nastring . \\
            -vcfinput >> 13.annovar/${label}_annovar.log 2>&1
    else
        echo "[ERROR] Unsupported label: ${label}" >> 13.annovar/${label}_annovar.log
        exit 1
    fi

    mv ${label}.hg19_multianno.txt 13.annovar/${label}.multianno.txt

    end_time=\$(date +%s)
    echo "Annotation completed in \$((end_time - start_time)) seconds" >> 13.annovar/${label}_annovar.log
    """
}
