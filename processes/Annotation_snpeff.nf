process Annotation_snpeff {
    tag "snpeff_annotate"

    input:
    tuple val(label), path(input_vcf), path(gff_file)

    output:
    path "11.snpeff/${label}.eff.vcf"
    path "11.snpeff/${label}.csv"
    path "11.snpeff/${label}.html"

    script:
    """
    mkdir -p 11.snpeff
    start_time=\$(date +%s)
    echo "[\$(date)] Running snpEff annotation for ${label}" > 11.snpeff/${label}_snpeff.log

    # Prepare Database Imput
    cp ${gff_file} 11.snpeff/data/${label}/genes.gff
    cp 2.index/ref.fa 11.snpeff/data/${label}/sequences.fa

    cd 11.snpeff

    # Generate Database 
    java -jar ~/snpEff.jar build -gff3 -v ${label} -noCheckCds -noCheckProtein >> ${label}_snpeff.log 2>&1

    # Annotate
    java -Xmx20G -jar ~/snpEff.jar -c ~/snpEff.config -v ${label} ${input_vcf} \\
        > ${label}.eff.vcf \\
        -csvStats ${label}.csv -stats ${label}.html 2>> ${label}_snpeff.log

    end_time=\$(date +%s)
    echo "Annotation completed in \$((end_time - start_time)) seconds" >> ${label}_snpeff.log
    """

}
