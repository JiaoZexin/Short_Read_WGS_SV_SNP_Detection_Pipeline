process sv_count_summary {
    tag { sample_id }

    input:
    tuple val(sample_id), path(original_vcf), path(filtered_vcf)

    output:
    path "8.duphold/${sample_id}_sv_summary.tsv"

    script:
    """
    mkdir -p 8.duphold
    echo -e "sample	status	SVTYPE	count" > 8.duphold/${sample_id}_sv_summary.tsv

    for file in ${original_vcf} ${filtered_vcf}; do
        status=\$(basename \$file | grep -q 'filter' && echo "filtered" || echo "unfiltered")
        bcftools query -f '%INFO/SVTYPE
' \$file |
        sort | uniq -c |
        awk -v sname=${sample_id} -v stat=\$status '{print sname"	"stat"	"\$2"	"\$1}' >> 8.duphold/${sample_id}_sv_summary.tsv
    done
    """
}