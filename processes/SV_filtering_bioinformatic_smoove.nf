process SV_filtering_bioinformatic_smoove {
    tag "final_filter"
    label 'standard'
    conda 'envs/smoove.yml'
    publishDir "./Results/", mode: 'copy'

    input:
    path merged_vcf

    output:
    path "9.final/*_filter.vcf"
    path "9.final/final_filter.log"

    script:
    """
    mkdir -p 9.final
    start_time=\$(date +%s)
    merge_name=\$(basename "${merged_vcf}" .smoove.square.vcf.gz)
    echo "[\$(date)] Starting final filtering for \$merge_name" > 9.final/final_filter.log

    # Step 1: split vcf
    gunzip -c ${merged_vcf} > 9.final/\${merge_name}.vcf
    bcftools view -h 9.final/\${merge_name}.vcf > 9.final/head.txt
    bcftools view -H 9.final/\${merge_name}.vcf > 9.final/body.txt

    # Step 2 & 3: remove scafold and BND
    cut -f 1 2.index/ref.fa.fai > 9.final/scaffold.id.txt
    if [ -f 9.final/scaffold.filter.id.txt ]; then
        grep -v -F -f 9.final/scaffold.filter.id.txt 9.final/body.txt | grep -v 'SVTYPE=BND' > 9.final/body.filtered
    else
        cp 9.final/body.txt 9.final/body.filtered
        echo "[WARN] scaffold.filter.id.txt not found, skipping scaffold filtering" >> 9.final/final_filter.log
    fi

    # Step 4: gap region 
    bash ./script/gap_bed.sh 2.index/ref.fa 9.final/ref.gap.txt
    awk '{print \$1"\t"\$2"\t"\$3}' 9.final/ref.gap.txt > 9.final/gap.bed

    # Step 5: mosdepth for high depth region
    bedtools makewindows -g 2.index/ref.fa.fai -w 1000 > 9.final/ref.1000.bed
    
    mkdir -p 9.final/mosdepth
    for i in \$(cat sample.txt); do
        mkdir -p 9.final/mosdepth/\${i}
        mosdepth -t 2 -b 9.final/ref.1000.bed 9.final/mosdepth/\${i}/\${i} 5.sorted/\${i}.sorted.bam
    done
    
    mkdir -p 9.final/mosdepth.final
    cp 9.final/mosdepth/*/*.regions.bed.gz 9.final/mosdepth.final/
    
    for i in \$(cat sample.txt); do
        gunzip -c 9.final/mosdepth.final/\${i}.regions.bed.gz \
        | awk '\$4 >= 100 {print \$1 "\t" \$2 "\t" \$3 "\t" \$4}' \
        > 9.final/\${i}.filter.bed
    done
    
    cat 9.final/*.filter.bed | awk '{print \$1 "\t" \$2 "\t" \$3}' | sort -k1,1 -k2,2n | uniq > 9.final/high_coverage.bed


    # Step 6: for bed
    awk '\$3>\$2{OFS="	"; print \$1"\t"\$2"\t"\$3"\t"\$1"_"\$2}' 9.final/body.filtered > 9.final/body.bed.withid

    # Step 7: remove gap & high coverage
    if [ -f 9.final/high_coverage.bed ]; then
        bedtools intersect -v -a 9.final/body.bed.withid -b 9.final/gap.bed 9.final/high_coverage.bed > 9.final/body.bed.clean
    else
        bedtools intersect -v -a 9.final/body.bed.withid -b 9.final/gap.bed > 9.final/body.bed.clean
        echo "[WARN] high_coverage.bed not found, only gap filtering applied" >> 9.final/final_filter.log
    fi

    # Step 8: GT filtering 
    cut -f4 9.final/body.bed.clean > 9.final/clean.id
    
    bash ./script/gt_count.sh 9.final/body.bed.clean 9.final/total_missing.txt 9.final/total_00.txt 9.final/total_01.txt 9.final/total_11.txt
    
    paste 9.final/total_missing.txt 9.final/total_00.txt 9.final/total_01.txt 9.final/total_11.txt 9.final/clean.id > 9.final/gt_stats.txt
    
    sample_count=\$(wc -l < sample.txt)
    
    awk -v N=\$sample_count '{
        het_hom = \$3 + \$4
        total = \$1 + \$2 + \$3 + \$4
        if (( (\$1+\$2)==N ) ||
        (\$3==1 && \$4==0) ||
        (\$3==2 && \$4==0) ||
        (\$3==0 && \$4==1) ||
        (\$3+\$4==N) ||
        (\$4==N-1) ||
        (\$3==2 && \$4==N-2) ||
        (het_hom + \$2 >= 0.9 * N))
        print \$5
    }' 9.final/gt_stats.txt > 9.final/gt_filtered_ids.txt
    
    grep -F -f 9.final/gt_filtered_ids.txt 9.final/body.bed.clean > 9.final/body.bed.gt_pass
    cut -f1-3 9.final/body.bed.gt_pass > 9.final/body.vcf.filtered


    # Step 9: vcf file
    cat 9.final/head.txt 9.final/body.vcf.filtered > 9.final/\${merge_name}_filter.vcf

    end_time=\$(date +%s)
    echo "Filtering completed in \$((end_time - start_time)) seconds" >> 9.final/final_filter.log
    """
}
