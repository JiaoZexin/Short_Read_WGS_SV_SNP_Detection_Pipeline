process gwas_analysis {
    tag "gwas"

    output:
    path "14.gwas/*"

    script:
    """
    mkdir -p 14.gwas
    cd 14.gwas

    start_time=$(date +%s)
    echo "[`date`] Starting GWAS analysis" > gwas.log

    # Step 1: Plink prepare 
    plink --vcf *.vcf --make-bed --out snp_SV_data >> gwas.log 2>&1

    # Step 2: GCTA for grm file
    gcta64 --bfile snp_SV_data --make-grm --out gcta_grm >> gwas.log 2>&1

    # Step 3: GWAS analysis
    gcta64 --grm gcta_grm --pheno phenotype.txt --reml --out gwas_result >> gwas.log 2>&1

    end_time=$(date +%s)
    echo "GWAS analysis completed in $((end_time - start_time)) seconds" >> gwas.log
    """
}