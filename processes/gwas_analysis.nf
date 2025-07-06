process GWAS_analysis {
    tag "gwas"
    label 'bigmem'
    conda 'envs/gwas.yml'
    publishDir "./Results/14.gwas", mode: 'copy'

    input:
    tuple path(vcf_file), path(pheno_file)

    output:
    path "14.gwas/*"

    script:
    """
    mkdir -p 14.gwas
    cd 14.gwas

    start_time=\$(date +%s)
    echo "[\$(date)] Starting GWAS analysis" > gwas.log

    # Step 1: VCF → PLINK
    plink --vcf ${vcf_file} --make-bed --out snp_SV_data >> gwas.log 2>&1

    # Step 2: GCTA GRM
    gcta64 --bfile snp_SV_data --make-grm --out gcta_grm >> gwas.log 2>&1

    # Step 3: GWAS
    gcta64 --grm gcta_grm --pheno ${pheno_file} --reml --out gwas_result >> gwas.log 2>&1

    end_time=\$(date +%s)
    echo "GWAS analysis completed in \$((end_time - start_time)) seconds" >> gwas.log
    """
}
