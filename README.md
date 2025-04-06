# SV & SNP Analysis Pipeline using Nextflow

This pipeline is implemented in Nextflow and modularized to support structural variation (SV), SNP detection, filtering, annotation, and GWAS analysis. It is suitable for large-scale genome studies in aquaculture species and other organisms.

---

## 📁 Directory Structure

```
.
├── main.nf                  # Main workflow
├── nextflow.config          # Nextflow configuration
├── channels.nf              # Channel definitions
├── processes/               # All modular processes
│   ├── alignment.nf
│   ├── build_index.nf
│   ├── call_smoove.nf
│   ├── duphold_filter.nf
│   ├── extract_samples.nf
│   ├── final_filter.nf
│   ├── fimpute_imputation.nf
│   ├── genotyping.nf
│   ├── go_kegg_annotation.nf
│   ├── gwas_analysis.nf
│   ├── init_setup.nf
│   ├── merge_smoove_vcfs.nf
│   ├── sam_to_bam.nf
│   ├── snp_call.nf
│   ├── snp_merge_filter.nf
│   ├── snpeff_annotate.nf
│   ├── sort_index_bam.nf
│   ├── sv_count_summary.nf
│   ├── sv_paste_merge.nf
└── README.md
└── script
```

---

## Module Functions

| Step | Module | Function |
|------|--------|----------|
| 0    | `init_setup` | Create folder structure and remind user to prepare input files |
| 1    | `build_index` | Index the reference genome using BWA and Samtools |
| 2    | `extract_samples` | Extract sample IDs automatically |
| 3    | `align_bwa` | Align reads to reference genome using BWA MEM |
| 4    | `sam_to_bam` + `sort_index_bam` | Convert SAM to BAM, sort, and index |
| 5    | `call_smoove` | Call SVs using smoove |
| 6    | `merge_smoove_vcfs` | Merge VCFs from all samples |
| 7    | `genotype_smoove` | Genotype merged SVs per sample |
| 8    | `duphold_filter` | Filter SVs by DHFFC and SVTYPE |
| 9    | `sv_count_summary` | Count SV types per sample |
| 10   | `sv_paste_merge` | Merge SVs with smoove paste |
| 11   | `final_filter` | Multi-step filtering: scaffold, BND, gap, depth, GT etc. |
| 12   | `snp_call` + `snp_merge_filter` | SNP calling, merging, filtering |
| 13   | `snpeff_annotate` | Annotate SVs and SNPs using snpEff |
| 14   | `go_kegg_annotation` | Build GO/KEGG annotation DB using eggNOG & GO obo |
| 15   | `fimpute_imputation` | Impute missing genotypes using FImpute3 |
| 16   | `gwas_analysis` | Perform GWAS using Plink and GCTA64 |

---

## How to Run

```bash
nextflow run main.nf
```

You can modify `nextflow.config` to customize resource usage, parameters, or profiles.

---

## Script

All bash and python script will upload in `script` documentary.

## Note

This pipeline was tested in a group workstation, please check the enviroment before use.
Contact: [jiaozexin@outlook.com](jiaozexin@outlook.com)