Here is a complete and professional README draft tailored for your GitHub repository:

---

# Short-Read WGS SV/SNP Detection and Analysis Pipeline

This repository hosts the **Short-Read Whole Genome Sequencing (WGS) SV/SNP Detection Pipeline**, a modular, scalable, and reproducible workflow designed for structural variation (SV) and single nucleotide polymorphism (SNP) detection, filtering, genotyping, annotation, imputation, and downstream analyses including GWAS and enrichment studies. This pipeline is especially optimized for large-scale aquaculture genome projects but is adaptable to other species and research contexts.


---

## 📌 Key Features

* **End-to-End Analysis**: From raw short-read WGS data to high-confidence SVs, SNPs, annotations, imputation, and GWAS.
* **Highly Modular Design**: Composed of 5 core Nextflow modules and multiple sub-modules.
* **Reproducibility and Scalability**: Implemented using **Nextflow DSL2 (v25.04.6)** with **Conda (v4.12.0)** and **Docker (v26.1.3)** for environment and dependency management.
* **Multi-Caller Support**: Enables SV calling via **Lumpy (Smoove framework)**, **Delly**, **CNVpytor**; SNP calling via **Bcftools**, **Freebayes**, **DeepVariant**.
* **Optional Functional Analysis**: Includes visualization, annotation, feature overlap, KEGG/GO enrichment, and GWAS.
* **Parallel Processing & HPC-Friendly**: Supports SLURM, Docker, Conda, and local environments.

---

## 🗂 Pipeline Modules

| Module                                      | Purpose                                                                      |
| ------------------------------------------- | ---------------------------------------------------------------------------- |
| **Main.nf**                                 | Core WGS data preparation, alignment, SV/SNP calling                         |
| **SV\_genotyping\_smoove.nf**               | SV genotyping using Smoove and SVtyper                                       |
| **SV\_bioinformatic\_filtering\_smoove.nf** | Multi-step filtering of SVs to reduce false positives                        |
| **SNP\_filtering\_from\_bcftools.nf**       | SNP merging and quality filtering (Bcftools outputs)                         |
| **Optional.nf**                             | Downstream analyses: visualization, annotation, GWAS, imputation, enrichment |

---

## 🔍 Pipeline Overview

![Pipeline Diagram](./docs/figure2.2.png)

**Figure:** Workflow structure detailing each module's role. Core modules (green), optional analyses (brown), inputs (blue), and outputs (red) are visualized.

---

## 🚀 Quick Start

### Requirements

* Nextflow >=25.04.6
* Conda >=4.12.0 or Docker >=26.1.3
* (Optional) SLURM for HPC environments

### Clone Repository

```bash
git clone https://github.com/JiaoZexin/Short_Read_WGS_SV_SNP_Detection_Pipeline.git
cd Short_Read_WGS_SV_SNP_Detection_Pipeline
```

### Basic Test Run Example

```bash
nextflow run main.nf -profile conda 
```

For a full test including all modules:

```bash
nextflow run main.nf -profile conda \
  --use_delly=true \
  --use_cnvpytor=true \
  --use_bcftools=true \
  --use_deepvariant=true \
  --use_freebayes=true
```

### Output Directory Structure Example

```
Results/
│── 1.trim/
│── 2.index/
│── 3.sam/
│── ...
│── 9.final/          # Final SV dataset after filtering
│── 10.bcftools_final.snp/  # Final filtered SNP dataset
```

---

## 📊 Example Execution

![Execution Example](./docs/figure2.3.png)

Detailed visualizations of inputs, outputs, and intermediate results are provided in `docs/`.

---

## 🔧 Configuration

Customize `nextflow.config` to define:

* Input files (e.g., FASTQ, reference genome, GFF3)
* Software paths or versions
* Resource allocations
* Optional modules via Boolean parameters (e.g., `use_delly`, `use_gwas`)

Example:

```bash
params {
  use_delly = true
  use_sv_visualize = true
  use_gwas = true
}
```

---

## 📖 Documentation

### Core Modules

* **Main.nf**: Preprocessing, alignment, variant calling.
* **SV\_genotyping\_smoove.nf**: Joint genotyping of SVs using Smoove + SVtyper.
* **SV\_bioinformatic\_filtering\_smoove.nf**: Depth-based and bioinformatic SV filtering using Duphold, Bedtools.
* **SNP\_filtering\_from\_bcftools.nf**: SNP VCF merging + VCFtools filtering.
* **Optional.nf**: SV visualization (Samplot), annotation (SnpEff, ANNOVAR), overlap analysis, GO/KEGG enrichment, GWAS (PLINK+GCTA), imputation (FImpute3).

For detailed descriptions and citations, refer to the **Methods** section of our \[associated manuscript (if public)].

---

## 🖼 Visual Guides

* **Figure 2.1**: Schematic overview of SV/SNP analysis workflow
* **Figure 2.2**: Nextflow pipeline architecture
* **Figures 2.3-2.4**: Execution examples per module
* **Figure 2.5**: SV quality improvements post-filtering

All figures are provided in the `/docs/` folder.


---

## 🛠 Future Work

* Additional support for long-read data integration
* Full workflow containerization for cloud-native execution
* Dynamic parameterization for species-specific defaults

---

## 🤝 Contributing

Contributions are welcome! Please open issues or pull requests for bug fixes, improvements, or new feature suggestions.

---

## 📬 Contact

For questions, collaborations, or feedback:
**Contact**: \[[jiaozexin@outlook.com](mailto:jiaozexin@outlook.com)]
**Lead Developer**: Jiao Zexin

---

> *If you use this pipeline in your research, please cite our associated publication (coming soon) and acknowledge this repository.*

---

Let me know if you want me to generate the **figures in markdown preview format**, **add license details**, or help with the **CITATION file**!
