// ==============================
// channels.nf
// Define all channels used by main.nf
// ==============================

// Read each sample name from sample.txt, and build paired FASTQ paths
fastq_pairs = Channel.fromPath("sample.txt")
    .splitText()
    .map { id ->
        def r1 = file("0.data/${id}_R1.fastq.gz")
        def r2 = file("0.data/${id}_R2.fastq.gz")
        tuple(id, r1, r2)
    }

// Reference genome path (used by align_bwa and other modules)
ref_fasta_ch = Channel.of(file("2.index/ref.fa"))

// Construct each sample's SAM path for samtools view
sam_files_ch = Channel.fromPath("sample.txt")
    .splitText()
    .map { id ->
        tuple(id, file("3.sam/${id}.sam"))
    }

// Construct each sample's BAM path for samtools sort
bam_files_ch = Channel.fromPath("sample.txt")
    .splitText()
    .map { id ->
        tuple(id, file("4.bam/${id}.bam"))
    }

// Provide sorted BAM, BAI, and reference genome for smoove call
smoove_input_ch = Channel.fromPath("sample.txt")
    .splitText()
    .map { id ->
        def bam     = file("5.sorted/${id}.sorted.bam")
        def bai     = file("5.sorted/${id}.sorted.bam.bai")
        def ref_fa  = file("2.index/ref.fa")
        tuple(id, bam, bai, ref_fa)
    }

// Collect all sample-level .vcf.gz files (for merging), send once as list
smoove_vcf_ch = Channel.fromPath("6.smoove/*/*.vcf.gz")
    .collect()

// Get the single merged sites.vcf.gz file (for genotyping)
sites_vcf = Channel.fromPath("7.merge/*.sites.vcf.gz")
    .first()

// Genotyping input channel: each sample's BAM + reference + merged sites.vcf
genotype_input_ch = Channel.fromPath("sample.txt")
    .splitText()
    .map { id ->
        def bam     = file("5.sorted/${id}.sorted.bam")
        def ref_fa  = file("2.index/ref.fa")
        tuple(id, bam, ref_fa, sites_vcf)
    }

// Each sample processed individually for duphold filtering
duphold_input_ch = Channel.fromPath("8.genotype/*/*.vcf.gz")
    .map { vcf ->
        def sample = vcf.getBaseName().replaceAll("-joint-smoove\.genotyped", "")
        tuple(sample, vcf)
    }

// Count SV types before and after filtering
sv_count_input_ch = Channel.fromPath("8.duphold/*.vcf.gz")
    .groupTuple(by: { file ->
        def name = file.getBaseName()
        name.replaceAll(/-filter-name-joint-smoove\.genotyped|original_summary_/, "")
    })
    .map { sample, files ->
        def original = files.find { it.name.startsWith("original_summary_") }
        def filtered = files.find { it.name.contains("filter-name") }
        tuple(sample, original, filtered)
    }

// Provide sorted BAM and reference genome for SNP calling
snp_call_input_ch = Channel.fromPath("sample.txt")
    .splitText()
    .map { id ->
        def bam = file("5.sorted/${id}.sorted.bam")
        def ref = file("2.index/ref.fa")
        tuple(id, bam, ref)
    }

// Channel for snpEff annotations: annotate both SV and SNP
snpeff_input_ch = Channel.of(
    tuple("sv", file("9.final/*_filter.vcf"), file("2.index/ref.gff3")),
    tuple("snp", file("10.snp/snp.final.vcf.recode.vcf"), file("2.index/ref.gff3"))
)
// Protein FASTA required for GO/KEGG annotation
kegg_go_input_ch = Channel.of(file("0,data/pep.fa"))

// Control file used by FImpute imputation
fimpute_input_ch = Channel.of(file("13.fimpute/control_file_final.txt"))




