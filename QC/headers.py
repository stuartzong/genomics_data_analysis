INPUT_FILE_HEADER = [
    "patient",
    "status",
    "DNA_lib",
    "DNA_bam",
]

SUM_SNV_HEADER = [
    "gene",
    "num_patients_gene_level",
    "num_SNVs_gene_level",
    "chromosome",
    "position",
    "ref_base",
    "alt_base",
    "num_patients_SNV_level",
    "patient_ID",
    "snp_ID",
    "gmaf",
    "cosmic64_ID",
    "snpeff_details",
    "in_DNA_single_mpileup",
    "in_paired_mpileup",
    "in_mutseq",
    "in_strelka",
    "in_RNA_single_mpileup"
]

SUM_SNV_HEADER_NON = (
    SUM_SNV_HEADER +
    ["t_DNA_cov",
     "t_DNA_RefC",
     "t_DNA_AltC",
     "t_DNA_AF",
     "DNA_tc",
     "t_RNA_cov",
     "t_RNA_RefC",
     "t_RNA_AltC",
     "t_RNA_AF",
     "RNA_tc"]
)

SUM_SNV_HEADER_WN = (
    SUM_SNV_HEADER +
    ["n_DNA_cov",
     "n_DNA_RefC",
     "n_DNA_AltC",
     "n_DNA_AF",
     "t_DNA_cov",
     "t_DNA_RefC",
     "t_DNA_AltC",
     "t_DNA_AF",
     "adj_t_DNA_cov",
     "adj_t_DNA_RefC",
     "adj_t_DNA_AltC",
     "adj_t_DNA_AF",
     "DNA_fisher_pvalue",
     "DNA_tc",
     "n_RNA_cov",
     "n_RNA_RefC",
     "n_RNA_AltC",
     "n_RNA_AF",
     "t_RNA_cov",
     "t_RNA_RefC",
     "t_RNA_AltC",
     "t_RNA_AF",
     "adj_t_RNA_cov",
     "adj_t_RNA_RefC",
     "adj_t_RNA_AltC",
     "adj_t_RNA_AF",
     "RNA_fisher_pvalue",
     "RNA_tc"]
)

SUM_INDEL_HEADER = [
    "gene",
    "num_patients_gene_level",
    "num_INDELs_gene_level",
    "chromosome",
    "position",
    "ref_base",
    "alt_base",
    "num_patients_INDEL_level",
    "patient_ID",
    "snp_ID", "gmaf",
    "cosmic64_ID",
    "snpeff_details",
    "in_DNA_single_mpileup",
    "in_paired_mpileup",
    "in_strelka",
    "pileup_cov",
    "pileup_RefC",
    "pileup_AltC",
    "pileup_AF",
    "strelka_n_Cov",
    "strelka_n_RefC",
    "strelka_n_AltC",
    "strelka_n_AF",
    "strelka_t_Cov",
    "strelka_t_RefC",
    "strelka_t_AltC",
    "strelka_t_AF"
]
