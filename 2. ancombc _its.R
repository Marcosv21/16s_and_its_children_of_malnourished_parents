library(ANCOMBC)
library(phyloseq)
library(tidyverse)

# 1. Define paths
matrix_path <- "/home/marcos/colaboracao/paulo_de_melo/RICHARDT-GAMA-LANDGRAF_results/RICHARDT-GAMA-LANDGRAF_results/downstream/4_kma_taxa/ncbi/otu_abundance_matrix_CLEAN.csv"
meta_path <- "/home/marcos/colaboracao/paulo_de_melo/RICHARDT-GAMA-LANDGRAF_results/RICHARDT-GAMA-LANDGRAF_results/downstream/5_statistics/kda/metadata_filtered.csv"
out_dir <- "/home/marcos/colaboracao/paulo_de_melo/RICHARDT-GAMA-LANDGRAF_results/RICHARDT-GAMA-LANDGRAF_results/downstream/5_statistics/ancombc"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# 2. Load Matrix and Metadata
otu_df <- read.csv(matrix_path, check.names = FALSE)
meta_df <- read.csv(meta_path, row.names = 1)
meta_df$Group <- sapply(strsplit(rownames(meta_df), "-"), `[`, 1)

# 3. Process RAW Counts Matrix - Genus Level
tax_cols <- c("tax_id", "superkingdom", "phylum", "class", "order", "family", "genus", "species")

otu_genus <- otu_df %>%
  select(-any_of(tax_cols[!tax_cols %in% "genus"])) %>%
  group_by(genus) %>%
  summarise(across(where(is.numeric), sum)) %>%
  filter(genus != "Unknown" & genus != "Unclassified") %>%
  column_to_rownames("genus")

# 4. Create Phyloseq Object
common_samples <- intersect(colnames(otu_genus), rownames(meta_df))
otu_genus <- otu_genus[, common_samples]
meta_df <- meta_df[common_samples, ]

OTU <- otu_table(as.matrix(otu_genus), taxa_are_rows = TRUE)
META <- sample_data(meta_df)
physeq <- phyloseq(OTU, META)

print(paste("Phyloseq created with", ntaxa(physeq), "genera and", nsamples(physeq), "samples."))

# 5. Execute ANCOM-BC2 with PAIRWISE enabled
print("Running ANCOM-BC2 (pairwise comparisons)...")
set.seed(123) 

out <- ancombc2(
  data = physeq, 
  assay_name = "counts",
  tax_level = NULL,
  fix_formula = "Group", 
  p_adj_method = "fdr", 
  prv_cut = 0.10,     
  lib_cut = 1000,     
  group = "Group", 
  struc_zero = TRUE,  
  neg_lb = TRUE, 
  alpha = 0.05, 
  global = TRUE,
  pairwise = TRUE
)

# 6. Extract and Save Results
write.csv(out$res, file.path(out_dir, "ancombc_vs_CTL.csv"), row.names = FALSE)
write.csv(out$res_pair, file.path(out_dir, "ancombc_pairwise_tratamentos.csv"), row.names = FALSE)

print("Results successfully saved!")