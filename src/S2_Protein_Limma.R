library(tidyverse)
library(limma)


df <- read_tsv("txt/proteinGroups.txt", guess_max = 20000) %>%
{set_names(., gsub(" ", "_", names(.)))}

sample_df <- read_tsv("Sample_metadata.tsv") %>%
  mutate(Sample_name = factor(Sample_name, levels = unique(Sample_name)),
         Model_group = factor(Model_group, levels = unique(Model_group))) %>%
  #For protein
  filter(grepl("Lysate", Enrichment))

contrast_df <- read_tsv("contrast_matrix.tsv")

df1 <- df %>%
  filter(is.na(Reverse),
         is.na(Potential_contaminant),
         is.na(Only_identified_by_site)) %>%
  select(id, one_of(sample_df$Sample)) %>%
  gather(Sample, Intensity, one_of(sample_df$Sample))


df2 <- df1 %>%
  filter(Intensity > 0) %>%
  mutate(Intensity = log2(Intensity))

df2n <- df2 %>%
  right_join(sample_df) %>%
  group_by(id, Batch) %>%
  mutate(Intensity = Intensity - mean(Intensity[which(Pool == "POOL")])) %>%
  group_by(Sample) %>%
  mutate(Intensity = Intensity - median(Intensity, na.rm = TRUE)) %>%
  ungroup()

df2n %>%
  select(-Sample) %>%
  write_tsv("data/Normalized_proteingroup_intensities.tsv")

# Limma -------------------------------------------------------------------

# Make model matrix and contrasts

model_df <- sample_df %>%
  filter(!is.na(Sample_name))

design <- model.matrix(~ 0 + model_df$Model_group)
colnames(design) <- unique(model_df$Model_group)


cont_table <- makeContrasts(contrasts = as.list(contrast_df$Contrast_name), levels = model_df$Model_group)
Comparisons <- dimnames(cont_table)$Contrasts

# Format for Limma

df3 <- df2n %>%
  filter(!is.na(Sample_name)) %>%
  select(id, Sample_name, Intensity) %>%
  #Imputation occurs here, but this is to deal with batch-related missing values
  spread(Sample_name, Intensity, fill = 0) %>% 
  as.data.frame() %>%
  write_tsv("data/Protein_limma_input.tsv") %>%
  column_to_rownames("id")



fit <- lmFit(df3, design)

cont_fit <- contrasts.fit(fit, cont_table)

fit2 <- eBayes(cont_fit)


#Extract topTable for each contrast
f1 <- function(x1){
  a1_name <- colnames(fit2$coefficients)[[x1]]
  a1 <- topTable(fit2, x1, number = Inf) %>%
    rownames_to_column("id") %>%
    mutate(Comparison = a1_name) %>%
    mutate(id = as.integer(id)) %>%
    as_tibble() %>%
    left_join(df %>% select(id),
              by = "id")
  
  a1
}

df4 <- map_df(seq_along(Comparisons), f1) %>%
  filter(!is.na(adj.P.Val))

df4 %>%
  write_tsv("data/Protein_limma_output.tsv")


spread_protein_limma <- df4 %>%
  select(id, logFC, adj.P.Val, Comparison) %>%
  gather(Type, Value, logFC, adj.P.Val) %>%
  unite(Type, Comparison, Type, sep = " ") %>%
  spread(Type, Value) %>%
  mutate(id = as.integer(id))

protein_meta <- read_tsv("data/Protein_metadata.tsv")

protein_quantitative <- df3 %>%
  rownames_to_column("id") %>%
  mutate(id = as.integer(id))


#For filtering, only include proteins identified in at least 1 sample
samples <- as.character(na.omit(sample_df$Sample_name))

summarized_protein <- protein_meta %>%
  left_join(protein_quantitative) %>%
  left_join(spread_protein_limma) %>%
  filter_at(vars(one_of(samples)), any_vars(!is.na(.))) %>%
  write_tsv("data/Protein_summarized_data.tsv")

