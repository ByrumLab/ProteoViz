library(tidyverse)
library(limma)

#Run Protein_Limma.R first.

df <- read_tsv("txt/Phospho (STY)Sites.txt", guess_max = 20000) %>%
{set_names(., gsub(" ", "_", names(.)))} %>%
  filter(
    is.na(Reverse),
    is.na(Potential_contaminant)
  )

sample_df <- read_tsv("Sample_metadata.tsv") %>%
  mutate(Sample_name = factor(Sample_name, levels = unique(Sample_name)),
         Model_group = factor(Model_group, levels = unique(Model_group))) %>%
  #For protein
  filter(grepl("Phos", Enrichment))

protein_df <- read_tsv("data/Normalized_proteingroup_intensities.tsv", guess_max = 100000) %>%
  rename(Protein_intensity = Intensity,
         Protein_group_IDs = id) %>%
  mutate(Protein_group_IDs = as.character(Protein_group_IDs),
         Sample_name = factor(Sample_name, levels = unique(Sample_name))) %>%
  select(Protein_group_IDs, Protein_intensity, Sample_name)

contrast_df <- read_tsv("contrast_matrix.tsv")


#Filters for class I sites, selects relevant columns, separates by phosphosite number
df1 <- df %>%
  filter(Localization_prob > .75) %>%
  select(id, Protein_group_IDs, matches("corrected.*[[:alpha:]].*\\_{3}[[:digit:]]$")) %>%
  gather(Sample, Intensity, everything(), -id, -Protein_group_IDs) %>%
  separate(Sample, into = c("Sample", "Phos_number"), sep = "___") %>%
  right_join(sample_df)

#Normalizes to pool, then centers distribution around median = 0, then normalizes to protein fold change
df2 <- df1 %>%
  filter(Intensity > 0) %>%
  mutate(Intensity = log2(Intensity)) %>%
  group_by(id, Batch, Phos_number) %>%
  mutate(Intensity = Intensity - mean(Intensity[which(Pool == "POOL")])) %>%
  group_by(Sample) %>%
  mutate(Intensity = Intensity - median(Intensity, na.rm = TRUE)) %>%
  ungroup() %>%
  left_join(protein_df) %>%
  mutate(Protein_group_match = ifelse(is.na(Protein_intensity), "0", "1"),
         Protein_group_match = ifelse(grepl("\\;", Protein_group_IDs), ">1", Protein_group_match))


#Keeps only phosphosites with 1:1 correspondence to protein
df3 <- df2 %>%
  filter(Protein_group_match == "1") %>%
  mutate(Intensity = Intensity - Protein_intensity)



# Limma -------------------------------------------------------------------


model_df <- sample_df %>%
  filter(!is.na(Sample_name))

design <- model.matrix(~ 0 + model_df$Model_group)
colnames(design) <- unique(model_df$Model_group)

cont_table <- makeContrasts(contrasts = as.list(contrast_df$Contrast_name), levels = model_df$Model_group)

#Format for Limma
df4 <- df3 %>%
  filter(!is.na(Sample_name)) %>%
  mutate(Intensity = ifelse(is.nan(Intensity), NA_real_, Intensity)) %>%
  unite(id_phos, id, Phos_number) %>%
  select(id_phos, Sample_name, Intensity) %>%
  spread(Sample_name, Intensity, fill = 0) %>%
  write_tsv("data/Phospho_limma_input.tsv") %>%
  as.data.frame() %>%
  column_to_rownames("id_phos")

fit <- lmFit(df4, design)


Comparisons <- dimnames(cont_table)$Contrasts

cont_fit <- contrasts.fit(fit, cont_table)

fit2 <- eBayes(cont_fit)

ncol(fit2$contrasts)
#Extract topTable for each contrast
f1 <- function(x1){
  a1_name <- colnames(fit2$coefficients)[[x1]]
  a1 <- topTable(fit2, x1, number = Inf) %>%
    rownames_to_column("id") %>%
    separate(id, into = c("id", "Phos_number")) %>%
    mutate(Comparison = a1_name) %>%
    mutate(id = as.integer(id)) %>%
    as.tibble() %>%
    left_join(df %>% select(id),
              by = "id")
  
  a1
}

df5 <- map_df(seq_along(Comparisons), f1) %>%
  filter(!is.na(adj.P.Val))

df5 %>%
  write_tsv("data/Phospho_limma_output.tsv")



spread_phospho_limma <- df5 %>%
  unite(id, id, Phos_number) %>%
  select(id, logFC, adj.P.Val, Comparison) %>%
  gather(Type, Value, logFC, adj.P.Val) %>%
  unite(Type, Comparison, Type, sep = " ") %>%
  spread(Type, Value) %>%
  separate(id, into = c("id", "Phos_number"), sep = "_") %>%
  mutate(id = as.integer(id))

phospho_meta <- read_tsv("data/Phospho_metadata.tsv")

phos_quantitative <- df4 %>%
  rownames_to_column("id") %>%
  separate(id, into = c("id", "Phos_number"), sep = "_") %>%
  mutate(id = as.integer(id))

phospho_meta %>%
  left_join(phos_quantitative) %>%
  left_join(spread_phospho_limma) %>%
  filter(!is.na(Phos_number)) %>%
  write_tsv("data/Phospho_summarized_data.tsv")

