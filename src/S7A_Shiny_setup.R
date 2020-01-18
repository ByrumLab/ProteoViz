library(tidyverse)
library(heatmaply)
library(reshape2)
library(shiny)
library(cowplot)
library(shinydashboard)
library(ggseqlogo)
library(rmotifx)
source("src/gct-io.R")

theme_set(theme_cowplot())

sample_df <- read_tsv("Sample_metadata.tsv") %>%
  filter(!is.na(Sample_name)) # %>%
#mutate(Sample = gsub("^.*corrected\\_", "", Sample))

contrast_df <- read_tsv("contrast_matrix.tsv")

protein_limma <- read_tsv("data/Protein_limma_output.tsv")
protein_heat <- read_tsv("data/Protein_limma_input.tsv")
protein_df <- read_tsv("data/Protein_metadata.tsv") %>%
  mutate(Protein_IDs = str_extract(Fasta_headers, "(?<=sp\\|)[^\\|]+(?=\\|)"),
         Protein_IDs = ifelse(is.na(Protein_IDs),
                              str_extract(Fasta_headers, "(?<=\\|)[^\\|]+(?=\\|)"),
                              Protein_IDs),
         Gene_name = ifelse(is.na(Gene_name), "", Gene_name),
         Description = ifelse(is.na(Description), "", Description),
         Protein_name = str_c(str_sub(Description, 1, 50), Gene_name, sep = " "),
         Protein_name2 = str_c(Protein_name, id, sep = "_")) %>%
  filter(!grepl("^.\\_", Protein_name2))

protein_samples <- sample_df %>%
  filter(Enrichment == "Lysate") %>%
  mutate(Sample_name = factor(Sample_name, levels = Sample_name))


phospho_samples <- read_tsv("Sample_metadata.tsv") %>%
  filter(!is.na(Sample_name),
         Enrichment == "Phos") %>%
  mutate(Sample_name = factor(Sample_name, levels = Sample_name))

phospho_limma <- read_tsv("data/Phospho_limma_output.tsv") %>%
  unite(id_phos, id, Phos_number, remove = FALSE)
phospho_meta <- read_tsv("data/Phospho_metadata.tsv") %>%
  mutate(Gene_name = ifelse(is.na(Gene_name), "", Gene_name),
         Phos_name = str_c(str_sub(Description, 1, 50), " ", Gene_name, " ", Amino_acid, Position))

phospho_heat <- read_tsv("data/Phospho_limma_input.tsv")


#PTMsig
PTMsig_output <- read_tsv("data/PTMsig/output/output-combined.gct", skip = 2)
PTMsig_path_types <- str_extract(PTMsig_output$id, "^[^\\-]+(?=\\-)") %>%
  unique()

PTM_sigs <- parse.gmt("databases/ptm.sig.db.all.flanking.human.v1.8.1.gmt")
PTM_sigs <- tibble(sig_id = map_chr(PTM_sigs, "head"),
                        sig_entries = purrr::map(PTM_sigs, "entry") %>%
                          purrr::map(~gsub("\\;.*$", "", .x)))

PTMsig_input <- read_tsv("data/PTMsig/phospho_PTMsig_input.gct", skip = 2)

#EGSEA

EGSEA_input <- read_tsv("data/Protein_EGSEA_input.tsv") %>%
  mutate(GeneID = as.character(GeneID))


GSEA_sigs <- parse.gmt("databases/msigdb.v6.2.entrez.gmt")
GSEA_sigs <- tibble(STANDARD_NAME = map_chr(GSEA_sigs, "head"),
                    sig_entries = purrr::map(GSEA_sigs, "entry"),
                    Broad_URL = map_chr(GSEA_sigs, "desc"))
GSEA_sig_meta <- read_csv("databases/parsed_msigdb_v6.2.csv", guess_max = 20000)
tidy_annots <- GSEA_sigs %>%
  left_join(GSEA_sig_meta) %>%
  select(-c(HISTORICAL_NAMES, ORGANISM, EXACT_SOURCE, CONTRIBUTOR, CONTRIBUTOR_ORG, CHIP,
            FOUNDER_NAMES, REFINEMENT_DATASETS, VALIDATION_DATASETS))

GSEA_results <- read_tsv("data/EGSEA/EGSEA_test_results.tsv", guess_max = 10000) %>%
  mutate(STANDARD_NAME = Pathway) %>%
  left_join(select(tidy_annots, STANDARD_NAME, CATEGORY_CODE))
GSEA_comparison <- read_tsv("data/EGSEA/EGSEA_comparison.tsv", guess_max = 10000)


Comparisons <- contrast_df$Contrast_name


bg_seq1 <- phospho_meta %>%
  select(id, Sequence_window) %>%
  mutate(Sequence_window = gsub("\\;.*$", "", Sequence_window),
         Window_size = str_length(Sequence_window),
         Center_position = (Window_size+1)/2,
         Center_AA = str_sub(Sequence_window, Center_position, Center_position)) %>%
  filter(Window_size > 0) %>%
  distinct()


