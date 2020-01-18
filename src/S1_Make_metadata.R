library(tidyverse)

if(!dir.exists("data")){dir.create("data")}


df <- read_tsv("txt/proteinGroups.txt", guess_max = 10000) %>%
{set_names(., gsub(" ", "_", names(.)))}


df %>%
  select(Majority_protein_IDs, Fasta_headers, Score, id, `Phospho_(STY)_site_IDs`) %>%
  mutate(Description = str_extract(Fasta_headers, "(?<= )[^\\|]+(?= OS\\=)"),
         Gene_name = str_extract(Fasta_headers, "(?<=GN\\=)[^\\|]+(?= PE\\=)"),
         Uniprot_ID = str_extract(Fasta_headers, "(?<=\\|)[^\\|]+(?=\\|)")) %>%
  write_tsv("data/Protein_metadata.tsv")

phos_df <- read_tsv("txt/Phospho (STY)Sites.txt", guess_max = 20000) %>%
{set_names(., gsub(" ", "_", names(.)))}

phos_df %>%
  select(Proteins:Score, Amino_acid, Sequence_window, `Phospho_(STY)_Probabilities`, Charge, id:Evidence_IDs) %>%
  mutate(Description = str_extract(Fasta_headers, "(?<= )[^\\|]+(?= OS\\=)"),
         Gene_name = str_extract(Fasta_headers, "(?<=GN\\=)[^\\|]+(?= PE\\=)"),
         Uniprot_ID = str_extract(Fasta_headers, "(?<=\\|)[^\\|]+(?=\\|)"),
         Flanking = gsub("\\;.*$", "", Sequence_window) %>%
           str_sub(9,23) %>%
           paste0("-p")) %>%
  write_tsv("data/Phospho_metadata.tsv")
