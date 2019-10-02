library(tidyverse)

#Run Protein_Limma.R first.

if(!dir.exists("data/PTMsig")){dir.create("data/PTMsig")}

df <- read_tsv("txt/Phospho (STY)Sites.txt", guess_max = 20000) %>%
{set_names(., gsub(" ", "_", names(.)))} %>%
  filter(
    is.na(Reverse),
    is.na(Potential_contaminant)
  )


sample_df <- read_tsv("Sample_metadata.tsv") %>%
  mutate(Sample_name = factor(Sample_name, levels = unique(Sample_name)),
         Model_group = factor(Model_group, levels = unique(Model_group))) %>%
  #For phospho
  filter(grepl("Phos", Enrichment))


df1 <- df %>%
  filter(Localization_prob > .75) %>%
  select(id, Sequence_window, matches("corrected.*[[:alpha:]].*\\_{3}[[:digit:]]$")) %>%
  gather(Sample, Intensity, everything(), -id, -Sequence_window) %>%
  separate(Sample, into = c("Sample", "Phos_number"), sep = "___") %>%
  right_join(sample_df) %>%
  filter(Intensity > 0)

df2 <- df1 %>%
  mutate(Flanking = gsub("\\;.*$", "", Sequence_window) %>%
           str_sub(9,23) %>%
           paste0("-p")) %>%
  group_by(Flanking, Sample_name, Batch, Pool) %>%
  summarize(Intensity = sum(Intensity)) %>% #Sums intensity for all phos numbers to yield single phosphosite intensity
  ungroup()


df3 <- df2 %>%
  filter(Intensity > 0) %>%
  mutate(Intensity = log2(Intensity)) %>%
  group_by(Flanking, Batch) %>%
  mutate(Intensity = Intensity - mean(Intensity[which(Pool == "POOL")])) %>% #Normalize to pool
  group_by(Sample_name) %>%
  mutate(Intensity = Intensity - median(Intensity, na.rm = TRUE)) %>% #This probably isn't even necessary for a competetive GSEA
  ungroup() %>%
  filter(!is.na(Sample_name)) %>% #Removes pool
  select(Flanking, Sample_name, Intensity)


to_gct <- function(x1, path){
  df3 <- x1 %>%
    spread(Sample_name, Intensity)
  
  n_sample <- n_distinct(x1$Sample_name)
  n_row <- nrow(df3)
  n_metacol <- ncol(df3) - n_sample - 1
  n_metarow <- 0
  write_file("#1.3\n", path)
  write_file(paste(n_row, n_sample, n_metacol, n_metarow, sep = "\t"), path, append = TRUE)
  write_file("\n", path, append = TRUE)
  write_tsv(df3, path, append = TRUE, col_names = TRUE)
}



to_gct(df3, "data/PTMsig/phospho_PTMsig_input.gct")
