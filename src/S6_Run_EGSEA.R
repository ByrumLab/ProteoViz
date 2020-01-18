library(tidyverse)

# library(EGSEA) #Loaded later, this takes a while and has namespace collisions.

#Notes:
#Serious conflict with dplyr::select and bioconductor::select

#Current imputation method: Fill NA values with fold change = 0.

if(!dir.exists("data/EGSEA")){dir.create("data/EGSEA")}

df <- read_tsv("txt/proteinGroups.txt", guess_max = 6350) %>%
  {set_names(., gsub(" ", "_", names(.)))}

entrez_df <- read_tsv("databases/Human_entrez_map.tsv.gz") %>%
  set_names(c("Protein_IDs", "Uniprot", "GeneID")) %>%
  mutate(GeneID = gsub("\\;.*$", "", GeneID)) %>%
  dplyr::select(Protein_IDs, GeneID, Uniprot)

sample_df <- read_tsv("Sample_metadata.tsv") %>%
  mutate(Sample_name = factor(Sample_name, levels = unique(Sample_name)),
         Model_group = factor(Model_group, levels = unique(Model_group))) %>%
  filter(Enrichment == "Lysate")

contrast_df <- read_tsv("contrast_matrix.tsv")

# Model matrix and contrast matrix ----------------------------------------

model_df <- sample_df %>%
  filter(!is.na(Sample_name))

design <- model.matrix(~ 0 + model_df$Model_group)
colnames(design) <- unique(model_df$Model_group)

cont_table <- limma::makeContrasts(contrasts = as.list(contrast_df$Contrast_name), levels = model_df$Model_group)

group1 <- model_df$Model_group


# Tidy and filter ---------------------------------------------------------


df1 <- df %>%
  filter(is.na(Reverse),
         is.na(Potential_contaminant),
         is.na(Only_identified_by_site)) %>%
  dplyr::select(id, Fasta_headers, matches("corrected\\_.+\\_TMT[[:digit:]]")) %>%
  gather(Sample, Intensity, matches("corrected\\_.+\\_TMT[[:digit:]]"))



df2 <- df1 %>%
  filter(Intensity > 0,
         !grepl("phos", Sample)) %>%
  mutate(Intensity = log2(Intensity))



# Normalize and impute ----------------------------------------------------

df2n <- df2 %>%
  left_join(sample_df) %>%
  group_by(id, Batch) %>%
  mutate(Intensity = Intensity - mean(Intensity[which(Pool == "POOL")])) %>%
  group_by(Sample) %>%
  mutate(Intensity = Intensity - median(Intensity, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(!is.na(Sample_name))


#Join with Entrez IDs
df4 <- df2n %>%
  mutate(Protein_IDs = str_extract(Fasta_headers, "(?<=sp\\|)[^\\|]+(?=\\|)"),
         Protein_IDs = ifelse(is.na(Protein_IDs),
                              str_extract(Fasta_headers, "(?<=\\|)[^\\|]+(?=\\|)"),
                              Protein_IDs)) %>%
  left_join(entrez_df) %>%
  filter(!is.na(GeneID))


df4n <- df4 %>%
  dplyr::select(Protein_IDs, GeneID, Sample_name, Intensity) %>%
  group_by(Sample_name, GeneID) %>%
  filter(length(GeneID) == 1) %>%
  ungroup() %>%
  spread(Sample_name, Intensity, fill = 0) %>%                  #This is where imputation occurs
  write_tsv("data/Protein_EGSEA_input.tsv")

# Data frame for EGSEA ----------------------------------------------------
# Be careful with rownames from here out

df5 <- df4n %>%
  na.omit()

#First column is possibly irrelevant, second column needs to be EntrezID.
anno <- df5 %>%
  mutate(GeneID2 = GeneID) %>%
  dplyr::select(GeneID, GeneID2, Protein_IDs) %>%
  as.data.frame()

#Make matrix with rownames = EntrezIDs
test1 <- df5 %>%
  dplyr::select(GeneID, one_of(as.character(model_df$Sample_name))) %>%
  as.data.frame() %>%
  column_to_rownames("GeneID") %>%
  as.matrix()

#Rownames of test1 need to match the second column of anno

library(EGSEA)

egsea.base()
baseMethods <- egsea.base()[c(1,3,4,8)]


gs.annots <- buildGMTIdx(geneIDs = rownames(test1), gmt.file = "databases/msigdb.v6.2.entrez.gmt")


# run with FDR cutoff < 0.05
gsam1 = egsea.ma(
  expr = test1,
  group = group1,
  probe.annot = as.matrix(anno),
  design = design,
  contrasts = cont_table,
  probeMap.method = "avg",
  gs.annots = gs.annots,
  baseGSEAs = baseMethods,
  minSize = 10,
  display.top = 30,
  combineMethod = "wilkinson",
  combineWeights = NULL,
  sort.by = "p.adj",
  report.dir = NULL,
  kegg.dir = NULL,
  logFC.cutoff = 0,
  fdr.cutoff = 0.05,
  sum.plot.axis = "p.adj",
  sum.plot.cutoff = NULL,
  vote.bin.width = 5,
  num.threads = 6,
  report = FALSE,
  interactive = FALSE,
  keep.base = FALSE,
  verbose = TRUE,
  keep.limma = TRUE,
  keep.set.scores = TRUE
)


out1 <- gsam1@results

#Prints a bar plot for each contrast
for(i in seq_along(colnames(cont_table))){
  plotBars(gsam1, contrast = i, file.name = paste0("bars_plot", colnames(cont_table)[i]))
}

# Admittedly, this is convoluted/unreadable code. I'm pulling the test results into a tidy data frame
# from the list of list of data frames.  The trick is to pull the rownames out during the process.

# If you see a map.poly error, this is a namespace collision problem. The loaded package uses the function map instead of purrr::map
# You can fix it by calling purrr::map instead of map

a3 <- purrr::map(out1, "test.results") %>%
  map2_df(., names(.), 
          ~map2_df(.x, names(.x),
                   ~mutate(.x, Comparison = .y, Pathway = row.names(.x))) %>%
            mutate(Big_pathway = .y)
  )


a6 <- purrr::map(out1, "comparison") %>%
  map2_df(., names(.), 
          ~as.data.frame(.x) %>%
            rownames_to_column("Pathway") %>%
            mutate(Big_pathway = .y))


a3 %>%
  write_tsv("data/EGSEA/EGSEA_test_results.tsv")

a6 %>%
  write_tsv("data/EGSEA/EGSEA_comparison.tsv")




