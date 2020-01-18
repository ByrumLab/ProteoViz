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



# UI ----------------------------------------------------------------------
# Phospho tab -------------------------------------------------------------
# Tab chunks
Volcano_box <- {
  box(
    width = 12,
    fluidRow(
      width = 12,
      column(3, selectInput("volcanoComparison", label = "Comparison", 
                            choices = Comparisons,
                            multiple = TRUE,
                            selected = Comparisons[[1]])),
      column(1,
             textInput("vp1FC", "FC cutoff:",
                       value = 0)),
      column(1,
             textInput("vp1pvalue", "p-value cutoff:",
                       value = 1)),
      column(1,
             textInput("motifText", "Motif filter:",
                       value = "")),
      column(1,
             checkboxInput("motifRemove", "Remove motif:")),
      column(1,
             actionButton("vp1button", "Update Volcano"))
    ),
    fluidRow(
      box(
        width = 6,
        plotlyOutput("volcanoplot1")
        #    ,
        # fluidRow(
        #   column(6,
        #          plotlyOutput("volcanoplot1"),
        #          selectInput("PhosphoSearch", label = "Highlight a phosphosite:",
        #                      choices = character(0), multiple = TRUE)
      )
    )
  )
}
Volcano_heatmap <- {
  box(
    title = "Heat map of selected phosphosites",
    collapsible = TRUE,
    collapsed = FALSE,
    width = 12,
    fluidRow(
      column(2, checkboxInput("vhcheck", label = "Plot heatmap")),
      column(3,
             actionButton("vhbutton", "Update heatmap"))
    ),
    fluidRow(
      style = "height:900px",
      tabBox(title = "Selected phosphosites", id = "tab1", width = 12,
             tabPanel("Select phosphosites", 
                      plotlyOutput("volcanoHeat", height = "900px"))
      )
    )
  )
}
GGseq_row <- {
  fluidRow(
    box(title = "Sequence motif of selected phosphosites",
        collapsible = TRUE,
        width = 12,
        checkboxInput("ggseqcheck", "Run motifx and plot enrichment", value = TRUE),
        box(selectInput("GGplot_option",
                        "Plot type",
                        choices = list("bits", "prob")),
            plotOutput("ggseq"),
            column(2,downloadButton("downloadGGseqlogo", "Save image")),
            column(2,textInput("GGseqWidth",  "Width (in):", value = 8)),
            column(2,textInput("GGseqHeight", "Height (in):", value = 4)),
            column(2,textInput("GGseqDpi", "Dpi:", value = 600)),
            column(2,textInput("GGseqFilename", "Filename", value = "Ggseq"))
        ),
        box(
          title = "Motif enrichment parameters",
          textInput("motif_text", "Enter motif:",
                    value = "ST"),
          sliderInput(
            "min_seq",
            "Minimum number of sequences:",
            min = 3,
            max = 30,
            value = 5
          ),
          textInput("p_cutoff", "P.value cutoff:",
                    value = .05),
          sliderInput(
            inputId = "Window_size",
            label = "Window size:",
            min = 3,
            max = 30,
            value = 10
          ),
          "Table of significantly enrichmed motifs, if present:",
          dataTableOutput("ggseq2")
        )
    )
  )
}
Select_table <- {
  fluidRow(
    box(
      title = "Maxquant data for selected phosphosites",
      width = 12,
      collapsible = TRUE,
      collapsed = TRUE,
      fluidRow(
        column(3, downloadButton("downloadPhosphoData", "Download"))
      ),
      dataTableOutput("selected_phosphosites")
    )
  )
}

# Tab
Phospho_tab <- {
  tabItem("Phospho",
          Volcano_box,
          GGseq_row,
          Volcano_heatmap,
          Select_table
  )
}


# Protein tab -------------------------------------------------------------
# Tab chunks
Protein_volcano <- {
  box(
    width = 12,
    fluidRow(
      column(3, selectInput("ProteinVolcanoComparison", label = "Comparison", 
                            choices = Comparisons,
                            multiple = TRUE,
                            selected = Comparisons[[1]])),
      column(1,
             textInput("pv1FC", "FC cutoff:",
                       value = 0)),
      column(1,
             textInput("pv1pvalue", "p-value cutoff:",
                       value = 1)),
      column(2,
             actionButton("pv1button", "Update volcano plot"))
      
    ),
    fluidRow(
      column(6,
             plotlyOutput("ProteinVolcanoplot"),
             selectInput("ProteinSearch", label = "Highlight a protein:",
                         choices = character(0), multiple = TRUE)
      ))
  )
  
}

Protein_heatmap <- {
  box(
    title = "Heat map of selected proteins",
    collapsible = TRUE,
    collapsed = FALSE,
    width = 12,
    fluidRow(
      column(2, checkboxInput("phcheck", label = "Plot heatmap"))
    ),
    fluidRow(
      style = "height:900px",
      tabBox(title = "Selected proteins", id = "tab1", width = 12,
             tabPanel("Select proteins", 
                      plotlyOutput("proteinHeat", height = "900px")),
             tabPanel("Summarized", 
                      plotlyOutput("ProteinHeatAverage", height = "900px"))
      )
    )
  )
}

# Tab
Protein_tab <- {
  tabItem("Protein",
          Protein_volcano,
          Protein_heatmap
  )
}

# PTMSig Tab --------------------------------------------------------------


PTMSig_tab <- {
  tabItem("PTMSig",
          fluidRow(
            box(
              width = 12,
              column(4, 
                     selectInput("selectSig", label = "Signature type",
                                 choices = PTMsig_path_types,
                                 selected = PTMsig_path_types[1],
                                 multiple = TRUE)
              ),
              column(4,
                     actionButton("plot1button", "Update heatmap"))
            )
          ),
          fluidRow(style = "height:900px",
                   box(
                     width = 12, plotlyOutput("PTMSig1", height = "900px")
                   )),
          fluidRow(style = "height:900px",
                   box(
                     width = 12, plotlyOutput("PTMSig2", height = "900px")
                   ))
  )
}


# EGSEA Tab ---------------------------------------------------------------

EGSEA1_inputs <- 
  fluidRow(
    column(2, 
           selectInput("EGSEAselectSig", label = "Signature type",
                       choices = sort(unique(tidy_annots$CATEGORY_CODE)),
                       selected = sort(unique(tidy_annots$CATEGORY_CODE)),
                       multiple = TRUE)
    ),
    column(3, 
           selectInput("EGSEAselectComparisons", label = "Comparisons",
                       choices = sort(unique(GSEA_results$Comparison)),
                       selected = sort(unique(GSEA_results$Comparison)),
                       multiple = TRUE)
    ),
    column(3,
           sliderInput("EGSEASigSlider", "Significance range",
                       min = 0, max = 100,
                       value = c(60, 100))),
    column(3,
           actionButton("EGSEAplot1button", "Update heatmap"))
  )


EGSEA_tab <- {
  tabItem("EGSEA",
          fluidRow(
            box(
              width = 12,
              EGSEA1_inputs,
              fluidRow(
                column(1, textInput("EGSEA1_rowname_size", "Row name size: ",
                                    value = 10)),
                column(1, textInput("EGSEA1_plot_height", "Plot height: ",
                                    value = 900))
              )) 
          ),
          fluidRow(
            box(
              width = 12,
              collapsible = TRUE,
              collapsed = FALSE,
              plotlyOutput("EGSEA1", height = "900px")
            )
          ),
          
          fluidRow(
            box(
              width = 12, dataTableOutput("EGSEA_pathway_table")
            )
          ),
          
          fluidRow(
            collapsible = TRUE,
            collapsed = FALSE,
            tabBox(
              title = "EGSEA selected pathway",
              width = 12,
              tabPanel(
                "OnClick",
                plotlyOutput("EGSEA2", height = "900px")
              ),
              tabPanel(
                "Select a pathway",
                selectizeInput("EGSEA3select", "Select3",
                               choices = character(0), multiple = TRUE),
                plotlyOutput("EGSEA3", height = "900px")
              )
            )
            
          ),
          fluidRow(
            box(
              collapsible = TRUE,
              collapsed = FALSE,
              width = 12, dataTableOutput("EGSEA2_table")
            )
          ),
          includeScript("custom.js")
  )
}



# UI setup ----------------------------------------------------------------


ui <- dashboardPage(
  dashboardHeader(),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Phospho", tabName = "Phospho", icon = icon("th")),
      menuItem("Protein", tabName = "Protein", icon = icon("th")),
      menuItem("EGSEA results", tabName = "EGSEA", icon = icon("atom")),
      menuItem("PTM signature", tabName = "PTMSig", icon = icon("atom"))
    )
  ),
  dashboardBody(
    tabItems(PTMSig_tab,
             Phospho_tab,
             EGSEA_tab,
             Protein_tab
    )
  )
)

# Server ------------------------------------------------------------------

server <- function(input, output, session) {
  select_phospho_df <- reactive({
    d <- event_data("plotly_selected", source = "volcano")
    req(!is.null(d))
    phospho_meta %>%
      filter(id_phos %in% d$key)
  })
  dselect_test <- reactive({
    data1 <- event_data("plotly_selected", source = "volcano")
    req(!is.null(data1) && nrow(data1) > 0)
    data1
  })
  
  
  output$volcanoplot1 <- renderPlotly({
    input$vp1button
    volcanoComparison <- isolate(input$volcanoComparison)
    FCcutoff <- isolate(as.numeric(input$vp1FC))
    Pcutoff <- isolate(as.numeric(input$vp1pvalue))
    motif_filter <- isolate(as.character(input$motifText))
    motif_remove <- isolate(input$motifRemove)
    
    
    p1 <- phospho_limma %>%
      filter(Comparison %in% volcanoComparison,
             abs(logFC) > FCcutoff,
             adj.P.Val < Pcutoff) %>%
      left_join(phospho_meta, by = "id")
    
    
    if(str_length(motif_filter) > 1){
      query <- as.character(motif_filter) %>%
        sub("^\\.+", ".", .) %>%
        gsub("\\t", "", .) %>%
        gsub(" ", "", .)
      
      #Motif length
      x3 <- sub("\\[.+\\]", "C", query)
      
      #Motif length
      x4 <- str_length(x3)
      
      #Motif center position
      x5 <- str_length(sub("C.*", "C", x3))
      
      
      query_matches <- bg_seq1 %>%
        mutate(z1 = (Center_position - x5) + 1,
               z2 = Center_position + (x4 - x5),
               short_window = str_sub(Sequence_window,
                                      z1,
                                      z2)) %>%
        filter(grepl(query, short_window)) %>%
        .$id
      if(motif_remove){
        p1 <- filter(p1, !id %in% query_matches)
      } else {
        p1 <- filter(p1, id %in% query_matches)
      }
    }
    
    
    p2 <- p1 %>%
      ggplot(aes(x = logFC, y = -log10(adj.P.Val), key = id_phos,
                 text = paste(
                   " Log2 fold change:",
                   signif(logFC, 5),
                   "\n",
                   "p.value",
                   signif(adj.P.Val, 5),
                   "\n",
                   "id:",
                   id,
                   "\n",
                   "Protein:",
                   Description,
                   "\n",
                   "Gene name:",
                   Gene_name,
                   "\n",
                   "Uniprot ID:",
                   Uniprot_ID,
                   "\n",
                   "Position:",
                   Position,
                   "\n",
                   "Score:",
                   Score,
                   "\n",
                   "Localization prob:",
                   Localization_prob,
                   "\n",
                   "Score diff:",
                   Score_diff,
                   "\n",
                   "Sequence:",
                   `Phospho_(STY)_Probabilities`,
                   "\n"
                 )
      )) +
      geom_point() +
      geom_hline(yintercept = -log10(.05), linetype = "dotted") +
      geom_vline(xintercept = 1, linetype = "dotted") +
      geom_vline(xintercept = -1, linetype = "dotted") +
      labs(x = "log2 FC", y= "-log10 adj.p.value")
    
    ggplotly(p2, source = "volcano", tooltip = "text") %>%
      layout(dragmode = "select")
  })
  output$volcanoplot2 <- renderPlotly({
    
    d <- event_data("plotly_click", source = "volcano")
    req(!is.null(d))
    p <- phos_input2 %>%
      filter(id == d$key)
    
    p2 <- p %>%
      ggplot(aes(x = Sample, y = Value, color = Treatment)) +
      geom_point() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    ggplotly(p2, source = "volcano2") %>%
      layout(dragmode = "zoom")
    
  })
  
  output$volcanoHeat <- renderPlotly({
    req(input$vhcheck)
    input$vhbutton
    d <- event_data("plotly_selected", source = "volcano")
    req(!is.null(d))
    
    p <- phospho_heat %>%
      filter(id_phos %in% d$key) %>% #Comment this to debug
      separate(id_phos, into = c("id", "phos_number"), sep = "_") %>%
      mutate(id = as.integer(id)) %>%
      left_join(select(phospho_meta, id, Phos_name), by = "id") %>%
      unite(Phos_name, Phos_name, phos_number) %>%
      select(-id)
    
    p %>%
      gather(Sample_name, Value, everything(), -Phos_name) %>%
      mutate(Sample_name = factor(Sample_name, levels = phospho_samples$Sample_name)) %>%
      left_join(phospho_samples) %>%
      # group_by(id_phos, Batch) %>%
      # mutate(Value = Value - mean(Value[which(Drug == "DMSO")], na.rm = TRUE)) %>%
      # ungroup() %>%
      select(Phos_name, Value, Sample_name) %>%
      spread(Sample_name, Value, fill = 0) %>%
      as.data.frame() %>%
      column_to_rownames("Phos_name") %>%
      heatmaply(
        scale_fill_gradient_fun = scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                                                       midpoint = 0, name = "log2 FC"),
        Colv = FALSE
      )
    
    # p %>%
    #   gather(Sample, Value, everything(), -id_phos) %>%
    #   spread(Sample, Value, fill = 0) %>%
    #   as.data.frame() %>%
    #   column_to_rownames("id_phos") %>%
    #   heatmaply(
    #     scale_fill_gradient_fun = scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)
    #   )
  })
  
  # GGseq -------------------------------------------------------------------
  bg_seqs <- reactive({
    req(input$ggseqcheck)
    req(input$volcanoComparison)
    volcanoComparison <- input$volcanoComparison
    bg_comparison <- phospho_limma %>%
      filter(Comparison %in% volcanoComparison)
    half_window <- round(input$Window_size / 2, 0)
    bg_seq1 %>%
      filter(!id %in% gsub("\\_.*$", "", dselect_test()$key),
             id %in% bg_comparison$id) %>%
      mutate(
        Sequence_window2 = str_sub(
          Sequence_window,
          Center_position - half_window,
          Center_position + half_window
        )
      ) %>%
      distinct(Sequence_window2) %>%
      .$Sequence_window2})
  fg_seqs <- reactive({
    req(input$ggseqcheck)
    req(input$volcanoComparison)
    volcanoComparison <- input$volcanoComparison
    bg_comparison <- phospho_limma %>%
      filter(Comparison %in% volcanoComparison)
    half_window <- round(input$Window_size / 2, 0)
    bg_seq1 %>%
      filter(id %in% gsub("\\_.*$", "", dselect_test()$key),
             id %in% bg_comparison$id) %>%
      mutate(
        Sequence_window2 = str_sub(
          Sequence_window,
          Center_position - half_window,
          Center_position + half_window
        )
      ) %>%
      distinct(Sequence_window2) %>%
      .$Sequence_window2})
  output$ggseq <- renderPlot({
    req(input$ggseqcheck)
    d <- event_data("plotly_selected", source = "volcano")
    if (!is.null(d) && length(d$key) > 3){
      #GGseqlogo
      half_window <- round(input$Window_size / 2, 0)
      ggseq1 <- fg_seqs() %>%
        ggseqlogo(method = as.character(input$GGplot_option))
      ggseq1
    }
  })
  output$ggseq2 <- renderDataTable({
    req(input$ggseqcheck)
    req(dselect_test())
    half_window <- round(input$Window_size / 2, 0)
    if (!is.null(bg_seqs()) &&
        !is.null(fg_seqs()) &&
        length(bg_seqs()) > 3 &&
        length(fg_seqs()) > 3) {
      mot <- motifx(
        fg_seqs(),
        bg_seqs(),
        central.res = as.character(input$motif_text),
        min.seqs = input$min_seq,
        pval.cutoff = as.numeric(input$p_cutoff)
      )
      as.tibble(mot) %>%
        mutate_if(is.numeric, funs(round(., 2))) %>%
        arrange(desc(score))
    }
  }, options = list(pageLength = 5),
  
  escape = FALSE)
  output$downloadGGseqlogo <- downloadHandler(
    
    filename = function() {
      str_c(input$GGseqFilename, ".tiff")
    },
    
    content = function(file){
      req(input$ggseqcheck)
      d <- event_data("plotly_selected", source = "volcano")
      if (!is.null(d) && length(d$key) > 3){
        #GGseqlogo
        half_window <- round(input$Window_size / 2, 0)
        ggseq1 <- fg_seqs() %>%
          ggseqlogo(method = as.character(input$GGplot_option))
        ggseq1
      }
      ggsave(file,
             width = as.numeric(input$GGseqWidth),
             height = as.numeric(input$GGseqHeight),
             dpi = as.numeric(input$GGseqDpi),
             compression = "lzw")
      
    }
  )
  
  
  # Select phosphosites table -----------------------------------------------
  output$selected_phosphosites <- renderDataTable({
    req(dselect_test())
    phospho_meta %>%
      filter(id %in% gsub("\\_.*$", "", dselect_test()$key)) %>%
      select(
        Gene_name,
        Protein,
        Position = Phos_name,
        Localization_prob,
        Flanking,
        Sequence_window,
        `Phospho_(STY)_Probabilities`
      ) %>%
      mutate(Position = gsub("^.* ", "", Position),
             Sequence_window = gsub("\\;.*$", "", Sequence_window))
  }, escape = FALSE)
  output$downloadPhosphoData <- downloadHandler(
    filename = "Selected_phospho_data.tsv",
    content = function(file){
      req(dselect_test())
      df1 <- phospho_meta %>%
        filter(id %in% gsub("\\_.*$", "", dselect_test()$key)) %>%
        select(
          Gene_name,
          Protein,
          Position = Phos_name,
          Localization_prob,
          Flanking,
          Sequence_window,
          `Phospho_(STY)_Probabilities`
        ) %>%
        mutate(Position = gsub("^.* ", "", Position),
               Protein_names = gsub("\\;.*$", "", Protein_names),
               Sequence_window = gsub("\\;.*$", "", Sequence_window))
      write_tsv(df1, file)
    }
  )
  
  
  
  
  
  
  
  # PTMsig ------------------------------------------------------------------
  
  output$PTMSig1 <- renderPlotly({
    input$plot1button
    select_sig <- isolate(input$selectSig)
    
    df2_score <- PTMsig_output %>%
      select(id, one_of(sample_df$Sample_name)) %>%
      as.data.frame() %>%
      column_to_rownames("id")
    
    paths <- str_extract(rownames(df2_score), "^[^\\-]+(?=\\-)")
    
    df_plot <- as.matrix(df2_score) %>%
      .[paths %in% select_sig,] %>%
      na.omit()
    
    
    clust <- df_plot %>%
      dist() %>%
      hclust()
    # Get order
    ord <- clust$order
    
    #Clusters columns
    # ord2 <- df_plot %>%
    #   t() %>%
    #   dist() %>%
    #   hclust() %>%
    #   .$order
    
    # Re-arrange based on order
    df1 <- df_plot[ord, ]
    
    # Plot
    p1 <- df1 %>%
      melt()
    
    key <- rownames(p1)
    
    p2 <- p1 %>%
      ggplot(aes(x = Var2, y = Var1, key = Var1)) +
      geom_tile(aes(fill = value)) +
      geom_point(alpha = 0) +
      #geom_text(data = fdr_matrix, aes(label = value)) +
      scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                           name = "Enrichment") +
      labs(x = "", y= "") +
      theme(axis.text.y = element_text(size = 10),
            axis.text.x = element_text(angle = 45, hjust = 1))
    #theme(axis.text.y = element_blank())
    
    ggplotly(p2, source = "phos_GSEA") %>%
      layout(dragmode = "zoom")
  })
  output$PTMSig2 <- renderPlotly({
    d <- event_data("plotly_click", source = "phos_GSEA")
    req(!is.null(d))
    selected_path <- PTM_sigs %>%
      filter(sig_id %in% d$key)
    
    #For debugging
    # selected_path <- PTM_sigs %>%
    #   filter(sig_id == "PERT-PSP_NOCODAZOLE")
    
    path_entries <- unlist(selected_path$sig_entries)
    
    
    
    df2 <- PTMsig_input %>%
      filter(Flanking %in% path_entries) %>%
      .[apply(., 1, function(x){any(!is.na(x))}),]
    
    a <- phospho_meta %>%
      filter(Flanking %in% df2$Flanking) %>%
      mutate(Gene_name = ifelse(is.na(Gene_name), "", Gene_name)) %>%
      mutate(Phos_name = str_c(Description, " ", Gene_name, " ", Amino_acid, Position)) %>%
      select(Flanking, Phos_name)
    
    
    df3 <- df2 %>%
      left_join(a, by = "Flanking") %>%
      select(-Flanking) %>%
      gather(Sample, Value, everything(), -Phos_name) %>%
      mutate(Sample = factor(Sample, levels = unique(sample_df$Sample_name))) %>%
      spread(Sample, Value, fill = 0) %>%
      as.data.frame() %>%
      column_to_rownames("Phos_name") %>%
      as.matrix()
    
    
    heatmaply(df3,
              Colv = FALSE,
              margins = c(50, 400,75,0),
              main = paste0(d$key, "\n", nrow(df2), "/", length(path_entries), " phosphosites in signature"),
              scale_fill_gradient_fun = scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0))
  })
  
  
  
  
  # EGSEA ######
  
  output$EGSEA1 <- renderPlotly({
    
    select_sig <- isolate(input$EGSEAselectSig)
    sig_range <- isolate(input$EGSEASigSlider)
    select_comparisons <- isolate(input$EGSEAselectComparisons)
    rowname_size <- isolate(ifelse(isTruthy(input$EGSEA1_rowname_size),
                                   as.integer(input$EGSEA1_rowname_size),
                                   10))
    input$EGSEAplot1button
    
    # #For debugging:
    # select_sig <- unique(tidy_annots$CATEGORY_CODE)
    # sig_range <- c(60, 100)
    # select_comparisons <- unique(GSEA_results$Comparison)[1]
    
    
    
    df2 <- GSEA_results %>%
      filter(CATEGORY_CODE %in% select_sig,
             Comparison %in% select_comparisons) %>%
      group_by(Pathway) %>%
      filter(any(p.adj < .05),
             any(between(significance, sig_range[1], sig_range[2]))) %>%
      ungroup() %>%
      select(Pathway, significance, Comparison) %>%
      spread(Comparison, significance) %>%
      as.data.frame() %>%
      column_to_rownames("Pathway")
    
    
    df_plot <- as.matrix(df2) %>%
      #.[paths %in% select_sig,] %>%
      na.omit()
    
    
    clust <- df_plot %>%
      dist() %>%
      hclust()
    # Get order
    ord <- clust$order
    
    
    
    if(ncol(df_plot) >1){
      #Clusters columns
      ord2 <- df_plot %>%
        t() %>%
        dist() %>%
        hclust() %>%
        .$order
      
      # Re-arrange based on order
      df1 <- df_plot[ord, ord2]
      
      # Plot
      p1 <- df1 %>%
        melt()
      
    } else {
      df1 <- data.frame(df_plot) %>%
        rownames_to_column("Var1") %>%
        mutate(Var2 = colnames(df_plot)) %>%
        set_names(c("Var1", "value", "Var2")) %>%
        arrange(value) %>%
        mutate(Var1 = factor(Var1, levels = Var1))
      
      p1 <- df1
    }
    
    
    key <- rownames(p1)
    
    p2 <- p1 %>%
      ggplot(aes(x = Var2, y = Var1)) +
      geom_tile(aes(fill = value)) +
      geom_point(alpha = 0,
                 aes(key = Var1,
                     text = paste0("Pathway: ", Var1, "\n",
                                   "Significance: ", signif(value, 4)))) +
      #geom_text(data = fdr_matrix, aes(label = value)) +
      # scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
      #                      name = "Enrichment") +
      scale_fill_viridis_c(name = "Significance \nscore", limits = c(0,100))+
      labs(x = "", y= "") +
      theme(axis.text.y = element_text(size = rowname_size),
            axis.text.x = element_text(angle = 45, hjust = 1))
    #theme(axis.text.y = element_blank())
    
    ggplotly(p2, source = "EGSEA", tooltip = "text") %>%
      layout(dragmode = "zoom")
  })
  
  EGSEA_selected_path <- reactive({
    d <- event_data("plotly_click", source = "EGSEA")
    req(!is.null(d))
    selected_path <- tidy_annots %>%
      filter(STANDARD_NAME %in% d$key)
  })
  
  output$EGSEA_pathway_table <- renderDataTable({
    
    selected_path <- EGSEA_selected_path()
    # For debugging
    # selected_path <- tidy_annots %>%
    #   filter(STANDARD_NAME == "HALLMARK_HYPOXIA")
    
    selected_path %>%
      select(-sig_entries) %>%
      mutate(Broad_URL =  paste0("<a href='",Broad_URL,"'>","Broad URL","</a>"),
             PMID = ifelse(!is.na(PMID),
                           paste0("<a href=https://www.ncbi.nlm.nih.gov/pubmed/",PMID,">", PMID,"</a>"),
                           PMID)) %>%
      select_if(function(x) !is.na(x))
    
    
    
  }, escape = FALSE)
  
  output$EGSEA2 <- renderPlotly({
    
    selected_path <- EGSEA_selected_path()
    #For debugging
    # selected_path <- tidy_annots %>%
    #   filter(STANDARD_NAME == "HALLMARK_G2M_CHECKPOINT")
    
    
    path_entries <- unlist(selected_path$sig_entries)
    
    df3 <- EGSEA_input %>%
      filter(GeneID %in% path_entries) %>%
      select(-GeneID) %>%
      left_join(select(protein_df, Protein_IDs, Protein_name),
                by = "Protein_IDs") %>%
      select(Protein_name, one_of(sample_df$Sample_name)) %>%
      distinct() %>%
      as.data.frame() %>%
      column_to_rownames("Protein_name")
    
    heatmaply(df3,
              Colv = FALSE,
              margins = c(50, 400,75,0),
              main = paste0(selected_path$STANDARD_NAME[[1]], "\n",
                            nrow(df3), "/", length(path_entries), " proteins in gene set"),
              
              scale_fill_gradient_fun = scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                                                             midpoint = 0, name = "log2 FC"))
    
  })
  
  updateSelectizeInput(session, "EGSEA3select", choices = tidy_annots$STANDARD_NAME, server = TRUE)
  output$EGSEA3 <- renderPlotly({
    req(input$EGSEA3select)
    selected_pathname <- as.character(input$EGSEA3select)
    selected_path <- tidy_annots %>%
      filter(STANDARD_NAME %in% selected_pathname)
    
    #For debugging
    # selected_path <- tidy_annots %>%
    #   filter(STANDARD_NAME == "HALLMARK_G2M_CHECKPOINT")
    
    
    path_entries <- unlist(selected_path$sig_entries)
    
    df3 <- EGSEA_input %>%
      filter(GeneID %in% path_entries) %>%
      select(-GeneID) %>%
      left_join(select(protein_df, Protein_IDs, Protein_name),
                by = "Protein_IDs") %>%
      select(Protein_name, one_of(sample_df$Sample_name)) %>%
      distinct() %>%
      as.data.frame() %>%
      column_to_rownames("Protein_name")
    
    heatmaply(df3,
              Colv = FALSE,
              margins = c(50, 400,75,0),
              main = paste0(selected_path$STANDARD_NAME[[1]], "\n",
                            nrow(df3), "/", length(path_entries), " proteins in gene set"),
              
              scale_fill_gradient_fun = scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                                                             midpoint = 0, name = "log2 FC"))
    
  })
  
  output$EGSEA2_table <- renderDataTable({
    
    selected_path <- EGSEA_selected_path()
    # For debugging
    # selected_path <- tidy_annots %>%
    #   filter(STANDARD_NAME == "HALLMARK_HYPOXIA")
    
    path_entries <- unlist(selected_path$sig_entries)
    
    df3 <- EGSEA_input %>%
      filter(GeneID %in% path_entries) %>%
      .$Protein_IDs
    
    protein_df %>%
      filter(Protein_IDs %in% df3) %>%
      select(Description, Gene_name, Protein_IDs, Score) %>%
      mutate(Gene_name = paste0("<a href=https://www.genecards.org/cgi-bin/carddisp.pl?gene=", Gene_name, ">",
                                Gene_name, "</a>"),
             Protein_IDs = paste0("<a href=https://www.uniprot.org/uniprot/", Protein_IDs, ">",
                                  Protein_IDs, "</a>"))
    
  }, escape = FALSE)
  
  # Protein ####
  updateSelectizeInput(session, 'ProteinSearch', choices = protein_df$Protein_name2, server = TRUE)
  
  output$ProteinVolcanoplot <- renderPlotly({
    input$pv1button
    volcanoComparison <- isolate(input$ProteinVolcanoComparison)
    FCcutoff <- isolate(as.numeric(input$pv1FC))
    Pcutoff <- isolate(as.numeric(input$pv1pvalue))
    
    
    p1 <- protein_limma %>%
      filter(Comparison %in% volcanoComparison,
             abs(logFC) > FCcutoff,
             adj.P.Val < Pcutoff) %>%
      left_join(protein_df, by = "id")
    
    p2 <- p1 %>%
      ggplot(aes(x = logFC, y = -log10(adj.P.Val), key = id,
                 text = paste0(
                   "Protein: ",
                   Description,
                   "\n",
                   "Gene name: ",
                   Gene_name,
                   "\n",
                   "Uniprot ID: ",
                   Uniprot_ID,
                   "\n",
                   "Score: ",
                   Score,
                   "\n",
                   "Log2 fold change: ",
                   signif(logFC, 4),
                   "\n",
                   "p.value: ",
                   signif(adj.P.Val, 4),
                   "\n",
                   "id: ",
                   id,
                   "\n"
                 )
      )) +
      geom_point() +
      geom_hline(yintercept = -log10(.05), linetype = "dotted") +
      geom_vline(xintercept = 1, linetype = "dotted") +
      geom_vline(xintercept = -1, linetype = "dotted") +
      labs(x = "log2 FC", y= "-log10 adj.p.value")
    
    if(isTruthy(input$ProteinSearch)){
      p_search <- protein_df %>%
        filter(Protein_name2 %in% input$ProteinSearch,
        ) %>%
        left_join(protein_limma, by = "id") %>%
        filter(Comparison %in% volcanoComparison)
      
      if(nrow(p_search) > 0){
        p2 <- p2 +
          geom_point(data = p_search, color = "red", size = 5)
      }
      
    }
    
    ggplotly(p2, source = "proteinVolcano", tooltip = "text") %>%
      layout(dragmode = "select")
  })
  
  
  #
  output$ProteinHeatAverage <- renderPlotly({
    
    req(input$phcheck)
    d <- event_data("plotly_selected", source = "proteinVolcano")
    req(!is.null(d))
    
    p <- protein_heat %>%
      filter(id %in% d$key) %>% #Comment this to debug
      left_join(select(protein_df, id, Description, Gene_name), by = "id") %>%
      unite(Protein_name, Description, Gene_name, id, sep = " ")
    
    p %>%
      gather(Sample_name, Value, everything(), -Protein_name) %>%
      left_join(select(sample_df, Sample_name, Model_group)) %>%
      group_by(Protein_name, Model_group) %>%
      summarize(Value = mean(Value, na.rm = TRUE)) %>%
      mutate(Model_group = factor(Model_group, levels = unique(sample_df$Model_group))) %>%
      spread(Model_group, Value, fill = 0) %>%
      as.data.frame() %>%
      column_to_rownames("Protein_name") %>%
      heatmaply(
        scale_fill_gradient_fun = scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                                                       midpoint = 0, name = "log2 FC"),
        margins = c(50, 400,75,0),
        Colv = FALSE
      )
  })
  
  output$proteinHeat <- renderPlotly({
    req(input$phcheck)
    d <- event_data("plotly_selected", source = "proteinVolcano")
    req(!is.null(d))
    
    p <- protein_heat %>%
      filter(id %in% d$key) %>% #Comment this to debug
      left_join(select(protein_df, id, Description, Gene_name), by = "id") %>%
      mutate(Description = str_sub(Description, 1, 75)) %>%
      unite(Protein_name, Description, Gene_name, id, sep = " ")
    
    p %>%
      gather(Sample_name, Value, everything(), -Protein_name) %>%
      mutate(Sample_name = factor(Sample_name, levels = protein_samples$Sample_name)) %>%
      left_join(protein_samples) %>%
      # group_by(id_phos, Batch) %>%
      # mutate(Value = Value - mean(Value[which(Drug == "DMSO")], na.rm = TRUE)) %>%
      # ungroup() %>%
      select(Protein_name, Value, Sample_name) %>%
      spread(Sample_name, Value, fill = 0) %>%
      as.data.frame() %>%
      column_to_rownames("Protein_name") %>%
      heatmaply(
        scale_fill_gradient_fun = scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                                                       midpoint = 0, name = "log2 FC"),
        Colv = FALSE
      )
  })
  
  
  #######
}




shinyApp(ui, server)
