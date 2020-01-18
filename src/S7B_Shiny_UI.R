
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

