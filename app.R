library(shiny)
library(Rsamtools)
library(DT)
library(dplyr)
library(shinycssloaders)
library(shinydisconnect)

# Load taxonomy file
taxonomy <- readRDS("gtdb_taxonomy.rds")

ui <- fluidPage(
  disconnectMessage(text = "😵 Wi-Fi faceplanted. Refresh to revive!", refresh = "Wake it up"),
  actionButton("disconnect", "Leave"),
  tags$head(
    tags$style(HTML("
      table.dataTable tbody td {
        font-size: 12px;
      }
      table.dataTable thead th {
        font-size: 13px;
      }
      #headerImage {
        position: absolute;
        top: 10px;
        right: 10px;
        width: 250px;
        height: auto;
      }
    "))
  ),

  img(src = "chu.png", id = "headerImage"),

  titlePanel("Multi-BAM Paired-End Summary"),
  fileInput("bamfiles", "Upload BAM Files", accept = ".bam", multiple = TRUE),
  uiOutput("tabs_ui")
)

server <- function(input, output, session) {
  bam_data <- reactiveValues()

  observeEvent(input$bamfiles, {
    req(input$bamfiles)
    bam_data_list <- list()

    for (i in seq_len(nrow(input$bamfiles))) {
      file <- input$bamfiles[i, ]
      bam_path <- file$datapath
      bam_name <- tools::file_path_sans_ext(file$name)

      if (!file.exists(paste0(bam_path, ".bai"))) {
        indexBam(bam_path)
      }

      param <- ScanBamParam(what = c("qname", "rname", "pos", "flag", "mapq"))
      aln <- scanBam(bam_path, param = param)[[1]]

      df <- data.frame(
        Read = aln$qname,
        Reference = aln$rname,
        Position = aln$pos,
        Flag = aln$flag,
        Quality = aln$mapq,
        stringsAsFactors = FALSE
      ) %>%
        filter(Flag == 99, !is.na(Reference), !is.na(Quality))

      bam_data_list[[bam_name]] <- df
    }

    bam_data$list <- bam_data_list
  })
output$tabs_ui <- renderUI({
  req(bam_data$list)

  # Add summary tab first
  tabs <- list(
    tabPanel("Summary",
      sidebarLayout(
        sidebarPanel(
          width = 2,
          sliderInput(
            inputId = "summary_mapq",
            label = "Mapping quality",
            min = 10,
            max = 40,
            value = 20,
            step = 1,
            width = "100%"
          )
        ),
        mainPanel(
          width = 10,
          DTOutput("summary_table") %>%
            withSpinner(type = 4, color = "#0dc5c1")
        )
      )
    )
  )

  # Build per-BAM tabs
  bam_tabs <- lapply(seq_along(bam_data$list), function(i) {
    bam_name <- names(bam_data$list)[i]
    df <- bam_data$list[[bam_name]]
    bam_path <- input$bamfiles$datapath[which(tools::file_path_sans_ext(input$bamfiles$name) == bam_name)]

    tabPanel(
      bam_name,
      sidebarLayout(
        sidebarPanel(
          width = 2,
          sliderInput(
            inputId = paste0("mapq_", bam_name),
            label = "Mapping quality",
            min = 10,
            max = 40,
            value = 20,
            step = 1,
            width = "100%"
          )
        ),
        mainPanel(
          width = 10,
          DTOutput(paste0("table_", bam_name)) %>%
            withSpinner(type = 4, color = "#0dc5c1"),
          tags$hr(),
          tags$h5("IGV Viewer"),
          tags$iframe(
            src = paste0("https://igv.org/app/?file=", bam_path),
            width = "100%", height = "600px", frameborder = "0"
          )
        )
      )
    )
  })

  do.call(tabsetPanel, c(tabs, bam_tabs))
})


  observe({
    req(bam_data$list)

    # Per-BAM table renderers
    for (bam_name in names(bam_data$list)) {
      local({
        name <- bam_name
        df <- bam_data$list[[name]]

        output[[paste0("table_", name)]] <- renderDT({
          req(input[[paste0("mapq_", name)]])
          df_filtered <- df %>%
            filter(Quality >= input[[paste0("mapq_", name)]]) %>% select(-Flag) %>%
            merge(taxonomy, ., by = "Reference") %>% select(Read,Taxonomy,Reference, Position, Quality)
          datatable(df_filtered, options = list(pageLength = 10), rownames = FALSE)
        })
      })
    }

    # Summary table with MAPQ filter
    output$summary_table <- renderDT({
      req(length(bam_data$list) > 0)

      min_mapq <- input$summary_mapq

      count_list <- lapply(names(bam_data$list), function(name) {
        df <- bam_data$list[[name]]
        df <- df[df$Quality >= min_mapq, ]
        ref_counts <- as.data.frame(table(df$Reference))
        colnames(ref_counts) <- c("Reference", name)
        ref_counts
      })

      if (length(count_list) == 0) return(NULL)

      merged_counts <- Reduce(function(x, y) merge(x, y, by = "Reference", all = TRUE), count_list)
      merged_counts[is.na(merged_counts)] <- 0

      summary_with_tax <- merge(taxonomy, merged_counts, by = "Reference")
      summary_with_tax$Total <- rowSums(summary_with_tax[, -(1:2)])
      summary_with_tax <- summary_with_tax[summary_with_tax$Total > 0, ]
      summary_with_tax$Total <- NULL
      summary_with_tax <- summary_with_tax[order(-rowSums(summary_with_tax[, -(1:2)])), ]

      datatable(
        summary_with_tax,
        options = list(
          pageLength = 20,
          scrollX = TRUE,
          autoWidth = TRUE,
          dom = "Bfrtip"
        ),
        rownames = FALSE,
        class = "display nowrap"
      )
    })
  })
}

shinyApp(ui = ui, server = server, options = list(port = 5000, host = "0.0.0.0"))
