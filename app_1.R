library(shiny)
library(Rsamtools)

ui <- fluidPage(
  titlePanel("Multi-BAM Viewer"),
  fileInput("bamfiles", "Upload BAM Files", accept = ".bam", multiple = TRUE),
  uiOutput("tabs_ui")
)

server <- function(input, output, session) {
  # Create one reactive per uploaded BAM
  bam_data <- reactive({
    req(input$bamfiles)
    bam_files <- input$bamfiles

    bam_list <- list()

    for (i in seq_len(nrow(bam_files))) {
      bam_path <- bam_files$datapath[i]
      bam_name <- bam_files$name[i]

      # Index BAM if not present
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
      )

      df <- df[!is.na(df$Quality), ]  # Remove NA mapping qualities
      bam_list[[bam_name]] <- df
    }

    bam_list
  })

  # Dynamically generate tabs for each BAM
  output$tabs_ui <- renderUI({
    req(bam_data())
    bam_list <- bam_data()

    tabs <- lapply(names(bam_list), function(bam_name) {
      df <- bam_list[[bam_name]]

      tabPanel(
        title = bam_name,
        sidebarLayout(
          sidebarPanel(
            sliderInput(
              inputId = paste0("mapq_", bam_name),
              label = "Minimum Mapping Quality:",
              min = min(df$Quality, na.rm = TRUE),
              max = max(df$Quality, na.rm = TRUE),
              value = min(df$Quality, na.rm = TRUE)
            )
          ),
          mainPanel(
            tableOutput(outputId = paste0("table_", bam_name))
          )
        )
      )
    })

    do.call(tabsetPanel, tabs)
  })

  # Dynamically create table renderers for each BAM
  observe({
    req(bam_data())
    bam_list <- bam_data()

    for (bam_name in names(bam_list)) {
      local({
        bam_id <- bam_name
        df <- bam_list[[bam_id]]
        output[[paste0("table_", bam_id)]] <- renderTable({
          req(input[[paste0("mapq_", bam_id)]])
          filtered <- df[df$Quality >= input[[paste0("mapq_", bam_id)]], ]
          head(filtered, 100)
        })
      })
    }
  })
}

shinyApp(ui = ui, server = server, options = list(port = 5000, host = "0.0.0.0"))