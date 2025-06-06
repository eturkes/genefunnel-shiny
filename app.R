#    This file is part of genefunnel-shiny.
#    Copyright (C) 2025  Emir Turkes, UK DRI at UCL
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#    Emir Turkes can be contacted at emir.turkes@eturkes.com

library(shiny)
library(GSEABase)
library(genefunnel)
library(BiocParallel)
library(future)
library(tools)

source("R/helpers.R")

example_matrix <- read.csv("data/sample_mat/human_symbol.csv", row.names = 1)

select_default_geneset <- function(mat) {
  row_ids <- rownames(mat)

  if (any(grepl("^ENSG", row_ids))) {
    return("data/gene_sets/gprofiler_full_hsapiens.ENSG.gmt")
  } else if (any(grepl("^ENSMUS", row_ids))) {
    return("data/gene_sets/gprofiler_full_mmusculus.ENSG.gmt")
  } else {
    prefixes <- gsub("[0-9]+$", "", row_ids)
    cleaned <- prefixes[nchar(prefixes) >= 3 & grepl("^[A-Za-z]+$", prefixes)]

    all_caps <- sum(grepl("^[A-Z]{2,}$", cleaned))
    capitalised <- sum(grepl("^[A-Z][a-z]+$", cleaned))

    if (all_caps >= capitalised) {
      return("data/gene_sets/gprofiler_full_hsapiens.name.gmt")
    } else {
      return("data/gene_sets/gprofiler_full_mmusculus.name.gmt")
    }
  }
}

ui <- fluidPage(
  titlePanel("GeneFunnel App"),

  sidebarLayout(
    sidebarPanel(
      fileInput(
        "matrix_file", "Upload Gene × Sample Matrix (CSV)", accept = ".csv"
      ),
      fileInput("geneset_file", "Upload Gene Sets (GMT)", accept = ".gmt"),
      actionButton("run", "Run GeneFunnel"),
      br(), br(),
      textOutput("auto_geneset_path"),
      verbatimTextOutput("error_text"),
      downloadButton("download", "Download Result")
    ),
    mainPanel(
      h4("Output Preview"),
      tableOutput("result_preview")
    )
  )
)

server <- function(input, output, session) {
  matrix_data <- reactive({
    req(input$matrix_file)
    read.csv(input$matrix_file$datapath, row.names = 1)
  })

  selected_geneset_path <- reactive({
    mat <- if (is.null(input$matrix_file)) example_matrix else matrix_data()
    if (!is.null(input$geneset_file)) {
      return(NULL)
    }
    select_default_geneset(mat)
  })

  gene_sets <- reactive({
    if (!is.null(input$geneset_file)) {
      getGmt(input$geneset_file$datapath)
    } else {
      gmt_path <- selected_geneset_path()
      getGmt(gmt_path)
    }
  })

  error_message <- reactiveVal(NULL)
  result_data <- reactiveVal(NULL)

  observeEvent(input$run, {
    mat <- if (is.null(input$matrix_file)) example_matrix else matrix_data()
    geneset_list <- geneIds(gene_sets())
    param <- MulticoreParam(availableCores())

    tryCatch({
      res <- genefunnel(mat, geneset_list, BPPARAM = param)
      result_data(res)
      error_message(NULL)
    }, error = function(e) {
      result_data(NULL)
      error_message(conditionMessage(e))
    })
  })

  output$result_preview <- renderTable({
    req(result_data())
    head(result_data(), 10)
  }, rownames = TRUE)

  output$auto_geneset_path <- renderText({
    path <- selected_geneset_path()
    if (is.null(path)) {
      "User-supplied gene set file"
    } else {
      paste("Auto-selected gene set:", basename(path))
    }
  })

  output$error_text <- renderText({
    error_message()
  })

  output$download <- downloadHandler(
    filename = function() {
      paste0("genefunnel_outputs_", Sys.Date(), ".zip")
    },
    content = function(file) {
      mat <- if (is.null(input$matrix_file)) example_matrix else matrix_data()
      gmt_path <- if (!is.null(input$geneset_file)) {
        input$geneset_file$datapath
      } else {
        selected_geneset_path()
      }

      zip_path <- prepare_result_archive(
        result_matrix = result_data(),
        base_gmt_path = gmt_path,
        gene_set_dir = "data/gene_sets"
      )
      file.copy(zip_path, file, overwrite = TRUE)
    }
  )
}

shinyApp(ui, server)
