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
library(shinyjs)
source("R/helpers.R")

options(shiny.maxRequestSize = 5 * 1024^3)

example_matrix <- read.csv(
  "data/sample_mat/otero-garcia-pseudobulk.csv", row.names = 1
)

select_default_geneset <- function(mat) {
  row_ids <- rownames(mat)

  if (any(grepl("^ENSG", row_ids))) {
    return("data/gene_sets/gprofiler_full_hsapiens.ENSG.qs")
  } else if (any(grepl("^ENSMUS", row_ids))) {
    return("data/gene_sets/gprofiler_full_mmusculus.ENSG.qs")
  } else {
    prefixes <- gsub("[0-9]+$", "", row_ids)
    cleaned <- prefixes[nchar(prefixes) >= 3 & grepl("^[A-Za-z]+$", prefixes)]

    all_caps <- sum(grepl("^[A-Z]{2,}$", cleaned))
    capitalised <- sum(grepl("^[A-Z][a-z]+$", cleaned))

    if (all_caps >= capitalised) {
      return("data/gene_sets/gprofiler_full_hsapiens.name.qs")
    } else {
      return("data/gene_sets/gprofiler_full_mmusculus.name.qs")
    }
  }
}

ui <- fluidPage(
  useShinyjs(),

  tags$head(
    tags$link(
      href = "https://fonts.googleapis.com/css2?family=Permanent+Marker&display=swap",
      rel = "stylesheet"
    ),
    tags$style(HTML("
      .custom-title {
        font-family: 'Permanent Marker', normal;
        font-size: 32px;
        margin-top: 20px;
        margin-bottom: 10px;
      }
    "))
  ),
  div(class = "custom-title", "GeneFunnel"),

  sidebarLayout(
    sidebarPanel(
      width = 5,
      fileInput(
        "matrix_file", "Upload Gene × Sample Matrix (CSV)", accept = ".csv"
      ),
      fileInput("geneset_file", "Upload Gene Set File (GMT)", accept = ".gmt"),
      uiOutput("auto_matrix_path"),
      uiOutput("auto_geneset_path"),
      br(),
      fluidRow(
        column(
          width = 4,
          div(
            style = "display: flex; align-items: center; gap: 10px;",
            actionButton(
              "run",
              tagList(
                tags$span(
                  icon("play"), style = "margin-right: 2px; margin-left: 2px"
                ),
                "Run GeneFunnel"
              ),
              class = "btn-primary"
            ),
            div(
              id = "run_spinner",
              style = "display: none;",
              tags$i(
                class = "fa fa-spinner fa-spin",
                style = "font-size: 18px; color: #007bff;"
              )
            )
          )
        ),
        column(
          width = 8,
          div(
            style = "display: flex; justify-content: flex-end; align-items: center; gap: 10px;",
            div(
              id = "example_spinner",
              style = "display: none;",
              tags$i(
                class = "fa fa-spinner fa-spin",
                style = "font-size: 18px; color: #007bff;"
              )
            ),
            downloadButton("download_example_matrix", "Download Example Matrix")
          )
        )
      ),
      div(style = "margin-top: 5px;"),
      fluidRow(
        column(
          width = 4,
          div(
            style = "display: flex; align-items: center; gap: 10px;",
            downloadButton("download", "Download Results"),
            div(
              id = "download_spinner",
              style = "display: none;",
              tags$i(
                class = "fa fa-spinner fa-spin",
                style = "font-size: 18px; color: #007bff;"
              )
            )
          )
        ),
        column(
          width = 8,
          div(
            style = "display: flex; justify-content: flex-end; align-items: center; gap: 10px;",
            div(
              id = "default_spinner",
              style = "display: none;",
              tags$i(
                class = "fa fa-spinner fa-spin",
                style = "font-size: 18px; color: #007bff;"
              )
            ),
            downloadButton(
              "download_default_genesets", "Download Default Gene Sets"
            )
          )
        )
      ),
      verbatimTextOutput("error_text")
    ),
    mainPanel(
      width = 7,
      h4("Output Preview"),
      tableOutput("result_preview")
    )
  ),
  tags$div(
    style = "display: flex; justify-content: space-between; align-items: center; padding: 40px 20px 20px 20px;",

    tags$div(
      style = "display: flex; gap: 20px;",
      tags$img(src = "ucl-logo.png", height = "60px"),
      tags$img(src = "ukdri-logo.jpg", height = "60px")
    ),

    tags$div(
      style = "text-align: right;",
      tags$p(
        tags$a(
          href = "https://forms.gle/Nzy5KyQoC9F5KVFA9",
          target = "_blank", "Submit Feedback"
        ),
        style = "margin-bottom: 6px; display: block;"
      ),
      tags$p(
        tags$a(
          href = "https://drive.google.com/file/d/1Q7wc8Ct7OkSdv6tB0HqjLDsQQ3NVCJnw/view?usp=sharing",
          target = "_blank", "Temporary Citation Info (Thesis)"
        ),
        style = "margin-bottom: 6px; display: block;"
      ),
      tags$p(
        tags$a(
          href = "https://github.com/eturkes/genefunnel-shiny",
          target = "_blank", "Source Code & Documentation"
        ),
        style = "margin: 0;"
      )
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
      gmt <- GSEABase::getGmt(input$geneset_file$datapath)
      out <- GSEABase::geneIds(gmt)
      for (i in seq_along(gmt)) {
        desc <- gmt[[i]]@shortDescription
        name <- gmt[[i]]@setName
        if (!is.null(desc) && nzchar(desc)) {
          names(out)[i] <- paste0(name, ": ", desc)
        }
      }
      return(out)
    } else {
      gmt_path <- selected_geneset_path()
      qs::qread(gmt_path)
    }
  })

  error_message <- reactiveVal(NULL)
  result_data <- reactiveVal(NULL)

  observeEvent(input$run, {
    shinyjs::show("run_spinner")
    on.exit(shinyjs::hide("run_spinner"), add = TRUE)

    mat <- if (is.null(input$matrix_file)) example_matrix else matrix_data()
    geneset_list <- gene_sets()
    param <- BiocParallel::MulticoreParam(future::availableCores())

    tryCatch({
      res <- genefunnel::genefunnel(mat, geneset_list, BPPARAM = param)
      result_data(res)
      error_message(NULL)
    }, error = function(e) {
      result_data(NULL)
      error_message(conditionMessage(e))
    })
  })

  output$result_preview <- renderTable({
    req(result_data())
    sorted <- result_data()[order(rowMeans(result_data()), decreasing = TRUE), ]
    head(sorted[, 1:5, drop = FALSE], 15)
  }, rownames = TRUE)

  output$auto_matrix_path <- renderUI({
    header <- if (is.null(input$matrix_file)) {
      "<strong>Using example matrix:</strong> otero-garcia-pseudobulk.csv"
    } else {
      paste0(
        "<strong>User-supplied matrix file:</strong> ", input$matrix_file$name
      )
    }

    details <- paste(
      "<ul style='font-style: italic; padding-left: 0; list-style-position: inside;'>",
      "<li>CSV uploads must be non-negative and all numerical (but zeros, missing values, NA values are OK).</li>",
      "<li>5GB file size limit, will be increasing this to 30GB or so soon.</li>",
      "<li>QC not required (can be applied after obtaining all results, does not affect accuracy).</li>",
      "<li>Preferably unprocessed data (i.e. no log-transform, normalisation, etc., but not a hard requirement).</li>",
      "<li>Output resembles the distribution of the input data and similar pipelines can be applied for analysis.</li>",
      "<li>Proteomics, metabolomics, etc. all supported, but needs gene set file if below criteria not met.</li>",
      "</ul>",
      sep = "\n"
    )

    HTML(paste(header, details, sep = "<br/>"))
  })

  output$auto_geneset_path <- renderUI({
    path <- selected_geneset_path()
    if (is.null(path)) {
      HTML(
        paste0(
          "<strong>User-supplied gene set file:</strong> ",
          input$geneset_file$name
        )
      )
    } else {
      info <- paste(
        "<ul style='font-style: italic; padding-left: 0; list-style-position: inside;'>",
        "<li>Includes all gene sets available at ",
        "<a href='https://biit.cs.ut.ee/gprofiler/gost' target='_blank'>g:Profiler</a> (concatenated and separate).</li>",
        "<li>Human and mouse species supported.</li>",
        "<li>Standard gene symbols and ENSEMBL IDs supported.</li>",
        "<li>If these criteria not met, or interested in custom sets, try submitting your own in the GMT format.</li>",
        "</ul>",
        sep = "\n"
      )

      HTML(
        paste0(
          "<strong>Auto-selected gene set file:</strong> ",
          sub("\\.qs$", ".gmt", basename(path)),
          "<br>",
          info
        )
      )
    }
  })

  output$error_text <- renderText({
    error_message()
  })

  output$download <- downloadHandler(
    filename = function() {
      if (is.null(input$geneset_file)) {
        paste0(
          "genefunnel_default_genesets_", format(Sys.Date(), "%Y_%m_%d"),
          ".zip"
        )
      } else {
        paste0(
          "genefunnel_custom_genesets_", format(Sys.Date(), "%Y_%m_%d"),
          ".zip"
        )
      }
    },
    content = function(file) {
      shinyjs::show("download_spinner")
      on.exit(shinyjs::hide("download_spinner"), add = TRUE)

      mat <- if (is.null(input$matrix_file)) example_matrix else matrix_data()
      gmt_path <- if (!is.null(input$geneset_file)) {
        input$geneset_file$datapath
      } else {
        selected_geneset_path()
      }

      original_gmt_name <- if (!is.null(input$geneset_file)) {
        input$geneset_file$name
      } else {
        basename(gmt_path)
      }

      zip_name <- if (is.null(input$geneset_file)) {
        paste0(
          "genefunnel_default_genesets_", format(Sys.Date(), "%Y_%m_%d"),
          ".zip"
        )
      } else {
        paste0(
          "genefunnel_custom_genesets_", format(Sys.Date(), "%Y_%m_%d"),
          ".zip"
        )
      }

      zip_path <- prepare_result_archive(
        result_matrix = result_data(),
        base_gmt_path = gmt_path,
        gene_set_dir = "data/gene_sets",
        zip_name = zip_name,
        original_filename = original_gmt_name
      )

      file.copy(zip_path, file, overwrite = TRUE)
    }
  )

  output$download_example_matrix <- downloadHandler(
    filename = function() {
      paste0(
        "genefunnel_example_matrix_", format(Sys.Date(), "%Y_%m_%d"), ".zip"
      )
    },
    content = function(file) {
      shinyjs::show("example_spinner")
      on.exit(shinyjs::hide("example_spinner"), add = TRUE)

      temp <- tempfile(fileext = ".csv")
      file.copy("data/sample_mat/otero-garcia-pseudobulk.csv", temp)
      utils::zip(zipfile = file, files = temp, flags = "-j")
    }
  )

  output$download_default_genesets <- downloadHandler(
    filename = function() {
      paste0(
        "genefunnel_default_genesets_", format(Sys.Date(), "%Y_%m_%d"), ".zip"
      )
    },
    content = function(file) {
      shinyjs::show("default_spinner")
      on.exit(shinyjs::hide("default_spinner"), add = TRUE)

      tmpdir <- tempfile()
      dir.create(tmpdir)
      file.copy(
        list.files("data/gene_sets_download", full.names = TRUE), tmpdir
      )
      utils::zip(
        zipfile = file, files = list.files(tmpdir, full.names = TRUE),
        flags = "-j"
      )
    }
  )
}

shinyApp(ui, server)
