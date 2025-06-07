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

prepare_result_archive <- function(result_matrix, base_qs_path, gene_set_dir) {
  stopifnot(file.exists(base_qs_path), dir.exists(gene_set_dir))

  zip_file <- tempfile(fileext = ".zip")
  organism_prefix <- sub(
    "^gprofiler_full_", "", sub("\\..*", "", basename(base_qs_path))
  )
  full_csv_path <- file.path(tempdir(), "genefunnel_gprofiler_full.csv")

  write.csv(result_matrix, full_csv_path)
  files_to_zip <- c(full_csv_path)

  all_qs <- list.files(
    gene_set_dir,
    pattern = "\\.qs$",
    full.names = TRUE
  )

  subset_files <- setdiff(
    all_qs,
    file.path(gene_set_dir, paste0("gprofiler_full_", organism_prefix, ".qs"))
  )

  subset_files <- subset_files[
    grepl(paste0("^", organism_prefix, "\\."), basename(subset_files))
  ]

  for (subset_path in subset_files) {
    gene_set <- tryCatch({
      qs::qread(subset_path)
    }, error = function(e) {
      warning("Could not load ", subset_path, ": ", conditionMessage(e))
      return(NULL)
    })

    if (!is.null(gene_set)) {
      ids <- names(gene_set)
      overlapping_ids <- intersect(ids, rownames(result_matrix))

      if (length(overlapping_ids) >= 1) {
        subset_matrix <- result_matrix[overlapping_ids, , drop = FALSE]
        shortname <- basename(subset_path)
        shortname <- sub("^hsapiens\\.|^mmusculus\\.", "", shortname)
        shortname <- sub("\\.qs$", "", shortname)
        shortname <- gsub(":", "_", shortname)
        subset_csv <- file.path(
          tempdir(), paste0("genefunnel_", shortname, ".csv")
        )
        write.csv(subset_matrix, subset_csv)
        files_to_zip <- c(files_to_zip, subset_csv)
      }
    }
  }

  zip(zipfile = zip_file, files = files_to_zip, flags = "-j")
  return(zip_file)
}
