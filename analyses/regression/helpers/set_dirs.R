#!/usr/bin/env Rscript
# ============================================================================ #
## set_dirs.R
## MGNGUncertainty set input and output directories for data processing.
# Johannes Algermissen, 2023.

set_dirs <- function(rootDir){
  
  ## Root directory:
  
  dirs <- list()
  
  # -------------------------------------------------------------------------- #
  ## Root directories:
  
  dirs$root <- "<insert your root directory here>"

  # -------------------------------------------------------------------------- #
  ## Code:
  dirs$codeDir    <- paste0(dirs$root, "analyses/regression/")
  
  # -------------------------------------------------------------------------- #
  ## Data:
  
  ## Data:
  dirs$dataDir <- paste0(dirs$root, "data/")
  
  ## Raw data:
  dirs$rawDataDir <- paste0(dirs$dataDir, "rawData/")
  dirs$questDir <- paste0(dirs$rawDataDir, "questionnaires/")
  
  # -------------------------------------------------------------------------- #
  ## Data sets:
  dirs$processedDataDir <- paste0(dirs$dataDir, "processedData/")
  dir.create(dirs$processedDataDir, recursive = TRUE, showWarnings = FALSE) # recursive = TRUE)
  
  # -------------------------------------------------------------------------- #
  ## Results:
  dirs$resultDir <- paste0(dirs$root, "results/regression/")
  dir.create(dirs$resultDir, showWarnings = FALSE)
  
  ## Models:
  dirs$modelDir <- paste0(dirs$resultDir, "models/")
  dir.create(dirs$modelDir, showWarnings = FALSE)
  
  ## Plots:
  dirs$plotDir <- paste0(dirs$resultDir, "plots/")
  dir.create(dirs$plotDir, showWarnings = FALSE)
  
  return(dirs)
  
}
