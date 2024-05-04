#!/usr/bin/env Rscript
# ============================================================================ #
## 01_mgngunc_preprocess_pupil.R
## MGNGUncertainty pre-processing of pupillometry data.
# Johannes Algermissen, 2023.

### Clear workspace and console.
cat("\014\n")
rm(list = ls()); 

### Retrieve inputs from bash:
args <- commandArgs(trailingOnly = TRUE)

### Extract info:
segMessage <- as.character(args[1])
if (is.na(segMessage)){
  segMessage <- "StartMask"
  cat("Errorneous input argument, set segMessage to", segMessage, "\n")
}

beforeSamples <- as.numeric(args[2])
if (is.na(beforeSamples)){
  beforeSamples <- 1000
  cat("Not beforeSamples provided, set to", beforeSamples, "\n")
}

afterSamples <- as.numeric(args[3])
if (is.na(afterSamples)){
  afterSamples <- 4000
  cat("Not afterSamples provided, set to", afterSamples, "\n")
}

# ============================================================================ #
#### Set directories, load custom functions: ####

## Set codeDir:
codeDir    <- dirname(rstudioapi::getSourceEditorContext()$path)
helperDir <- paste0(codeDir, "helpers/")

## Load packages:
require(stringr)

## Load directories:
source(paste0(helperDir, "set_dirs.R")) # Load packages and options settings
dirs <- set_dirs(rootDir)

# ------------------------------------------------- #
## Load custom functions:

source(paste0(dirs$codeDir, "functions/00_mgngunc_functions_preprocess_pupil.R")) # load functions for pupil preprocessing

## Use eyelinker:
# https://cran.r-project.org/web/packages/eyelinker/vignettes/basics.html
# https://cran.r-project.org/web/packages/eyelinker/eyelinker.pdf

# ============================================================================ #
#### Detect files: ####

## Detect files:
allFileNames <- list.files(paste0(dirs$rawDataDir, "pupil"), pattern = ".asc", full = TRUE)
nFile <- length(allFileNames)
cat(paste0("Found ", nFile, " raw data files\n"))

## Initialize count of trials with NA dilation:
nDeleted <- rep(NA, nFile)

# ============================================================================ #
#### Initialize target directory: ####

## Specify message for segmentation and # samples after message:

# segMessage <- "StartMask"; beforeSamples <- 1000; afterSamples <- 2966
# segMessage <- "StartCue"; beforeSamples <- 1000; afterSamples <- 2600
# segMessage <- "StartResponse"; beforeSamples <- 2000; afterSamples <- 3500
# segMessage <- "StartOutcome"; beforeSamples <- 1000; afterSamples <- 2000

# beforeSamples <- 1000
# afterSamples <- 4000 # 2000 2500 3000
lengthEdges <- 250

cat(paste0("Segment from ", beforeSamples, " before until ", afterSamples, " after ", segMessage, " messages\n"))

## Create target directories:
beforeSamplesStr <- str_pad(beforeSamples, width = 4, side = "left", pad = "0")
afterSamplesStr <- str_pad(afterSamples, width = 4, side = "left", pad = "0")
dirs$pupilTimeCourseDir <- paste0(dirs$processedDataDir, "pupil/", "pupil_timecourse_", segMessage, "_", beforeSamplesStr, "_", afterSamplesStr, "ms/")
dir.create(dirs$pupilTimeCourseDir, showWarnings = FALSE, recursive = TRUE)
dirs$pupilTrialDir <- paste0(dirs$processedDataDir, "pupil/", "pupil_trial_", segMessage, "_", beforeSamplesStr, "_", afterSamplesStr, "ms/")
dir.create(dirs$pupilTrialDir, showWarnings = FALSE, recursive = TRUE)

# ============================================================================ #
#### Short pre-processing: ####

## Define target directory including segmentation message and afterEvent duration:

cat("Start short pre-processing\n")

for (iSub in 1:nFile){ # iSub <- 1 
  
  cat("# -------------------------------------------------------------------------------------- #\n")
  cat(paste0(">>> Start subject ", iSub, "\n"))
  
  ## Retrieve file name:
  fileName <- allFileNames[iSub]
  
  ## Load data:
  dat <- readEyelink(fileName)
  # unique(dat$msg$text)
  
  ## Segmentation message (changed if response locked):
  segmentMessage <- segMessage
  if (segmentMessage == "StartResponse"){segmentMessage <- "StartCue"}
  
  
  ## Epoch trials:
  trialData <- readTrials(rawData = dat, segmentMessage = segmentMessage, 
                          beforeEvent = beforeSamples + lengthEdges, 
                          afterEvent = afterSamples + lengthEdges,
                          varMessages = c("Block", "Stimulus", "Condition", "Valence", "Required Action", "Manipulation", "Response", "ACC", "RT", "Outcome"),
                          recode2Num = c(T, F, T, T, T, T, T, T, T, T),
                          varNames = c("block", "stimulus", "condition", "valence", "reqAction", "manipulation", "response", "ACC", "RT", "outcome")
                          )
  
  if (segMessage == "StartResponse"){
    trialData <- respLock(trialData)
  }
  
  ## Detect and delete blinks and saccades:
  blinkData <- deleteBlkSac(rawData = dat, trialData = trialData, selVar = "pupil",
                            removeBlinks = TRUE, removeSaccades = TRUE,
                            zeropadBlinks = 150, zeropadSaccades = 25, percentReject = 50)
  # Sys.sleep(10)
  
  ## Compute derivative:
  derivData <- deleteHighDerivative(data = blinkData, selVar = "pupil",
                                    derivThresh = 2, clusterThresh = 10, isInteractive = F)
  
  ## Interpolate:
  interpData <- applyInterpolation(data = derivData, selVar = "pupil", 
                                   padEdges = FALSE, isInteractive = FALSE)
  
  ## Filter:
  filtData <- applyFilter(data = interpData, selVar = "pupil",
                          filtType = "Butterworth", order = 3, filtPos = "low", cutoff = 6, 
                          deleteEdges = T, lengthEdges = lengthEdges, padEdges = F, 
                          isInteractive = FALSE)
  
  ## Standardize per run:
  # names(filtData)
  # filtData <- defineBlocks(filtData, nBlock = 4)
  stdData <- applyRunStandardization(filtData, stdType = "percentSignalChange")
  
  # --------------------------------------------------------------------------------- #
  ## Save entire pupil time course per trial as matrix with trials in rows and samples in columns:

  completeTimeVec <- seq(-1*beforeSamples, afterSamples, 1) # all possibly expected time segments
  nTime <- length(completeTimeVec) # number time segments
  timeMat <- base::matrix(NA, nrow = nTrial, ncol = nTime) # initialize to NA
  cat("Reformat time series data into trial x time matrix \n")
  for (iTrial in 1:nTrial){ # iTrial <- 1
    rowIdx <- which(stdData$trialnr == iTrial) # rows for this trial
    extractTimeIdx <- which(stdData$trialTime[rowIdx] %in% completeTimeVec) # which time segments from this trial to extract--relative to rowIdx
    insertTimeIdx <- which(completeTimeVec %in% stdData$trialTime[rowIdx]) # which time segments available for this trial to be filled
    stopifnot(length(extractTimeIdx) == length(insertTimeIdx))
    timeMat[iTrial, insertTimeIdx] <- stdData$pupil[rowIdx[extractTimeIdx]] # insert
  }

  
  #--------------------------------------------------------------------------------- #
  ### Compute baselines and dilation per trial:
  
  aggrData <- applyAggregation(stdData, 
                               baselineLim = c(-500, 0), trialLim = c(0, afterSamples), # StartMask
                               # baselineLim = c(-1500, -1000), trialLim = c(0, afterSamples), # StartResponse
                               iSub = iSub)
  head(aggrData)
  
  ## Count trials with NAs:
  nDeleted[iSub] <- sum(is.na(aggrData$dilation))
  cat(paste0("Subject ", iSub, ": ", nDeleted[iSub], " trials without dilations\n"))
    
  ## Retrieve subject number:
  # subStr <- str_pad(iSub, width = 2, side = "left", pad = "0")
  subStr <- str_sub(fileName, -9, -8) # extract subject number from file name
  cat(paste0("Save subject ", subStr, "\n"))
  
  ## Save data:
  write.csv(timeMat, paste0(dirs$pupilTimeCourseDir, "MGNGUNC_Sub_", subStr, "_pupil.csv"), row.names = F)
  write.csv(aggrData, paste0(dirs$pupilTrialDir, "MGNGUNC_Sub_", subStr, "_pupil.csv"), row.names = F)
  
  ## Clean data:
  rm(dat, trialData, blinkData, derivData, interpData, filterData, stdData, aggrData)
  if (iSub == nFile){
    cat("# --------------------------------------------------------- #\n")
    cat("Finished pre-processing all subjects! :-)\n")
  }
}

# END OF FILE.
