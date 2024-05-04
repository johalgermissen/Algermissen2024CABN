#!/usr/bin/env Rscript
# ============================================================================ #
## 01_mgngunc_preprocess_gaze.R
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
  segMessage <- "StartCue"
  cat("Errorneous input argument, set segMessage to", segMessage, "\n")
}

beforeSamples <- as.numeric(args[2])
if (is.na(beforeSamples)){
  beforeSamples <- 866
  cat("Not beforeSamples provided, set to", beforeSamples, "\n")
}

afterSamples <- as.numeric(args[3])
if (is.na(afterSamples)){
  afterSamples <- 1300
  cat("Not afterSamples provided, set to", afterSamples, "\n")
}

# ============================================================================ #
#### Load packages: ####

library(devtools)
library(eyelinker)
library(plyr)
library(stringr) # for str_sub

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

source(paste0(dirs$codeDir, "functions/00_mgngunc_functions_regression.R")) # Load functions
source(paste0(dirs$codeDir, "functions/00_mgngunc_functions_preprocess_gaze.R")) # load functions for pupil preprocessing

## Use eyelinker:
# https://cran.r-project.org/web/packages/eyelinker/vignettes/basics.html
# https://cran.r-project.org/web/packages/eyelinker/eyelinker.pdf

# ============================================================================ #
#### ROI Settings: ####

scrx <- 516*2 # screen size x dimension
scry <- 392*2 # screen size y dimension
eccentricity <- 0.30 # eccentricty of potential outcome presentation
tolFix <- 150 # originally 200 tolerance for whether gaze is in ROI
plotFactor <- 3 # how many tolFix to deviate from midline in plots

selVariables <-  c("Block", "TrialNr", "Stimulus", "Condition", "Valence", "Required Action", "Manipulation", "GoValidity", "NoGoValidity", "ISI", "ITI", 
                   "Response", "ACC", "RT", "Outcome")
selVariablesNew <- selVariables
behavVariables <- c("TrialNr", "Stimulus", "Condition", "Valence", "ReqAction", "Manipulation", "Validity", 
                    "Response", "ACC", "RT", "Outcome")

# ============================================================================ #
#### Read in behavioral data: ####

# ----------------------------------- #
## Read in behavioral data:

behavData <- read_behavior()
behavData$TrialNr <- behavData$Trialnr

# ============================================================================ #
#### Detect eye gaze raw data files: ####

## Detect files:
allFileNames <- list.files(paste0(dirs$rawDataDir, "pupil"), pattern = ".asc", full = TRUE)
nFile <- length(allFileNames)
cat(paste0("Found ", nFile, " raw data files\n"))

# ============================================================================ #
#### Create target directories: ####

lengthEdges <- 250

beforeSamplesStr <- str_pad(beforeSamples, width = 4, side = "left", pad = "0")
afterSamplesStr <- str_pad(afterSamples, width = 4, side = "left", pad = "0")
newDir <- paste0("timecourse_", segMessage, "_", beforeSamplesStr, "_", afterSamplesStr, "ms")
cat(paste0("Create new directory ", newDir, "\n"))

dirs$xpDir <- paste0(dirs$processedDataDir, "gaze/xp_", newDir, "/")
dirs$ypDir <- paste0(dirs$processedDataDir, "gaze/yp_", newDir, "/")
dirs$euclDistDir <- paste0(dirs$processedDataDir, "gaze/distanceGroupMean_", newDir, "/")
dir.create(dirs$xpDir, showWarnings = FALSE, recursive = TRUE)
dir.create(dirs$ypDir, showWarnings = FALSE, recursive = TRUE)
dir.create(dirs$euclDistDir, showWarnings = FALSE, recursive = TRUE)

# ============================================================================ #
#### Loop over subjects, save aggregated data: ####

for (iSub in 1:nFile){ # iSub <- 1 
  
  cat("# --------------------------------------------------------- #\n")
  ### Retrieve file name:
  fileName <- allFileNames[iSub]
  
  ### Retrieve subject number:
  
  subStr <- str_sub(fileName, -9, -8) # extract subject number from file name
  cat(paste0(">>> Subject ", subStr, ": Start pre-processing\n"))
  
  ## Load data:
  dat <- readEyelink(fileName)
  # unique(dat$msg$text)[1:40]
  # StartTrial, Block, TrialNr, Stimulus, Condition, Valence, Required Action, Manipulation, GoValidity, NoGoValidity, ISI, ITI, 
  # StartMask, EndMask, StartCue, StopCue, StartOutcome, StopOutcome, EndTrial,
  # Response, ACC, RT, Outcome
  
  ## Extract behavioral data for this subject:
  # names(behavData)
  subBehavData <- behavData[which(behavData$Subject == iSub), behavVariables]
  
  # -------------------------------------------------------------- #
  ### Epoch trials:
  
  cat(paste0(">>> Subject ", subStr, ": Epoch data into trials\n"))
  
  trialData <- readTrials(dat = dat, segmentMessage = segMessage, 
                          beforeEvent = beforeSamples, 
                          afterEvent = afterSamples,
                          variablesMessages = selVariables,
                          recode2Num = rep(F, length(selVariables)),
                          variablesNames = selVariablesNew
  )
  
  ## Convert behavioral variables to numeric:
  trialData <- recode_behavioral_variables(trialData)

  # -------------------------------------------------------------- #
  ### Test if behavioral data correctly saved or trials missing:
  
  cat(paste0(">>> Subject ", subStr, ": Check if behavioral data correctly saved or trials misssing\n"))
  
  ## Select behavioral data put into eye-tracking data (one row per trial):
  subEyeData <- trialData[which(trialData$trialTime == 1), c(names(subBehavData))]
  
  ## Check whether behavioral variables match with those from behavioral data set:
  # Behavioral variables first (use to select non-NA rows), eye data second
  isSame <- compare_data_sets(subBehavData, subEyeData)
  # cbind(subEyeData$respSide, subBehavData$respSide)
  
  ## If not: overwrite
  if (isSame == F){
    cat("Mismatch in behavioral data, overwrite eyetracking data with behavioral variables from behavioral data")
    
    trialData <- trialData[, c("trialnr", "absTime", "trialTime", "xp", "yp", "pupil")] # only keep eye-tracking variables
    subBehavData$trialnr <- subBehavData$TrialNr # rename trial number variable in behavioral data
    
    ## Merge with behavioral data:
    trialData <- merge(trialData, subBehavData, by = "trialnr")
    # head(trialData)
  } else {
    cat("Behavioral and eye-tracking data contain the same trial-by-trial behavioral data :-)\n")
  }
  
  # -------------------------------------------------------------- #
  ### Detect and delete blinks and saccades:

  ## Delete data identified by Eyelink as blinks or saccades (and zeropad around them):
  cat(paste0(">>> Subject ", subStr, ": Delete blinks and saccades\n"))
  blinkData <- trialData # copy over
  blinkData <- deleteBlkSac(rawData = dat, trialData = blinkData, selVar = "xp",
                            removeBlinks = TRUE, removeSaccades = TRUE,
                            zeropadBlinks = 150, zeropadSaccades = 25, percentReject = 50)
  blinkData <- deleteBlkSac(rawData = dat, trialData = blinkData, selVar = "yp",
                            removeBlinks = TRUE, removeSaccades = TRUE,
                            zeropadBlinks = 150, zeropadSaccades = 25, percentReject = 50)
  # Sys.sleep(10)
  
  ## Compute derivative, delete high derivative:
  cat(paste0(">>> Subject ", subStr, ": Compute derivative, delete outliers\n"))
  derivData <- blinkData
  derivData <- deleteHighDerivative(data = blinkData, selVar = "xp",
                                    derivThresh = 2, clusterThresh = 10, isInteractive = F)
  derivData <- deleteHighDerivative(data = blinkData, selVar = "yp",
                                    derivThresh = 2, clusterThresh = 10, isInteractive = F)
  
  # --------------------------------------------------------------------------------- #
  ### Compute Euclidean distance to center of screen:
  
  cat(paste0(">>> Subject ", subStr, ": Epoch data into trials\n"))
  cat(paste0(">>> Compute Euclidean distance to center of screen for subject ", subStr, "\n"))
  derivData$euclDist <- sqrt((derivData$xp - scrx/2)^2 + (derivData$yp - scry/2)^2)
  
  # --------------------------------------------------------------------------------- #
  ## Save entire pupil time course per trial as matrix with trials in rows and samples in columns:
  
  ## Data dimensions:
  nTrial <- length(unique(derivData$trialnr))
  nTime <- length(unique(derivData$trialTime))
  
  ## Convert into trial x time matrices:
  xpMat <- base::matrix(derivData$xp, nrow = nTrial, ncol = nTime, byrow = TRUE)
  ypMat <- base::matrix(derivData$yp, nrow = nTrial, ncol = nTime, byrow = TRUE)
  euclDistMat <- base::matrix(derivData$euclDist, nrow = nTrial, ncol = nTime, byrow = TRUE)
  euclDistMat <- round(euclDistMat, 2) # round to save space and computation time when loading in again
  
  ## Save data:
  cat(paste0(">>> Subject ", subStr, ": Save data ...\n"))
  write.csv(xpMat, paste0(dirs$xpDir, "MGNGUNC_Sub_", subStr, "_gaze.csv"), row.names = F)
  write.csv(ypMat, paste0(dirs$ypDir, "MGNGUNC_Sub_", subStr, "_gaze.csv"), row.names = F)
  write.csv(euclDistMat, paste0(dirs$euclDistDir, "MGNGUNC_Sub_", subStr, "_gaze.csv"), row.names = F)
  cat(paste0(">>> Subject ", subStr, ": Finished saving data! :-)\n"))
   
  if(iSub == nFile){cat("Finished all subjects!!!! :-)\n")}
   
}
cat("Finished preprocessing :-) \n")

# ---------------------------------------------------------------------------- # 
warnings()


# ============================================================================ #
# ============================================================================ #
# ============================================================================ #
# ============================================================================ #
#### Compute Euclidean distance from trial-by-trial baseline: ####

aggrMethod <- "distanceBaselineTrial"

source(paste0(dirs$codeDir, "functions/00_mgngunc_functions_preprocess_gaze.R")) # load functions for pupil preprocessing
nFile <- 35
for (iSub in 1:nFile){
   aggregate_gaze(subjectID = iSub, aggrMethod = aggrMethod, 
                  segMessage = "StartCue", beforeSamples = 866, afterSamples = 1300,
                          baselineTime = c(-866, -367))
}

# END OF FILE.
