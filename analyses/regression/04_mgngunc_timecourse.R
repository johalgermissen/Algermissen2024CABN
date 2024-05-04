#!/usr/bin/env Rscript
# ============================================================================ #
## 04_mgngunc_timecourse.R
## MGNGUncertainty plots and permutation test on ms-level time courses.
# Johannes Algermissen, 2023.

rm(list = ls())

# ============================================================================ #
#### Set directories, load packages and custom functions: ####

## Set codeDir:
codeDir    <- dirname(rstudioapi::getSourceEditorContext()$path)
helperDir <- paste0(codeDir, "helpers/")

## Load directories:
source(paste0(helperDir, "set_dirs.R")) # Load packages and options settings
dirs <- set_dirs(rootDir)

## Load colours:
source(paste0(helperDir, "load_colours.R")) # Load colours
colours <- load_colours()

## Load packages:
source(paste0(helperDir, "package_manager.R")) # Load packages and options settings

# ------------------------------------------------- #
## Load custom functions:

source(paste0(dirs$codeDir, "functions/00_mgngunc_functions_timecourse.R")) # Load functions
source(paste0(dirs$codeDir, "functions/00_mgngunc_functions_regression.R")) # Load functions

# ============================================================================ #
# ============================================================================ #
# ============================================================================ #
# ============================================================================ #
#### READ IN BEHAVIORAL DATA: ####

# ----------------------------------- #
## Read in behavioral data:
data <- read_behavior()

# ----------------------------------- #
## Preprocessing: 
data <- wrapper_preprocessing(data)
names(data)
table(data$subject_n)

# ----------------------------------- #
## Add pupil data:
segMessage <- "StartMask"; beforeSamples <- 1000; afterSamples <- 2966 # prime-locked/ stimulus-locked time courses
# segMessage <- "StartResponse"; beforeSamples <- 2000; afterSamples <- 3500 # response-locked time courses
data <- add_pupil(data, segMessage = segMessage, afterSamples = afterSamples)

# ============================================================================ #
# ============================================================================ #
# ============================================================================ #
# ============================================================================ #
#### Mask-locked plots: ####

## Settings for loading time courses:
segMessage <- "StartMask"
beforeSamples <- 1000
afterSamples <- 2966

## Cue-locked:
zVar <- "all_f"
p <- plot_timecourse(data = data, zVar = zVar, beforeSamples = beforeSamples, afterSamples = afterSamples, 
                           isBaselineCorGroup = F, addLegend = F, axesLWD = T)
p <- plot_timecourse(data = data, zVar = zVar, beforeSamples = beforeSamples, afterSamples = afterSamples, 
                           isBaselineCorGroup = T, addLegend = F)

zVar <- "response_f"
p <- plot_timecourse(data = data, zVar = zVar, beforeSamples = beforeSamples, afterSamples = afterSamples,
                           isBaselineCorGroup = F, addLegend = F)
p <- plot_timecourse(data = data, zVar = zVar, beforeSamples = beforeSamples, afterSamples = afterSamples,
                           isBaselineCorGroup = T, addLegend = F)

zVar <- "valence_f"
p <- plot_timecourse(data = data, zVar = zVar, beforeSamples = beforeSamples, afterSamples = afterSamples,
                           isBaselineCorGroup = F, addLegend = F)
p <- plot_timecourse(data = data, zVar = zVar, beforeSamples = beforeSamples, afterSamples = afterSamples,
                           isBaselineCorGroup = T, addLegend = F)

zVar <- "goValence_f"
p <- plot_timecourse(data = data, zVar = zVar, beforeSamples = beforeSamples, afterSamples = afterSamples,
                           isBaselineCorGroup = F, addLegend = F)
p <- plot_timecourse(data = data, zVar = zVar, beforeSamples = beforeSamples, afterSamples = afterSamples,
                           isBaselineCorGroup = T, addLegend = F)

zVar <- "arousal_f"
p <- plot_timecourse(data = data, zVar = zVar, beforeSamples = beforeSamples, afterSamples = afterSamples,
                           isBaselineCorGroup = F, addLegend = F)
p <- plot_timecourse(data = data, zVar = zVar, beforeSamples = beforeSamples, afterSamples = afterSamples,
                           isBaselineCorGroup = T, addLegend = F)

# ------------------------------------------------------------------------ #
## Condition involving performed response:

zVar <- "respCond_f"
p <- plot_timecourse(data = data, zVar = zVar, beforeSamples = beforeSamples, afterSamples = afterSamples,
                           selLineType = c(1, 1, 2, 2), 
                           isBaselineCorGroup = F, addLegend = F)
p <- plot_timecourse(data = data, zVar = zVar, beforeSamples = beforeSamples, afterSamples = afterSamples,
                           selLineType = c(1, 1, 2, 2), 
                           isBaselineCorGroup = T, addLegend = F)

## First half of cue repetitions:
p <- plot_timecourse(data = data, zVar = zVar, beforeSamples = beforeSamples, afterSamples = afterSamples,
                           selTrialVar = "firstHalfCueRep_f", selTrialVal = "first",
                           selLineType = c(1, 1, 2, 2), 
                           isBaselineCorGroup = F, addLegend = F)
p <- plot_timecourse(data = data, zVar = zVar, beforeSamples = beforeSamples, afterSamples = afterSamples,
                           selTrialVar = "firstHalfCueRep_f", selTrialVal = "first",
                           selLineType = c(1, 1, 2, 2), 
                           isBaselineCorGroup = T, addLegend = F)

## Second half of cue repetitions:
p <- plot_timecourse(data = data, zVar = zVar, beforeSamples = beforeSamples, afterSamples = afterSamples,
                           selTrialVar = "firstHalfCueRep_f", selTrialVal = "second",
                           selLineType = c(1, 1, 2, 2), 
                           isBaselineCorGroup = F, addLegend = F)
p <- plot_timecourse(data = data, zVar = zVar, beforeSamples = beforeSamples, afterSamples = afterSamples,
                           selTrialVar = "firstHalfCueRep_f", selTrialVal = "second",
                           selLineType = c(1, 1, 2, 2), 
                           isBaselineCorGroup = T, addLegend = F)

# ============================================================================ #
#### Response-locked plots: ####

source(paste0(dirs$codeDir, "functions/00_mgngunc_functions_timecourse.R")) # Load functions

segMessage <- "StartResponse"; beforeSamples <- 2000; afterSamples <- 3500

xLim <- c(-1000, 2500)
baselineTime <- -1000

zVar <- "all_f"
p <- plot_timecourse(data = data, zVar = zVar, segMessage = segMessage, beforeSamples = beforeSamples, afterSamples = afterSamples, 
                     xLim = c(-2000, 2000),
                     isBaselineCorGroup = F, isBaselineCorSub = F, baselineTime = -1000,
                     addLegend = F, axesLWD = T)

zVar <- "respCond_f"
p <- plot_timecourse(data = data, zVar = zVar, segMessage = segMessage, beforeSamples = beforeSamples, afterSamples = afterSamples,
                     selLineType = c(1, 1, 2, 2), 
                     xLim = xLim,
                     isBaselineCorGroup = T, isBaselineCorSub = F, baselineTime = baselineTime, # baseline-correct per group
                     addLegend = F)


# ============================================================================ #
#### Outcome-locked plots: ####

## Settings for loading time courses:
segMessage <- "StartOutcome"
beforeSamples <- 1000
afterSamples <- 2500

names(data)
table(data[, zVar])

zVar <- "outcome_f"
p <- plot_timecourse(data = data, zVar = zVar, segMessage = segMessage, afterSamples = afterSamples,
                           isBaselineCorGroup = F, addLegend = F)
p <- plot_timecourse(data = data, zVar = zVar, segMessage = segMessage, afterSamples = afterSamples,
                           isBaselineCorGroup = T, addLegend = F)
p <- plot_timecourse(data = data, zVar = zVar, segMessage = segMessage, afterSamples = afterSamples,
                           isBaselineCorGroup = F, addLegend = T)
p <- plot_timecourse(data = data, zVar = zVar, segMessage = segMessage, afterSamples = afterSamples,
                           isBaselineCorGroup = T, addLegend = T)

zVar <- "outcome_rel_f"
p <- plot_timecourse(data = data, zVar = zVar, segMessage = segMessage, afterSamples = afterSamples,
                           isBaselineCorGroup = F, addLegend = F)
p <- plot_timecourse(data = data, zVar = zVar, segMessage = segMessage, afterSamples = afterSamples,
                           isBaselineCorGroup = T, addLegend = F)
p <- plot_timecourse(data = data, zVar = zVar, segMessage = segMessage, afterSamples = afterSamples,
                           isBaselineCorGroup = F, addLegend = T)
p <- plot_timecourse(data = data, zVar = zVar, segMessage = segMessage, afterSamples = afterSamples,
                           isBaselineCorGroup = T, addLegend = T)

zVar <- "outcome_all_f"
p <- plot_timecourse(data = data, zVar = zVar, segMessage = segMessage, afterSamples = afterSamples,
                           isBaselineCorGroup = F, addLegend = F)
p <- plot_timecourse(data = data, zVar = zVar, segMessage = segMessage, afterSamples = afterSamples,
                           isBaselineCorGroup = T, addLegend = F)
p <- plot_timecourse(data = data, zVar = zVar, segMessage = segMessage, afterSamples = afterSamples,
                           isBaselineCorGroup = F, addLegend = T)
p <- plot_timecourse(data = data, zVar = zVar, segMessage = segMessage, afterSamples = afterSamples,
                           isBaselineCorGroup = T, addLegend = T)

zVar <- "outcome_response_f"
p <- plot_timecourse(data = data, zVar = zVar, segMessage = segMessage, afterSamples = afterSamples,
                           selLineType = c(1, 2, 1, 2, 1, 2), 
                           isBaselineCorGroup = F, addLegend = F)
p <- plot_timecourse(data = data, zVar = zVar, segMessage = segMessage, afterSamples = afterSamples,
                           selLineType = c(1, 2, 1, 2, 1, 2), 
                           isBaselineCorGroup = T, addLegend = F)
p <- plot_timecourse(data = data, zVar = zVar, segMessage = segMessage, afterSamples = afterSamples,
                           selLineType = c(1, 2, 1, 2, 1, 2), 
                           isBaselineCorGroup = F, addLegend = T)
p <- plot_timecourse(data = data, zVar = zVar, segMessage = segMessage, afterSamples = afterSamples,
                           selLineType = c(1, 2, 1, 2, 1, 2), 
                           isBaselineCorGroup = T, addLegend = T)

## Enforce y-axis limit (baseline-corrected):
# rename with suffix manually
zVar <- "outcome_rel_f"
p <- plot_timecourse(data = data, zVar = zVar, segMessage = segMessage, afterSamples = afterSamples,
                           yLim = c(-30, 60), suffix = "_yLim-30-60",
                           isBaselineCorGroup = T, addLegend = F)
zVar <- "outcome_f"
p <- plot_timecourse(data = data, zVar = zVar, segMessage = segMessage, afterSamples = afterSamples,
                           yLim = c(-30, 60), suffix = "_yLim-30-60",
                           isBaselineCorGroup = T, addLegend = F)
zVar <- "outcome_all_f"
p <- plot_timecourse(data = data, zVar = zVar, segMessage = segMessage, afterSamples = afterSamples,
                           yLim = c(-30, 60), suffix = "_yLim-30-60",
                           isBaselineCorGroup = T, addLegend = F)

# ============================================================================ #
# ============================================================================ #
# ============================================================================ #
# ============================================================================ #
#### Perform permutation test: ####

# ---------------------------------------------------------------------------- #
### Data settings:

segMessage <- "StartMask"; beforeSamples <- 1000; afterSamples <- 2966
segMessage <- "StartResponse"; beforeSamples <- 2000; afterSamples <- 3500

# ---------------------------------------------------------------------------- #
### Perform permutation test:

names(data)

## Select variable to permute:
permVar <- "response_f"
permVar <- "goValence_f"
permVar <- "arousal_f"

## Perform permutation test:
source(paste0(dirs$codeDir, "functions/00_mgngunc_functions_timecourse.R")) # Load functions
set.seed(70)
p <- test_permutation_timecourse(data = data, permVar = permVar, 
                                 selTrialVar = "firstHalfCueRep_f", # subselect data given selected variable 
				 selTrialVal = "second", # value of selected variable: first second
                                 segMessage = segMessage, beforeSamples = beforeSamples, afterSamples = afterSamples, 
                                 # timeLim = c(-1000, afterSamples),
                                 timeLim = c(1137, 2966), # significant difference in Go/NoGo (second half, baseline-uncorrected)
                                 # decimate = 10, # downsample or not
                                 isBaselineCorGroupGroup = F, isBaselineCorGroupSub = F, 
                                 tThresh = 2, nIter = 10000)

# END OF FILE.
