#!/usr/bin/env Rscript
# ============================================================================ #
## 06_mgngunc_questionnaires.R.
## MGNGUncertainty correlations between regression coefficients from behaviour and questionnaires.
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

source(paste0(dirs$codeDir, "functions/00_mgngunc_functions_regression.R")) # Load functions

# ============================================================================ #
# ============================================================================ #
# ============================================================================ #
#### Read in behavioral and pupil data: ####

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
segMessage <- "StartMask"
beforeSamples <- 1000
afterSamples <- 2966
data <- add_pupil(data, segMessage = segMessage, beforeSamples = beforeSamples, afterSamples = afterSamples)

# ----------------------------------- #
## Select data, standardize variables:

modData <- select_standardize(data)

# ============================================================================ #
# ============================================================================ #
#### 02) Read in questionnaires: ####

questData <- load_preprocess_questionnaires()

# ============================================================================ #
#### 03) Fit models and plot correlations: ####

# ----------------------------------------- #
## Specify formula:
formula <- "response_n ~ reqAction_f * valence_f + (reqAction_f * valence_f|subject_f)"; yLab <- "Pavlovian bias (choices)"; coefIdx <- 3
formula <- "RTcleaned_z ~ reqAction_f * valence_f + (reqAction_f * valence_f|subject_f)"; yLab <- "Pavlovian bias (RTs)"; coefIdx <- 3

formula <- "response_n ~ reqAction_f * valence_f * arousal_f + (reqAction_f * valence_f * arousal_f|subject_f)"; yLab <- "Effect of prime on choice"; coefIdx <- 4
formula <- "RTcleaned_z ~ reqAction_f * valence_f * arousal_f + (reqAction_f * valence_f * arousal_f|subject_f)"; yLab <- "Effect of prime on RTs"; coefIdx <- 4

formula <- "response_n ~ reqAction_f * valence_f * dilation_z + (reqAction_f * valence_f * dilation_z|subject_f)"; yLab <- "Association dilations and Go"; coefIdx <- 4
formula <- "RTcleaned_z ~ reqAction_f * valence_f * dilation_z + (reqAction_f * valence_f * dilation_z|subject_f)"; yLab <- "Association dilations and RTs"; coefIdx <- 4

## Load model:
mod <- fit_lmem(formula)

## Extract coefficients, add to questData:
cat(paste0("Extract term ", names(coef(mod)$subject_f)[coefIdx], "\n"))
questData$regCoef_n <- coef(mod)$subject_f[, coefIdx]

# ----------------------------------------- #
## Select questionnaire data:

names(questData)
xVar <- "stai_mean_n"; xLab <- "Trait anxiety (STAI)"
xVar <- "UPPSP_mean_n"; xLab <- "Impulsivity (UPPS-P)"
xVar <- "UPPSP_negative_urgency_n"; xLab <- "Negative urgency (UPPS-P)"
xVar <- "UPPSP_lack_perseverance_n"; xLab <- "Lack of perseverance (UPPS-P)"
xVar <- "UPPSP_lack_premeditation_n"; xLab <- "Lack of premeditation (UPPS-P)"
xVar <- "UPPSP_sensation_seeking_n"; xLab <- "Sensation seeking (UPPS-P)"
xVar <- "UPPSP_positive_urgency_n"; xLab <- "Positive urgency (UPPS-P)"

# ----------------------------------------- #
## Plot:

p <- plot_correlation(data = questData, xVar = xVar, yVar = "regCoef_n", xLab = xLab, yLab = yLab)

## Save:
plotName <- paste0("correlation_", dvName, "~", yLab, "_", xVar)
cat(paste0("Save as ", plotName, "\n"))
png(paste0(dirs$plotDir, plotName, ".png"), width = 480, height = 480)
print(p)
dev.off()

# END OF FILE.
