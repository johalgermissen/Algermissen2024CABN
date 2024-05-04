#!/usr/bin/env Rscript
# ============================================================================ #
## 03_mgngunc_gaze.R
## MGNGUncertainty plots of gaze.
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
# ============================================================================ #
#### Read in behavioral data: ####

# ----------------------------------- #
## Read in behavioral data:
data <- read_behavior()

# ----------------------------------- #
## Preprocessing: 
data <- wrapper_preprocessing(data)
names(data)
table(data$subject_n)

# ============================================================================ #
# ============================================================================ #
# ============================================================================ #
# ============================================================================ #
#### Time course plots of gaze: ####

source(paste0(dirs$codeDir, "functions/00_mgngunc_functions_timecourse.R")) # Load functions

## Settings:
segMessage <- "StartCue"
beforeSamples <- 866
afterSamples <- 1300

aggrMethod <- "distanceBaselineTrial"
isBaselineCor <- F

## Valence:
names(data)
zVar <- "valence_f"
zVar <- "valence_trial01_05_f" # valence only defined for cue repetitions 1-5
zVar <- "valence_trial06_10_f" # valence only defined for cue repetitions 6-10
zVar <- "valence_trial11_15_f" # valence only defined for cue repetitions 11-15
zVar <- "valence_trial11_16_f" # valence only defined for cue repetitions 11-16
table(data[, zVar], data$cueRep_n)

## Enforce x-axis limit:
xLim <- c(-100, 300)
xLim <- c(-100, 400)

p <- plot_timecourse(data = data, zVar = zVar, 
                     eyeMeasure = "gaze", aggrMetric = aggrMethod,
                     segMessage = segMessage, beforeSamples = beforeSamples, afterSamples = afterSamples,
                     xLim = xLim, suffix = paste0("_xLim", xLim[1], "-", xLim[2]), # zoom in on x-axis
                     # yLim = c(9.5, 17.5), # enforce y-axis limit
                     # yLim = c(0, 30), suffix = "_yLim00-30",
                     # yLim = c(10, 20), suffix = "_yLim10-20",
                     isBaselineCor = isBaselineCor, addLegend = F)

# ============================================================================ #
#### Permutation test: ####

source(paste0(dirs$codeDir, "functions/00_mgngunc_functions_timecourse.R")) # Load functions

## Settings:
segMessage <- "StartCue"
beforeSamples <- 866
afterSamples <- 1300

aggrMethod <- "distanceBaselineTrial"

permVar <- "valence_f"
permVar <- "valence_trial01_05_f" # valence only defined for cue repetitions 1-5
permVar <- "valence_trial06_10_f" # valence only defined for cue repetitions 6-10
permVar <- "valence_trial11_15_f" # valence only defined for cue repetitions 11-15

set.seed(20231202)
p <- test_permutation_timecourse(data = data, permVar = permVar, 
                                 eyeMeasure = "gaze", aggrMetric = aggrMethod,
                                 segMessage = segMessage, beforeSamples = beforeSamples, afterSamples = afterSamples, 
                                 # timeLim = c(200, 300),
                                 # decimate = 10, # downsample
                                 isBaselineCorGroup = F, isBaselineCorSub = F, 
                                 tThresh = 2, nIter = 1000)


# ============================================================================ #
# ============================================================================ #
# ============================================================================ #
# ============================================================================ #
#### Add trial-by-trial gaze data to behavioural data, plot bar plots: ####

segMessage <- "StartCue"
beforeSamples <- 866
afterSamples <- 1300

gazeMetric <- "distanceBaselineTrial"

aggrMethod <- "max"

timeWindow <- c(202, 278)

## Names of new variables:
newVarShort <- paste0(gazeMetric, "_", aggrMethod); newVarShort
newVarLong <- paste0(gazeMetric, "_", aggrMethod, "_", str_pad(timeWindow[1], width = 4, side = "left", pad = "0"), "-",
                     str_pad(timeWindow[2], width = 4, side = "left", pad = "0"), "ms");newVarLong

## Add gaze data:
source(paste0(dirs$codeDir, "functions/00_mgngunc_functions_regression.R")) # Load functions
data <- add_gaze_summary(data, gazeMetric = gazeMetric, aggrMethod = aggrMethod, timeWindow = timeWindow,
                         segMessage = segMessage, beforeSamples = beforeSamples, afterSamples = afterSamples)

# ----------------------------------- #
## Select data, standardize variables:

modData <- select_standardize(data)
plotData <- modData

# ============================================================================ #
#### Pretty plots: ####

yLimGaze <- c(10, 32); yName <- "Max" # Max 202-278

custom_barplot1(plotData, yVar = newVarLong, xVar = "valence_f", subVar = "subject_f", yLim = yLimGaze,
                yLab = paste0(yName, " dist. center (", timeWindow[1], "-", timeWindow[2], " ms)"),
                suffix = "pretty")

# ============================================================================ #
#### Fit linear regression to gaze data: ####

## Main effects:
formula <- "gaze_z ~ valence_f + (valence_f|subject_f)"

## Fit linear regression:
mod <- lmer(formula = formula, data = modData,
            control = lmerControl(optCtrl = list(maxfun = 1e+9), calc.derivs = F, optimizer = c("bobyqa")))
summary(mod)

## Compute p-values:
mod_LRT <- mixed(formula = formula, data = modData, method = "LRT", type = "III", # all_fit = T,
                 control = lmerControl(optCtrl = list(maxfun = 1e+9), calc.derivs = F, optimizer = c("bobyqa"))) # nloptwrap
anova(mod_LRT)

## Effect plots:
plot(effect("valence_f", mod))

## Coef plots:
library(MetBrewer)
nCoef <- length(fixef(mod)) - 1; nCoef
selCol <- met.brewer("Demuth", n = nCoef, type = "continuous"); colName <- "Demuth"
custom_coefplot(mod, plotSub = T, plotText = T, dropIntercept = T, revOrder = T, selCol = selCol, FTS = 10)

# END OF FILE.
