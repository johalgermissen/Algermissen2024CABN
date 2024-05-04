#!/usr/bin/env Rscript
# ============================================================================ #
## 02_mgngunc_plot.R
## MGNGUncertainty plots of behaviour.
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
segMessage <- "StartMask"; beforeSamples <- 1000; afterSamples <- 2966
# segMessage <- "StartResponse"; beforeSamples <- 1000; afterSamples <- 2600
data <- add_pupil(data, segMessage = segMessage, beforeSamples = beforeSamples, afterSamples = afterSamples)

# ----------------------------------- #
## Select data, standardize variables:

modData <- select_standardize(data)

# ============================================================================ #
#### DENSITY OF RTs AND DILATIONS: ####

## Select data:
plotData <- modData

t(stat.desc(data$RT_n))
p <- customplot_density2(plotData, xVar = "RT_n", zVar = "reqAction_f") 
p <- customplot_density2(plotData, xVar = "RT_n", zVar = "valence_f") 
p <- customplot_density2(plotData, xVar = "RT_n", zVar = "arousal_f") 

p <- customplot_density2(plotData, xVar = "RT_n", zVar = "reqAction_f", addLegend = T) 
p <- customplot_density2(plotData, xVar = "RT_n", zVar = "valence_f", addLegend = T) 
p <- customplot_density2(plotData, xVar = "RT_n", zVar = "arousal_f", addLegend = T) 

t(stat.desc(data$dilation_n))
p <- customplot_density2(plotData, xVar = "dilation_n", zVar = "response_f") 
p <- customplot_density2(plotData, xVar = "dilation_n", zVar = "reqAction_f") 
p <- customplot_density2(plotData, xVar = "dilation_n", zVar = "valence_f") 
p <- customplot_density2(plotData, xVar = "dilation_n", zVar = "arousal_f") 


# ============================================================================ #
#### CUSTOM BAR PLOTS averaged over time: ####

# ============================================================================ #
#### Bar plots RESPONSE: ####

plotData <- modData

## 1D bar plot: Response per block:
custom_barplot1(plotData, xVar = "block_f", yVar = "response_n", subVar = "subject_f")
custom_barplot1(plotData, xVar = "cueRep_f", yVar = "response_n", subVar = "subject_f")

## 2D bar plot: Response per required action per required action/ valence/ arousal:
custom_barplot2(plotData, xVar = "block_f", yVar = "response_n", zVar = "reqAction_f", subVar = "subject_f")
custom_barplot2(plotData, xVar = "block_f", yVar = "response_n", zVar = "valence_f", subVar = "subject_f")
custom_barplot2(plotData, xVar = "block_f", yVar = "response_n", zVar = "arousal_f", subVar = "subject_f")

custom_barplot2(plotData, xVar = "cueRep_f", yVar = "response_n", zVar = "reqAction_f", subVar = "subject_f", isPoint = F, isConnect = F)
custom_barplot2(plotData, xVar = "cueRep_f", yVar = "response_n", zVar = "valence_f", subVar = "subject_f", isPoint = F, isConnect = F)
custom_barplot2(plotData, xVar = "cueRep_f", yVar = "response_n", zVar = "arousal_f", subVar = "subject_f", isPoint = F, isConnect = F)

## 1D bar plot: Response per required action:
custom_barplot1(plotData, xVar = "reqAction_f", yVar = "response_n", subVar = "subject_f")

## 1D bar plot: Response per valence:
custom_barplot1(plotData, xVar = "valence_f", yVar = "response_n", subVar = "subject_f")

## 2D bar plot: Response per required action per valence:
custom_barplot2(plotData, xVar = "valence_f", yVar = "response_n", zVar = "reqAction_f", subVar = "subject_f")
custom_barplot2(plotData, xVar = "reqAction_f", yVar = "response_n", zVar = "valence_f", subVar = "subject_f")

## 1D bar plot: Response per condition:
custom_barplot1(plotData, xVar = "condition_f", yVar = "response_n", subVar = "subject_f")

## 1D bar plot: Response per arousal:
custom_barplot1(modData, xVar = "arousal_f", yVar = "response_n", subVar = "subject_f")

## 2D bar plot: Response per arousal per required action/valence:
custom_barplot2(plotData, xVar = "arousal_f", yVar = "response_n", zVar = "reqAction_f", subVar = "subject_f")
custom_barplot2(plotData, xVar = "reqAction_f", yVar = "response_n", zVar = "arousal_f", subVar = "subject_f")
custom_barplot2(plotData, xVar = "arousal_f", yVar = "response_n", zVar = "valence_f", subVar = "subject_f")
custom_barplot2(plotData, xVar = "valence_f", yVar = "response_n", zVar = "arousal_f", subVar = "subject_f")
custom_barplot2(plotData, xVar = "condition_f", yVar = "response_n", zVar = "arousal_f", subVar = "subject_f")
custom_barplot2(plotData, xVar = "condition_short1_f", yVar = "response_n", zVar = "arousal_f", subVar = "subject_f")

## 1D bar plot: Response per congruency (required action):
custom_barplot1(plotData, xVar = "reqCongruency_f", yVar = "response_n", subVar = "subject_f")

## 2D bar plot: Response per congruency (required action) per arousal:
custom_barplot2(plotData, xVar = "reqCongruency_short1_f", yVar = "response_n", zVar = "arousal_f", subVar = "subject_f")

## 1D bar plot: Response per outcome last trial:
custom_barplot1(plotData, xVar = "outcome_last_f", yVar = "response_n", subVar = "subject_f")

## 2D bar plot: Response per required action per outcome last trial:
custom_barplot2(plotData, xVar = "outcome_last_short_f", yVar = "response_n", zVar = "reqAction_f", subVar = "subject_f")
custom_barplot2(plotData, xVar = "outcome_last_short_f", yVar = "response_n", zVar = "valence_f", subVar = "subject_f")
custom_barplot2(plotData, xVar = "reqAction_f", yVar = "response_n", zVar = "outcome_last_short_f", subVar = "subject_f")
custom_barplot2(plotData, xVar = "valence_f", yVar = "response_n", zVar = "outcome_last_short_f", subVar = "subject_f")

## 3-way interaction: response per required action x valence x arousal:
plotData <- subset(modData, !is.na(response_n) & !is.na(reqAction_f) & !is.na(valence_f) & !is.na(arousal_f))
ggplot(plotData, aes(y = response_n, x = reqAction_f, fill = valence_f)) + 
  stat_summary(fun = mean, geom = "bar", aes(fill = valence_f), position = "dodge") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", position = position_dodge(width=0.9), width = .1) +
  facet_wrap(vars(arousal_f))
custom_barplot3(data = modData, yVar = "response_n", xVar = "reqAction_f", zVar = "valence_f", splitVar = "arousal_f", subVar = "subject_f")
custom_barplot3(data = modData, yVar = "response_n", xVar = "reqAction_f", zVar = "arousal_f", splitVar = "valence_f", subVar = "subject_f")
custom_barplot3(data = modData, yVar = "response_n", xVar = "valence_f", zVar = "arousal_f", splitVar = "reqAction_f", subVar = "subject_f")

# ============================================================================ #
#### Bar plots ACCURACY: ####

plotData <- modData

## 1D bar plot: Accuracy per block:
custom_barplot1(plotData, xVar = "block_f", yVar = "ACC_n", subVar = "subject_f")
custom_barplot1(plotData, xVar = "cueRep_f", yVar = "ACC_n", subVar = "subject_f")

## 2D bar plot: Response per required action per required action/ valence/ arousal:
custom_barplot2(plotData, xVar = "block_f", yVar = "ACC_n", zVar = "reqAction_f", subVar = "subject_f")
custom_barplot2(plotData, xVar = "block_f", yVar = "ACC_n", zVar = "valence_f", subVar = "subject_f")
custom_barplot2(plotData, xVar = "block_f", yVar = "ACC_n", zVar = "arousal_f", subVar = "subject_f")
custom_barplot2(plotData, xVar = "block_f", yVar = "ACC_n", zVar = "reqCongruency_short2_f", subVar = "subject_f")

custom_barplot2(plotData, xVar = "cueRep_f", yVar = "ACC_n", zVar = "reqAction_f", subVar = "subject_f", isPoint = F, isConnect = F)
custom_barplot2(plotData, xVar = "cueRep_f", yVar = "ACC_n", zVar = "valence_f", subVar = "subject_f", isPoint = F, isConnect = F)
custom_barplot2(plotData, xVar = "cueRep_f", yVar = "ACC_n", zVar = "arousal_f", subVar = "subject_f", isPoint = F, isConnect = F)

## 1D bar plot: Accuracy per required action:
custom_barplot1(plotData, xVar = "reqAction_f", yVar = "ACC_n", subVar = "subject_f")

## 1D bar plot: Accuracy per valence:
custom_barplot1(plotData, xVar = "valence_f", yVar = "ACC_n", subVar = "subject_f")

## 2D bar plot: Accuracy per required action per valence:
custom_barplot2(plotData, xVar = "reqAction_f", yVar = "ACC_n", zVar = "valence_f", subVar = "subject_f")
custom_barplot2(plotData, xVar = "valence_f", yVar = "ACC_n", zVar = "reqAction_f", subVar = "subject_f")

## 1D bar plot: Accuracy per condition:
custom_barplot1(plotData, xVar = "condition_f", yVar = "response_n", subVar = "subject_f")

## 1D bar plot: Accuracy per arousal:
custom_barplot1(plotData, xVar = "arousal_f", yVar = "ACC_n", subVar = "subject_f")

## 2D bar plot: Accuracy per arousal per required action/ valence/ congruency/ condition:
custom_barplot2(plotData, xVar = "arousal_f", yVar = "ACC_n", zVar = "reqAction_f", subVar = "subject_f")
custom_barplot2(plotData, xVar = "reqAction_f", yVar = "ACC_n", zVar = "arousal_f", subVar = "subject_f")
custom_barplot2(plotData, xVar = "arousal_f", yVar = "ACC_n", zVar = "valence_f", subVar = "subject_f")
custom_barplot2(plotData, xVar = "valence_f", yVar = "ACC_n", zVar = "arousal_f", subVar = "subject_f")

custom_barplot2(plotData, xVar = "arousal_f", yVar = "ACC_n", zVar = "reqCongruency_short2_f", subVar = "subject_f")
custom_barplot2(plotData, xVar = "reqCongruency_short1_f", yVar = "ACC_n", zVar = "arousal_f", subVar = "subject_f")
custom_barplot2(plotData, xVar = "condition_short1_f", yVar = "ACC_n", zVar = "arousal_f", subVar = "subject_f") # , FTS = 28)

## 1D bar plot: Accuracy per congruency (required action):
custom_barplot1(plotData, xVar = "reqCongruency_f", yVar = "ACC_n", subVar = "subject_f")

## 3-way interaction: Accuracy per required action x valence x arousal:
plotData <- subset(modData, !is.na(ACC_n) & !is.na(reqAction_f) & !is.na(valence_f) & !is.na(arousal_f))
ggplot(plotData, aes(y = ACC_n, x = reqAction_f, fill = valence_f)) + 
  stat_summary(fun = mean, geom = "bar", aes(fill = valence_f), position = "dodge") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", position = position_dodge(width=0.9), width = .1) +
  facet_wrap(vars(arousal_f))
custom_barplot3(data = modData, yVar = "ACC_n", xVar = "reqAction_f", zVar = "valence_f", splitVar = "arousal_f", subVar = "subject_f")

# ============================================================================ #
#### Bar plots RTs (raw): ####

plotData <- modData

# yLimRT <- c(0, 1.0)
yLimRT <- c(0.4, 1.0)

## 1D bar plot: RT per block:
custom_barplot1(plotData, xVar = "block_f", yVar = "RTcleaned_n", subVar = "subject_f", yLim = yLimRT)
custom_barplot1(plotData, xVar = "cueRep_f", yVar = "RTcleaned_n", subVar = "subject_f", yLim = yLimRT)

## 2D bar plot: RT per congruency (required action) per block:
custom_barplot2(plotData, xVar = "block_f", yVar = "RTcleaned_n", zVar = "reqAction_f", subVar = "subject_f", yLim = yLimRT)
custom_barplot2(plotData, xVar = "block_f", yVar = "RTcleaned_n", zVar = "valence_f", subVar = "subject_f", yLim = yLimRT)
custom_barplot2(plotData, xVar = "block_f", yVar = "RTcleaned_n", zVar = "reqCongruency_short2_f", subVar = "subject_f", yLim = yLimRT)

custom_barplot2(plotData, xVar = "cueRep_f", yVar = "RTcleaned_n", zVar = "reqAction_f", subVar = "subject_f", isPoint = F, isConnect = F)
custom_barplot2(plotData, xVar = "cueRep_f", yVar = "RTcleaned_n", zVar = "valence_f", subVar = "subject_f", isPoint = F, isConnect = F)
custom_barplot2(plotData, xVar = "cueRep_f", yVar = "RTcleaned_n", zVar = "arousal_f", subVar = "subject_f", isPoint = F, isConnect = F)

## 1D bar plot: RT per required action:
custom_barplot1(plotData, xVar = "reqAction_f", yVar = "RTcleaned_n", subVar = "subject_f", yLim = yLimRT)

## 1D bar plot: RT per valence:
custom_barplot1(plotData, xVar = "valence_f", yVar = "RTcleaned_n", subVar = "subject_f", yLim = yLimRT)

## 2D bar plot: RT per valence per required action:
custom_barplot2(plotData, xVar = "valence_f", yVar = "RTcleaned_n", zVar = "reqAction_f", subVar = "subject_f", yLim = yLimRT)
custom_barplot2(plotData, xVar = "reqAction_f", yVar = "RTcleaned_n", zVar = "valence_f", subVar = "subject_f", yLim = yLimRT)

## 1D bar plot: RT per cue condition:
custom_barplot1(modData, xVar = "condition_f", yVar = "RTcleaned_n", subVar = "subject_f", yLim = yLimRT)

## 1D bar plot: RT per arousal:
custom_barplot1(plotData, xVar = "arousal_f", yVar = "RTcleaned_n", subVar = "subject_f", yLim = yLimRT)

## 2D bar plot: RT per required action per arousal:
custom_barplot2(plotData, xVar = "reqAction_f", yVar = "RTcleaned_n", zVar = "arousal_f", subVar = "subject_f", yLim = yLimRT)
custom_barplot2(plotData, xVar = "arousal_f", yVar = "RTcleaned_n", zVar = "reqAction_f", subVar = "subject_f", yLim = yLimRT)
custom_barplot2(plotData, xVar = "valence_f", yVar = "RTcleaned_n", zVar = "arousal_f", subVar = "subject_f", yLim = yLimRT)
custom_barplot2(plotData, xVar = "arousal_f", yVar = "RTcleaned_n", zVar = "valence_f", subVar = "subject_f", yLim = yLimRT)

custom_barplot2(plotData, xVar = "reqCongruency_short1_f", yVar = "RTcleaned_n", zVar = "arousal_f", subVar = "subject_f", yLim = yLimRT)
custom_barplot2(plotData, xVar = "arousal_f", yVar = "RTcleaned_n", zVar = "reqCongruency_short2_f", subVar = "subject_f", yLim = yLimRT)
custom_barplot2(plotData, xVar = "condition_short1_f", yVar = "RTcleaned_n", zVar = "arousal_f", subVar = "subject_f", yLim = yLimRT) # , FTS = 28)

## 1D bar plot: RT per congruency (required action):
custom_barplot1(modData, xVar = "reqCongruency_f", yVar = "RTcleaned_n", subVar = "subject_f", yLim = yLimRT)

# -------------------------------------------- #
## 1D bar plot: RT per outcome last trial:

custom_barplot1(modData, xVar = "outcome_last_f", yVar = "RTcleaned_n", subVar = "subject_f", yLim = yLimRT)
custom_barplot1(modData, xVar = "outcome_last_rel_f", yVar = "RTcleaned_n", subVar = "subject_f", yLim = yLimRT)

## 2D bar plot: RT per required action per outcome last trial:
custom_barplot2(modData, xVar = "outcome_last_f", yVar = "RTcleaned_n", zVar = "reqAction_f", subVar = "subject_f", yLim = yLimRT)
custom_barplot2(modData, xVar = "reqAction_f", yVar = "RTcleaned_n", zVar = "outcome_last_short_f", subVar = "subject_f", yLim = yLimRT)

## 2D bar plot: RT per valence per outcome last trial:
# custom_barplot2(modData, xVar = "outcome_last_short_f", yVar = "RTcleaned_n", zVar = "valence_f", subVar = "subject_f", yLim = yLimRT)
# custom_barplot2(modData, xVar = "valence_f", yVar = "RTcleaned_n", zVar = "outcome_last_short_f", subVar = "subject_f", yLim = yLimRT)

## 2D bar plot: RT per congruency (required action) per outcome last trial:
custom_barplot2(modData, xVar = "outcome_last_short_f", yVar = "RTcleaned_n", zVar = "reqCongruency_short2_f", subVar = "subject_f", yLim = yLimRT)
custom_barplot2(modData, xVar = "reqCongruency_short1_f", yVar = "RTcleaned_n", zVar = "outcome_last_short_f", subVar = "subject_f", yLim = yLimRT)

custom_barplot2(modData, xVar = "outcome_last_rel_short_f", yVar = "RTcleaned_n", zVar = "reqCongruency_short2_f", subVar = "subject_f", yLim = yLimRT)
custom_barplot2(modData, xVar = "reqCongruency_short1_f", yVar = "RTcleaned_n", zVar = "outcome_last_rel_short_f", subVar = "subject_f", yLim = yLimRT)

## 2D bar plot: RT per repeat/switch per outcome last trial:
custom_barplot2(modData, xVar = "outcome_last_short_f", yVar = "RTcleaned_n", zVar = "repeat_f", subVar = "subject_f", yLim = yLimRT)
custom_barplot2(modData, xVar = "repeat_f", yVar = "RTcleaned_n", zVar = "outcome_last_short_f", subVar = "subject_f", yLim = yLimRT)

## 3-way interaction: RT per required action x valence x arousal:
plotData <- subset(modData, !is.na(RTcleaned_n) & !is.na(reqAction_f) & !is.na(valence_f) & !is.na(arousal_f))
ggplot(plotData, aes(y = RTcleaned_n, x = reqAction_f, fill = valence_f)) + 
  stat_summary(fun = mean, geom = "bar", aes(fill = valence_f), position = "dodge") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", position = position_dodge(width=0.9), width = .1) +
  facet_wrap(vars(arousal_f))
custom_barplot3(data = modData, yVar = "RTcleaned_n", xVar = "reqAction_f", zVar = "valence_f", splitVar = "arousal_f", subVar = "subject_f") # , yLim = yLimRT)
custom_barplot3(data = modData, yVar = "RTcleaned_n", xVar = "reqAction_f", zVar = "arousal_f", splitVar = "valence_f", subVar = "subject_f") # , yLim = yLimRT)

# ============================================================================ #
#### Bar plots log(RTs): ####

yLimRTLog <- c(4, 7)

## 1D bar plot: log(RT) per emotion:
custom_barplot1(modData, xVar = "valence_f", yVar = "RT_ms_log_n", subVar = "subject_f", yLim = yLimRTLog)

## 1D bar plot: log(RT) per valence:
custom_barplot1(modData, xVar = "valence_f", yVar = "RT_ms_log_n", subVar = "subject_f", yLim = yLimRTLog)

## 1D bar plot: log(RT) per required action:
custom_barplot1(modData, xVar = "reqAction_f", yVar = "RT_ms_log_n", subVar = "subject_f", yLim = yLimRTLog)

## 1D bar plot: log(RT) per congruency (required action):
custom_barplot1(modData, xVar = "reqCongruency_f", yVar = "RT_ms_log_n", subVar = "subject_f", yLim = yLimRTLog)

## 1D bar plot: log(RT) per congruency (actual response):
custom_barplot1(modData, xVar = "actCongruency_f", yVar = "RT_ms_log_n", subVar = "subject_f", yLim = yLimRTLog)

## 2D bar plot: log(RT) per required action per emotion:
custom_barplot2(modData, xVar = "valence_f", yVar = "RT_ms_log_n", zVar = "reqAction_f", subVar = "subject_f", yLim = yLimRTLog)

# ============================================================================ #
#### Bar plots BASELINE: ####

plotData <- modData

yLimBase <- c(-10, 10)

## S06:
## 1D bar plot: Baseline per response/ Accuracy/ repeat/ RT:
yLimBase <- c(0, 29)
plotData$baseline12_n <- plotData$baseline_n + 12
custom_barplot1(plotData, yVar = "baseline12_n", xVar = "response_f", subVar = "subject_f", yLim = yLimBase)
custom_barplot1(plotData, yVar = "baseline12_n", xVar = "ACC_f", subVar = "subject_f", yLim = yLimBase)
custom_barplot1(plotData, yVar = "baseline12_n", xVar = "repeat_f", subVar = "subject_f", yLim = yLimBase)
custom_barplot1(plotData, yVar = "baseline12_n", xVar = "RTcleaned_fast_f", subVar = "subject_f", yLim = yLimBase)
custom_barplot2(plotData, xVar = "cueRep_f", yVar = "baseline_n", zVar = "response_f", subVar = "subject_f", isPoint = F, isConnect = F, yLim = c(-25, 25))
custom_barplot2(plotData, xVar = "cueRep_f", yVar = "baseline_n", zVar = "ACC_f", subVar = "subject_f", isPoint = F, isConnect = F, yLim = c(-25, 25))
custom_barplot2(plotData, xVar = "cueRep_f", yVar = "baseline_n", zVar = "RTcleaned_fast_f", subVar = "subject_f", isPoint = F, isConnect = F, yLim = c(-25, 25))

## 1D bar plot: Baseline per block:
custom_barplot1(plotData, xVar = "block_f", yVar = "baseline_n", subVar = "subject_f", yLim = yLimBase)
custom_barplot1(plotData, xVar = "cueRep_f", yVar = "baseline_n", subVar = "subject_f", yLim = c(-25, 25))
custom_barplot1(plotData, xVar = "trialnr_f", yVar = "baseline_n", subVar = "subject_f", isPoint = F, isConnect = F, yLim = c(-25, 25))
custom_barplot1(plotData, xVar = "trialnr_block_f", yVar = "baseline_n", subVar = "subject_f", isPoint = F, isConnect = F, yLim = c(-25, 25))

## 2D bar plot: Response per required action per required action/ valence/ arousal:
custom_barplot2(plotData, xVar = "block_f", yVar = "baseline_n", zVar = "reqAction_f", subVar = "subject_f", yLim = yLimBase)
custom_barplot2(plotData, xVar = "block_f", yVar = "baseline_n", zVar = "valence_f", subVar = "subject_f", yLim = yLimBase)
custom_barplot2(plotData, xVar = "block_f", yVar = "baseline_n", zVar = "arousal_f", subVar = "subject_f", yLim = yLimBase)
custom_barplot2(plotData, xVar = "block_f", yVar = "baseline_n", zVar = "reqCongruency_short2_f", subVar = "subject_f", yLim = yLimBase)
custom_barplot2(plotData, xVar = "block_f", yVar = "baseline_n", zVar = "response_f", subVar = "subject_f", yLim = yLimBase)
custom_barplot2(plotData, xVar = "block_f", yVar = "baseline_n", zVar = "ACC_f", subVar = "subject_f", yLim = yLimBase)

custom_barplot2(plotData, xVar = "cueRep_f", yVar = "baseline_n", zVar = "reqAction_f", subVar = "subject_f", isPoint = F, isConnect = F, yLim = c(-25, 25))
custom_barplot2(plotData, xVar = "cueRep_f", yVar = "baseline_n", zVar = "valence_f", subVar = "subject_f", isPoint = F, isConnect = F, yLim = c(-25, 25))
custom_barplot2(plotData, xVar = "cueRep_f", yVar = "baseline_n", zVar = "arousal_f", subVar = "subject_f", isPoint = F, isConnect = F, yLim = c(-25, 25))
custom_barplot2(plotData, xVar = "cueRep_f", yVar = "baseline_n", zVar = "response_f", subVar = "subject_f", isPoint = F, isConnect = F, yLim = c(-25, 25))
custom_barplot2(plotData, xVar = "cueRep_f", yVar = "baseline_n", zVar = "ACC_f", subVar = "subject_f", isPoint = F, isConnect = F, yLim = c(-25, 25))
custom_barplot2(plotData, xVar = "cueRep_f", yVar = "baseline_n", zVar = "RTcleaned_fast_f", subVar = "subject_f", isPoint = F, isConnect = F, yLim = c(-25, 25))

## 1D bar plot: Baseline per required action:
custom_barplot1(plotData, xVar = "reqAction_f", yVar = "baseline_n", subVar = "subject_f", yLim = yLimBase)
custom_barplot1(plotData, xVar = "response_f", yVar = "baseline_n", subVar = "subject_f", yLim = yLimBase)
custom_barplot1(plotData, xVar = "ACC_f", yVar = "baseline_n", subVar = "subject_f", yLim = yLimBase)
custom_barplot1(plotData, xVar = "repeat_f", yVar = "baseline_n", subVar = "subject_f", yLim = yLimBase)

## 1D bar plot: Baseline per valence:
custom_barplot1(plotData, xVar = "valence_f", yVar = "baseline_n", subVar = "subject_f", yLim = yLimBase)

## 2D bar plot: Baseline per required action per valence:
custom_barplot2(plotData, xVar = "reqAction_f", yVar = "baseline_n", zVar = "valence_f", subVar = "subject_f", yLim = yLimBase)
custom_barplot2(plotData, xVar = "valence_f", yVar = "baseline_n", zVar = "reqAction_f", subVar = "subject_f", yLim = yLimBase)

## 1D bar plot: Baseline per condition:
custom_barplot1(plotData, xVar = "condition_f", yVar = "baseline_n", subVar = "subject_f", yLim = yLimBase)
custom_barplot1(plotData, xVar = "respCond_f", yVar = "baseline_n", subVar = "subject_f", yLim = yLimBase)

## 1D bar plot: Baseline per arousal:
custom_barplot1(plotData, xVar = "arousal_f", yVar = "baseline_n", subVar = "subject_f", yLim = yLimBase)

## 2D bar plot: Baseline per arousal per required action/ valence/ congruency/ condition:
custom_barplot2(plotData, xVar = "arousal_f", yVar = "baseline_n", zVar = "reqAction_f", subVar = "subject_f", yLim = yLimBase)
custom_barplot2(plotData, xVar = "reqAction_f", yVar = "baseline_n", zVar = "arousal_f", subVar = "subject_f", yLim = yLimBase)
custom_barplot2(plotData, xVar = "arousal_f", yVar = "baseline_n", zVar = "valence_f", subVar = "subject_f", yLim = yLimBase)
custom_barplot2(plotData, xVar = "valence_f", yVar = "baseline_n", zVar = "arousal_f", subVar = "subject_f", yLim = yLimBase)

custom_barplot2(plotData, xVar = "arousal_f", yVar = "baseline_n", zVar = "reqCongruency_short2_f", subVar = "subject_f", yLim = yLimBase)
custom_barplot2(plotData, xVar = "reqCongruency_short1_f", yVar = "baseline_n", zVar = "arousal_f", subVar = "subject_f", yLim = yLimBase)
custom_barplot2(plotData, xVar = "condition_short1_f", yVar = "baseline_n", zVar = "arousal_f", subVar = "subject_f", yLim = yLimBase) # , FTS = 28)

## 1D bar plot: Baseline per congruency (required action):
custom_barplot1(plotData, xVar = "reqCongruency_f", yVar = "baseline_n", subVar = "subject_f", yLim = yLimBase)

## 1D bar plot: Baseline per response/ Accuracy/ repeat/ RT:
yLimBase <- c(0, 29)
plotData$baseline12_n <- plotData$baseline_n + 12
custom_barplot1(plotData, yVar = "baseline12_n", xVar = "response_f", subVar = "subject_f", yLim = yLimBase)
custom_barplot1(plotData, yVar = "baseline12_n", xVar = "ACC_f", subVar = "subject_f", yLim = yLimBase)
custom_barplot1(plotData, yVar = "baseline12_n", xVar = "repeat_f", subVar = "subject_f", yLim = yLimBase)
custom_barplot1(plotData, yVar = "baseline12_n", xVar = "RTcleaned_fast_f", subVar = "subject_f", yLim = yLimBase)

## 1D bar plot: Baseline per outcome last trial:
custom_barplot1(plotData, xVar = "outcome_last_f", yVar = "baseline_n", subVar = "subject_f", yLim = yLimBase)
custom_barplot1(plotData, xVar = "outcome_last_rel_f", yVar = "baseline_n", subVar = "subject_f", yLim = yLimBase)
# custom_barplot1(plotData, xVar = "outcome_last_all_f", yVar = "baseline_n", subVar = "subject_f", yLim = yLimBase)

## 3-way interaction: Baseline per required action x valence x arousal:
plotData <- subset(modData, !is.na(baseline_n) & !is.na(reqAction_f) & !is.na(valence_f) & !is.na(arousal_f))
ggplot(plotData, aes(y = baseline_n, x = reqAction_f, fill = valence_f)) + 
  stat_summary(fun = mean, geom = "bar", aes(fill = valence_f), position = "dodge") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", position = position_dodge(width=0.9), width = .1) +
  facet_wrap(vars(arousal_f))
custom_barplot3(data = plotData, yVar = "baseline_n", xVar = "reqAction_f", zVar = "valence_f", splitVar = "arousal_f", subVar = "subject_f", yLim = yLimBase)

## 3-way interaction: Baseline per response x valence x arousal:
plotData <- subset(modData, !is.na(baseline_n) & !is.na(response_f) & !is.na(valence_f) & !is.na(arousal_f))
ggplot(plotData, aes(y = baseline_n, x = response_f, fill = valence_f)) + 
  stat_summary(fun = mean, geom = "bar", aes(fill = valence_f), position = "dodge") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", position = position_dodge(width=0.9), width = .1) +
  facet_wrap(vars(arousal_f))
custom_barplot3(data = plotData, yVar = "baseline_n", xVar = "response_f", zVar = "valence_f", splitVar = "arousal_f", subVar = "subject_f", yLim = yLimBase)

# ============================================================================ #
#### Bar plots DILATION: ####

plotData <- modData

# yLimDil <- c(5, 25)
yLimDil <- c(6, 28)

# ------------------------------------------- #
### S05:

## 1D plots:
custom_barplot1(plotData, xVar = "arousal_f", yVar = "dilation_n", subVar = "subject_f", yLim = yLimDil)
custom_barplot1(plotData, xVar = "arousal_f", yVar = "dilation_n", subVar = "subject_f", yLim = yLimDil)
custom_barplot1(plotData, xVar = "ACC_f", yVar = "dilation_n", subVar = "subject_f", yLim = yLimDil)
custom_barplot1(plotData, xVar = "RTcleaned_fast_f", yVar = "dilation_n", subVar = "subject_f", yLim = yLimDil)
custom_barplot1(plotData, xVar = "repeat_f", yVar = "dilation_n", subVar = "subject_f", yLim = yLimDil)

## 2D plots:
custom_barplot2(plotData, xVar = "ACC_f", yVar = "dilation_n", zVar = "response_f", subVar = "subject_f", yLim = yLimDil)
custom_barplot2(plotData, xVar = "repeat_f", yVar = "dilation_n", zVar = "response_f", subVar = "subject_f", yLim = yLimDil)
custom_barplot2(plotData, xVar = "response_f", yVar = "dilation_n", zVar = "ACC_f", subVar = "subject_f", yLim = yLimDil)
custom_barplot2(plotData, xVar = "response_f", yVar = "dilation_n", zVar = "ACC_short_f", subVar = "subject_f", yLim = yLimDil)
custom_barplot2(plotData, xVar = "response_f", yVar = "dilation_n", zVar = "repeat_f", subVar = "subject_f", yLim = yLimDil)
custom_barplot2(plotData, xVar = "response_f", yVar = "dilation_n", zVar = "repeat_short_f", subVar = "subject_f", yLim = yLimDil)

custom_barplot2(plotData, xVar = "RTcleaned_fast_f", yVar = "dilation_n", zVar = "ACC_f", subVar = "subject_f", yLim = yLimDil)
custom_barplot2(plotData, xVar = "RTcleaned_fast_f", yVar = "dilation_n", zVar = "ACC_short_f", subVar = "subject_f", yLim = yLimDil)
custom_barplot2(plotData, xVar = "ACC_f", yVar = "dilation_n", zVar = "valence_f", subVar = "subject_f", yLim = yLimDil)
custom_barplot2(plotData, xVar = "RTcleaned_fast_f", yVar = "dilation_n", zVar = "valence_f", subVar = "subject_f", yLim = yLimDil)
custom_barplot2(plotData, xVar = "repeat_f", yVar = "dilation_n", zVar = "valence_f", subVar = "subject_f", yLim = yLimDil)

## 2D plots on Go responses only:
# yLimDil <- c(6, 28)
yLimDil <- c(4, 36)
custom_barplot2(subset(plotData, response_f == "Go"), xVar = "ACC_f", yVar = "dilation_n", zVar = "valence_f", subVar = "subject_f", yLim = yLimDil,
                suffix = paste0("_respGo_yLim", yLimDil[1], "-", yLimDil[2]))
custom_barplot2(subset(plotData, response_f == "Go"), xVar = "RTcleaned_fast_f", yVar = "dilation_n", zVar = "valence_f", subVar = "subject_f", yLim = yLimDil,
                suffix = paste0("_respGo_yLim", yLimDil[1], "-", yLimDil[2]))
custom_barplot2(subset(plotData, response_f == "Go"), xVar = "repeat_f", yVar = "dilation_n", zVar = "valence_f", subVar = "subject_f", yLim = yLimDil,
                suffix = paste0("_respGo_yLim", yLimDil[1], "-", yLimDil[2]))
custom_barplot2(subset(plotData, response_f == "Go"), xVar = "RTcleaned_fast_f", yVar = "dilation_n", zVar = "ACC_short_f", subVar = "subject_f", yLim = yLimDil,
                suffix = paste0("_respGo_yLim", yLimDil[1], "-", yLimDil[2]))

## 3D plots:
custom_barplot3(data = plotData, yVar = "dilation_n", xVar = "response_f", zVar = "valence_f", splitVar = "ACC_f", subVar = "subject_f", yLim = yLimDil)
custom_barplot3(data = plotData, yVar = "dilation_n", xVar = "response_f", zVar = "valence_f", splitVar = "RTcleaned_fast_f", subVar = "subject_f", yLim = yLimDil)
custom_barplot3(data = plotData, yVar = "dilation_n", xVar = "response_f", zVar = "valence_f", splitVar = "repeat_f", subVar = "subject_f", yLim = yLimDil)

# ------------------------------------------- #
## 1D bar plot: Dilation per block:
custom_barplot1(plotData, xVar = "block_f", yVar = "dilation_n", subVar = "subject_f", yLim = yLimDil)
custom_barplot1(plotData, xVar = "cueRep_f", yVar = "dilation_n", subVar = "subject_f", yLim = yLimDil)
custom_barplot1(plotData, xVar = "trialnr_block_f", yVar = "dilation_n", subVar = "subject_f", isPoint = F, isConnect = F, yLim = yLimDil)
custom_barplot1(plotData, xVar = "trialnr_f", yVar = "dilation_n", subVar = "subject_f", isPoint = F, isConnect = F, yLim = yLimDil)

## 2D bar plot: Dilation per block per required action/ valence/ arousal:
custom_barplot2(plotData, xVar = "block_f", yVar = "dilation_n", zVar = "reqAction_f", subVar = "subject_f", yLim = yLimDil)
custom_barplot2(plotData, xVar = "block_f", yVar = "dilation_n", zVar = "valence_f", subVar = "subject_f", yLim = yLimDil)
custom_barplot2(plotData, xVar = "block_f", yVar = "dilation_n", zVar = "arousal_f", subVar = "subject_f", yLim = yLimDil)
custom_barplot2(plotData, xVar = "block_f", yVar = "dilation_n", zVar = "reqCongruency_short2_f", subVar = "subject_f", yLim = yLimDil)
custom_barplot2(plotData, xVar = "block_f", yVar = "dilation_n", zVar = "response_f", subVar = "subject_f", yLim = yLimDil)
custom_barplot2(plotData, xVar = "block_f", yVar = "dilation_n", zVar = "ACC_f", subVar = "subject_f", yLim = yLimDil)

## 2D bar plot: Dilation per cue repetition per required action/ valence/ arousal:
custom_barplot2(plotData, xVar = "cueRep_f", yVar = "dilation_n", zVar = "reqAction_f", subVar = "subject_f", isPoint = F, isConnect = F, yLim = yLimDil)
custom_barplot2(plotData, xVar = "cueRep_f", yVar = "dilation_n", zVar = "valence_f", subVar = "subject_f", isPoint = F, isConnect = F, yLim = yLimDil)
custom_barplot2(plotData, xVar = "cueRep_f", yVar = "dilation_n", zVar = "arousal_f", subVar = "subject_f", isPoint = F, isConnect = F, yLim = yLimDil)
custom_barplot2(plotData, xVar = "cueRep_f", yVar = "dilation_n", zVar = "response_f", subVar = "subject_f", isPoint = F, isConnect = F, yLim = yLimDil)
custom_barplot2(plotData, xVar = "cueRep_f", yVar = "dilation_n", zVar = "ACC_f", subVar = "subject_f", isPoint = F, isConnect = F, yLim = yLimDil)
custom_barplot2(plotData, xVar = "cueRep_f", yVar = "dilation_n", zVar = "repeat_f", subVar = "subject_f", isPoint = F, isConnect = F, yLim = yLimDil)

## 1D bar plot: Dilation per required action/ response/ ACC/ repeat:
custom_barplot1(plotData, xVar = "reqAction_f", yVar = "dilation_n", subVar = "subject_f", yLim = yLimDil)
custom_barplot1(plotData, xVar = "response_f", yVar = "dilation_n", subVar = "subject_f", yLim = yLimDil)
custom_barplot1(plotData, xVar = "ACC_f", yVar = "dilation_n", subVar = "subject_f", yLim = yLimDil)
custom_barplot1(plotData, xVar = "repeat_f", yVar = "dilation_n", subVar = "subject_f", yLim = yLimDil)

## 1D bar plot: Dilation per valence:
custom_barplot1(plotData, xVar = "valence_f", yVar = "dilation_n", subVar = "subject_f", yLim = yLimDil)

## 2D bar plot: Dilation per required action per valence:
custom_barplot2(plotData, xVar = "reqAction_f", yVar = "dilation_n", zVar = "valence_f", subVar = "subject_f", yLim = yLimDil)
custom_barplot2(plotData, xVar = "valence_f", yVar = "dilation_n", zVar = "reqAction_f", subVar = "subject_f", yLim = yLimDil)

## 2D bar plot: Dilation per response per valence:
custom_barplot2(plotData, xVar = "response_f", yVar = "dilation_n", zVar = "valence_f", subVar = "subject_f", yLim = yLimDil)
custom_barplot2(plotData, xVar = "valence_f", yVar = "dilation_n", zVar = "response_f", subVar = "subject_f", yLim = yLimDil)

## 1D bar plot: Dilation per condition:
custom_barplot1(plotData, xVar = "condition_f", yVar = "dilation_n", subVar = "subject_f", yLim = yLimDil)
custom_barplot1(plotData, xVar = "respCond_f", yVar = "dilation_n", subVar = "subject_f", yLim = yLimDil)

## 1D bar plot: Dilation per arousal:
custom_barplot1(plotData, xVar = "arousal_f", yVar = "dilation_n", subVar = "subject_f", yLim = yLimDil)

## 2D bar plot: Dilation per arousal per required action/ valence/ congruency/ condition:
custom_barplot2(plotData, xVar = "arousal_f", yVar = "dilation_n", zVar = "reqAction_f", subVar = "subject_f", yLim = yLimDil)
custom_barplot2(plotData, xVar = "arousal_f", yVar = "dilation_n", zVar = "reqAction_f", subVar = "subject_f", yLim = yLimDil, addLegend = T)
custom_barplot2(plotData, xVar = "reqAction_f", yVar = "dilation_n", zVar = "arousal_f", subVar = "subject_f", yLim = yLimDil, addLegend = T)
custom_barplot2(plotData, xVar = "arousal_f", yVar = "dilation_n", zVar = "valence_f", subVar = "subject_f", yLim = yLimDil)
custom_barplot2(plotData, xVar = "valence_f", yVar = "dilation_n", zVar = "arousal_f", subVar = "subject_f", yLim = yLimDil)

custom_barplot2(plotData, xVar = "arousal_f", yVar = "dilation_n", zVar = "reqCongruency_short2_f", subVar = "subject_f", yLim = yLimDil)
custom_barplot2(plotData, xVar = "reqCongruency_short1_f", yVar = "dilation_n", zVar = "arousal_f", subVar = "subject_f", yLim = yLimDil)
custom_barplot2(plotData, xVar = "condition_short1_f", yVar = "dilation_n", zVar = "arousal_f", subVar = "subject_f", yLim = yLimDil) # , FTS = 28)
custom_barplot2(plotData, xVar = "respCond_f", yVar = "dilation_n", zVar = "arousal_f", subVar = "subject_f", yLim = yLimDil) # , FTS = 28)

## 1D bar plot: Dilation per congruency (required action):
custom_barplot1(plotData, xVar = "reqCongruency_f", yVar = "dilation_n", subVar = "subject_f", yLim = yLimDil)
custom_barplot1(plotData, xVar = "actCongruency_f", yVar = "dilation_n", subVar = "subject_f", yLim = yLimDil)

## 1D bar plot: Dilation per accuracy/ repeat/ RT fast/ slow:
custom_barplot1(plotData, xVar = "ACC_f", yVar = "dilation_n", subVar = "subject_f", yLim = yLimDil)
custom_barplot1(plotData, xVar = "repeat_f", yVar = "dilation_n", subVar = "subject_f", yLim = yLimDil)
custom_barplot1(plotData, xVar = "RTcleaned_fast_f", yVar = "dilation_n", subVar = "subject_f", yLim = yLimDil)

yLimDil <- c(5, 28)
custom_barplot2(plotData, yVar = "dilation_n", xVar = "ACC_f", zVar = "valence_f", subVar = "subject_f", yLim = yLimDil)
custom_barplot2(plotData, yVar = "dilation_n", xVar = "valence_f", zVar = "ACC_short_f", subVar = "subject_f", yLim = yLimDil)
custom_barplot2(plotData, yVar = "dilation_n", xVar = "RTcleaned_fast_f", zVar = "valence_f", subVar = "subject_f", yLim = yLimDil)
custom_barplot2(plotData, yVar = "dilation_n", xVar = "valence_f", zVar = "RTcleaned_fast_f", subVar = "subject_f", yLim = yLimDil)

custom_barplot2(plotData, yVar = "dilation_n", xVar = "ACC_f", zVar = "RTcleaned_fast_f", subVar = "subject_f", yLim = yLimDil)
custom_barplot2(plotData, yVar = "dilation_n", xVar = "RTcleaned_fast_f", zVar = "ACC_short_f", subVar = "subject_f", yLim = yLimDil)
custom_barplot2(plotData, yVar = "dilation_n", xVar = "ACC_f", zVar = "ACC_short_f", subVar = "subject_f", yLim = yLimDil)
custom_barplot2(plotData, yVar = "dilation_n", xVar = "valence_f", zVar = "valence_f", subVar = "subject_f", yLim = yLimDil)

## 1D bar plot: Dilation per outcome last trial:
custom_barplot1(plotData, yVar = "dilation_n", xVar = "outcome_last_f",subVar = "subject_f", yLim = yLimDil)
custom_barplot1(plotData, yVar = "dilation_n", xVar = "outcome_last_rel_f", subVar = "subject_f", yLim = yLimDil)
# custom_barplot1(plotData, xVar = "outcome_last_all_f", yVar = "dilation_n", subVar = "subject_f", yLim = yLimDil)

## 3-way interaction: Dilation per required action x valence x arousal:
yLimDil <- c(0, 26)
plotData <- subset(modData, !is.na(dilation_n) & !is.na(reqAction_f) & !is.na(valence_f) & !is.na(arousal_f))
ggplot(plotData, aes(y = dilation_n, x = reqAction_f, fill = valence_f)) + 
  stat_summary(fun = mean, geom = "bar", aes(fill = valence_f), position = "dodge") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", position = position_dodge(width=0.9), width = .1) +
  facet_wrap(vars(arousal_f))
custom_barplot3(data = plotData, yVar = "dilation_n", xVar = "reqAction_f", zVar = "valence_f", splitVar = "arousal_f", subVar = "subject_f", yLim = yLimDil)

## 3-way interaction: Dilation per response x valence x arousal:
yLimDil <- c(0, 26)
plotData <- subset(modData, !is.na(dilation_n) & !is.na(response_f) & !is.na(valence_f) & !is.na(arousal_f))
ggplot(plotData, aes(y = dilation_n, x = response_f, fill = valence_f)) + 
  stat_summary(fun = mean, geom = "bar", aes(fill = valence_f), position = "dodge") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", position = position_dodge(width=0.9), width = .1) +
  facet_wrap(vars(arousal_f))
custom_barplot3(data = plotData, yVar = "dilation_n", xVar = "response_f", zVar = "valence_f", splitVar = "arousal_f", subVar = "subject_f", yLim = yLimDil)

yLimDil <- c(0, 26)
custom_barplot3(data = plotData, yVar = "dilation_n", xVar = "response_f", zVar = "valence_f", splitVar = "firstHalfTask_f", subVar = "subject_f", yLim = yLimDil)
custom_barplot3(data = plotData, yVar = "dilation_n", xVar = "response_f", zVar = "valence_f", splitVar = "firstHalfBlock_f", subVar = "subject_f", yLim = yLimDil)
custom_barplot3(data = plotData, yVar = "dilation_n", xVar = "response_f", zVar = "valence_f", splitVar = "firstHalfCueRep_f", subVar = "subject_f", yLim = yLimDil)

yLimDil <- c(0, 28)
custom_barplot3(data = plotData, yVar = "dilation_n", xVar = "ACC_short_f", zVar = "valence_f", splitVar = "RTcleaned_fast_f", subVar = "subject_f", yLim = yLimDil)
custom_barplot3(data = plotData, yVar = "dilation_n", xVar = "valence_f", zVar = "ACC_short_f", splitVar = "RTcleaned_fast_f", subVar = "subject_f", yLim = yLimDil)
custom_barplot3(data = plotData, yVar = "dilation_n", xVar = "valence_f", zVar = "RTcleaned_fast_f", splitVar = "ACC_short_f", subVar = "subject_f", yLim = yLimDil)
custom_barplot3(data = plotData, yVar = "dilation_n", xVar = "RTcleaned_fast_f", zVar = "valence_f", splitVar = "ACC_short_f", subVar = "subject_f", yLim = yLimDil)
custom_barplot3(data = plotData, yVar = "dilation_n", xVar = "ACC_short_f", zVar = "RTcleaned_fast_f", splitVar = "valence_f", subVar = "subject_f", yLim = yLimDil)
custom_barplot3(data = plotData, yVar = "dilation_n", xVar = "RTcleaned_fast_f", zVar = "ACC_short_f", splitVar = "valence_f", subVar = "subject_f", yLim = yLimDil)

# ============================================================================ #
#### Bar plots DILATION at OUTCOME: ####

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
segMessage <- "StartOutcome"
beforeSamples <- 1000
afterSamples <- 2000
data <- add_pupil(data, segMessage = segMessage, beforeSamples = beforeSamples, afterSamples = afterSamples)

modData <- select_standardize(data)
plotData <- modData

## 1D plots:
yLimDil <- c(0, 20)
custom_barplot1(plotData, yVar = "dilation_n", xVar = "outcome_rel_f", subVar = "subject_f", yLim = yLimDil)
custom_barplot1(plotData, yVar = "dilation_n", xVar = "outcome_f", subVar = "subject_f", yLim = yLimDil)
custom_barplot1(plotData, yVar = "dilation_n", xVar = "outcome_all_f", subVar = "subject_f", yLim = yLimDil)
custom_barplot1(plotData, yVar = "dilation_n", xVar = "outcome_all_short_f", subVar = "subject_f", yLim = yLimDil)

## 2D plots:
yLimDil <- c(0, 20)
custom_barplot2(plotData, yVar = "dilation_n", xVar = "outcome_rel_f", zVar = "response_f", subVar = "subject_f", yLim = yLimDil)
custom_barplot2(plotData, yVar = "dilation_n", xVar = "outcome_f", zVar = "response_f", subVar = "subject_f", yLim = yLimDil)
custom_barplot2(plotData, yVar = "dilation_n", xVar = "outcome_short_f", zVar = "response_f", subVar = "subject_f", yLim = yLimDil)
custom_barplot2(plotData, yVar = "dilation_n", xVar = "outcome_short_f", zVar = "response_f", subVar = "subject_f", yLim = yLimDil)
custom_barplot2(plotData, yVar = "dilation_n", xVar = "outcome_all_short_f", zVar = "response_f", subVar = "subject_f", yLim = yLimDil)

custom_barplot2(plotData, yVar = "dilation_n", xVar = "outcome_all_short_f", zVar = "ACC_f", subVar = "subject_f", yLim = yLimDil)
custom_barplot2(plotData, yVar = "dilation_n", xVar = "outcome_all_short_f", zVar = "RTcleaned_fast_f", subVar = "subject_f", yLim = yLimDil)
custom_barplot2(plotData, yVar = "dilation_n", xVar = "outcome_all_short_f", zVar = "repeat_f", subVar = "subject_f", yLim = yLimDil)

cat("Finished! :-)\n")

# END OF FILE.
