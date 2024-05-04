#!/usr/bin/env Rscript
# ============================================================================ #
## 08_mgngunc_plots_final.R
## MGNGUncertainty plots of behavior, gaze, and pupil data for main manuscript and supplementary material.
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
source(paste0(dirs$codeDir, "functions/00_mgngunc_functions_timecourse.R")) # Load functions

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
segMessage <- "StartMask"
beforeSamples <- 1000
afterSamples <- 2966
data <- add_pupil(data, segMessage = segMessage, beforeSamples = beforeSamples, afterSamples = afterSamples)

# ----------------------------------- #
## Select data, standardize variables:

modData <- select_standardize(data)
plotData <- modData

# ============================================================================ #
# ============================================================================ #
# ============================================================================ #
# ============================================================================ #
# ============================================================================ #
# ============================================================================ #
#### SELECTION OF PLOTS FOR MAIN MANUSCRIPT: ####

## Global plotting settings:
plotHeight <- 480
plotWidth <- 480
plotWidthTime <- 600
plotWidthGAMM <- 600
plotWidthCoef <- 720 # 900

# ============================================================================ #
#### Figure01: grand mean and hypotheses: ####

# -------------------------------------------------------------- #
### Figure01B: Time course: dilation_n ~ all_f:

## Settings for loading time courses:
segMessage <- "StartMask"
beforeSamples <- 1000
afterSamples <- 2966
zVar <- "all_f"
p <- plot_timecourse(data = data, zVar = zVar, beforeSamples = beforeSamples, afterSamples = afterSamples, 
                     isBaselineCorGroup = F, addLegend = F)
png(paste0(dirs$plot, "final/Figure01B.png"), width = 600, height = plotHeight) # plotWidth
print(p)
dev.off()

# -------------------------------------------------------------- #
### Figure01E & F: Hypotheses Dilation ~ reqAction x valence:

## Define data:
hypData <- data.frame(condition_n <- 1:4)
hypData$reqAction_f <- factor(c("Go", "Go", "NoGo", "NoGo"))
hypData$valence_f <- factor(c("1Win", "2Avoid", "1Win", "2Avoid"))
hypData$congruency_f <- factor(c("cong", "incong", "incong", "cong"))
hypData$Cog <- c(0.30, 0.70, 0.70, 0.30)
hypData$Phys <- c(0.70, 0.80, 0.30, 0.30)

## Plotting settings:
LWD <- retrieve_plot_defaults("LWD") # 1.3
FTS <- 32 # retrieve_plot_defaults("FTS")
dodgeVal <- 0.6
colAlpha <- 1
yLim <- c(0, 1)
selCol <- retrieve_colour("condition_f")

## Select variable:
hypData$y <- hypData$Cog; yName <- "01E"
hypData$y <- hypData$Phys; yName <- "01F"

## Plot:
p <- ggplot(hypData, aes(x = reqAction_f, fill = valence_f, y = y, pattern = congruency_f)) + 
  stat_summary(fun = mean, geom = "bar", position = "dodge", width = dodgeVal,
               lwd = LWD, color = "black") + 
  scale_fill_manual(values = selCol, limits = levels(hypData$valence_f)) + 
  # scale_pattern_manual(values = c(incong = "stripe", cong = "none")) +
  labs(x = "Response", y = "Pupil dilation (a.u.)", fill = "Valence") +
  coord_cartesian(ylim = yLim) + 
  theme_classic() + 
  theme(axis.text = element_text(size = FTS),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_text(size = FTS), 
        plot.title = element_text(size = FTS, hjust = 0.5), 
        legend.title = element_blank(), legend.position = "none",
        axis.line = element_line(colour = 'black')) # , linewidth = LWD)) # fixed font sizes
print(p)  

## Save:
plotName <- paste0("Figure", yName)
png(paste0(dirs$plot, "final/", plotName, ".png"), width = 240, height = 480)
print(p)
dev.off()

# ============================================================================ #
#### Figure02: responses/ RTs ~ reqAction x valence ####

# -------------------------------------------------------------- #
### Figure02A: Line plot: response_n ~ reqAction_f * valence_f:
p <- custom_lineplot_gg(plotData, xVar = "cueRep_n", yVar = "response_n", zVar = "condition_f", subVar = "subject_f", 
                   selLineType = c(1, 1, 2, 2), breakVec = c(1, 5, 10, 15, 20), addLegend = F, savePNG = F)
png(paste0(dirs$plot, "final/Figure02A.png"), width = plotWidth, height = plotHeight)
print(p)
dev.off()

# -------------------------------------------------------------- #
### Figure 02B: Bar plot: response_n ~ reqAction_f * valence_f:
p <- custom_barplot2(plotData, xVar = "reqAction_f", yVar = "response_n", zVar = "valence_f", subVar = "subject_f")
png(paste0(dirs$plot, "final/Figure02B.png"), width = plotWidth, height = plotHeight)
print(p)
dev.off()

# -------------------------------------------------------------- #
## Figure 02C: Coefplot: response_n ~ reqAction_f * valence_f:
formula <- "response_n ~ reqAction_f * valence_f + (reqAction_f * valence_f|subject_f)"
mod <- fit_lmem(formula)
nCoef <- length(fixef(mod)) - 1; nCoef
selCol <- met.brewer("Demuth", n = nCoef, type = "continuous"); colName <- "Demuth"
p <- custom_coefplot(mod, plotSub = T, plotText = T, dropIntercept = T, revOrder = T, selCol = selCol)
png(paste0(dirs$plot, "final/Figure02C.png"), width = plotWidthCoef, height = plotHeight)
print(p)
dev.off()

# -------------------------------------------------------------- #
## Figure 02D: Density plot: RTcleaned_n ~ valence_f:
p <- customplot_density2(plotData, xVar = "RT_n", zVar = "valence_f", addLegend = F) 
png(paste0(dirs$plot, "final/Figure02D.png"), width = plotWidth, height = plotHeight)
print(p)
dev.off()

# -------------------------------------------------------------- #
## Figure 02E: Bar plot: RTcleaned_n ~ reqAction_f * valence_f:
yLimRT <- c(0.4, 1.0)
p <- custom_barplot2(plotData, xVar = "reqAction_f", yVar = "RTcleaned_n", zVar = "valence_f", subVar = "subject_f", yLim = yLimRT)
png(paste0(dirs$plot, "final/Figure02E.png"), width = plotWidth, height = plotHeight)
print(p)
dev.off()

# -------------------------------------------------------------- #
## Figure 02F: Coef plot: RTcleaned_z ~ reqAction_f * valence_f:
formula <- "RTcleaned_z ~ reqAction_f * valence_f + (reqAction_f * valence_f|subject_f)"
mod <- fit_lmem(formula)
nCoef <- length(fixef(mod)) - 1; nCoef
selCol <- met.brewer("Demuth", n = nCoef, type = "continuous"); colName <- "Demuth"
p <- custom_coefplot(mod, plotSub = T, plotText = T, dropIntercept = T, revOrder = T, selCol = selCol, xLim = c(-0.5, 0.2))
png(paste0(dirs$plot, "final/Figure02F.png"), width = plotWidthCoef, height = plotHeight)
print(p)
dev.off()

# ============================================================================ #
#### Figure03: Freezing of gaze: ####

# -------------------------------------------------------------- #
## Settings:
segMessage <- "StartCue"
beforeSamples <- 866
afterSamples <- 1300
aggrMethod <- "distanceBaselineTrial"

# -------------------------------------------------------------- #
### Figure 03A: Gaze distance ~ valence, zoomed out:

zVar <- "valence_f"
isBaselineCorGroup <- F
p <- plot_timecourse(data = data, zVar = zVar, 
                     eyeMeasure = "gaze", aggrMetric = aggrMethod,
                     segMessage = segMessage, beforeSamples = beforeSamples, afterSamples = afterSamples,
                     isBaselineCorGroup = isBaselineCorGroup, addLegend = F)
png(paste0(dirs$plot, "final/Figure03A.png"), width = plotWidthTime, height = plotHeight)
print(p)
dev.off()

# -------------------------------------------------------------- #
### Figure 03B: Gaze distance ~ valence, zoomed in:

xLim <- c(-100, 400)
zVar <- "valence_f"
isBaselineCorGroup <- F
p <- plot_timecourse(data = data, zVar = zVar, 
                     eyeMeasure = "gaze", aggrMetric = aggrMethod,
                     segMessage = segMessage, beforeSamples = beforeSamples, afterSamples = afterSamples,
                     xLim = xLim, suffix = paste0("_xLim", xLim[1], "-", xLim[2]),
                     yLim = c(9.5, 17.5),
                     isBaselineCorGroup = isBaselineCorGroup, addLegend = F)
png(paste0(dirs$plot, "final/Figure03B.png"), width = plotWidthTime, height = plotHeight)
print(p)
dev.off()

# -------------------------------------------------------------- #
### Figure 03C: Max distance from center ~ valence:

## Load data:
segMessage <- "StartCue"
beforeSamples <- 866
afterSamples <- 1300
gazeMetric <- "distanceBaselineTrial" 
aggrMethod <- "max"
timeWindow <- c(202, 278)
newVarLong <- paste0(gazeMetric, "_", aggrMethod, "_", str_pad(timeWindow[1], width = 4, side = "left", pad = "0"), "-",
                     str_pad(timeWindow[2], width = 4, side = "left", pad = "0"), "ms");newVarLong
data <- add_gaze_summary(data, gazeMetric = gazeMetric, aggrMethod = aggrMethod, timeWindow = timeWindow,
                         segMessage = segMessage, beforeSamples = beforeSamples, afterSamples = afterSamples)
modData <- select_standardize(data)
plotData <- modData

## Plot:
yLimGaze <- c(10, 32); yName <- "Max" # Max 202-278
p <- custom_barplot1(plotData, yVar = newVarLong, xVar = "valence_f", subVar = "subject_f", yLim = yLimGaze,
                     yLab = paste0(yName, " dist. center (", timeWindow[1], "-", timeWindow[2], " ms)"),
                     suffix = "pretty")
png(paste0(dirs$plot, "final/Figure03C.png"), width = plotWidth, height = plotHeight)
print(p)
dev.off()

# -------------------------------------------------------------- #
### Figure 03D: Gaze distance ~ valence, trials 1-5:

segMessage <- "StartCue"
beforeSamples <- 866
afterSamples <- 1300
aggrMethod <- "distanceBaselineTrial"

xLim <- c(-100, 400)
zVar <- "valence_trial01_05_f"
isBaselineCorGroup <- F
p <- plot_timecourse(data = data, zVar = zVar, 
                     eyeMeasure = "gaze", aggrMetric = aggrMethod,
                     segMessage = segMessage, beforeSamples = beforeSamples, afterSamples = afterSamples,
                     xLim = xLim, suffix = paste0("_xLim", xLim[1], "-", xLim[2]),
                     yLim = c(9.5, 17.5),
                     isBaselineCorGroup = isBaselineCorGroup, addLegend = F)
png(paste0(dirs$plot, "final/Figure03D.png"), width = plotWidthTime, height = plotHeight)
print(p)
dev.off()

# -------------------------------------------------------------- #
### Figure 03E: Gaze distance ~ valence, trials 1-5:

segMessage <- "StartCue"
beforeSamples <- 866
afterSamples <- 1300
aggrMethod <- "distanceBaselineTrial"

xLim <- c(-100, 400)
zVar <- "valence_trial06_10_f"
table(data$valence_trial06_10_f, data$cueRep_n)
isBaselineCorGroup <- F
p <- plot_timecourse(data = data, zVar = zVar, 
                     eyeMeasure = "gaze", aggrMetric = aggrMethod,
                     segMessage = segMessage, beforeSamples = beforeSamples, afterSamples = afterSamples,
                     xLim = xLim, suffix = paste0("_xLim", xLim[1], "-", xLim[2]),
                     yLim = c(9.5, 17.5),
                     isBaselineCorGroup = isBaselineCorGroup, addLegend = F)
png(paste0(dirs$plot, "final/Figure03E.png"), width = plotWidthTime, height = plotHeight)
print(p)
dev.off()

# -------------------------------------------------------------- #
### Figure 03F: Gaze distance ~ valence, trials 1-5:

segMessage <- "StartCue"
beforeSamples <- 866
afterSamples <- 1300
aggrMethod <- "distanceBaselineTrial"

xLim <- c(-100, 400)
zVar <- "valence_trial11_15_f"
isBaselineCorGroup <- F
p <- plot_timecourse(data = data, zVar = zVar, 
                     eyeMeasure = "gaze", aggrMetric = aggrMethod,
                     segMessage = segMessage, beforeSamples = beforeSamples, afterSamples = afterSamples,
                     xLim = xLim, suffix = paste0("_xLim", xLim[1], "-", xLim[2]),
                     yLim = c(9.5, 17.5),
                     isBaselineCorGroup = isBaselineCorGroup, addLegend = F)
png(paste0(dirs$plot, "final/Figure03F.png"), width = plotWidthTime, height = plotHeight)
print(p)
dev.off()

# ============================================================================ #
#### Figure04: dilations ~ reqAction x valence: ####

# -------------------------------------------------------------- #
### Figure04A: Bar plot: dilation ~ reqAction_f * valence_f:

yLimDil <- c(6, 26)
p <- custom_barplot2(plotData, xVar = "reqAction_f", yVar = "dilation_n", zVar = "valence_f", subVar = "subject_f", yLim = yLimDil)
png(paste0(dirs$plot, "final/Figure04A.png"), width = plotWidth, height = plotHeight)
print(p)
dev.off()

# -------------------------------------------------------------- #
### Figure04B: Coef plot: dilation_z ~ reqAction_f * valence_f:

formula <- "dilation_z ~ response_f * valence_f + (response_f * valence_f|subject_f)"
mod <- fit_lmem(formula)
nCoef <- length(fixef(mod)) - 1; nCoef
selCol <- met.brewer("Demuth", n = nCoef, type = "continuous"); colName <- "Demuth"
p <- custom_coefplot(mod, plotSub = T, plotText = T, dropIntercept = T, revOrder = T, selCol = selCol)
png(paste0(dirs$plot, "final/Figure04B.png"), width = plotWidthCoef, height = plotHeight)
print(p)
dev.off()

# -------------------------------------------------------------- #
### Figure04C: Time course: dilation_z ~ reqAction_f * valence_f, mask-locked:

## Settings for loading time courses:
segMessage <- "StartMask"
beforeSamples <- 1000
afterSamples <- 2966
zVar <- "respCond_f"
p <- plot_timecourse(data = data, zVar = zVar, segMessage = segMessage, beforeSamples = beforeSamples, afterSamples = afterSamples,
                     selLineType = c(1, 1, 2, 2), xLim = c(-500, 2966), yLim = c(-25, 100),
                     isBaselineCorGroup  = T, baselineTime = 0,
                     addLegend = F)
png(paste0(dirs$plot, "final/Figure04C.png"), width = 720, height = plotHeight) # plotWidthTime
print(p)
dev.off()

# -------------------------------------------------------------- #
### Figure04D: GAMM plot: dilation_z ~ reqAction_f * valence_f:

## Read in behavioral data:
data <- read_behavior()
## Preprocessing: 
data <- wrapper_preprocessing(data)
## Add pupil data:
segMessage <- "StartMask"
beforeSamples <- 1000
afterSamples <- 2966
data <- add_pupil(data, segMessage = segMessage, beforeSamples = beforeSamples, afterSamples = afterSamples)
modData <- select_standardize(data)
modData$baseline12_n <- modData$baseline_n + 12
plotData <- modData

## Start plot:
yVar <- "dilation_n"; timeVar <- "cueRep_n"; groupVar <- "subject_block_f"
splitVar <- "respCond_f";
mod <- fit_gamm(data = modData, yVar = yVar, timeVar = timeVar, splitVar = splitVar, groupVar = groupVar,
                bsType = "fs", addARIMA = F, useScat = T)
png(paste0(dirs$plot, "final/Figure04D.png"), width = plotWidthGAMM, height = plotHeight)
custom_plot_smooth(mod, ltyVec = c(1, 1, 2, 2), isPNG = F)
dev.off()

# -------------------------------------------------------------- #
### Figure04E: GAMM diff plot: dilation_z ~ valence_f for Go responses:

yVar <- "dilation_n"; timeVar <- "cueRep_n"; groupVar <- "subject_block_f"
splitVar <- "respCond_f"; splitLevels <- c("G2A", "G2W")
mod <- fit_gamm_diff(data = modData, yVar = yVar, timeVar = timeVar, splitVar = splitVar, groupVar = groupVar, splitLevels = splitLevels,
                     bsType = "fs", addARIMA = F, useScat = T)
png(paste0(dirs$plot, "final/Figure04E.png"), width = plotWidthGAMM, height = plotHeight)
custom_plot_diff(mod, isPNG = F)
dev.off()

# -------------------------------------------------------------- #
### Figure04F: Time course: dilation_z ~ reqAction_f * valence_f, response-locked:

## Settings for loading time courses:
segMessage <- "StartResponse"
beforeSamples <- 2000
afterSamples <- 3500
zVar <- "respCond_f"
p <- plot_timecourse(data = data, zVar = zVar, segMessage = segMessage, beforeSamples = beforeSamples, afterSamples = afterSamples,
                     selLineType = c(1, 1, 2, 2), 
                     xLim = c(-1000, 2500), yLim = c(-15, 100),
                     isBaselineCorGroup = T, isBaselineCorSub = F, baselineTime = -1433, # baseline-correct per group
                     addLegend = F)
png(paste0(dirs$plot, "final/Figure04F.png"), width = 720, height = plotHeight) # plotWidthTime
print(p)
dev.off()


# ============================================================================ #
# ============================================================================ #
# ============================================================================ #
# ============================================================================ #
# ============================================================================ #
#### FigureS01: Effect of required action/ valence on RTs over time: ####

source(paste0(dirs$codeDir, "functions/00_mgngunc_functions_regression.R")) # Load functions

# -------------------------------------------------------------- #
### FigureS01A: GAMM plot: RT_n ~ reqAction_f * valence_f:

## Read in behavioral data:
data <- read_behavior()
## Preprocessing: 
data <- wrapper_preprocessing(data)

## Start plot:
yVar <- "RTcleaned_n"; timeVar <- "cueRep_n"; groupVar <- "subject_block_f"
splitVar <- "condition_f";
mod <- fit_gamm(data = modData, yVar = yVar, timeVar = timeVar, splitVar = splitVar, groupVar = groupVar,
                bsType = "fs", addARIMA = F, useScat = T)
png(paste0(dirs$plot, "final/FigureS01A.png"), width = plotWidthGAMM, height = plotHeight)
custom_plot_smooth(mod, ltyVec = c(1, 1, 2, 2), isPNG = F)
dev.off()

# -------------------------------------------------------------- #
### Figure01B: GAMM diff plot: RT_n ~ reqAction_f (i.e. required action = accuracy):

yVar <- "RTcleaned_n"; timeVar <- "cueRep_n"; groupVar <- "subject_block_f"
splitVar <- "reqAction_f"; splitLevels <- c("NoGo", "Go")
mod <- fit_gamm_diff(data = modData, yVar = yVar, timeVar = timeVar, splitVar = splitVar, groupVar = groupVar, splitLevels = splitLevels,
                     bsType = "fs", addARIMA = F, useScat = T)
png(paste0(dirs$plot, "final/FigureS01B.png"), width = plotWidthGAMM, height = plotHeight)
custom_plot_diff(mod, isPNG = F)
dev.off()

# -------------------------------------------------------------- #
### Figure01C: GAMM diff plot: RT_n ~ valence_f:

yVar <- "RTcleaned_n"; timeVar <- "cueRep_n"; groupVar <- "subject_block_f"
splitVar <- "valence_f"; splitLevels <- c("Avoid", "Win")
mod <- fit_gamm_diff(data = modData, yVar = yVar, timeVar = timeVar, splitVar = splitVar, groupVar = groupVar, splitLevels = splitLevels,
                     bsType = "fs", addARIMA = F, useScat = T)
png(paste0(dirs$plot, "final/FigureS01C.png"), width = plotWidthGAMM, height = plotHeight)
custom_plot_diff(mod, isPNG = F)
dev.off()

# ============================================================================ #
#### FigureS02: Correlations Pav bias choices/ RTs ~ questionnaires: ####

# ----------------------------------- #
## Read in behavioral data:
data <- read_behavior()

# ----------------------------------- #
## Pre-processing: 
data <- wrapper_preprocessing(data)
names(data)
table(data$subject_n)

# ----------------------------------- #
## Select data, standardize variables:
modData <- select_standardize(data)

# ----------------------------------- #
## Read in questionnaires:
questData <- load_preprocess_questionnaires()

# -------------------------------------------------------------- #
### FigureS02A: Correlation: response ~ valence with STAI:
xVar <- "stai_mean_n"; xLab <- "Trait anxiety (STAI)"
formula <- "response_n ~ reqAction_f * valence_f + (reqAction_f * valence_f|subject_f)"; yLab <- "Pavlovian bias (choices)"; coefIdx <- 3
mod <- fit_lmem(formula)
questData$regCoef_n <- coef(mod)$subject_f[, coefIdx]
p <- plot_correlation(data = questData, xVar = xVar, yVar = "regCoef_n", xLab = xLab, yLab = yLab)
png(paste0(dirs$plot, "final/FigureS02A.png"), width = 480, height = 480)
print(p)
dev.off()

# -------------------------------------------------------------- #
### FigureS02B: Correlation: response ~ valence with UPPS-P negative urgency:
xVar <- "UPPSP_negative_urgency_n"; xLab <- "Negative urgency (UPPS-P)"
formula <- "response_n ~ reqAction_f * valence_f + (reqAction_f * valence_f|subject_f)"; yLab <- "Pavlovian bias (choices)"; coefIdx <- 3
mod <- fit_lmem(formula)
questData$regCoef_n <- coef(mod)$subject_f[, coefIdx]
p <- plot_correlation(data = questData, xVar = xVar, yVar = "regCoef_n", xLab = xLab, yLab = yLab)
png(paste0(dirs$plot, "final/FigureS02B.png"), width = 480, height = 480)
print(p)
dev.off()

# -------------------------------------------------------------- #
### FigureS02C: Correlation: response ~ valence with UPPS-P lack of perseverance:
xVar <- "UPPSP_lack_perseverance_n"; xLab <- "Lack of perseverance (UPPS-P)"
formula <- "response_n ~ reqAction_f * valence_f + (reqAction_f * valence_f|subject_f)"; yLab <- "Pavlovian bias (choices)"; coefIdx <- 3
mod <- fit_lmem(formula)
questData$regCoef_n <- coef(mod)$subject_f[, coefIdx]
p <- plot_correlation(data = questData, xVar = xVar, yVar = "regCoef_n", xLab = xLab, yLab = yLab)
png(paste0(dirs$plot, "final/FigureS02C.png"), width = 480, height = 480)
print(p)
dev.off()

# -------------------------------------------------------------- #
### FigureS02D: Correlation: response ~ valence with UPPS-P lack of premeditation:
xVar <- "UPPSP_lack_premeditation_n"; xLab <- "Lack of premeditation (UPPS-P)"
formula <- "response_n ~ reqAction_f * valence_f + (reqAction_f * valence_f|subject_f)"; yLab <- "Pavlovian bias (choices)"; coefIdx <- 3
mod <- fit_lmem(formula)
questData$regCoef_n <- coef(mod)$subject_f[, coefIdx]
p <- plot_correlation(data = questData, xVar = xVar, yVar = "regCoef_n", xLab = xLab, yLab = yLab)
png(paste0(dirs$plot, "final/FigureS02D.png"), width = 480, height = 480)
print(p)
dev.off()

# -------------------------------------------------------------- #
### FigureS02E: Correlation: response ~ valence with UPPS-P sensation seeking:
xVar <- "UPPSP_sensation_seeking_n"; xLab <- "Sensation seeking (UPPS-P)"
formula <- "response_n ~ reqAction_f * valence_f + (reqAction_f * valence_f|subject_f)"; yLab <- "Pavlovian bias (choices)"; coefIdx <- 3
mod <- fit_lmem(formula)
questData$regCoef_n <- coef(mod)$subject_f[, coefIdx]
p <- plot_correlation(data = questData, xVar = xVar, yVar = "regCoef_n", xLab = xLab, yLab = yLab)
png(paste0(dirs$plot, "final/FigureS02E.png"), width = 480, height = 480)
print(p)
dev.off()

# -------------------------------------------------------------- #
### FigureS02F: Correlation: response ~ valence with UPPS-P positive urgency:
xVar <- "UPPSP_positive_urgency_n"; xLab <- "Positive urgency (UPPS-P)"
formula <- "response_n ~ reqAction_f * valence_f + (reqAction_f * valence_f|subject_f)"; yLab <- "Pavlovian bias (choices)"; coefIdx <- 3
mod <- fit_lmem(formula)
questData$regCoef_n <- coef(mod)$subject_f[, coefIdx]
p <- plot_correlation(data = questData, xVar = xVar, yVar = "regCoef_n", xLab = xLab, yLab = yLab)
png(paste0(dirs$plot, "final/FigureS02F.png"), width = 480, height = 480)
print(p)
dev.off()

# -------------------------------------------------------------- #
### FigureS02G: Correlation: RTs ~ valence with STAI:
xVar <- "stai_mean_n"; xLab <- "Trait anxiety (STAI)"
formula <- "RTcleaned_z ~ reqAction_f * valence_f + (reqAction_f * valence_f|subject_f)"; yLab <- "Pavlovian bias (RTs)"; coefIdx <- 3
mod <- fit_lmem(formula)
questData$regCoef_n <- coef(mod)$subject_f[, coefIdx]
p <- plot_correlation(data = questData, xVar = xVar, yVar = "regCoef_n", xLab = xLab, yLab = yLab)
png(paste0(dirs$plot, "final/FigureS02G.png"), width = 480, height = 480)
print(p)
dev.off()

# -------------------------------------------------------------- #
### FigureS02H: Correlation: RTs ~ valence with UPPS-P negative urgency:
xVar <- "UPPSP_negative_urgency_n"; xLab <- "Negative urgency (UPPS-P)"
formula <- "RTcleaned_z ~ reqAction_f * valence_f + (reqAction_f * valence_f|subject_f)"; yLab <- "Pavlovian bias (RTs)"; coefIdx <- 3
mod <- fit_lmem(formula)
questData$regCoef_n <- coef(mod)$subject_f[, coefIdx]
p <- plot_correlation(data = questData, xVar = xVar, yVar = "regCoef_n", xLab = xLab, yLab = yLab)
png(paste0(dirs$plot, "final/FigureS02H.png"), width = 480, height = 480)
print(p)
dev.off()

# -------------------------------------------------------------- #
### FigureS02I: Correlation: RTs ~ valence with UPPS-P lack of perseverance:
xVar <- "UPPSP_lack_perseverance_n"; xLab <- "Lack of perseverance (UPPS-P)"
formula <- "RTcleaned_z ~ reqAction_f * valence_f + (reqAction_f * valence_f|subject_f)"; yLab <- "Pavlovian bias (RTs)"; coefIdx <- 3
mod <- fit_lmem(formula)
questData$regCoef_n <- coef(mod)$subject_f[, coefIdx]
p <- plot_correlation(data = questData, xVar = xVar, yVar = "regCoef_n", xLab = xLab, yLab = yLab)
png(paste0(dirs$plot, "final/FigureS02I.png"), width = 480, height = 480)
print(p)
dev.off()

# -------------------------------------------------------------- #
### FigureS02J: Correlation: RTs ~ valence with UPPS-P lack of premeditation:
xVar <- "UPPSP_lack_premeditation_n"; xLab <- "Lack of premeditation (UPPS-P)"
formula <- "RTcleaned_z ~ reqAction_f * valence_f + (reqAction_f * valence_f|subject_f)"; yLab <- "Pavlovian bias (RTs)"; coefIdx <- 3
mod <- fit_lmem(formula)
questData$regCoef_n <- coef(mod)$subject_f[, coefIdx]
p <- plot_correlation(data = questData, xVar = xVar, yVar = "regCoef_n", xLab = xLab, yLab = yLab)
png(paste0(dirs$plot, "final/FigureS02J.png"), width = 480, height = 480)
print(p)
dev.off()

# -------------------------------------------------------------- #
### FigureS02K: Correlation: RTs ~ valence with UPPS-P sensation seeking:
xVar <- "UPPSP_sensation_seeking_n"; xLab <- "Sensation seeking (UPPS-P)"
formula <- "RTcleaned_z ~ reqAction_f * valence_f + (reqAction_f * valence_f|subject_f)"; yLab <- "Pavlovian bias (RTs)"; coefIdx <- 3
mod <- fit_lmem(formula)
questData$regCoef_n <- coef(mod)$subject_f[, coefIdx]
p <- plot_correlation(data = questData, xVar = xVar, yVar = "regCoef_n", xLab = xLab, yLab = yLab)
png(paste0(dirs$plot, "final/FigureS02K.png"), width = 480, height = 480)
print(p)
dev.off()

# -------------------------------------------------------------- #
### FigureS02L: Correlation: RTs ~ valence with UPPS-P positive urgency:
xVar <- "UPPSP_positive_urgency_n"; xLab <- "Positive urgency (UPPS-P)"
formula <- "RTcleaned_z ~ reqAction_f * valence_f + (reqAction_f * valence_f|subject_f)"; yLab <- "Pavlovian bias (RTs)"; coefIdx <- 3
mod <- fit_lmem(formula)
questData$regCoef_n <- coef(mod)$subject_f[, coefIdx]
p <- plot_correlation(data = questData, xVar = xVar, yVar = "regCoef_n", xLab = xLab, yLab = yLab)
png(paste0(dirs$plot, "final/FigureS02L.png"), width = 480, height = 480)
print(p)
dev.off()

# ============================================================================ #
#### FigureS03: Dilations ~ prime manipulation: ####

# -------------------------------------------------------------- #
### FigureS03A: Bar plot: dilation ~ manipulation_f:

yLimDil <- c(6, 26)
p <- custom_barplot1(plotData, xVar = "arousal_f", yVar = "dilation_n", subVar = "subject_f", yLim = yLimDil)
png(paste0(dirs$plot, "final/FigureS03A.png"), width = plotWidth, height = plotHeight)
print(p)
dev.off()

# -------------------------------------------------------------- #
### FigureS03B: Time course: dilation ~ arousal:

## Settings for loading time courses:
segMessage <- "StartMask"
beforeSamples <- 1000
afterSamples <- 2966
zVar <- "arousal_f"
p <- plot_timecourse(data = data, zVar = zVar, beforeSamples = beforeSamples, afterSamples = afterSamples,
                     selLineType = c(1, 1, 2, 2), yLim = c(-15, 75),
                     isBaselineCorGroup = T, addLegend = F)
png(paste0(dirs$plot, "final/FigureS03B.png"), width = plotWidthTime, height = plotHeight)
print(p)
dev.off()

# -------------------------------------------------------------- #
### FigureS03C: GAMM plot: dilation ~ arousal:

## Read in behavioral data:
data <- read_behavior()
## Preprocessing: 
data <- wrapper_preprocessing(data)
## Add pupil data:
segMessage <- "StartMask"
beforeSamples <- 1000
afterSamples <- 2966
data <- add_pupil(data, segMessage = segMessage, beforeSamples = beforeSamples, afterSamples = afterSamples)
modData <- select_standardize(data)
modData$baseline12_n <- modData$baseline_n + 12
plotData <- modData

## Start plot:
yVar <- "dilation_n"; timeVar <- "cueRep_n"; groupVar <- "subject_block_f"
splitVar <- "arousal_f";
mod <- fit_gamm(data = modData, yVar = yVar, timeVar = timeVar, splitVar = splitVar, groupVar = groupVar,
                bsType = "fs", addARIMA = F, useScat = T)
png(paste0(dirs$plot, "final/FigureS03C.png"), width = plotWidthGAMM, height = plotHeight)
custom_plot_smooth(mod, ltyVec = c(1, 1, 2, 2), isPNG = F)
dev.off()

# ============================================================================ #
#### FigureS04: responses/ RTs ~ prime: ####

# -------------------------------------------------------------- #
### Figure04A: Bar plot: responses ~ arousal:

p <- custom_barplot1(plotData, xVar = "arousal_f", yVar = "response_n", subVar = "subject_f")
png(paste0(dirs$plot, "final/Figure04A.png"), width = plotWidth, height = plotHeight)
print(p)
dev.off()

# -------------------------------------------------------------- #
### Figure04B: Bar plot: responses ~ arousal x condition:

p <- custom_barplot2(plotData, xVar = "condition_short1_f", yVar = "response_n", zVar = "arousal_f", subVar = "subject_f")
png(paste0(dirs$plot, "final/Figure04B.png"), width = plotWidth, height = plotHeight)
print(p)
dev.off()

# -------------------------------------------------------------- #
### Figure04C: Coef plot: responses ~ reqAction x valence x arousal:

formula <- "response_n ~ reqAction_f * valence_f * arousal_f + (reqAction_f * valence_f * arousal_f|subject_f)"
mod <- fit_lmem(formula)
nCoef <- length(fixef(mod)) - 1; nCoef
selCol <- met.brewer("Manet", n = nCoef, type = "continuous"); colName <- "Manet" # brown, brown, green
p <- custom_coefplot(mod, plotSub = T, plotText = T, dropIntercept = T, revOrder = T, selCol = selCol)
png(paste0(dirs$plot, "final/Figure04C.png"), width = plotWidthCoef, height = plotHeight)
print(p)
dev.off()

# -------------------------------------------------------------- #
### Figure04D: Density plot: RT ~ arousal:

p <- customplot_density2(plotData, xVar = "RT_n", zVar = "arousal_f", addLegend = F) 
png(paste0(dirs$plot, "final/Figure04D.png"), width = plotWidth, height = plotHeight)
print(p)
dev.off()

# -------------------------------------------------------------- #
### Figure04E: Bar plot: RT ~ arousal x condition:

yLimRT <- c(0.4, 1.0)
p <- custom_barplot2(plotData, xVar = "condition_short1_f", yVar = "RTcleaned_n", zVar = "arousal_f", subVar = "subject_f", yLim = yLimRT)
png(paste0(dirs$plot, "final/Figure04E.png"), width = plotWidth, height = plotHeight)
print(p)
dev.off()

# -------------------------------------------------------------- #
### Figure04F: Coef plot: RT ~ reqAction x valence x arousal:

formula <- "RTcleaned_z ~ reqAction_f * valence_f * arousal_f + (reqAction_f * valence_f * arousal_f|subject_f)"
mod <- fit_lmem(formula)
nCoef <- length(fixef(mod)) - 1; nCoef
selCol <- met.brewer("Manet", n = nCoef, type = "continuous"); colName <- "Manet" # brown, brown, green
p <- custom_coefplot(mod, plotSub = T, plotText = T, dropIntercept = T, revOrder = T, selCol = selCol, xLim = c(-0.5, 0.2))
png(paste0(dirs$plot, "final/Figure04F.png"), width = plotWidthCoef, height = plotHeight)
print(p)
dev.off()

# ============================================================================ #
#### FigureS05: Correlation prime effect on choices/ RTs with prime effect on dilations: ####

## Fit models:
formula <- "response_n ~ arousal_f + (arousal_f|subject_f)"
mod <- fit_lmem(formula)
respVec <- coef(mod)[[1]][, 2]
formula <- "RTcleaned_z ~ arousal_f + (arousal_f|subject_f)"
mod <- fit_lmem(formula)
rtVec <- coef(mod)[[1]][, 2]
formula <- "dilation_z ~ arousal_f + (arousal_f|subject_f)"
mod <- fit_lmem(formula)
dilVec <- coef(mod)[[1]][, 2]

## Data frame:
corData <- data.frame(cbind(respVec, rtVec, dilVec))
names(corData) <- c("response_n", "RTcleaned_z", "dilation_z")

# -------------------------------------------------------------- #
### FigureS05A: Correlation response ~ arousal with dilation ~ arousal:

## Plot:
xVar <- "response_n"; yVar <- "dilation_z"; xLab <- "responses"; yLab <- "dilations"
p <- plot_correlation(data = corData, xVar = xVar, yVar = yVar, 
                      xLab = paste0("Effect of arousal manipulation on ", xLab), 
                      yLab = paste0("Effect of arousal manipulation on ", yLab))
png(paste0(dirs$plot, "final/FigureS05A.png"), width = plotWidth, height = plotHeight)
print(p)
dev.off()

# -------------------------------------------------------------- #
### Figure S05B: Correlation RT ~ arousal with dilation ~ arousal:

## Plot:
xVar <- "RTcleaned_z"; yVar <- "dilation_z"; xLab <- "RTs"; yLab <- "dilations"
p <- plot_correlation(data = corData, xVar = xVar, yVar = yVar, 
                      xLab = paste0("Effect of arousal manipulation on ", xLab), 
                      yLab = paste0("Effect of arousal manipulation on ", yLab))
png(paste0(dirs$plot, "final/FigureS05B.png"), width = plotWidth, height = plotHeight)
print(p)
dev.off()

# ============================================================================ #
#### FigureS06: Correlations prime effect on choices/ RTs ~ questionnaires: ####

# ----------------------------------- #
## Read in behavioral data:
data <- read_behavior()

# ----------------------------------- #
## Pre-processing: 
data <- wrapper_preprocessing(data)
names(data)
table(data$subject_n)

# ----------------------------------- #
## Select data, standardize variables:
modData <- select_standardize(data)

# ----------------------------------- #
## Read in questionnaires:
questData <- load_preprocess_questionnaires()

# -------------------------------------------------------------- #
### FigureS06A: Correlation: response ~ prime with STAI:
xVar <- "stai_mean_n"; xLab <- "Trait anxiety (STAI)"
formula <- "response_n ~ reqAction_f * valence_f * arousal_f + (reqAction_f * valence_f * arousal_f|subject_f)"; yLab <- "Effect of prime on choice"; coefIdx <- 4
mod <- fit_lmem(formula)
questData$regCoef_n <- coef(mod)$subject_f[, coefIdx]
p <- plot_correlation(data = questData, xVar = xVar, yVar = "regCoef_n", xLab = xLab, yLab = yLab)
png(paste0(dirs$plot, "final/FigureS06A.png"), width = 480, height = 480)
print(p)
dev.off()

# -------------------------------------------------------------- #
### FigureS06B: Correlation: response ~ prime with UPPS-P negative urgency:
xVar <- "UPPSP_negative_urgency_n"; xLab <- "Negative urgency (UPPS-P)"
formula <- "response_n ~ reqAction_f * valence_f * arousal_f + (reqAction_f * valence_f * arousal_f|subject_f)"; yLab <- "Effect of prime on choice"; coefIdx <- 4
mod <- fit_lmem(formula)
questData$regCoef_n <- coef(mod)$subject_f[, coefIdx]
p <- plot_correlation(data = questData, xVar = xVar, yVar = "regCoef_n", xLab = xLab, yLab = yLab)
png(paste0(dirs$plot, "final/FigureS06B.png"), width = 480, height = 480)
print(p)
dev.off()

# -------------------------------------------------------------- #
### FigureS06C: Correlation: response ~ prime with UPPS-P lack of perseverance:
xVar <- "UPPSP_lack_perseverance_n"; xLab <- "Lack of perseverance (UPPS-P)"
formula <- "response_n ~ reqAction_f * valence_f * arousal_f + (reqAction_f * valence_f * arousal_f|subject_f)"; yLab <- "Effect of prime on choice"; coefIdx <- 4
mod <- fit_lmem(formula)
questData$regCoef_n <- coef(mod)$subject_f[, coefIdx]
p <- plot_correlation(data = questData, xVar = xVar, yVar = "regCoef_n", xLab = xLab, yLab = yLab)
png(paste0(dirs$plot, "final/FigureS06C.png"), width = 480, height = 480)
print(p)
dev.off()

# -------------------------------------------------------------- #
### FigureS06D: Correlation: response ~ prime with UPPS-P lack of premeditation:
xVar <- "UPPSP_lack_premeditation_n"; xLab <- "Lack of premeditation (UPPS-P)"
formula <- "response_n ~ reqAction_f * valence_f * arousal_f + (reqAction_f * valence_f * arousal_f|subject_f)"; yLab <- "Effect of prime on choice"; coefIdx <- 4
mod <- fit_lmem(formula)
questData$regCoef_n <- coef(mod)$subject_f[, coefIdx]
p <- plot_correlation(data = questData, xVar = xVar, yVar = "regCoef_n", xLab = xLab, yLab = yLab)
png(paste0(dirs$plot, "final/FigureS06D.png"), width = 480, height = 480)
print(p)
dev.off()

# -------------------------------------------------------------- #
### FigureS06E: Correlation: response ~ prime with UPPS-P sensation seeking:
xVar <- "UPPSP_sensation_seeking_n"; xLab <- "Sensation seeking (UPPS-P)"
formula <- "response_n ~ reqAction_f * valence_f * arousal_f + (reqAction_f * valence_f * arousal_f|subject_f)"; yLab <- "Effect of prime on choice"; coefIdx <- 4
mod <- fit_lmem(formula)
questData$regCoef_n <- coef(mod)$subject_f[, coefIdx]
p <- plot_correlation(data = questData, xVar = xVar, yVar = "regCoef_n", xLab = xLab, yLab = yLab)
png(paste0(dirs$plot, "final/FigureS06E.png"), width = 480, height = 480)
print(p)
dev.off()

# -------------------------------------------------------------- #
### FigureS06F: Correlation: response ~ prime with UPPS-P positive urgency:
xVar <- "UPPSP_positive_urgency_n"; xLab <- "Positive urgency (UPPS-P)"
formula <- "response_n ~ reqAction_f * valence_f * arousal_f + (reqAction_f * valence_f * arousal_f|subject_f)"; yLab <- "Effect of prime on choice"; coefIdx <- 4
mod <- fit_lmem(formula)
questData$regCoef_n <- coef(mod)$subject_f[, coefIdx]
p <- plot_correlation(data = questData, xVar = xVar, yVar = "regCoef_n", xLab = xLab, yLab = yLab)
png(paste0(dirs$plot, "final/FigureS06F.png"), width = 480, height = 480)
print(p)
dev.off()

# -------------------------------------------------------------- #
### FigureS06G: Correlation: RTs ~ prime with STAI:
xVar <- "stai_mean_n"; xLab <- "Trait anxiety (STAI)"
formula <- "RTcleaned_z ~ reqAction_f * valence_f * arousal_f + (reqAction_f * valence_f * arousal_f|subject_f)"; yLab <- "Effect of prime on RTs"; coefIdx <- 4
mod <- fit_lmem(formula)
questData$regCoef_n <- coef(mod)$subject_f[, coefIdx]
p <- plot_correlation(data = questData, xVar = xVar, yVar = "regCoef_n", xLab = xLab, yLab = yLab)
png(paste0(dirs$plot, "final/FigureS06G.png"), width = 480, height = 480)
print(p)
dev.off()

# -------------------------------------------------------------- #
### FigureS06H: Correlation: RTs ~ prime with UPPS-P negative urgency:
xVar <- "UPPSP_negative_urgency_n"; xLab <- "Negative urgency (UPPS-P)"
formula <- "RTcleaned_z ~ reqAction_f * valence_f * arousal_f + (reqAction_f * valence_f * arousal_f|subject_f)"; yLab <- "Effect of prime on RTs"; coefIdx <- 4
mod <- fit_lmem(formula)
questData$regCoef_n <- coef(mod)$subject_f[, coefIdx]
p <- plot_correlation(data = questData, xVar = xVar, yVar = "regCoef_n", xLab = xLab, yLab = yLab)
png(paste0(dirs$plot, "final/FigureS06H.png"), width = 480, height = 480)
print(p)
dev.off()

# -------------------------------------------------------------- #
### FigureS06I: Correlation: RTs ~ prime with UPPS-P lack of perseverance:
xVar <- "UPPSP_lack_perseverance_n"; xLab <- "Lack of perseverance (UPPS-P)"
formula <- "RTcleaned_z ~ reqAction_f * valence_f * arousal_f + (reqAction_f * valence_f * arousal_f|subject_f)"; yLab <- "Effect of prime on RTs"; coefIdx <- 4
mod <- fit_lmem(formula)
questData$regCoef_n <- coef(mod)$subject_f[, coefIdx]
p <- plot_correlation(data = questData, xVar = xVar, yVar = "regCoef_n", xLab = xLab, yLab = yLab)
png(paste0(dirs$plot, "final/FigureS06I.png"), width = 480, height = 480)
print(p)
dev.off()

# -------------------------------------------------------------- #
### FigureS06J: Correlation: RTs ~ prime with UPPS-P lack of premeditation:
xVar <- "UPPSP_lack_premeditation_n"; xLab <- "Lack of premeditation (UPPS-P)"
formula <- "RTcleaned_z ~ reqAction_f * valence_f * arousal_f + (reqAction_f * valence_f * arousal_f|subject_f)"; yLab <- "Effect of prime on RTs"; coefIdx <- 4
mod <- fit_lmem(formula)
questData$regCoef_n <- coef(mod)$subject_f[, coefIdx]
p <- plot_correlation(data = questData, xVar = xVar, yVar = "regCoef_n", xLab = xLab, yLab = yLab)
png(paste0(dirs$plot, "final/FigureS06J.png"), width = 480, height = 480)
print(p)
dev.off()

# -------------------------------------------------------------- #
### FigureS06K: Correlation: RTs ~ prime with UPPS-P sensation seeking:
xVar <- "UPPSP_sensation_seeking_n"; xLab <- "Sensation seeking (UPPS-P)"
formula <- "RTcleaned_z ~ reqAction_f * valence_f * arousal_f + (reqAction_f * valence_f * arousal_f|subject_f)"; yLab <- "Effect of prime on RTs"; coefIdx <- 4
mod <- fit_lmem(formula)
questData$regCoef_n <- coef(mod)$subject_f[, coefIdx]
p <- plot_correlation(data = questData, xVar = xVar, yVar = "regCoef_n", xLab = xLab, yLab = yLab)
png(paste0(dirs$plot, "final/FigureS06K.png"), width = 480, height = 480)
print(p)
dev.off()

# -------------------------------------------------------------- #
### FigureS06L: Correlation: RTs ~ prime with UPPS-P positive urgency:
xVar <- "UPPSP_positive_urgency_n"; xLab <- "Positive urgency (UPPS-P)"
formula <- "RTcleaned_z ~ reqAction_f * valence_f * arousal_f + (reqAction_f * valence_f * arousal_f|subject_f)"; yLab <- "Effect of prime on RTs"; coefIdx <- 4
mod <- fit_lmem(formula)
questData$regCoef_n <- coef(mod)$subject_f[, coefIdx]
p <- plot_correlation(data = questData, xVar = xVar, yVar = "regCoef_n", xLab = xLab, yLab = yLab)
png(paste0(dirs$plot, "final/FigureS06L.png"), width = 480, height = 480)
print(p)
dev.off()

# ============================================================================ #
#### FigureS07: responses/ RTs ~ dilations: ####

# -------------------------------------------------------------- #
### FigureS07A: Regression line plot: responses ~ dilation 

formula <- "response_n ~ dilation_z + (dilation_z|subject_f)"
mod <- fit_lmem(formula)
selEff <- "dilation_z"
yLim <- c(0, 1) # responses
p <- custom_regressionline1(mod, selEff = selEff, yLim = yLim, useEffect = F)
png(paste0(dirs$plot, "final/FigureS07A.png"), width = plotWidth, height = plotHeight)
print(p)
dev.off()

# -------------------------------------------------------------- #
### FigureS07B: Regression line plot: responses ~ dilation x required action

formula <- "response_n ~ reqAction_f * dilation_z  + (reqAction_f * dilation_z|subject_f)"
mod <- fit_lmem(formula)
xVar <- "dilation_z"; zVar <- "reqAction_f"; 
yLim <- c(0, 1) # responses
p <- custom_regressionline2(mod = mod, xVar = xVar, zVar = zVar, yLim = yLim)
png(paste0(dirs$plot, "final/FigureS07B.png"), width = plotWidth, height = plotHeight)
print(p)
dev.off()

# -------------------------------------------------------------- #
### FigureS07C: Coef plot: responses ~ reqAction x valence x dilation 

formula <- "response_n ~ reqAction_f * valence_f * dilation_z + (reqAction_f * valence_f * dilation_z|subject_f)"
mod <- fit_lmem(formula)
nCoef <- length(fixef(mod)) - 1; nCoef
selCol <- met.brewer("Monet", n = nCoef, type = "continuous"); colName <- "Manet" # brown, brown, green
p <- custom_coefplot(mod, plotSub = T, plotText = T, dropIntercept = T, revOrder = T, selCol = selCol)
png(paste0(dirs$plot, "final/FigureS07C.png"), width = plotWidthCoef, height = plotHeight)
print(p)
dev.off()

# -------------------------------------------------------------- #
### FigureS07D: Regression line plot: RT ~ dilation x valence

formula <- "RTcleaned_z ~ valence_f * dilation_z  + (valence_f * dilation_z|subject_f)"
mod <- fit_lmem(formula)
xVar <- "dilation_z"; zVar <- "valence_f"
yLim <- c(-1, 2) # RT_z
p <- custom_regressionline2(mod = mod, xVar = xVar, zVar = zVar, yLim = yLim)
png(paste0(dirs$plot, "final/FigureS07D.png"), width = plotWidth, height = plotHeight)
print(p)
dev.off()

# -------------------------------------------------------------- #
### FigureS07E: Regression line plot: RT ~ dilation x required action:

formula <- "RTcleaned_z ~ reqAction_f * dilation_z  + (reqAction_f * dilation_z|subject_f)"
mod <- fit_lmem(formula)
xVar <- "dilation_z"; zVar <- "reqAction_f"
yLim <- c(-1, 2) # RT_z
p <- custom_regressionline2(mod = mod, xVar = xVar, zVar = zVar, yLim = yLim)
png(paste0(dirs$plot, "final/FigureS07E.png"), width = plotWidth, height = plotHeight)
print(p)
dev.off()

# -------------------------------------------------------------- #
### FigureS07F: Coef plot: RT ~ reqAction x valence x dilation:

formula <- "RTcleaned_z ~ reqAction_f * valence_f * dilation_z + (reqAction_f * valence_f * dilation_z|subject_f)"
mod <- fit_lmem(formula)
nCoef <- length(fixef(mod)) - 1; nCoef
selCol <- met.brewer("Monet", n = nCoef, type = "continuous"); colName <- "Manet" # brown, brown, green
p <- custom_coefplot(mod, plotSub = T, plotText = T, dropIntercept = T, revOrder = T, selCol = selCol, xLim = c(-0.5, 0.2))
png(paste0(dirs$plot, "final/FigureS07F.png"), width = plotWidthCoef, height = plotHeight)
print(p)
dev.off()

# ============================================================================ #
#### FigureS08: Correlations dilation effect on responses/ RTs ~ questionnaires: ####

# ----------------------------------- #
## Read in behavioral data:
data <- read_behavior()

# ----------------------------------- #
## Pre-processing: 
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

# ----------------------------------- #
## Read in questionnaires:
questData <- load_preprocess_questionnaires()

# -------------------------------------------------------------- #
### FigureS08A: Correlation: response ~ dilation with STAI:
xVar <- "stai_mean_n"; xLab <- "Trait anxiety (STAI)"
formula <- "response_n ~ reqAction_f * valence_f * dilation_z + (reqAction_f * valence_f * dilation_z|subject_f)"; yLab <- "Association dilations and Go"; coefIdx <- 4
mod <- fit_lmem(formula)
questData$regCoef_n <- coef(mod)$subject_f[, coefIdx]
p <- plot_correlation(data = questData, xVar = xVar, yVar = "regCoef_n", xLab = xLab, yLab = yLab)
png(paste0(dirs$plot, "final/FigureS08A.png"), width = 480, height = 480)
print(p)
dev.off()

# -------------------------------------------------------------- #
### FigureS08B: Correlation: response ~ dilation with UPPS-P negative urgency:
xVar <- "UPPSP_negative_urgency_n"; xLab <- "Negative urgency (UPPS-P)"
formula <- "response_n ~ reqAction_f * valence_f * dilation_z + (reqAction_f * valence_f * dilation_z|subject_f)"; yLab <- "Association dilations and Go"; coefIdx <- 4
mod <- fit_lmem(formula)
questData$regCoef_n <- coef(mod)$subject_f[, coefIdx]
p <- plot_correlation(data = questData, xVar = xVar, yVar = "regCoef_n", xLab = xLab, yLab = yLab)
png(paste0(dirs$plot, "final/FigureS08B.png"), width = 480, height = 480)
print(p)
dev.off()

# -------------------------------------------------------------- #
### FigureS08C: Correlation: response ~ dilation with UPPS-P lack of perseverance:
xVar <- "UPPSP_lack_perseverance_n"; xLab <- "Lack of perseverance (UPPS-P)"
formula <- "response_n ~ reqAction_f * valence_f * dilation_z + (reqAction_f * valence_f * dilation_z|subject_f)"; yLab <- "Association dilations and Go"; coefIdx <- 4
mod <- fit_lmem(formula)
questData$regCoef_n <- coef(mod)$subject_f[, coefIdx]
p <- plot_correlation(data = questData, xVar = xVar, yVar = "regCoef_n", xLab = xLab, yLab = yLab)
png(paste0(dirs$plot, "final/FigureS08C.png"), width = 480, height = 480)
print(p)
dev.off()

# -------------------------------------------------------------- #
### FigureS08D: Correlation: response ~ dilation with UPPS-P lack of premeditation:
xVar <- "UPPSP_lack_premeditation_n"; xLab <- "Lack of premeditation (UPPS-P)"
formula <- "response_n ~ reqAction_f * valence_f * dilation_z + (reqAction_f * valence_f * dilation_z|subject_f)"; yLab <- "Association dilations and Go"; coefIdx <- 4
mod <- fit_lmem(formula)
questData$regCoef_n <- coef(mod)$subject_f[, coefIdx]
p <- plot_correlation(data = questData, xVar = xVar, yVar = "regCoef_n", xLab = xLab, yLab = yLab)
png(paste0(dirs$plot, "final/FigureS08D.png"), width = 480, height = 480)
print(p)
dev.off()

# -------------------------------------------------------------- #
### FigureS08E: Correlation: response ~ dilation with UPPS-P sensation seeking:
xVar <- "UPPSP_sensation_seeking_n"; xLab <- "Sensation seeking (UPPS-P)"
formula <- "response_n ~ reqAction_f * valence_f * dilation_z + (reqAction_f * valence_f * dilation_z|subject_f)"; yLab <- "Association dilations and Go"; coefIdx <- 4
mod <- fit_lmem(formula)
questData$regCoef_n <- coef(mod)$subject_f[, coefIdx]
p <- plot_correlation(data = questData, xVar = xVar, yVar = "regCoef_n", xLab = xLab, yLab = yLab)
png(paste0(dirs$plot, "final/FigureS08E.png"), width = 480, height = 480)
print(p)
dev.off()

# -------------------------------------------------------------- #
### FigureS08F: Correlation: response ~ dilation with UPPS-P positive urgency:
xVar <- "UPPSP_positive_urgency_n"; xLab <- "Positive urgency (UPPS-P)"
formula <- "response_n ~ reqAction_f * valence_f * dilation_z + (reqAction_f * valence_f * dilation_z|subject_f)"; yLab <- "Association dilations and Go"; coefIdx <- 4
mod <- fit_lmem(formula)
questData$regCoef_n <- coef(mod)$subject_f[, coefIdx]
p <- plot_correlation(data = questData, xVar = xVar, yVar = "regCoef_n", xLab = xLab, yLab = yLab)
png(paste0(dirs$plot, "final/FigureS08F.png"), width = 480, height = 480)
print(p)
dev.off()

# -------------------------------------------------------------- #
### FigureS08G: Correlation: RTs ~ dilation with STAI:
xVar <- "stai_mean_n"; xLab <- "Trait anxiety (STAI)"
formula <- "RTcleaned_z ~ reqAction_f * valence_f * dilation_z + (reqAction_f * valence_f * dilation_z|subject_f)"; yLab <- "Association dilations and RTs"; coefIdx <- 4
mod <- fit_lmem(formula)
questData$regCoef_n <- coef(mod)$subject_f[, coefIdx]
p <- plot_correlation(data = questData, xVar = xVar, yVar = "regCoef_n", xLab = xLab, yLab = yLab)
png(paste0(dirs$plot, "final/FigureS08G.png"), width = 480, height = 480)
print(p)
dev.off()

# -------------------------------------------------------------- #
### FigureS08H: Correlation: RTs ~ dilation with UPPS-P negative urgency:
xVar <- "UPPSP_negative_urgency_n"; xLab <- "Negative urgency (UPPS-P)"
formula <- "RTcleaned_z ~ reqAction_f * valence_f * dilation_z + (reqAction_f * valence_f * dilation_z|subject_f)"; yLab <- "Association dilations and RTs"; coefIdx <- 4
mod <- fit_lmem(formula)
questData$regCoef_n <- coef(mod)$subject_f[, coefIdx]
p <- plot_correlation(data = questData, xVar = xVar, yVar = "regCoef_n", xLab = xLab, yLab = yLab)
png(paste0(dirs$plot, "final/FigureS08H.png"), width = 480, height = 480)
print(p)
dev.off()

# -------------------------------------------------------------- #
### FigureS08I: Correlation: RTs ~ dilation with UPPS-P lack of perseverance:
xVar <- "UPPSP_lack_perseverance_n"; xLab <- "Lack of perseverance (UPPS-P)"
formula <- "RTcleaned_z ~ reqAction_f * valence_f * dilation_z + (reqAction_f * valence_f * dilation_z|subject_f)"; yLab <- "Association dilations and RTs"; coefIdx <- 4
mod <- fit_lmem(formula)
questData$regCoef_n <- coef(mod)$subject_f[, coefIdx]
p <- plot_correlation(data = questData, xVar = xVar, yVar = "regCoef_n", xLab = xLab, yLab = yLab)
png(paste0(dirs$plot, "final/FigureS08I.png"), width = 480, height = 480)
print(p)
dev.off()

# -------------------------------------------------------------- #
### FigureS08J: Correlation: RTs ~ dilation with UPPS-P lack of premeditation:
xVar <- "UPPSP_lack_premeditation_n"; xLab <- "Lack of premeditation (UPPS-P)"
formula <- "RTcleaned_z ~ reqAction_f * valence_f * dilation_z + (reqAction_f * valence_f * dilation_z|subject_f)"; yLab <- "Association dilations and RTs"; coefIdx <- 4
mod <- fit_lmem(formula)
questData$regCoef_n <- coef(mod)$subject_f[, coefIdx]
p <- plot_correlation(data = questData, xVar = xVar, yVar = "regCoef_n", xLab = xLab, yLab = yLab)
png(paste0(dirs$plot, "final/FigureS08J.png"), width = 480, height = 480)
print(p)
dev.off()

# -------------------------------------------------------------- #
### FigureS08K: Correlation: RTs ~ dilation with UPPS-P sensation seeking:
xVar <- "UPPSP_sensation_seeking_n"; xLab <- "Sensation seeking (UPPS-P)"
formula <- "RTcleaned_z ~ reqAction_f * valence_f * dilation_z + (reqAction_f * valence_f * dilation_z|subject_f)"; yLab <- "Association dilations and RTs"; coefIdx <- 4
mod <- fit_lmem(formula)
questData$regCoef_n <- coef(mod)$subject_f[, coefIdx]
p <- plot_correlation(data = questData, xVar = xVar, yVar = "regCoef_n", xLab = xLab, yLab = yLab)
png(paste0(dirs$plot, "final/FigureS08K.png"), width = 480, height = 480)
print(p)
dev.off()

# -------------------------------------------------------------- #
### FigureS08L: Correlation: RTs ~ dilation with UPPS-P positive urgency:
xVar <- "UPPSP_positive_urgency_n"; xLab <- "Positive urgency (UPPS-P)"
formula <- "RTcleaned_z ~ reqAction_f * valence_f * dilation_z + (reqAction_f * valence_f * dilation_z|subject_f)"; yLab <- "Association dilations and RTs"; coefIdx <- 4
mod <- fit_lmem(formula)
questData$regCoef_n <- coef(mod)$subject_f[, coefIdx]
p <- plot_correlation(data = questData, xVar = xVar, yVar = "regCoef_n", xLab = xLab, yLab = yLab)
png(paste0(dirs$plot, "final/FigureS08L.png"), width = 480, height = 480)
print(p)
dev.off()

# ============================================================================ #
#### FigureS09: Dilations ~ ACC/ RT/ repeat x response: ####

## Read in behavioral data:
data <- read_behavior()
## Preprocessing: 
data <- wrapper_preprocessing(data)
## Add pupil data:
segMessage <- "StartMask"
beforeSamples <- 1000
afterSamples <- 2966
data <- add_pupil(data, segMessage = segMessage, beforeSamples = beforeSamples, afterSamples = afterSamples)
modData <- select_standardize(data)
modData$baseline12_n <- modData$baseline_n + 12
plotData <- modData

# -------------------------------------------------------------- #
### FigureS09A: bar plot: dilation ~ response x ACC:

yLimDil <- c(6, 28)
p <- custom_barplot2(plotData, xVar = "response_f", yVar = "dilation_n", zVar = "ACC_f", subVar = "subject_f", yLim = yLimDil)
png(paste0(dirs$plot, "final/FigureS09A.png"), width = plotWidth, height = plotHeight)
print(p)
dev.off()

# -------------------------------------------------------------- #
### FigureS09B: bar plot: dilation ~ RT

yLimDil <- c(6, 28)
p <- custom_barplot1(plotData, xVar = "RTcleaned_fast_f", yVar = "dilation_n", subVar = "subject_f", yLim = yLimDil)
png(paste0(dirs$plot, "final/FigureS09B.png"), width = plotWidth, height = plotHeight)
print(p)
dev.off()

# -------------------------------------------------------------- #
### FigureS09C: bar plot: dilation ~ response x repeat:

yLimDil <- c(6, 28)
p <- custom_barplot2(plotData, xVar = "response_f", yVar = "dilation_n", zVar = "repeat_f", subVar = "subject_f", yLim = yLimDil)
png(paste0(dirs$plot, "final/FigureS09C.png"), width = plotWidth, height = plotHeight)
print(p)
dev.off()

# -------------------------------------------------------------- #
### FigureS09D: GAMM plot: dilation ~ response x ACC:

yVar <- "dilation_n"; timeVar <- "cueRep_n"; groupVar <- "subject_block_f"
splitVar <- "ACC_response_f";
mod <- fit_gamm(data = modData, yVar = yVar, timeVar = timeVar, splitVar = splitVar, groupVar = groupVar,
                bsType = "fs", addARIMA = F, useScat = T)
png(paste0(dirs$plot, "final/FigureS09D.png"), width = plotWidthGAMM, height = plotHeight)
custom_plot_smooth(mod, ltyVec = c(1, 2, 1, 2), isPNG = F)
dev.off()

# -------------------------------------------------------------- #
### FigureS09E: GAMM plot: dilation ~ response x RT:

yVar <- "dilation_n"; timeVar <- "cueRep_n"; groupVar <- "subject_block_f"
splitVar <- "RTcleaned_fast_f";
mod <- fit_gamm(data = modData, yVar = yVar, timeVar = timeVar, splitVar = splitVar, groupVar = groupVar,
                bsType = "fs", addARIMA = F, useScat = T)
png(paste0(dirs$plot, "final/FigureS09E.png"), width = plotWidthGAMM, height = plotHeight)
custom_plot_smooth(mod, isPNG = F)
dev.off()

# -------------------------------------------------------------- #
### FigureS09F: GAMM plot: dilation ~ response x repeat:

yVar <- "dilation_n"; timeVar <- "cueRep_n"; groupVar <- "subject_block_f"
splitVar <- "repeat_response_f";
mod <- fit_gamm(data = modData, yVar = yVar, timeVar = timeVar, splitVar = splitVar, groupVar = groupVar,
                bsType = "fs", addARIMA = F, useScat = T)
png(paste0(dirs$plot, "final/FigureS09F.png"), width = plotWidthGAMM, height = plotHeight)
custom_plot_smooth(mod, ltyVec = c(1, 2, 1, 2), isPNG = F)
dev.off()

# ============================================================================ #
#### FigureS10: Dilations ~ ACC x RT ####

# -------------------------------------------------------------- #
### FigureS10A: bar plot: dilation ~ response x repeat:

yLimDil <- c(6, 28)
p <- custom_barplot2(plotData, xVar = "RTcleaned_fast_f", yVar = "dilation_n", zVar = "ACC_f", subVar = "subject_f", yLim = yLimDil)
png(paste0(dirs$plot, "final/FigureS10A.png"), width = plotWidth, height = plotHeight)
print(p)
dev.off()

# -------------------------------------------------------------- #
### FigureS10B: GAMM plot: dilation ~ ACC x RT:

yVar <- "dilation_n"; timeVar <- "cueRep_n"; groupVar <- "subject_block_f"
splitVar <- "RTcleaned_ACC_f";
mod <- fit_gamm(data = modData, yVar = yVar, timeVar = timeVar, splitVar = splitVar, groupVar = groupVar,
                bsType = "fs", addARIMA = F, useScat = T)
png(paste0(dirs$plot, "final/FigureS10B.png"), width = plotWidthGAMM, height = plotHeight)
custom_plot_smooth(mod, ltyVec = c(1, 1, 2, 2), isPNG = F)
dev.off()

# ============================================================================ #
#### FigureS11: Dilations ~ ACC/ RT/ repeat x valence ####

yLimDil <- c(4, 36)

# -------------------------------------------------------------- #
### FigureS11A: bar plot: dilation ~ ACC x valence:

p <- custom_barplot2(subset(plotData, response_f == "Go"), xVar = "ACC_f", yVar = "dilation_n", zVar = "valence_f", subVar = "subject_f", yLim = yLimDil,
                     suffix = paste0("_respGo_yLim", yLimDil[1], "-", yLimDil[2]))
png(paste0(dirs$plot, "final/FigureS11A.png"), width = plotWidth, height = plotHeight)
print(p)
dev.off()

# -------------------------------------------------------------- #
### FigureS11B: bar plot: dilation ~ RT x valence:

p <- custom_barplot2(subset(plotData, response_f == "Go"), xVar = "RTcleaned_fast_f", yVar = "dilation_n", zVar = "valence_f", subVar = "subject_f", yLim = yLimDil,
                     suffix = paste0("_respGo_yLim", yLimDil[1], "-", yLimDil[2]))
png(paste0(dirs$plot, "final/FigureS11B.png"), width = plotWidth, height = plotHeight)
print(p)
dev.off()

# -------------------------------------------------------------- #
### FigureS11C: bar plot: dilation ~ repeat x valence:

p <- custom_barplot2(subset(plotData, response_f == "Go"), xVar = "repeat_f", yVar = "dilation_n", zVar = "valence_f", subVar = "subject_f", yLim = yLimDil,
                     suffix = paste0("_respGo_yLim", yLimDil[1], "-", yLimDil[2]))
png(paste0(dirs$plot, "final/FigureS11C.png"), width = plotWidth, height = plotHeight)
print(p)
dev.off()

# -------------------------------------------------------------- #
### FigureS11D: GAMM plot: dilation ~ ACC x valence:

yVar <- "dilation_n"; timeVar <- "cueRep_n"; groupVar <- "subject_block_f"
splitVar <- "ACC_valence_Go_f";
mod <- fit_gamm(data = modData, yVar = yVar, timeVar = timeVar, splitVar = splitVar, groupVar = groupVar,
                bsType = "fs", addARIMA = F, useScat = T)
png(paste0(dirs$plot, "final/FigureS11D.png"), width = plotWidthGAMM, height = plotHeight)
custom_plot_smooth(mod, ltyVec = c(1, 1, 2, 2), isPNG = F)
dev.off()

# -------------------------------------------------------------- #
### FigureS11E: GAMM plot: dilation ~ RT x valence:

yVar <- "dilation_n"; timeVar <- "cueRep_n"; groupVar <- "subject_block_f"
splitVar <- "RTcleaned_valence_f";
mod <- fit_gamm(data = modData, yVar = yVar, timeVar = timeVar, splitVar = splitVar, groupVar = groupVar,
                bsType = "fs", addARIMA = F, useScat = T)
png(paste0(dirs$plot, "final/FigureS11D.png"), width = plotWidthGAMM, height = plotHeight)
custom_plot_smooth(mod, ltyVec = c(1, 1, 2, 2), isPNG = F)
dev.off()

# -------------------------------------------------------------- #
### FigureS11F: GAMM plot: dilation ~ repeat x valence:

yVar <- "dilation_n"; timeVar <- "cueRep_n"; groupVar <- "subject_block_f"
splitVar <- "repeat_valence_Go_f";
mod <- fit_gamm(data = modData, yVar = yVar, timeVar = timeVar, splitVar = splitVar, groupVar = groupVar,
                bsType = "fs", addARIMA = F, useScat = T)
png(paste0(dirs$plot, "final/FigureS11F.png"), width = plotWidthGAMM, height = plotHeight)
custom_plot_smooth(mod, ltyVec = c(1, 1, 2, 2), isPNG = F)
dev.off()

# ============================================================================ #
#### FigureS12: Baselines: ####

yLimBase <- c(0, 29)
plotData$baseline12_n <- plotData$baseline_n + 12

# -------------------------------------------------------------- #
### FigureS12A: Bar plot: baseline ~ response:

p <- custom_barplot1(plotData, yVar = "baseline12_n", xVar = "response_f", subVar = "subject_f", yLim = yLimBase)
png(paste0(dirs$plot, "final/FigureS12A.png"), width = plotWidth, height = plotHeight)
print(p)
dev.off()

# -------------------------------------------------------------- #
### FigureS12B: Bar plot: baseline ~ ACC:

p <- custom_barplot1(plotData, yVar = "baseline12_n", xVar = "ACC_f", subVar = "subject_f", yLim = yLimBase)
png(paste0(dirs$plot, "final/FigureS12B.png"), width = plotWidth, height = plotHeight)
print(p)
dev.off()

# -------------------------------------------------------------- #
### FigureS12C: Bar plot: baseline ~ RTcleaned_fast_f:

p <- custom_barplot1(plotData, yVar = "baseline12_n", xVar = "RTcleaned_fast_f", subVar = "subject_f", yLim = yLimBase)
png(paste0(dirs$plot, "final/FigureS12C.png"), width = plotWidth, height = plotHeight)
print(p)
dev.off()

# -------------------------------------------------------------- #
### FigureS12D: GAMM plot: baseline ~ response:

## Read in behavioral data:
data <- read_behavior()
## Preprocessing: 
data <- wrapper_preprocessing(data)
## Add pupil data:
segMessage <- "StartMask"
beforeSamples <- 1000
afterSamples <- 2966
data <- add_pupil(data, segMessage = segMessage, beforeSamples = beforeSamples, afterSamples = afterSamples)
modData <- select_standardize(data)
modData$baseline12_n <- modData$baseline_n + 12
plotData <- modData

## Start plot:
yVar <- "baseline12_n"; timeVar <- "cueRep_n"; groupVar <- "subject_block_f"
splitVar <- "response_f";
mod <- fit_gamm(data = modData, yVar = yVar, timeVar = timeVar, splitVar = splitVar, groupVar = groupVar,
                bsType = "fs", addARIMA = F, useScat = T)
png(paste0(dirs$plot, "final/FigureS12D.png"), width = plotWidthGAMM, height = plotHeight)
custom_plot_smooth(mod, isPNG = F)
dev.off()

# -------------------------------------------------------------- #
### FigureS12E: GAMM plot: baseline ~ ACC:

yVar <- "baseline12_n"; timeVar <- "cueRep_n"; groupVar <- "subject_block_f"
splitVar <- "ACC_f";
mod <- fit_gamm(data = modData, yVar = yVar, timeVar = timeVar, splitVar = splitVar, groupVar = groupVar,
                bsType = "fs", addARIMA = F, useScat = T)
png(paste0(dirs$plot, "final/FigureS12E.png"), width = plotWidthGAMM, height = plotHeight)
custom_plot_smooth(mod, isPNG = F)
dev.off()

# -------------------------------------------------------------- #
### FigureS12F: GAMM plot: baseline ~ RTcleaned_fast_f:

yVar <- "baseline12_n"; timeVar <- "cueRep_n"; groupVar <- "subject_block_f"
splitVar <- "RTcleaned_fast_f";
mod <- fit_gamm(data = modData, yVar = yVar, timeVar = timeVar, splitVar = splitVar, groupVar = groupVar,
                bsType = "fs", addARIMA = F, useScat = T)
png(paste0(dirs$plot, "final/FigureS12F.png"), width = plotWidthGAMM, height = plotHeight)
custom_plot_smooth(mod, isPNG = F)
dev.off()

# ============================================================================ #
#### FigureS13: Pupil time course non-baseline corrected: ####

# -------------------------------------------------------------- #
### FigureS13A: Time course: dilation_z ~ reqAction_f * valence_f uncorrected:

## Read in behavioral data:
data <- read_behavior()
## Preprocessing: 
data <- wrapper_preprocessing(data)
## Settings for loading time courses:
segMessage <- "StartMask"
beforeSamples <- 1000
afterSamples <- 2966

zVar <- "respCond_f"
p <- plot_timecourse(data = data, zVar = zVar, beforeSamples = beforeSamples, afterSamples = afterSamples,
                     selLineType = c(1, 1, 2, 2), xLim = c(-500, 2966), yLim = c(920, 1060),
                     isBaselineCorGroup = F, addLegend = F)
png(paste0(dirs$plot, "final/FigureS13A.png"), width = 700, height = plotHeight)
print(p)
dev.off()

# -------------------------------------------------------------- #
### FigureS13B: Time course: dilation_z ~ reqAction_f * valence_f uncorrected, first half of blocks:

zVar <- "respCond_f"
p <- plot_timecourse(data = data, zVar = zVar, beforeSamples = beforeSamples, afterSamples = afterSamples,
                     selTrialVar = "firstHalfCueRep_f", selTrialVal = "first", # first second
                     selLineType = c(1, 1, 2, 2), xLim = c(-500, 2966), yLim = c(970, 1120),
                     isBaselineCorGroup = F, addLegend = F)
png(paste0(dirs$plot, "final/FigureS13B.png"), width = 700, height = plotHeight)
print(p)
dev.off()

# -------------------------------------------------------------- #
### FigureS13C: Time course: dilation_z ~ reqAction_f * valence_f uncorrected, second half of blocks:

zVar <- "respCond_f"
p <- plot_timecourse(data = data, zVar = zVar, beforeSamples = beforeSamples, afterSamples = afterSamples,
                     selTrialVar = "firstHalfCueRep_f", selTrialVal = "second", # first second
                     selLineType = c(1, 1, 2, 2), xLim = c(-500, 2966), yLim = c(880, 1010),
                     isBaselineCorGroup = F, addLegend = F)
png(paste0(dirs$plot, "final/FigureS13C.png"), width = 700, height = plotHeight)
print(p)
dev.off()

# -------------------------------------------------------------- #
### FigureS13D: GAMM plot: baseline ~ response x valence:

## Read in behavioral data:
data <- read_behavior()
## Preprocessing: 
data <- wrapper_preprocessing(data)
## Add pupil data:
segMessage <- "StartMask"
beforeSamples <- 1000
afterSamples <- 2966
data <- add_pupil(data, segMessage = segMessage, beforeSamples = beforeSamples, afterSamples = afterSamples)
modData <- select_standardize(data)
modData$baseline12_n <- modData$baseline_n + 12
plotData <- modData

## Start plot:
yVar <- "baseline12_n"; timeVar <- "cueRep_n"; groupVar <- "subject_block_f"
splitVar <- "respCond_f";
mod <- fit_gamm(data = modData, yVar = yVar, timeVar = timeVar, splitVar = splitVar, groupVar = groupVar,
                bsType = "fs", addARIMA = F, useScat = T)
png(paste0(dirs$plot, "final/FigureS13D.png"), width = plotWidthGAMM, height = plotHeight)
custom_plot_smooth(mod, ltyVec = c(1, 1, 2, 2), isPNG = F)
dev.off()

# -------------------------------------------------------------- #
### FigureS13E: GAMM plot: difference baseline ~ goValence:

## Start plot:
yVar <- "baseline12_n"; timeVar <- "cueRep_n"; groupVar <- "subject_block_f"
splitVar <- "respCond_f"; splitLevels <- c("G2A", "G2W")
mod <- fit_gamm_diff(data = modData, yVar = yVar, timeVar = timeVar, splitVar = splitVar, groupVar = groupVar, splitLevels = splitLevels,
                     bsType = "fs", addARIMA = F, useScat = T)
png(paste0(dirs$plot, "final/FigureS13E.png"), width = plotWidthGAMM, height = plotHeight)
custom_plot_diff(mod, isPNG = F)
dev.off()

# ============================================================================ #
#### FigureS14: Outcomes: ####

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
# afterSamples <- 2500
data <- add_pupil(data, segMessage = segMessage, beforeSamples = beforeSamples, afterSamples = afterSamples)

modData <- select_standardize(data)
plotData <- modData

## 1D plots:
yLimDil <- c(0, 20)

# -------------------------------------------------------------- #
### FigureS14A: bar plot: dilation ~ outcome valence:

p <- custom_barplot1(plotData, yVar = "dilation_n", xVar = "outcome_rel_f", subVar = "subject_f", yLim = yLimDil)
png(paste0(dirs$plot, "final/FigureS14A.png"), width = plotWidth, height = plotHeight)
print(p)
dev.off()

# -------------------------------------------------------------- #
### FigureS14B: bar plot: dilation ~ outcome displayed:

p <- custom_barplot1(plotData, yVar = "dilation_n", xVar = "outcome_f", subVar = "subject_f", yLim = yLimDil)
png(paste0(dirs$plot, "final/FigureS14B.png"), width = plotWidth, height = plotHeight)
print(p)
dev.off()

# -------------------------------------------------------------- #
### FigureS14C: bar plot: dilation ~ outcome interpreted:

p <- custom_barplot1(plotData, yVar = "dilation_n", xVar = "outcome_all_short_f", subVar = "subject_f", yLim = yLimDil)
png(paste0(dirs$plot, "final/FigureS14C.png"), width = plotWidth, height = plotHeight)
print(p)
dev.off()

# -------------------------------------------------------------- #
## Settings for loading time courses:

segMessage <- "StartOutcome"
beforeSamples <- 1000
afterSamples <- 2500

# -------------------------------------------------------------- #
### FigureS14D: time course: dilation ~ outcome valence:

zVar <- "outcome_rel_f"
p <- plot_timecourse(data = data, zVar = zVar, segMessage = segMessage, afterSamples = afterSamples,
                     isBaselineCorGroup = T, addLegend = F, yLim = c(-30, 60))
png(paste0(dirs$plot, "final/FigureS14D.png"), width = plotWidthTime, height = plotHeight)
print(p)
dev.off()

# -------------------------------------------------------------- #
### FigureS14E: time course: dilation ~ outcome valence:

zVar <- "outcome_f"
p <- plot_timecourse(data = data, zVar = zVar, segMessage = segMessage, afterSamples = afterSamples,
                     isBaselineCorGroup = T, addLegend = F, yLim = c(-30, 60))
png(paste0(dirs$plot, "final/FigureS14E.png"), width = plotWidthTime, height = plotHeight)
print(p)
dev.off()

# -------------------------------------------------------------- #
### FigureS14F: time course: dilation ~ outcome valence:

zVar <- "outcome_all_f"
p <- plot_timecourse(data = data, zVar = zVar, segMessage = segMessage, afterSamples = afterSamples,
                     isBaselineCorGroup = T, addLegend = F, yLim = c(-30, 60))
png(paste0(dirs$plot, "final/FigureS14F.png"), width = plotWidthTime, height = plotHeight)
print(p)
dev.off()

# ============================================================================ #
#### FigureS15: Outcomes x action: ####

yLimDil <- c(0, 21)

# -------------------------------------------------------------- #
### FigureS15A: bar plot: dilation ~ outcome displayed x response:

p <- custom_barplot2(plotData, yVar = "dilation_n", xVar = "outcome_f", zVar = "response_f", subVar = "subject_f", yLim = yLimDil)
png(paste0(dirs$plot, "final/FigureS15A.png"), width = plotWidth, height = plotHeight)
print(p)
dev.off()

# -------------------------------------------------------------- #
## Settings for loading time courses:

segMessage <- "StartOutcome"
beforeSamples <- 1000
afterSamples <- 2500

# -------------------------------------------------------------- #
### FigureS15B: bar plot: dilation ~ outcome displayed x response, baseline-corrected:

zVar <- "outcome_response_f"
p <- plot_timecourse(data = data, zVar = zVar, segMessage = segMessage, afterSamples = afterSamples,
                     selLineType = c(1, 2, 1, 2, 1, 2), 
                     isBaselineCorGroup = T, addLegend = F)
png(paste0(dirs$plot, "final/FigureS15B.png"), width = plotWidthTime, height = plotHeight)
print(p)
dev.off()

# -------------------------------------------------------------- #
### FigureS15C: bar plot: dilation ~ outcome displayed x response, uncorrected:

zVar <- "outcome_response_f"
p <- plot_timecourse(data = data, zVar = zVar, segMessage = segMessage, afterSamples = afterSamples,
                     selLineType = c(1, 2, 1, 2, 1, 2), 
                     isBaselineCorGroup = F, addLegend = F)
png(paste0(dirs$plot, "final/FigureS15C.png"), width = plotWidthTime, height = plotHeight)
print(p)
dev.off()

# END OF FILE.
