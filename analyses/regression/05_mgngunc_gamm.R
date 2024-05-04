#!/usr/bin/env Rscript
# ============================================================================ #
## 05_mgngunc_gamm.R
## MGNGUncertainty fot generalized additive mixed-effects models on pupil data.
# Johannes Algermissen, 2023.

rm(list = ls())

# ============================================================================ #
#### 00) Set directories, load packages and custom functions: ####

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

## Set factors to treatment by default:
options(contrasts = c("contr.treatment", "contr.poly")) # need treatment coding for factors in GAMMs.

# ------------------------------------------------- #
## Load custom functions:

source(paste0(dirs$codeDir, "functions/00_mgngunc_functions_regression.R")) # Load functions

# ============================================================================ #
# ============================================================================ #
# ============================================================================ #
# ============================================================================ #
#### 01) Read in behavioral and pupil data: ####

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
names(data)[grepl("dilation", names(data))]

# ----------------------------------- #
## Select data, standardize:

modData <- select_standardize(data)
modData$baseline12_n <- modData$baseline_n + 12 # remove negative values for baseline by adding intercept
plotData <- modData
cat("Finished pre-processing! :-)\n")

# ============================================================================ #
# ============================================================================ #
# ============================================================================ #
#### 02) Fit time course for each CONDITION using GAMMs: ####
# http://jacolienvanrij.com/Tutorials/GAMM.html
# http://jacolienvanrij.com/PupilAnalysis/SupplementaryMaterials-2.html

# ---------------------------------------------------------------------------- #
#### 02a) Select variables: ####

## Dependent variable (must be numeric):

yVar <- "baseline12_n"

yVar <- "dilation_n"

## Time variable (must be numeric):
timeVar <- "cueRep_n"

## Split (modulation) variable (must be factor):
splitVar <- "respCond_f"
splitVar <- "condition_f"

splitVar <- "response_f"
splitVar <- "valence_f"
splitVar <- "arousal_f"

splitVar <- "ACC_f"
splitVar <- "ACC_response_f"
splitVar <- "ACC_valence_f"
splitVar <- "ACC_valence_Go_f"

splitVar <- "repeat_f"
splitVar <- "repeat_response_f"
splitVar <- "repeat_valence_f"
splitVar <- "repeat_valence_Go_f"

splitVar <- "RTcleaned_fast_f"
splitVar <- "RTcleaned_valence_f"
splitVar <- "RTcleaned_valence_uncor_f"
splitVar <- "RTcleaned_ACC_f"
splitVar <- "RTcleaned_ACC_uncor_f"

## Grouping variable (must be factor):
groupVar <- "subject_block_f"

# ============================================================================ #
#### 02b) Fit models: ####

## Fit 4 models with/ without ARIMA and/or scat:
mod <- fit_gamm(data = modData, yVar = yVar, timeVar = timeVar, splitVar = splitVar, groupVar = groupVar,
                bsType = "fs", addARIMA = F, useScat = F); mod4A <- mod
custom_plot_smooth(mod4A, ltyVec = c(1, 1, 2, 2), isPNG = F)
mod <- fit_gamm(data = modData, yVar = yVar, timeVar = timeVar, splitVar = splitVar, groupVar = groupVar,
                bsType = "fs", addARIMA = T, useScat = F); mod4B <- mod
custom_plot_smooth(mod4B, ltyVec = c(1, 1, 2, 2), isPNG = F)
mod <- fit_gamm(data = modData, yVar = yVar, timeVar = timeVar, splitVar = splitVar, groupVar = groupVar,
                bsType = "fs", addARIMA = F, useScat = T); mod4C <- mod
custom_plot_smooth(mod4C, ltyVec = c(1, 1, 2, 2), isPNG = F)
mod <- fit_gamm(data = modData, yVar = yVar, timeVar = timeVar, splitVar = splitVar, groupVar = groupVar,
                bsType = "fs", addARIMA = T, useScat = T); mod4D <- mod
custom_plot_smooth(mod4D, ltyVec = c(1, 1, 2, 2), isPNG = F)
cat("Finished fitting all models! :-)\n")

# ============================================================================ #
#### 02c) Model comparison: #####

# ----------------------------------- #
## a) Identify best model based on fREML:
c(mod4A$gcv.ubre, mod4B$gcv.ubre, mod4C$gcv.ubre, mod4D$gcv.ubre)
minIdx <- which.min(c(mod4A$gcv.ubre, mod4B$gcv.ubre, mod4C$gcv.ubre, mod4D$gcv.ubre)); minIdx
bestModName <- paste0("mod4", LETTERS[minIdx]); bestModName
# bestModName <- "mod4C"

## Plot:
custom_plot_smooth(eval(parse(text = bestModName)), ltyVec = c(1, 1, 2, 2), isPNG = T)
custom_plot_smooth(eval(parse(text = bestModName)), ltyVec = c(1, 1, 2, 2), isPNG = F)

## Print:
summary(eval(parse(text = bestModName)))

# ----------------------------------- #
## b) Identify best model based on AIC:
AIC(mod4A, mod4B, mod4C, mod4D)
AIC(mod4A, mod4B, mod4C, mod4D)[, 2]
minIdx <- which.min(AIC(mod4A, mod4B, mod4C, mod4D)[, 2]); minIdx
bestModName <- paste0("mod4", LETTERS[minIdx]); bestModName

## Plot:
custom_plot_smooth(eval(parse(text = bestModName)), ltyVec = c(1, 1, 2, 2), isPNG = T)
custom_plot_smooth(eval(parse(text = bestModName)), ltyVec = c(1, 1, 2, 2), isPNG = F)

## Print:
summary(eval(parse(text = bestModName)))
          
# ============================================================================ #
#### 02d) Model criticism: ####
# http://jacolienvanrij.com/Tutorials/GAMM.html
# http://jacolienvanrij.com/PupilAnalysis/SupplementaryMaterials-2.html
# https://journals.sagepub.com/doi/10.1177/2331216519832483
# https://osf.io/wgc4f/wiki/mgcv:%20model%20selection/
# https://stat.ethz.ch/R-manual/R-devel/library/mgcv/html/family.mgcv.html

## Summary:
summary(mod4C)
AIC(mod4C)
mod4C$gcv.ubre# fREML score
compareML(mod4C, mod4D)

## Model comparison:
compareML(mod4A, mod4B)
compareML(mod4B, mod4C)
compareML(mod4C, mod4D)
quickCImgcv(mod4)

## Re-select model for criticism:
mod <- mod4A
mod <- mod4B
mod <- mod4C
mod <- mod4D

## Plot:
custom_plot_smooth(mod, ltyVec = c(1, 1, 2, 2), isPNG = F)
custom_plot_diff(mod, isPNG = F)

# ----------------------------------- #
## Summary:

summary(mod)
mod$gcv.ubre

## Check model convergence:
gam.check(mod)

## Residuals:
check_resid(mod, split_pred = groupVar) # overview all major plots
par(ask = F) # switch of asking after check_resid

## QQ-plot:
par(ask = F) # switch of asking after check_resid
check_resid(mod, split_pred = groupVar, select = 1, ask = F) # distribution residuals against normal distribution
qq.gam(mod, cex = 5)

## Distribution of residuals:
par(ask = F) # switch of asking after check_resid
check_resid(mod, split_pred = groupVar, select = 2, ask = F) # distribution residuals against normal distribution
densityplot(mod$residuals)
sd(mod$residuals)
max(mod$residuals) - min(mod$residuals)

## Residuals vs. fitted values:
plot(mod$residuals, mod$fitted.values)

## Empirical vs. fitted values:
# yVec <- data[which(!(is.na(data[, yVar]))), yVar]
selData <- modData[which(data[, splitVar] %in% splitLevels), ]
yVec <- selData[which(!(is.na(selData[, yVar]))), yVar]
# length(yVec)
# length(mod$fitted.values)
plot(yVec, mod$fitted.values)
abline(a = 0, b = 1, col = "blue", lty = 2, lwd = 3) # identity line
abline(lm(mod$fitted.values ~ yVec), col = "red", lwd = 5) # model prediction
cor(yVec, mod$fitted.values)

## Auto-correlation in residuals:
acf(resid(mod)) # --> very low
acf(resid(mod), plot = F)
# acf_resid(mod, split_pred = "subject_block_f", main="ACF resid(mod)")

# ============================================================================ #
# ============================================================================ #
# ============================================================================ #
# ============================================================================ #
# ============================================================================ #
# ============================================================================ #
# ============================================================================ #
# ============================================================================ #
#### 03) Model DIFFERENCE between cue conditions as ordered factor using GAMMs: ####
# http://jacolienvanrij.com/Tutorials/GAMM.html#example-ordered-factor

# ============================================================================ #
#### 03a) Select variables and levels: ####

## Dependent variable (must be numeric):
yVar <- "baseline12_n"

yVar <- "dilation_n"

## Time variable (must be numeric):
timeVar <- "cueRep_n"

splitVar <- "respCond_f"; splitLevels <- c("G2A", "G2W")
splitVar <- "condition_f"; splitLevels <- c("G2A", "G2W")
splitVar <- "response_f"; splitLevels <- levels(data[, splitVar])
splitVar <- "valence_f"; splitLevels <- levels(data[, splitVar])

splitVar <- "ACC_f"; splitLevels <- levels(data[, splitVar])
splitVar <- "ACC_response_f"; splitLevels <- c("corGo", "incGo")
splitVar <- "ACC_response_f"; splitLevels <- c("corNoGo", "incNoGo")
splitVar <- "ACC_valence_f"; splitLevels <- c("cor2A", "cor2W")
splitVar <- "ACC_valence_f"; splitLevels <- c("inc2A", "inc2W")
splitVar <- "ACC_valence_Go_f"; splitLevels <- c("cor2A", "cor2W")
splitVar <- "ACC_valence_Go_f"; splitLevels <- c("inc2A", "inc2W")

splitVar <- "repeat_f"; splitLevels <- levels(data[, splitVar])
splitVar <- "repeat_response_f"; splitLevels <- c("repGo", "swiGo")
splitVar <- "repeat_response_f"; splitLevels <- c("repNoGo", "swiNoGo")
splitVar <- "repeat_valence_f"; splitLevels <- c("rep2A", "rep2W")
splitVar <- "repeat_valence_f"; splitLevels <- c("swi2A", "swi2W")
splitVar <- "repeat_valence_Go_f"; splitLevels <- c("rep2A", "rep2W")
splitVar <- "repeat_valence_Go_f"; splitLevels <- c("swi2A", "swi2W")

splitVar <- "RTcleaned_fast_f"; splitLevels <- levels(data[, splitVar])
splitVar <- "RTcleaned_valence_f"; splitLevels <- c("F2A", "F2W")
splitVar <- "RTcleaned_valence_f"; splitLevels <- c("S2A", "S2W")

splitVar <- "RTcleaned_ACC_f"; splitLevels <- c("Finc", "Fcor")
splitVar <- "RTcleaned_ACC_f"; splitLevels <- c("Sinc", "Scor")
splitVar <- "RTcleaned_ACC_f"; splitLevels <- c("Scor", "Fcor")
splitVar <- "RTcleaned_ACC_f"; splitLevels <- c("Sinc", "Finc")

splitVar <- "RTcleaned_ACC_uncor_f"; splitLevels <- c("Finc", "Fcor")
splitVar <- "RTcleaned_ACC_uncor_f"; splitLevels <- c("Sinc", "Scor")
splitVar <- "RTcleaned_ACC_uncor_f"; splitLevels <- c("Scor", "Fcor")
splitVar <- "RTcleaned_ACC_uncor_f"; splitLevels <- c("Sinc", "Finc")

## Grouping variable (must be factor):
groupVar <- "subject_block_f"

# ============================================================================ #
#### 03b) Fit difference models: ####

## Fit 4 models with/ without ARIMA and/or scat:
mod <- fit_gamm_diff(data = modData, yVar = yVar, timeVar = timeVar, splitVar = splitVar, groupVar = groupVar, splitLevels = splitLevels,
                     bsType = "fs", addARIMA = F, useScat = F); mod4A <- mod
custom_plot_diff(mod4A, isPNG = F)
mod <- fit_gamm_diff(data = modData, yVar = yVar, timeVar = timeVar, splitVar = splitVar, groupVar = groupVar, splitLevels = splitLevels,
                     bsType = "fs", addARIMA = T, useScat = F); mod4B <- mod
custom_plot_diff(mod4B, isPNG = F)
mod <- fit_gamm_diff(data = modData, yVar = yVar, timeVar = timeVar, splitVar = splitVar, groupVar = groupVar, splitLevels = splitLevels,
                     bsType = "fs", addARIMA = F, useScat = T); mod4C <- mod
custom_plot_diff(mod4C, isPNG = F)
mod <- fit_gamm_diff(data = modData, yVar = yVar, timeVar = timeVar, splitVar = splitVar, groupVar = groupVar, splitLevels = splitLevels,
                     bsType = "fs", addARIMA = T, useScat = T); mod4D <- mod
custom_plot_diff(mod4D, isPNG = F)
cat("Finished fitting all models! :-)\n")

# ============================================================================ #
#### 03c) Model comparison: #####
# See https://cran.r-project.org/web/packages/itsadug/vignettes/test.html:
# "An alternative test is AIC, but when an AR1 model is included, AIC does not provide a reliable test."

# ---------------------------------- #
## Identify best model based on fREML:
c(mod4A$gcv.ubre, mod4B$gcv.ubre, mod4C$gcv.ubre, mod4D$gcv.ubre)
minIdx <- which.min(c(mod4A$gcv.ubre, mod4B$gcv.ubre, mod4C$gcv.ubre, mod4D$gcv.ubre)); minIdx
bestModName <- paste0("mod4", LETTERS[minIdx]); bestModName

## Plot:
custom_plot_diff(eval(parse(text = bestModName)), isPNG = T)
custom_plot_diff(eval(parse(text = bestModName)), isPNG = F)

## Print:
quickCImgcv(eval(parse(text = bestModName)))
summary(eval(parse(text = bestModName)))

# ---------------------------------- #
## Identify best model based on AIC:
AIC(mod4A, mod4B, mod4C, mod4D)
AIC(mod4A, mod4B, mod4C, mod4D)[, 2]
minIdx <- which.min(AIC(mod4A, mod4B, mod4C, mod4D)[, 2]); minIdx
bestModName <- paste0("mod4", LETTERS[minIdx]); bestModName

## Plot:
custom_plot_diff(eval(parse(text = bestModName)), isPNG = T)
custom_plot_diff(eval(parse(text = bestModName)), isPNG = F)

# Print:
quickCImgcv(eval(parse(text = bestModName)))
summary(eval(parse(text = bestModName)))

# ============================================================================ #
#### 03d) Model criticism: ####
# http://jacolienvanrij.com/Tutorials/GAMM.html
# http://jacolienvanrij.com/PupilAnalysis/SupplementaryMaterials-2.html
# https://journals.sagepub.com/doi/10.1177/2331216519832483
# https://osf.io/wgc4f/wiki/mgcv:%20model%20selection/
# https://stat.ethz.ch/R-manual/R-devel/library/mgcv/html/family.mgcv.html

# ------------------------ #
## Effect of ARIMA:
compareML(mod4A, mod4B) # effect of ARIMA
acf(resid(mod4A)) # without ARIMA
acf(resid(mod4B)) # with ARIMA
custom_plot_diff(mod4A, isPNG = F)
custom_plot_diff(mod4B, isPNG = F)

## Effect of  SCAT:
compareML(mod4A, mod4C) # effect of scat
check_resid(mod4A, split_pred = groupVar, select = 2, ask = F) # overview all major plots
check_resid(mod4C, split_pred = groupVar, select = 2, ask = F) # overview all major plots
custom_plot_diff(mod4A, isPNG = F)
custom_plot_diff(mod4C, isPNG = F)

## Interactions:
compareML(mod4C, mod4D) # effect ARIMA on top of scat
check_resid(mod4C, split_pred = groupVar, select = 2, ask = F) # overview all major plots
check_resid(mod4D, split_pred = groupVar, select = 2, ask = F) # overview all major plots
custom_plot_diff(mod4C, isPNG = F)
custom_plot_diff(mod4D, isPNG = F)

## Consecutive comparisons:
compareML(mod4A, mod4B)
compareML(mod4A, mod4B)
compareML(mod4B, mod4C)
compareML(mod4C, mod4D)
AIC(mod4A, mod4B, mod4C, mod4D)
c(mod4A$gcv.ubre, mod4B$gcv.ubre, mod4C$gcv.ubre, mod4D$gcv.ubre)

# ---------------------------------------------- #
## Print b + 95%-CI for focal effect:

quickCImgcv(mod4A)
quickCImgcv(mod4B)
quickCImgcv(mod4C)
quickCImgcv(mod4D)

# ---------------------------------------------- #
### Plot:

## Automatic plots:
custom_plot_diff(mod, isPNG = F)
custom_plot_diff(mod, isPNG = T)
custom_plot_diff(mod, isPNG = F)

# END OF FILE.
