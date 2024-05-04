#!/usr/bin/env Rscript
# ============================================================================ #
## 07_mgngunc_power.R
## MGNGUncertainty estimate post-hoc power for responses, RTs, baselines, dilations
## using the approach described by Aarts, Verhage, Veenvliet, Dolan, & van der Sluis (2014). https://doi.org/10.1038/nn.3648.
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

## Load packages:
library(lme4)
library(pwr)

# ------------------------------------------------- #
## Load custom functions:

source(paste0(dirs$codeDir, "functions/00_mgngunc_functions_regression.R")) # Load function

# ============================================================================ #
# ============================================================================ #
# ============================================================================ #
# ============================================================================ #
#### 01) Read in data: ####

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
afterSamples <- 2966
data <- add_pupil(data, segMessage = segMessage, afterSamples = afterSamples)

# ----------------------------------- #
## Select data, standardize variables:

modData <- data
# modData <- subset(data, !(subject_n %in% c(52)))

# ============================================================================ #
#### 02) Estimate ICC: ####

# ----------------------------------------------------- #
### responses:

formula <- "response_n ~ 1 + (1|subject_f)"
mod <- glmer(formula = formula, data = modData, family = binomial(),
             control = glmerControl(optCtrl = list(maxfun = 1e+9), calc.derivs = F, optimizer = c("bobyqa")))
summary(mod)

## Extract variances and compute ICC:
varData <- as.data.frame(VarCorr(mod))
ICC <- varData$vcov[1] / (varData$vcov[1] + (pi^2)/3); ICC # 0.03501369

# Residual variance is fixed to pi^2/3, see 
# https://stats.stackexchange.com/questions/62770/calculating-icc-for-random-effects-logistic-regression
# https://stats.stackexchange.com/questions/128750/residual-variance-for-glmer

## response: 0.03501369

# ----------------------------------------------------- #
### RTs:

formula <- "RTcleaned_n ~ 1 + (1|subject_f)"
mod <- lmer(formula = formula, data = modData,
            control = lmerControl(optCtrl = list(maxfun = 1e+9), calc.derivs = F, optimizer = c("bobyqa")))
summary(mod)

varData <- as.data.frame(VarCorr(mod))
ICC <- varData$vcov[1] / (varData$vcov[1] + varData$vcov[2]); ICC # 0.1397013
# RTs: 0.1397013

# ----------------------------------------------------- #
### Baseline:

formula <- "baseline_n ~ 1 + (1|subject_f)"
mod <- lmer(formula = formula, data = modData,
            control = lmerControl(optCtrl = list(maxfun = 1e+9), calc.derivs = F, optimizer = c("bobyqa")))
summary(mod)

varData <- as.data.frame(VarCorr(mod))
ICC <- varData$vcov[1] / (varData$vcov[1] + varData$vcov[2]); ICC # 0.01835622
# Baseline: 0.01835622

# ----------------------------------------------------- #
### Dilation:

formula <- "dilation_n ~ 1 + (1|subject_f)"
mod <- lmer(formula = formula, data = modData,
            control = lmerControl(optCtrl = list(maxfun = 1e+9), calc.derivs = F, optimizer = c("bobyqa")))
summary(mod)

varData <- as.data.frame(VarCorr(mod))
ICC <- varData$vcov[1] / (varData$vcov[1] + varData$vcov[2]); ICC # 0.1689371
# Dilation: 0.1689371

# ============================================================================ #
#### 03) Compute r achieved with 80% given FIXED number of subjects: ####

## Fixed:
# ICC <- 0.03501369 # response
# ICC <- 0.1397013 # RT
# ICC <- 0.01835622 # baseline
# ICC <- 0.1689371 # dilation

# ------------------------------------ #
### 03a) Compute nTrials:
nCue <- 4
nRep <- 16
nTrialBlock <- nCue * nRep
nBlock <- 4
nTrial <- nBlock * nTrialBlock; nTrial

# ------------------------------------ #
### 03b) Compute effective sample size:
nSub <- 35
nData <- nSub*nTrial; nData
Neff <- nData/ (1 + (nSub - 1) * ICC); Neff

# responses: 4090.455
# RTs: 1558.303
# baselines: 5516.863
# dilations: 1328.616

# ------------------------------------ #
### 03c) Estimate smallest effect size (standardized regression coefficient) for which we have >=80% power:
power <- 0.80
pwr.r.test(n = Neff, sig.level = 0.05, power = power)

# responses: 0.04379679
# RTs: 0.07089884
# baselines: 0.03768582
# dilations: 0.07678902

# END OF FILE.
