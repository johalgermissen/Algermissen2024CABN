#!/usr/bin/env Rscript
# ============================================================================ #
## 02_mgngunc_regression.R
## MGNGUncertainty regression models on behaviour.
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
# names(data)
# table(data$subject_n)

# ----------------------------------- #
## Add pupil data:

segMessage <- "StartMask"; beforeSamples <- 1000; afterSamples <- 2966
data <- add_pupil(data, segMessage = segMessage, beforeSamples = beforeSamples, afterSamples = afterSamples)

# ----------------------------------- #
## Select data, standardize variables:

modData <- select_standardize(data)

# ============================================================================ #
#### 02a) Response: ####

## Main effects:
formula <- "response_n ~ reqAction_f + (reqAction_f|subject_f)"
formula <- "response_n ~ valence_f + (valence_f|subject_f)"
formula <- "response_n ~ dilation_z + (dilation_z|subject_f)"

## 2-way interactions:
formula <- "response_n ~ reqAction_f * valence_f + (reqAction_f * valence_f|subject_f)"
formula <- "response_n ~ reqAction_f * arousal_f + (reqAction_f * arousal_f|subject_f)"
formula <- "response_n ~ valence_f * arousal_f + (valence_f * arousal_f|subject_f)"

formula <- "response_n ~ condition_f * arousal_f + (condition_f * arousal_f|subject_f)"
formula <- "response_n ~ condition_f * dilation_z + (condition_f * dilation_z|subject_f)"

## 3-way interaction:
formula <- "response_n ~ reqAction_f * valence_f * arousal_f + (reqAction_f * valence_f * arousal_f|subject_f)"
formula <- "response_n ~ reqAction_f * valence_f * dilation_z + (reqAction_f * valence_f * dilation_z|subject_f)"

## 4-way interaction:
formula <- "response_n ~ reqAction_f * valence_f * arousal_f * dilation_z + (reqAction_f * valence_f * arousal_f * dilation_z|subject_f)"

formula <- "response_n ~ reqAction_f * valence_f * dilation_z + RTcleaned_z + (reqAction_f * valence_f * dilation_z + RTcleaned_z|subject_f)"
formula <- "response_n ~ reqAction_f * valence_f * dilation_z * RTcleaned_z + (reqAction_f * valence_f * dilation_z * RTcleaned_z|subject_f)"

### Fit logistic regression:
# mod <- glmer(formula = formula, data = modData, family = binomial(),
#              control = glmerControl(optCtrl = list(maxfun = 1e+9), calc.derivs = F, optimizer = c("bobyqa")))
# summary(mod)

## Read existing model back in:
mod <- fit_lmem(formula)
quickCI(mod, nRound = 3)
mod <- fit_lmem(formula, useLRT = T)

## Plots:
plot(effect("dilation_z", mod))
plot(effect("reqAction_f:dilation_z", mod))
plot(effect("valence_f:dilation_z", mod))
plot(effect("valence_f:dilation_z", mod, x.var = "valence_f"))

plot(effect("condition_f:dilation_z", mod))
plot(effect("condition_f:dilation_z", mod, x.var = "condition_f"))

## Read manually:
fitType <- "LRT"
modFamily <- "binom"
optimizer <- "nlminbwrap" # bobyqa nloptwrap nlminbwrap
modName <- paste0(fitType, "_", modFamily, "_", formula2handle(formula), "_", optimizer, ".rds")
mod <- readRDS(paste0(dirs$modelDir, modName))
anova(mod)

# ============================================================================ #
#### 02b) RTs: ####

## Main effects:
formula <- "RTcleaned_z ~ reqAction_f + (reqAction_f|subject_f)"
formula <- "RTcleaned_z ~ valence_f + (valence_f|subject_f)"
formula <- "RTcleaned_z ~ reqCongruency_f + (reqCongruency_f|subject_f)"
formula <- "RTcleaned_z ~ actCongruency_f + (actCongruency_f|subject_f)"
formula <- "RTcleaned_z ~ arousal_f + (arousal_f|subject_f)"

## 2-way interactions:
formula <- "RTcleaned_z ~ reqAction_f * valence_f + (reqAction_f * valence_f|subject_f)"
formula <- "RTcleaned_z ~ reqAction_f * arousal_f + (reqAction_f * arousal_f|subject_f)"
formula <- "RTcleaned_z ~ valence_f * arousal_f + (valence_f * arousal_f|subject_f)"

formula <- "RTcleaned_z ~ condition_f * arousal_f + (condition_f * arousal_f|subject_f)"
formula <- "RTcleaned_z ~ condition_f * dilation_z + (condition_f * dilation_z|subject_f)"

## 3-way interaction:
formula <- "RTcleaned_z ~ reqAction_f * valence_f * arousal_f + (reqAction_f * valence_f * arousal_f|subject_f)"
formula <- "RTcleaned_z ~ reqAction_f * valence_f * dilation_z + (reqAction_f * valence_f * dilation_z|subject_f)"

## 4-way interaction:
formula <- "RTcleaned_z ~ reqAction_f * valence_f * arousal_f * dilation_z + (reqAction_f * valence_f * arousal_f * dilation_z|subject_f)"

### Fit linear regression:
selData <- modData
selData <- subset(modData, cueRep_n %in% c(12:16))
formula <- "RTcleaned_z ~ ACC_f + (ACC_f|subject_f)"
mod <- lmer(formula = formula, data = selData,
            control = lmerControl(optCtrl = list(maxfun = 1e+9), calc.derivs = F, optimizer = c("bobyqa")))
summary(mod)
plot(effect("ACC_f", mod))

mod <- fit_lmem(formula)
quickCI(mod, nRound = 3)
mod <- fit_lmem(formula, useLRT = T)

plot(effect("ACC_f", mod))

plot(effect("reqAction_f:valence_f:arousal_f ", mod))
plot(effect("reqAction_f:valence_f:arousal_f ", mod, x.var = "arousal_f"))

plot(effect("dilation_z", mod))
plot(effect("reqAction_f:dilation_z", mod))
plot(effect("valence_f:dilation_z", mod))
plot(effect("valence_f:dilation_z", mod, x.var = "valence_f"))

plot(effect("condition_f:dilation_z", mod))
plot(effect("condition_f:dilation_z", mod, x.var = "condition_f"))

plot(effect("reqAction_f:valence_f:arousal_f", mod))
plot(effect("reqAction_f:valence_f:arousal_f", mod, x.var = "arousal_f"))
plot(effect("reqAction_f:valence_f:arousal_f:dilation_z", mod))
plot(effect("reqAction_f:valence_f:arousal_f:dilation_z", mod, x.var = "valence_f"))

# ============================================================================ #
#### 02c) Baselines: ####

## Main effects:
formula <- "baseline_z ~ reqAction_f + (reqAction_f|subject_f)"
formula <- "baseline_z ~ valence_f + (valence_f|subject_f)"
formula <- "baseline_z ~ arousal_f + (arousal_f|subject_f)"
formula <- "baseline_z ~ reqCongruency_f + (reqCongruency_f|subject_f)"
formula <- "baseline_z ~ actCongruency_f + (actCongruency_f|subject_f)"

formula <- "baseline_z ~ trialnr_z + (trialnr_z|subject_f)"
formula <- "baseline_z ~ trialnr_block_z + (trialnr_block_z|subject_f)"
formula <- "baseline_z ~ cueRep_z + (cueRep_z|subject_f)"

formula <- "baseline_z ~ response_f + (response_f|subject_f)"
formula <- "baseline_z ~ ACC_f + (ACC_f|subject_f)"
formula <- "baseline_z ~ RTcleaned_z + (RTcleaned_z|subject_f)"

## 2-way interactions:
formula <- "baseline_z ~ reqAction_f * valence_f + (reqAction_f * valence_f|subject_f)"
formula <- "baseline_z ~ response_f * valence_f + (response_f * valence_f|subject_f)"

formula <- "baseline_z ~ condition_f * arousal_f + (condition_f * arousal_f|subject_f)"

## 3-way interaction:
formula <- "baseline_z ~ reqAction_f * valence_f * arousal_f + (reqAction_f * valence_f * arousal_f|subject_f)"

### Fit linear regression:
# mod <- lmer(formula = formula, data = modData,
#             control = lmerControl(optCtrl = list(maxfun = 1e+9), calc.derivs = F, optimizer = c("bobyqa")))
# summary(mod)

## S06:
formula <- "baseline_z ~ response_f + (response_f|subject_f)"
formula <- "baseline_z ~ ACC_f + (ACC_f|subject_f)"
formula <- "baseline_z ~ RTcleaned_fast_f + (RTcleaned_fast_f|subject_f)"

formula <- "baseline_z ~ response_f * cueRep_z + (response_f * cueRep_z|subject_f)"
formula <- "baseline_z ~ ACC_f * cueRep_z + (ACC_f * cueRep_z|subject_f)"
formula <- "baseline_z ~ RTcleaned_fast_f * cueRep_z + (RTcleaned_fast_f * cueRep_z|subject_f)"

## Fit model:
mod <- fit_lmem(formula)
quickCI(mod, nRound = 3)
modLRT <- fit_lmem(formula, useLRT = T)

# ============================================================================ #
#### 02d) Dilations: ####

## Main effects:
formula <- "dilation_z ~ reqAction_f + (reqAction_f|subject_f)"
formula <- "dilation_z ~ response_f + (response_f|subject_f)"
formula <- "dilation_z ~ valence_f + (valence_f|subject_f)"
formula <- "dilation_z ~ arousal_f + (arousal_f|subject_f)"
formula <- "dilation_z ~ reqCongruency_f + (reqCongruency_f|subject_f)"
formula <- "dilation_z ~ actCongruency_f + (actCongruency_f|subject_f)"

## 2-way interactions:
formula <- "dilation_z ~ reqAction_f * valence_f + (reqAction_f * valence_f|subject_f)"
formula <- "dilation_z ~ reqAction_f * response_f + (reqAction_f * response_f|subject_f)"

formula <- "dilation_z ~ response_f * valence_f + (response_f * valence_f|subject_f)"
formula <- "dilation_z ~ valence_f * response_f + (valence_f * response_f|subject_f)"

formula <- "dilation_n ~ response_f * valence_f + (response_f * valence_f|subject_f)"
formula <- "dilation_n ~ valence_f * response_f + (valence_f * response_f|subject_f)"

formula <- "dilation_z ~ condition_f * arousal_f + (condition_f * arousal_f|subject_f)"

formula <- "dilation_z ~ RTcleaned_z * valence_f + (RTcleaned_z * valence_f|subject_f)"
formula <- "dilation_z ~ RTcleaned_z * arousal_f + (RTcleaned_z * arousal_f|subject_f)"
formula <- "dilation_z ~ RTcleaned_z * ACC_f + (RTcleaned_z * ACC_f|subject_f)"
formula <- "dilation_z ~ RTcleaned_z * repeat_f + (RTcleaned_z * repeat_f|subject_f)"

formula <- "dilation_z ~ ACC_f * RTcleaned_fast_f + (ACC_f * RTcleaned_fast_f|subject_f)"

## 3-way interaction:
formula <- "dilation_z ~ reqAction_f * valence_f * arousal_f + (reqAction_f * valence_f * arousal_f|subject_f)"
formula <- "dilation_z ~ response_f * valence_f * arousal_f + (response_f * valence_f * arousal_f|subject_f)"

formula <- "dilation_z ~ reqAction_f * valence_f * RTcleaned_z + (reqAction_f * valence_f * RTcleaned_z|subject_f)"

formula <- "dilation_z ~ response_f * valence_f * trialnr_block_z + (response_f * valence_f * trialnr_block_z|subject_f)"
formula <- "dilation_z ~ response_f * valence_f * cueRep_z + (response_f * valence_f * cueRep_z|subject_f)"

formula <- "dilation_z ~ RTcleaned_z * ACC_f * valence_f + (RTcleaned_z * ACC_f * valence_f|subject_f)"

## S04:
formula <- "dilation_z ~ ACC_f + (arousal_f|subject_f)"
formula <- "dilation_z ~ RTcleaned_fast_f + (RTcleaned_fast_f|subject_f)"
formula <- "dilation_z ~ repeat_f + (repeat_f|subject_f)"

formula <- "dilation_z ~ ACC_f * response_f + (ACC_f * response_f|subject_f)"
formula <- "dilation_z ~ repeat_f * response_f + (repeat_f * response_f|subject_f)"

formula <- "dilation_z ~ ACC_f * valence_f + (ACC_f * valence_f|subject_f)"
formula <- "dilation_z ~ RTcleaned_fast_f * valence_f + (RTcleaned_fast_f * valence_f|subject_f)"
formula <- "dilation_z ~ repeat_f * valence_f + (repeat_f * valence_f|subject_f)"

formula <- "dilation_z ~ ACC_f * RTcleaned_fast_f + (ACC_f * RTcleaned_fast_f|subject_f)"
formula <- "dilation_z ~ ACC_f * valence_f + (ACC_f * valence_f|subject_f)"
formula <- "dilation_z ~ RTcleaned_fast_f * valence_f + (RTcleaned_fast_f * valence_f|subject_f)"

## S07:
formula <- "dilation_z ~ outcome_rel_f + (outcome_rel_f|subject_f)"
formula <- "dilation_z ~ outcome_f + (outcome_f|subject_f)"
formula <- "dilation_z ~ outcome_all_f + (outcome_all_f|subject_f)"

formula <- "dilation_z ~ outcome_rel_f * response_f + (outcome_rel_f * response_f|subject_f)"
formula <- "dilation_z ~ outcome_f * response_f + (outcome_f * response_f|subject_f)"
formula <- "dilation_z ~ outcome_all_f * response_f + (outcome_all_f * response_f|subject_f)"

formula <- "dilation_z ~ outcome_f * ACC_f + (outcome_f * ACC_f|subject_f)"
formula <- "dilation_z ~ outcome_f * RTcleaned_fast_f + (outcome_f * RTcleaned_fast_f|subject_f)"
formula <- "dilation_z ~ outcome_f * repeat_f + (outcome_f * repeat_f|subject_f)"

## Fit model:
mod <- fit_lmem(formula)
# Anova(mod, type = 3)
quickCI(mod, nRound = 3)
modLRT <- fit_lmem(formula, useLRT = T)

### Fit linear regression:
selData <- modData # select all data
mod <- lmer(formula = formula, data = selData,
            control = lmerControl(optCtrl = list(maxfun = 1e+9), calc.derivs = F, optimizer = c("bobyqa")))
summary(mod)
Anova(mod, type = 3)
quickCI(mod, nRound = 3)
mod_LRT <- mixed(formula = formula, data = selData, method = "LRT", type = "III", # all_fit = T,
                 control = lmerControl(optCtrl = list(maxfun = 1e+9), calc.derivs = F, optimizer = c("bobyqa"))) # nloptwrap
anova(mod_LRT)

### Fit linear regression on Go data only:
formula <- "dilation_z ~ ACC_f * valence_f + (ACC_f * valence_f|subject_f)"
formula <- "dilation_z ~ RTcleaned_fast_f * valence_f + (RTcleaned_fast_f * valence_f|subject_f)"
formula <- "dilation_z ~ repeat_f * valence_f + (repeat_f * valence_f|subject_f)"

selData <- subset(modData, response_f == "Go")
mod <- lmer(formula = formula, data = selData,
            control = lmerControl(optCtrl = list(maxfun = 1e+9), calc.derivs = F, optimizer = c("bobyqa")))
summary(mod)
Anova(mod, type = 3)
quickCI(mod, nRound = 3)
mod_LRT <- mixed(formula = formula, data = selData, method = "LRT", type = "III", # all_fit = T,
                 control = lmerControl(optCtrl = list(maxfun = 1e+9), calc.derivs = F, optimizer = c("bobyqa"))) # nloptwrap
anova(mod_LRT)

## Plot:
plot(effect("response_f", mod))
plot(effect("valence_f", mod))
plot(effect("ACC_f", mod))
plot(effect("repeat_f", mod))
plot(effect("RTcleaned_z", mod))

plot(effect("ACC_f:response_f", mod))
plot(effect("repeat_f:response_f", mod))

plot(effect("reqAction_f:valence_f", mod, x.var = "valence_f"))
plot(effect("response_f:valence_f", mod))
plot(effect("response_f:valence_f", mod, x.var = "valence_f"))
plot(effect("response_f:trialnr_block_z", mod))
plot(effect("response_f:trialnr_block_z", mod, x.var = "response_f"))
plot(effect("response_f:cueRep_z", mod))
plot(effect("response_f:cueRep_z", mod, x.var = "response_f"))
plot(effect("response_f:trialnr_block_z", mod))

plot(effect("RTcleaned_z:ACC_f", mod))
plot(effect("RTcleaned_z:valence_f", mod))
plot(effect("RTcleaned_z:valence_f", mod, x.var = "valence_f"))
plot(effect("ACC_f:valence_f", mod))
plot(effect("ACC_f:valence_f", mod, x.var = "valence_f"))

# ============================================================================ #
# ============================================================================ #
# ============================================================================ #
#### 03) Inter-correlation effects of arousal manipulation on behavior and dilation: ####

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

## Correlations:
rcor.test(corData)
#           response_n RTcleaned_z dilation_z
# response_n   *****     -0.423      -0.202    
# RTcleaned_z  0.011      *****       0.121    
# dilation_z   0.243      0.487       *****    

## Select variables and labels:
xVar <- "response_n"; yVar <- "dilation_z"; xLab <- "responses"; yLab <- "dilations"
xVar <- "RTcleaned_z"; yVar <- "dilation_z"; xLab <- "RTs"; yLab <- "dilations"

## Plot:
p <- plot_correlation(data = corData, xVar = xVar, yVar = yVar, 
                      xLab = paste0("Effect of arousal manipulation on ", xLab), 
                      yLab = paste0("Effect of arousal manipulation on ", yLab))
print(p)

## Save:
plotName <- paste0("correlation_effect_manip_on_", xLab, "~", yLab)
cat(paste0("Save as ", plotName, "\n"))
png(paste0(dirs$plotDir, plotName, ".png"), width = 480, height = 480)
print(p)
dev.off()

# ============================================================================ #
#### 04) Follow-up with EMMEANS: ####
# https://aosmith.rbind.io/2019/03/25/getting-started-with-emmeans/

emmeans(mod, ~ arousal_f)
emmeans(mod, specs = pairwise ~ arousal_f, 
        regrid = "response")

## Arousal:
emmeans(mod, specs = pairwise ~ arousal_f | condition_f, 
        regrid = "response", interaction = "pairwise")
emmeans(mod, specs = pairwise ~ arousal_f | condition_f, 
        interaction = "pairwise")

emmeans(mod, specs = pairwise ~ condition_f | arousal_f, 
        regrid = "response", interaction = "pairwise")


## Dilation ~ response x valence:
emmeans(mod, specs = pairwise ~ valence_f | response_f, 
        interaction = "pairwise")
emmeans(mod, specs = pairwise ~ valence_f | response_f, 
        regrid = "response", interaction = "pairwise")

## Dilation as continuous IV:
# https://bookdown.org/mike/data_analysis/emmeans-package.html
emtrends(mod, ~ condition_f, var = "dilation_z")
emtrends(mod, pairwise ~ condition_f, var = "dilation_z")

test(emmeans(mod, pairwise ~condition_f * dilation_z))
test(emmeans(mod, pairwise ~condition_f * dilation_z, at = list(dilation_z = 0)))

emmeans(mod, specs = pairwise ~ dilation_z | condition_f, 
        regrid = "response", interaction = "pairwise")
emmeans(mod, specs = pairwise ~ dilation_z | condition_f, 
        interaction = "pairwise")

## Dilation as a function of outcome:
emmeans(mod, specs = pairwise ~ outcome_rel_f, 
        regrid = "response")
emmeans(mod, specs = pairwise ~ outcome_f, 
        regrid = "response")
emmeans(mod, specs = pairwise ~ outcome_all_f, 
        regrid = "response")

# ============================================================================ #
#### 05) Plot effects with plot(effect()): ####

## Main effect:
plot(effect("reqAction_f", mod)) # 
plot(effect("valence_f", mod)) # 
plot(effect("arousal_f", mod)) # 
plot(effect("reqCongruency_f", mod)) # 

## 2-way interactions:
plot(effect("reqAction_f:valence_f", mod)) # 
plot(effect("reqAction_f:valence_f", mod, x.var = "valence_f")) # 

plot(effect("valence_f:arousal_f", mod)) # 
plot(effect("valence_f:arousal_f", mod, x.var = "valence_f")) # 
plot(effect("valence_f:arousal_f", mod, x.var = "arousal_f")) # 

plot(effect("valence_f:dilation_z", mod)) # 
plot(effect("valence_f:dilation_z", mod, x.var = "valence_f")) # 
plot(effect("valence_f:dilation_z", mod, x.var = "dilation_z")) # 

plot(effect("condition_f:arousal_f", mod)) # 
plot(effect("condition_f:arousal_f", mod, x.var = "arousal_f")) # 

plot(effect("arousal_f:trialnr_z", mod)) # 
plot(effect("arousal_f:trialnr_z", mod, x.var = "arousal_f")) # 
plot(effect("arousal_f:cueRep_z", mod)) # 
plot(effect("arousal_f:cueRep_z", mod, x.var = "arousal_f")) # 
plot(effect("reqCongruency_f:arousal_f:trialnr_z", mod)) # 
plot(effect("reqCongruency_f:arousal_f:cueRep_z", mod)) # 

plot(effect("outcome_last_f:reqCongruency_f", mod)) # 
plot(effect("outcome_last_f:reqCongruency_f", mod, x.var = "reqCongruency_f")) # 

## 3-way interactions:
plot(effect("reqAction_f:valence_f:arousal_f", mod)) # 
plot(effect("reqAction_f:valence_f:arousal_f", mod, x.var = "arousal_f")) # 

plot(effect("valence_f:reqResp_f:outcome_lag1_f", mod, x.var = "outcome_lag1_f")) # 

# ============================================================================ #
#### 06) Plot effect with CUSTOM REGRESSION PLOTS: ####

plot(effect("reqAction_f", mod)) # 
plot(effect("valence_f", mod)) # 
plot(effect("reqCongruency_f", mod)) # 

# ---------------------------------------------- #
## Bar plot: 

selEff <- "reqAction_f"
selEff <- "valence_f"

yLim <- c(0, 1) # responses
yLim <- c(0.4, 1) # RTs

plot(effect(selEff, mod)) 
p <- custom_regressionbar1(mod, selEff = selEff, yLim = yLim)
png(paste0(dirs$plotDir, "regressionbar1_", all.vars(terms(mod))[1], "~", selEff, ".png"), width = 480, height = 480)
print(p)
dev.off()

# ---------------------------------------------- #
## Line plot main effect:

formula <- "response_n ~ dilation_z + (dilation_z|subject_f)"
mod <- fit_lmem(formula)

selEff <- "dilation_z"

yLim <- c(0, 1) # responses
yLim <- c(0.4, 1) # RTs
yLim <- c(-1, 2) # RT_z

plot(effect(selEff, mod)) 
p <- custom_regressionline1(mod, selEff = selEff, yLim = yLim, 
                            # xLim = c(-3, 6), 
                            useEffect = F)
png(paste0(dirs$plotDir, "regressionline1_", all.vars(terms(mod))[1], "~", selEff, ".png"), width = 480, height = 480)
print(p)
dev.off()

# ---------------------------------------------- #
## Line plot 2-way interaction between continuous and categorical regressor:

# reqAction_f * valence_f * dilation_z 

formula <- "response_n ~ reqAction_f * dilation_z  + (reqAction_f * dilation_z|subject_f)"
formula <- "RTcleaned_z ~ reqAction_f * dilation_z  + (reqAction_f * dilation_z|subject_f)"
formula <- "RTcleaned_z ~ valence_f * dilation_z  + (valence_f * dilation_z|subject_f)"

mod <- fit_lmem(formula)

xVar <- "dilation_z"; zVar <- "reqAction_f"; 
xVar <- "dilation_z"; zVar <- "valence_f"

yLim <- c(0, 1) # responses
yLim <- c(0.4, 1) # RTs
yLim <- c(-1, 2) # RT_z

plot(effect(paste0(zVar, ":", xVar), mod), multiline = T, ci.style = "bands")
plot(effect(paste0(xVar, ":", zVar), mod), multiline = T, ci.style = "bands") 
p <- custom_regressionline2(mod = mod, xVar = xVar, zVar = zVar, yLim = yLim)
                            # xVec = seq(-3, 6, 0.5))
png(paste0(dirs$plotDir, "regressionline1_", all.vars(terms(mod))[1], "~", xVar, "_x_", zVar, ".png"), width = 480, height = 480)
print(p)
dev.off()

# ============================================================================ #
#### 07) Plot coefficients with CUSTOM COEFPLOTS: ####

formula <- "response_n ~ reqAction_f * valence_f + (reqAction_f * valence_f|subject_f)"
formula <- "response_n ~ reqAction_f * valence_f * arousal_f + (reqAction_f * valence_f * arousal_f|subject_f)"
formula <- "response_n ~ reqAction_f * valence_f * dilation_z + (reqAction_f * valence_f * dilation_z|subject_f)"
formula <- "response_n ~ reqAction_f * valence_f * arousal_f * dilation_z + (reqAction_f * valence_f * arousal_f * dilation_z|subject_f)"

formula <- "RTcleaned_z ~ reqAction_f * valence_f + (reqAction_f * valence_f|subject_f)"
formula <- "RTcleaned_z ~ reqAction_f * valence_f * arousal_f + (reqAction_f * valence_f * arousal_f|subject_f)"
formula <- "RTcleaned_z ~ reqAction_f * valence_f * dilation_z + (reqAction_f * valence_f * dilation_z|subject_f)"
formula <- "RTcleaned_z ~ reqAction_f * valence_f * arousal_f * dilation_z + (reqAction_f * valence_f * arousal_f * dilation_z|subject_f)"

formula <- "dilation_z ~ response_f * valence_f + (response_f * valence_f|subject_f)"
formula <- "dilation_z ~ reqAction_f * valence_f + (reqAction_f * valence_f|subject_f)"

# --------------------------------------------- #
## Retrieve model:

mod <- fit_lmem(formula)

# --------------------------------------------- #
## Compute number of coefficients:

nCoef <- length(fixef(mod)) - 1; nCoef
# nCoef <- 3 # set manually

## Select colours:
library(MetBrewer)
selCol <- met.brewer("Demuth", n = nCoef, type = "continuous"); colName <- "Demuth"
selCol <- met.brewer("Manet", n = nCoef, type = "continuous"); colName <- "Manet" # brown, brown, green
selCol <- met.brewer("Monet", n = nCoef, type = "continuous"); colName <- "Monet" # dark green, green, light green

save_coefplots <- function(){
        
  ## Without subs:
  plotSub <- F
  p <- custom_coefplot(mod, plotSub = plotSub, plotText = T, dropIntercept = T, 
                       # xLim = c(-2, 1),
                       revOrder = T, selCol = selCol)
  plotHandle <- formula2handle(formula)
  plotName <- paste0("coef_glmer_", plotHandle, "_", colName, "_nCol", nCoef)
  png(paste0(dirs$plotDir, plotName, ".png"), width = 720, height = 480)
  print(p)
  dev.off()
  
  ## With subs:
  plotSub <- T
  p <- custom_coefplot(mod, plotSub = plotSub, plotText = T, dropIntercept = T, 
                       # xLim = c(-2, 1),
                       revOrder = T, selCol = selCol)
  plotHandle <- formula2handle(formula)
  plotName <- paste0("coef_glmer_", plotHandle, "_", colName, "_nCol", nCoef)
  if (plotSub){plotName <- paste0(plotName, "_subs")}
  png(paste0(dirs$plotDir, plotName, ".png"), width = 720, height = 480)
  print(p)
  dev.off()
  
}

save_coefplots() # defined above

# ---------------------------------------------------------------------------- #
#### Plot regressor inter-correlations: ####

corrplot_regressors(mod, perSub = F, savePNG = T) # intercorrelations overall
corrplot_regressors(mod, perSub = T, savePNG = T) # intercorrelations per subject averaged

# ---------------------------------------------------------------------------- #
#### Plot coefficient inter-correlations: ####

corrplot_coefficients(mod, savePNG = T) # plot inter-correlations between coefficients

# ---------------------------------------------------------------------------- #
#### Save coefficients per subject: ####

save_coefficients(mod) # save coefficients per subject
cat("Saved everything :-)\n")

# END OF FILE.
