#!/usr/bin/env Rscript
# ============================================================================ #
## 00_mgngunc_functions_timecourse.R
## MGNGUncertainty functions for analyzing gaze and pupil data on a ms-by-ms basis.
# Johannes Algermissen, 2023.

# =============================================================================================== #
#### 01) Define aggregation function: ####

aggregate_timecourse_condition <- function(data, zVar, selTrialVar = NULL, selTrialVal = NULL, 
                                           eyeMeasure = "pupil", aggrMetric = NULL,
                                           segMessage = "StartMask", beforeSamples = 1000, afterSamples = 3000){
  
  #' Aggregate pupil time courses per subject per condition.
  #' @param data data frame, trial-by-trial behavioural data.
  #' @param zVar scalar string, name of variable to split by. Variable needs to be a factor.
  #' @param selTrialVar string, variable based on which to select subset of trials (optional).
  #' @param selTrialVal integer or string, value to select for on selTrialVar (optional).
  #' @param eyeMeasure string, eye metric to load in ("pupil" or "gaze", default: "pupil").
  #' @param aggrMetric string, aggregation metric to use (default: same as eyeMeasure).
  #' @param segMessage scalar string, event relative to which time courses are epoched and to-be-loaded (default: StartMask).
  #' @param beforeSamples scalar integer, number of samples before epoching event (default: 1000).
  #' @param afterSamples scalar integer, number of samples after epoching event (default: 3000).
  #' @return return 3D-array with dimensions nSub, nSample, nCond. 
  
  require(stringr)
  
  ## Set default settings:
  if (eyeMeasure == "pupil"){aggrMetric <- eyeMeasure}
  
  cat(paste0("Use ", eyeMeasure, " data aggregated as ", aggrMetric, "\n"))
  cat(paste0("Aggregate data per condition per subject split per ", zVar, "\n"))
  cat(paste0("Load trial-by-trial pupil time courses with data epoched relative to ", 
             segMessage, " with ", beforeSamples, " ms before and ", afterSamples, " ms afterwards\n"))
  if (!(is.null(selTrialVar))){
    if(is.null(selTrialVal)){stop("selTrialVar provided, but not selTrialVal")}
    cat(paste0("Select only data for ", selTrialVar, " == ", selTrialVal, "\n"))
  }
  
  # --------------------------------------------------------------------------------------------- #
  ### Determine input directory:
  
  beforeSamplesStr <- str_pad(beforeSamples, width = 4, side = "left", pad = "0")
  afterSamplesStr <- str_pad(afterSamples, width = 4, side = "left", pad = "0")
  inputDir <- paste0(dirs$processedDataDir, eyeMeasure, "/", aggrMetric, "_timecourse_", 
                     segMessage, "_", beforeSamplesStr, "_", afterSamplesStr, "ms/")
  cat(paste0("Load time course data per trial from ", inputDir, "\n"))
  
  # --------------------------------------------------------------------------------------------- #
  ### Determine output directory and file name:
  
  ## Output directory:
  outputDir <- paste0(dirs$processedDataDir, eyeMeasure, "/", eyeMeasure, "_timecourse_condition/")
  dir.create(outputDir, showWarnings = FALSE, recursive = TRUE)
  
  ## Output file:
  outFileName <- paste0(aggrMetric, "_timecourse_per_", zVar, "_", segMessage, "_", beforeSamplesStr, "_", afterSamplesStr, "ms") 
  if(!(is.null(selTrialVar))){
    outFileName <- paste0(outFileName, "_selTrials_", selTrialVar, "_", selTrialVal)
  }
  fullFileName <- paste0(outputDir, outFileName, ".rds")
  
  # --------------------------------------------------------------------------------------------- #
  ### Check if file exists, if yes, then load:
  if(file.exists(fullFileName)){
    
    cat(paste0(fullFileName, " already exists, load...\n"))
    condArray <- readRDS(fullFileName)
    cat("Finished loading! :-)\n")
    
  } else {
    
    cat(paste0(fullFileName, " does not exist yet, create anew...\n"))
    
    # ------------------------------------------------------------------------------------------- #
    ### Determine number of conditions:
    
    condVec <- levels(data[, zVar])
    nCond <- length(condVec)
    cat(paste0("Split by variable \'", zVar, "\' with ", nCond, " conditions (", paste0(condVec, collapse = ", "), ")\n"))
    
    # ------------------------------------------------------------------------------------------- #
    ### Initialize output data:
    
    nSub <- length(unique(data$subject_n))
    nSample <- beforeSamples + afterSamples + 1 
    condArray <- array(NA, dim = c(nSub, nSample, nCond)) # via array
    
    # ------------------------------------------------------------------------------------------- #
    ### Loop over subjects:
    
    for (iSub in 1:nSub){ # iSub <- 1
      
      subStr <- str_pad(iSub, width = 2, side = "left", pad = "0")
      cat(paste0("Start subject ", subStr, ", aggregate per level of variable \'", zVar, "\'\n"))
      
      ## Extract subject data:
      subData <- subset(data, subject_n == iSub)
      
      ## Load time course data for this subject:
      subStr <- str_pad(iSub, width = 2, side = "left", pad = "0")
      timeData <- read.csv(paste0(inputDir, "MGNGUNC_Sub_", subStr, "_", eyeMeasure, ".csv"))
      
      ## Loop over levels:
      for (iCond in 1:nCond){ # iCond <- 1
        
        ## Trials for this condition:
        trlIdx <- which(subData[, zVar] == condVec[iCond])
        
        ## If restricted to certain trials in selected condition:
        if(!(is.null(selTrialVar))){ # only keep trials in selected condition
          selIdx <- which(subData[, selTrialVar] == selTrialVal)
          trlIdx <- trlIdx[trlIdx %in% selIdx]
        }
        condMean <- as.numeric(colMeans(timeData[trlIdx, ], na.rm = T)) # mean per condition
        # plot(condMean)
        condArray[iSub, , iCond] <- condMean # save
        
      } # end iCond
      
    } # end iSub
    
    ## Save:
    cat(paste0("Save data under", fullFileName, "\n"))
    require(plyr)
    saveRDS(condArray, fullFileName)
    
    cat("Finished!\n")
  }
  
  return(condArray)
}

# ============================================================================ #
#### 02) Plot time course per condition: ####

### Define plotting function:
plot_timecourse <- function(data, zVar, selTrialVar = NULL, selTrialVal = NULL,
                            eyeMeasure = "pupil", aggrMetric = NULL,
                            segMessage = "StartMask", beforeSamples = 1000, afterSamples = 3000, 
                            isBaselineCorGroup = FALSE, isBaselineCorSub = FALSE, baselineTime = 0,
                            zLab = NULL, main = NULL, FTS = NULL, LWD = NULL, axesLWD = FALSE,
                            selCol = NULL, selLineType = c(1), breakVec = NULL,
                            SEweight = 1, xLim = NULL, yLim = NULL, addLegend = F, isPNG = TRUE, suffix = NULL){
  #' Plot pupil or gaze time course as line plot per condition with group-level lines plus shades in ggplot.
  #' @param data data frame, trial-by-trial data.
  #' @param zVar scalar string, name of variable to split by (differently colored lines). Variable needs to be a factor.
  #' @param selTrialVar string, variable based on which to select subset of trials (optional).
  #' @param selTrialVal integer or string, value to select for on selTrialVar (optional).
  #' @param eyeMeasure string, eye metric to load in ("pupil" or "gaze", default: "pupil").
  #' @param aggrMetric string, aggregation metric to use (default: same as eyeMeasure).
  #' @param segMessage scalar string, event relative to which time courses are epoched and to-be-loaded (default: StartMask).
  #' @param beforeSamples scalar integer, number of samples before epoching event (default: 1000).
  #' @param afterSamples scalar integer, number of samples after epoching event (default: 3000).
  #' @param isBaselineCorGroup  Boolean, set average condition difference across participants to 0 at time point 0 (default: FALSE).
  #' @param isBaselineCorSub    Boolean, set condition difference for each participant to 0 at time point 0 (default: FALSE).
  #' @param baselineTime scalar integer (single time point), vector of two integers (start and stop), or nCond x 2 matrix with baseline window per condition in ms (default: 0).
  #' @param zLab scalar string, label to use as header for coloring legend (optional).
  #' @param main string, overall plot label (optional).
  #' @param FTS scalar numeric, font size, optional (default: retrieved via retrieve_plot_defaults()).
  #' @param LWD scalar numeric, line width, optional (default: retrieved via retrieve_plot_defaults()).
  #' @param axesLWD Boolean, apply LWD also to axis (default: FALSE).
  #' @param selCol vector of strings (HEX colors), colors for bars (default: retrieve via retrieve_colour()).
  #' @param selLineType vector of numerics, line types to use (default: c(1)).
  #' @param breakVec vector of numerics, x-axis tick marks (default: seq(-2, 3, 1)).
  #' @param SEweight scalar, weight to use for error shades (how many times SE; default: 1).
  #' @param xLim vector of two numerics, x-axis limits (optional; default: NULL).
  #' @param yLim vector of two numerics, y-axis limits (optional: default: NULL).
  #' @param addLegend Boolean, add legend at the top or not (default: FALSE).
  #' @param isPNG Boolean, save as .png file (default: FALSE).
  #' @param suffix string, appended to figure name (optional; default: NULL).
  #' @return creates (and saves) plot.
  
  # eyeMeasure = "pupil"; aggrMetric = "pupil"
  # selTrialVar = NULL; selTrialVal = NULL; isBaselineCorGroup = F; isBaselineCorSub = F; zLab = NULL; main = NULL; selCol = NULL; selLineType = c(1); breakVec = NULL; SEweight = 1; yLim = NULL; addLegend = F; isPNG = F
  require(ggplot2)
  
  ## Set default settings:
  if (eyeMeasure == "pupil"){aggrMetric <- eyeMeasure}
  
  # -------------------------------------------------------------------------- #
  ## Retrieve input dimensions:
  
  condArray <- aggregate_timecourse_condition(data = data, zVar = zVar, selTrialVar = selTrialVar, selTrialVal = selTrialVal,
                                              eyeMeasure = eyeMeasure, aggrMetric = aggrMetric,
                                              segMessage = segMessage, beforeSamples = beforeSamples, afterSamples = afterSamples)
  
  # -------------------------------------------------------------------------- #
  ## Retrieve input dimensions:
  
  nSub <- dim(condArray)[1]
  nSample <- dim(condArray)[2]
  nCond <- dim(condArray)[3]
  
  ## Time variable:
  timeVecMs <- seq(from = -1*beforeSamples, to = afterSamples, by = 1) # add 1 because of sample at 0
  
  if (is.null(xLim)){
    xLim <- c(min(timeVecMs), max(timeVecMs))
  }
  cat(paste0("Plot from ", xLim[1], " till ", xLim[2], " ms\n"))
  
  ## Select samples within specified x-limits:
  timeVecXlim <- seq(from = xLim[1], to = xLim[2], by = 1)
  timeIdx <- which(timeVecMs %in% timeVecXlim)
  
  ## Select valid time samples, convert to seconds:
  timeVecSec <- timeVecMs[timeIdx]/1000 # to sec.
  
  ## Condition names:
  condNames <- levels(data[, zVar])
  if(nCond != length(condNames)){stop("3rd dimension condArray and number levels zVar do not match")}
  
  # -------------------------------------------------------------------------- #
  ## Fixed plotting settings:
  
  if (is.null(FTS)){
    FTS <- 28 # retrieve_plot_defaults("FTS")
  }
  if (is.null(LWD)){
    LWD <- retrieve_plot_defaults("LWD") # 1.3
  }
  colAlpha <- 1
  
  # -------------------------------------------------------------------------- #
  ## Fill other plotting settings:
  
  ## Axis labels:
  xLab <- "Time (in sec.)"
  if (aggrMetric == "pupil"){yLab <- "Pupil diameter (a.u.)"}
  if (aggrMetric == "xp"){yLab <- "x-coordinate (pixels)"}
  if (aggrMetric == "yp"){yLab <- "y-coordinate (pixels)"}
  if (aggrMetric == "distanceGroupMean"){yLab <- "Distance from center (pixels)"}
  if (aggrMetric == "distanceBaselineTrial"){yLab <- "Distance from center (pixels)"}
  if (aggrMetric == "distanceAbsDeriv"){yLab <- "Absolute derivative (pixels)"}
  if (aggrMetric == "distanceAbsCumSum"){yLab <- "Cum sum of abs deriv (pixels)"}
  if(is.null(zLab)){zLab <- substitute_label(zVar)}
  
  ## Colours:
  if(is.null(selCol)){selCol <- retrieve_colour(zVar)}
  
  ## Line type:
  if (length(selLineType) == 1){selLineType <- rep(selLineType, nCond)}
  
  # -------------------------------------------------------------------------- #
  ## Initialize empty objects:
  
  meanMat <- matrix(NA, nrow = nCond, ncol = nSample)
  sdMat <- matrix(NA, nrow = nCond, ncol = nSample)
  
  ## Compute mean and sd across subjects (= aggregate) per condition:
  for (iCond in 1:nCond){ # iCond <- 1
    # https://stackoverflow.com/questions/13058978/mean-on-the-third-dimension-in-r
    
    ## Extract condition:
    selCondMat <- condArray[, , iCond]
    
    # baselineTime <- matrix(c(-2000, -1500, -1950, -1450, -2000, -1500, -1950, -1450), nrow = 4, ncol = 2, byrow = T)
    
    ## Apply baseline correction at subject-level:
    if(isBaselineCorSub){
      cat(paste0("Apply baseline correction at subject level for condition ", iCond, "\n"))
      if (length(baselineTime) == 1){baselineIdx <- which(timeVecMs == baselineTime)} # single number
      if (length(baselineTime) == 2){baselineIdx <- which(timeVecMs %in% seq(baselineTime[1], baselineTime[2], 1))} # vector of start/ end
      if (is.matrix(baselineTime)){baselineIdx <- which(timeVecMs %in% seq(baselineTime[iCond, 1], baselineTime[iCond, 2], 1))} # vector of start/end per condition
      if (length(baselineIdx) > 1){
        subBase <- rowMeans(selCondMat[, baselineIdx], na.rm = T) # average over selected time bins
      } else {
        subBase <- selCondMat[, baselineIdx] # average over selected time bins
      }
      baseMat <- replicate(nSample, subBase)
      selCondMat <- selCondMat - baseMat
      # stopifnot(all(selCondMat[, which(timeVecMs %in% baselineTime)] == 0))
    }
    
    ## Compute mean per condition across subjects:
    selCondVec <- as.numeric(colMeans(selCondMat, na.rm = T)) # average over subjects
    
    ## Perform baseline correction (subtract mean within baseline window for each condition separately):
    if(isBaselineCorGroup){
      cat("Apply baseline correction at group level\n")
      if (length(baselineTime) == 1){baselineIdx <- which(timeVecMs == baselineTime)} # single number
      if (length(baselineTime) == 2){baselineIdx <- which(timeVecMs %in% seq(baselineTime[1], baselineTime[2], 1))} # vector of start/ end
      if (is.matrix(baselineTime)){baselineIdx <- which(timeVecMs %in% seq(baselineTime[iCond, 1], baselineTime[iCond, 2], 1))} # vector of start/end per condition
      selCondVec <- selCondVec - mean(selCondVec[baselineIdx], na.rm = T)
    }
    meanMat[iCond, ] <- selCondVec # save
    
    ## Compute SE using Cousineau-Morey correction:
    subMean <- as.numeric(rowMeans(selCondMat,  na.rm = T)) # mean over time per subject per condition
    grandMean <- mean(subMean, na.rm = T) # grand mean per condition
    meanMatCor <- selCondMat - replicate(nSample, subMean) + grandMean # subtract subject mean, add grand mean
    sdMat[iCond, ] <- apply(meanMatCor, 2, sd, na.rm = T) # average over first dimension
  }
  
  ## Compute SE across subjects per condition:
  # SD divided by sqrt(nSub), correct for number conditions with Cousineau-Morey correction
  cmFactor <- sqrt(nCond / (nCond - 1))
  if(nCond == 1){cmFactor <- 1}
  seMat <- sdMat / sqrt(nSub) * cmFactor
  
  # -------------------------------------------------------------------------- #
  ### Prepare data sets:
  
  ## Set pro-forma data set to be able to initialize ggplot:
  summary_d <- data.frame(x = timeVecSec, y = meanMat[iCond, timeIdx])
  summary_d$z <- factor(1)
  
  ## Keep track of minimum and maximum for y-axis limits:
  yMinVec <- rep(NA, nCond)
  yMaxVec <- rep(NA, nCond)
  
  # -------------------------------------------------------------------------- #
  ### Plot with ggplot:
  p <- ggplot(data = summary_d, aes(x = x, y = y, col = z, fill = z))
  
  ## Loop over conditions, make shade and line on top:
  for (iCond in 1:nCond){ # iCond <- 1
    
    ## Select data, set upper and lower shade limits:
    condData <- data.frame(x = timeVecSec, y = meanMat[iCond, timeIdx], se = seMat[iCond, timeIdx]) # select data for this condition
    condData$ymin <- condData$y - SEweight * condData$se # lower edge of shade
    condData$ymax <- condData$y + SEweight * condData$se # upper edge of shade
    condData$z <- condNames[iCond]
    
    ## Y-axis limits:
    yMinVec[iCond] <- min(condData$ymin, na.rm = T)
    yMaxVec[iCond] <- max(condData$ymax, na.rm = T)
    
    ## Shade:
    p <- p + geom_ribbon(data = condData, aes(x = x, y = y, ymin = ymin, ymax = ymax, group = 1), 
                         # col = NA, # remove outer border of shades
                         linetype = 0, # remove outer border of shades
                         # fill = selCol[iCond], # comment out for legend in colours (not grey)
                         # show.legend = T, 
                         alpha = 0.2, lwd = 0)
    
    ## Line:
    p <- p + geom_path(data = condData, aes(x = x, y = y, group = 1),
                       # col = selCol[iCond], 
                       linetype = selLineType[iCond], size = LWD) 
  }
  
  # -------------------------------------------------------------------------- #
  ### Horizontal line at h = 0:
  
  
  if (segMessage == "StartMask"){
    p <- p + geom_vline(xintercept = 0, linetype = "dashed", colour = "#949494", lwd = 0.5) # forward mask onset
    p <- p + geom_vline(xintercept = 0.250, linetype = "dashed", colour = "#949494", lwd = 0.5) # prime onset
    p <- p + geom_vline(xintercept = 0.266, linetype = "dashed", colour = "#949494", lwd = 0.5) # backward mask onset
    p <- p + geom_vline(xintercept = 0.366, linetype = "dashed", colour = "#949494", lwd = 0.5) # cue onset
    p <- p + geom_vline(xintercept = 1.666, linetype = "dashed", colour = "#949494", lwd = 0.5) # cue offset
  }
  if (segMessage == "StartCue"){
    p <- p + geom_vline(xintercept = -0.366, linetype = "dashed", colour = "#949494", lwd = 0.5) # forward mask onset
    p <- p + geom_vline(xintercept = -0.255, linetype = "dashed", colour = "#949494", lwd = 0.5) # prime onset
    p <- p + geom_vline(xintercept = -0.250, linetype = "dashed", colour = "#949494", lwd = 0.5) # backward mask onset
    p <- p + geom_vline(xintercept = 0, linetype = "dashed", colour = "#949494", lwd = 0.5) # cue offset
    p <- p + geom_vline(xintercept = 1.300, linetype = "dashed", colour = "#949494", lwd = 0.5) # cue offset
  }
  if (segMessage == "StartResponse"){
    p <- p + geom_vline(xintercept = 0, linetype = "dashed", colour = "#949494", lwd = 0.5) # forward mask onset
  }
  if (segMessage == "StartOutcome"){
    p <- p + geom_vline(xintercept = 0, linetype = "dashed", colour = "#949494", lwd = 0.5) # outcome onset
    p <- p + geom_vline(xintercept = 0.700, linetype = "dashed", colour = "#949494", lwd = 0.5) # outcome offset
  }
  
  # -------------------------------------------------------------------------- #
  ### Axes:
  
  ## X-axis:
  xMin <- min(timeVecSec)
  xMax <- max(timeVecSec)
  if (is.null(breakVec)){
    # breakVec <- seq(-3, 5, 1)
    breakVec <- find_break_points(c(xMin, xMax), nTickTarget = 5)
  }
  p <- p + scale_x_continuous(limits = c(xMin, xMax), breaks = breakVec)
  
  ## Y-axis:
  if(is.null(yLim)){yLim <- c(min(yMinVec), max(yMaxVec))}
  cat(paste0("yMin = ", yLim[1], "; yMax = ", yLim[2], "\n"))
  
  if (yLim[1] == 0 & yLim[2] == 1){
    p <- p + scale_y_continuous(breaks = seq(0, 1, by = 0.5)) # only 0, 0.5, 1 as axis labels
  }
  if(!(is.null(yLim))){
    # p <- p + scale_y_continuous(breaks = seq(yLim[1], yLim[2], by = 10))
    p <- p + coord_cartesian(ylim = yLim)}
  
  # -------------------------------------------------------------------------- #
  ### Labels:
  
  ## Add labels:
  p <- p + labs(x = xLab, y = yLab, fill = condNames, col = condNames)
  
  ## Add title:
  if (!(is.null(main))){
    p <- p + ggtitle(main)  
  }
  
  
  # -------------------------------------------------------------------------- #
  ## Add line colors:
  
  p <- p + scale_fill_manual(values = selCol, limits = condNames) # set limits for correct order
  p <- p + scale_color_manual(values = selCol, limits = condNames, guide = "none") # set limits for correct order
  
  # -------------------------------------------------------------------------- #
  ### Add theme, font sizes:
  
  require(ggthemes)
  p <- p + theme_classic() + 
    theme(axis.text = element_text(size = FTS),
          axis.title = element_text(size = FTS), 
          plot.title = element_text(size = FTS, hjust = 0.5))
  
  ## Change line width of axes:
  if (axesLWD){
    p <- p + theme(
      axis.line = element_line(size = LWD)
    )
  }
  
  ## Add legend:
  if (addLegend){
    p <- p + theme(
      legend.position = c(0.13, 0.9), # "topleft",
      legend.box = "vertical",
      legend.title = element_blank(),
      # legend.title = element_text(size = FTS),
      legend.text = element_text(size = FTS)
    )
  } else {
    p <- p + theme(
      legend.title = element_blank(), legend.position = "none"
    )
  }
  
  # -------------------------------------------------------------------------- #
  ### Save:
  
  if (isPNG){
    ## Save:  
    figName <- paste0(eyeMeasure, "_", aggrMetric, "_timecourse_", zVar)
    if(!(is.null(selTrialVar))){figName <- paste0(figName, "_selTrials_", selTrialVar, "_", selTrialVal)}
    if(isBaselineCorGroup){figName <- paste0(figName, "_baselineCorGroup")}
    if(isBaselineCorSub){figName <- paste0(figName, "isBaselineCorSub")}
    if(addLegend){figName <- paste0(figName, "_addLegend")}
    if(!is.null(suffix)){figName <- paste0(figName, suffix)}
    cat(paste0("Save figure as file ", figName, ".png\n"))
    png(paste0(dirs$plotDir, figName, ".png"), width = 600, height = 480) # 480 x 480
    print(p)
    dev.off()
  }
  
  # -------------------------------------------------------------------------- #
  ### Return:
  
  print(p)
  cat("Finished! :-)\n")
  return(p)
  
}

# ============================================================================ #
#### 03) Find plausible axis break points: ####

find_break_points <- function(input, nTickTarget = 4){
  
  # iMag <- log10(abs(xMin)) # determine order of magnitude for rounding
  # iMag <- ifelse(iMag < 0, floor(iMag), ceil(iMag)) # round
  # xUnit <- 10 ^ iMag
  # xBreakMin <- floor(xMin/xUnit)*xUnit # remove post-digit part, round, transform back
  # 
  # ## Determine x-axis upper limit:
  # iMag <- log10(abs(xMax)) # determine order of magnitude for rounding
  # iMag <- ifelse(iMag < 0, floor(iMag), ceil(iMag)) # round
  # xUnit <- 10 ^ iMag
  # xBreakMax <- ceiling(xMax/xUnit)*xUnit # remove post-digit part, round, transform back
  
  ## Extract:
  xBreakMin <- input[1]
  xBreakMax <- input[2]
  
  ## xStep either 1 or 5; try out both:
  xStep <- find_step(c(xBreakMin, xBreakMax), nTickTarget = nTickTarget)
  
  ## Correct if one limit smaller than xStep:
  if (abs(xBreakMin) < xStep){xBreakMin <- 0}
  if (abs(xBreakMax) < xStep){xBreakMax <- 0}
  
  cat(paste0("xBreakMin = ", xBreakMin, ", xBreakMax = ", xBreakMax, ", xStep = ", xStep, "\n"))
  
  ## Distance between x-axis ticks:
  breakVec <- seq(xBreakMin, xBreakMax, xStep) # just very broad break points, aligned to magnitude
  if(!(0 %in% breakVec)){ # correct if zero not included
    closestToZero <- breakVec[abs(breakVec) == min(abs(breakVec))]
    breakVec <- breakVec - closestToZero # subtract so becomes zero
  }
  cat(paste0("Use axis break points from ", min(breakVec), " till ", max(breakVec), " in steps of ", breakVec[2] - breakVec[1], "\n"))
  
  return(breakVec)
  
}

# ============================================================================ #
#### 04) Progress bar: ####

ProgressBar <- function (x, max = 100){
  percent <- x / max * 100
  cat(sprintf("\r\t[%-50s] %d%%",
              paste(rep("=", percent / 2), collapse = ""),
              floor(percent)))
  if (x == max)
    cat("\n")}


# ============================================================================ #
#### 05) Compute cluster mass with for loop: ####

computeClusterMass <- function(tVec, tThresh = 2){
  
  ## Exclude NAs:  
  tVecValid <- tVec[!(is.na(tVec))]
  
  ## Detect clusters above threshold:
  aboveThreshVec <- as.numeric(abs(tVecValid) > tThresh)
  # plot(aboveThreshVec)
  
  # ------------------------------------------------------------ #
  ## Label clusters of connected components: 
  count <- 1 # initialize count of clusters
  clustIdx <- rep(NA, length(aboveThreshVec)) # initialize
  clustIdx[1] <- ifelse(aboveThreshVec[1] == 1, 1, NA) # first cluster or not
  for (iElem in 2:length(aboveThreshVec)){ # loop over other elements
    if(aboveThreshVec[iElem] == 1 & aboveThreshVec[iElem - 1] == 1){ # same cluster
      clustIdx[iElem] <- count
    } else if (aboveThreshVec[iElem] == 1 & aboveThreshVec[iElem - 1] == 0){ # different cluster
      count <- count + 1 # increment
      clustIdx[iElem] <- count
    } else { # no cluster
      clustIdx[iElem] <- NA
    }
  }
  # plot(clustIdx)
  
  # ------------------------------------------------------------ #
  ## Compute cluster mass:
  clustLevel <- unique(clustIdx[!is.na(clustIdx)]) # indices of all clusters
  nClust <- length(clustLevel) # number clusters
  clustMassVec <- rep(NA, nClust) # initialize
  
  ## Loop over cluster:
  for (iClust in 1:nClust){
    idx <- which(clustIdx == iClust) # elements belonging to this cluster
    clustMassVec[iClust] <- sum(tVecValid[idx])
  }
  
  return(max(abs(clustMassVec))) # return maximum value of any cluster
  
}

# ============================================================================ #
#### 06) Compute cluster mass with cumsum: ####

computeClusterMass2 <- function(tVec, tThresh = 2){
  
  ## Exclude NAs:  
  tVecValid <- tVec[!(is.na(tVec))]
  
  ## Detect clusters above threshold:
  aboveThreshVec <- as.numeric(abs(tVecValid) > tThresh)
  # plot(aboveThreshVec)
  
  ## Cumulative sum over elements:
  cumMassVec <- cumsum(abs(tVecValid)*aboveThreshVec)
  
  ## Detect end of clusters:
  difVec <- diff(c(aboveThreshVec, 0)) # pad by 0 at end; starts and stops of clusters detectable as +1/-1
  # which(difVec == -1)
  
  ## Extract cluster mass at end of clusters:
  clustMassVec <- cumMassVec[difVec == -1] # detect where jump from 1 to 0 (i.e. -1) --> total sum of that cluster
  clustMassVec <- diff(c(0, clustMassVec)) # correct for cluster mass of previous cluster
  
  ## Determine if valid output or not:
  output <- max(c(abs(clustMassVec), 0), na.rm = T)
  
  return(output)
  
}

## Just to measure maximal speed of clusterMass function:
computeClusterMass3 <- function(tVec, tThresh = 2){
  return(1)
}

# ============================================================================ #
#### 07) Permutation test function for difference between two conditions: ####

test_permutation_timecourse <- function(data, permVar, selTrialVar = NULL, selTrialVal = NULL,
                                        eyeMeasure = "pupil", aggrMetric = NULL,
                                        segMessage = "StartMask", beforeSamples = 1000, afterSamples = 3000, 
                                        timeLim = NULL, decimate = NULL, 
                                        isBaselineCorGroup = FALSE, isBaselineCorSub = FALSE, baselineTime = 0,
                                        tThresh = 2, nIter = 10000){
  #' Perform cluster-based permutation test on condition difference in pupil time course.
  #' First plots t-value time course of average condition difference across subjects.
  #' If no cluster above threshold, stop.
  #' If any cluster above threshold detected, perform permutations.
  #' At the end: plot permutation null distribution and empirical cluster mass value, output permutation p-value.
  #' @param data                data frame with trial-level behavioral data.
  #' @param permVar             scalar string, name of variable to permute over. Must be factor with two conditions.
  #' @param selTrialVar         string, variable based on which to select subset of trials (optional).
  #' @param selTrialVal         integer or string, value to select for on selTrialVar (optional).
  #' @param eyeMeasure          string, eye metric to load in ("pupil" or "gaze", default: "pupil").
  #' @param aggrMetric          string, aggregation metric to use (default: same as eyeMeasure).
  #' @param segMessage          scalar string, name of message relative to which data that should be loaded was epoched (default: "StartMask").
  #' @param beforeSamples       scalar integer, number of samples to be included before message relative to which data was epoched (default: 1000).
  #' @param afterSamples        scalar integer, number of samples to be included after message relative to which data was epoched (default: 1966).
  #' @param timeLim             vector of two integers, time window to select for permutation test (in units of ms; optional).
  #' @param decimate            scalar integer, factor by which to downsample data (optional).
  #' @param isBaselineCorGroup  Boolean, set average condition difference across participants to 0 at time point 0 (default: FALSE).
  #' @param isBaselineCorSub    Boolean, set condition difference for each participant to 0 at time point 0 (default: FALSE).
  #' @param baselineTime        scalar integer (single time point), vector of two integers (start and stop), or nCond x 2 matrix with baseline window per condition in ms (default: 0).
  #' @param tThresh             scalar numeric, absolute t-value to threshold by for creating clusters (default: 2).
  #' @param nIter               scalar integer, number of iterations (default: 10000).
  #' @return pVal               Return permutation p-value.
  
  cat(paste0("Perform cluster based permutation test given ", permVar, " for ", nIter, " iterations with threshold t = ", tThresh, "\n"))
  if(!(permVar %in% names(data))){stop("No variable ", permVar, " in data")}  
  if(is.numeric(data[, permVar])){stop("permVar must be a factor")}
  if(length(levels(data[, permVar])) != 2){stop("permVar must have 2 levels")}
  
  ## Set default settings:
  if (eyeMeasure == "pupil"){aggrMetric <- eyeMeasure}
  
  # -------------------------------------------------------------------------- #
  ### Load array of time course per subject per condition:
  
  condArray <- aggregate_timecourse_condition(data = data, zVar = permVar, selTrialVar = selTrialVar, selTrialVal = selTrialVal,
                                              eyeMeasure = eyeMeasure, aggrMetric = aggrMetric,
                                              segMessage = segMessage, beforeSamples = beforeSamples, afterSamples = afterSamples)
  
  # -------------------------------------------------------------------------- #
  ### Data settings:
  
  ## Retrieve input dimensions:
  nSub <- dim(condArray)[1]
  nSample <- dim(condArray)[2]
  nCond <- dim(condArray)[3]
  
  if(nCond != length(levels(data[, permVar]))){stop("Number conditions in condArray and levels of permVar in data do not match")}
  
  ## Create time variable:
  timeVecMs <- seq(from = -1*beforeSamples, to = afterSamples, by = 1) # add 1 because of sample at 0
  
  # -------------------------------------------------------------------------- #
  ### Compute difference between conditions in actual data implementation:
  
  difMat <- condArray[, , 1] - condArray[, , 2] # compute condition difference
  
  # -------------------------------------------------------------------------- #
  ### Perform baseline correction:
  
  ## Apply baseline correction per subject:
  if(isBaselineCorSub){
    cat("Apply baseline correction at subject level\n")
    if (length(baselineTime) == 1){baselineIdx <- which(timeVecMs == baselineTime)} # single number
    if (length(baselineTime) == 2){baselineIdx <- which(timeVecMs %in% seq(baselineTime[1], baselineTime[2], 1))} # vector of start/ end
    if (is.matrix(baselineTime)){baselineIdx <- which(timeVecMs %in% seq(baselineTime[iCond, 1], baselineTime[iCond, 2], 1))} # vector of start/end per condition
    subBase <- difMat[, baselineIdx]
    baseMat <- replicate(nSample, subBase)
    difMat <- difMat - baseMat
  }
  # difMat[, which(timeVecMs == 0)]
  # mean(difMat[, which(timeVecMs == 0)], na.rm = T)
  
  # Apply baseline correction on group-level:
  if(isBaselineCorGroup){
    cat("Apply baseline correction at group level\n")
    if (length(baselineTime) == 1){baselineIdx <- which(timeVecMs == baselineTime)} # single number
    if (length(baselineTime) == 2){baselineIdx <- which(timeVecMs %in% seq(baselineTime[1], baselineTime[2], 1))} # vector of start/ end
    # if (is.matrix(baselineTime)){baselineIdx <- which(timeVecMs %in% seq(baselineTime[iCond, 1], baselineTime[iCond, 2], 1))} # vector of start/end per condition
    groupBase <- mean(difMat[, baselineIdx], na.rm = T)
    baseMat <- replicate(nSample, rep(groupBase, nSub))
    difMat <- difMat - baseMat
  }
  # mean(difMat[, which(timeVecMs == 0)], na.rm = T)
  
  # -------------------------------------------------------------------------- #
  ### Sub-select time window:
  
  if(!is.null(timeLim)){
    if(length(timeLim) != 2){stop("timeLim must be vector of two integers")}
    cat(paste0("Select time window from ", timeLim[1], " - ", timeLim[2], " ms\n"))
    selTimeVec <- seq(timeLim[1], timeLim[2], 1)
    selTimeIdx <- timeVecMs %in% selTimeVec
    difMat <- difMat[, selTimeIdx]
    nSample <- dim(difMat)[2] # update
    timeVecMs <- selTimeVec # update
  }
  
  # -------------------------------------------------------------------------- #
  ### Down-sample data:
  
  if(!is.null(decimate)){
    cat(paste0("Downsample data by factor ", decimate, "\n"))
    if (min(timeVecMs) < 0 & max(timeVecMs) > 0){
      selTimeVec <- sort(unique(c(-1 * seq(0, abs(min(timeVecMs)), decimate), seq(0, abs(max(timeVecMs)), decimate))))
    } else {
      selTimeVec <- seq(ceil(min(timeVecMs)/decimate)*decimate, floor(max(timeVecMs)/decimate)*decimate, decimate)
    }
    selTimeIdx <- timeVecMs %in% selTimeVec
    difMat <- difMat[, selTimeIdx]
    nSample <- dim(difMat)[2] # update
    timeVecMs <- selTimeVec # update
  }  
  
  ## Convert to seconds:
  timeVecSec <- timeVecMs/1000 # to sec.
  
  # -------------------------------------------------------------------------- #
  ### Compute t-value mass for actual data implementation:
  
  meanVec <- colMeans(difMat, na.rm = T) # compute mean across subjects
  seVec <- apply(difMat, 2, sd, na.rm = T) / sqrt(nSub) # compute standard error across subjects 
  tVec <- meanVec / seVec # compute t-statistics
  # plot(tVec)
  
  # tVal <- computeClusterMass(tVec, tThresh)
  # tVal <- computeClusterMass2(tVec, tThresh)
  # tVal <- computeClusterMass3(tVec, tThresh)
  
  # -------------------------------------------------------------------------- #
  ### Detect cluster mass via cumsum:
  
  ## Exclude NAs:  
  valIdx <- which(!(is.na(tVec)))
  tVecValid <- tVec[valIdx]
  # plot(tVecValid)
  
  ## Detect clusters above threshold:
  aboveThreshVec <- as.numeric(abs(tVecValid) > tThresh)
  # plot(aboveThreshVec)
  
  ## Cumulative sum over elements:
  cumMassVec <- cumsum(abs(tVecValid)*aboveThreshVec)
  # plot(cumMassVec)
  
  ## Detect end of clusters:
  stopVec <- diff(c(aboveThreshVec, 0)) == -1 # pad by 0 at end; starts and stops of clusters detectable as +1/-1
  
  ## Extract cluster mass at end of clusters:
  clustMassVec <- cumMassVec[stopVec] # extract cluster mass at jumps from 1 to 0 (i.e. -1) --> total sum of that cluster
  clustMassVec <- diff(c(0, clustMassVec)) # correct for cluster mass of previous cluster
  tVal <- max(c(abs(clustMassVec), 0), na.rm = T)
  cat(paste0("Empirical t = ", tVal, "\n"))
  
  # https://cran.r-project.org/web/packages/BiocManager/vignettes/BiocManager.html
  # https://stackoverflow.com/questions/61777568/r-counting-the-number-of-objects-in-an-image-using-bwlabel
  # https://search.r-project.org/CRAN/refmans/mgc/html/ConnCompLabel.html#:~:text=Description,SpatialGridDataFrame'%20(sp%20package).
  # https://rdrr.io/cran/wvtool/man/cc.label.html
  # https://rdrr.io/cran/imager/man/label.html
  
  # -------------------------------------------------------------------------- #
  ### If no clusters above threshold: skip permutations:
  if (tVal == 0){
    
    pVal <- 1
    cat(paste0("No clusters above threshold"))
    
  } else {
    
    ## Start time:
    start.time <- Sys.time()
    
    # ------------------------------------------------------------------------ #
    ### Compute location of all clusters above threshold:
    
    nClust <- length(clustMassVec)
    difVec <- diff(c(0, aboveThreshVec)) 
    if(nClust > 1){
      for (iClust in 1:nClust){
        startIdx <- which(difVec == 1)[iClust] # begin of cluster
        stopIdx <- which(difVec == -1)[iClust] # end of cluster
        clustTime <- c(timeVecSec[valIdx[startIdx]], timeVecSec[valIdx[stopIdx]]) # save
        cat(paste0("Cluster above threshold no. ", iClust, ": t = ", round(clustMassVec[iClust], 3), 
                   " from ", clustTime[1], " - ", clustTime[2], " sec.\n"))
      }
    }
    
    # ------------------------------------------------------------------------ #
    ### Compute location of largest cluster:
    
    maxIdx <- which(clustMassVec == max(clustMassVec)) # which one is biggest
    difVec <- diff(c(0, aboveThreshVec)) 
    startIdx <- which(difVec == 1)[maxIdx] # begin of cluster
    difVec <- diff(c(aboveThreshVec, 0)) 
    stopIdx <- which(difVec == -1)[maxIdx] # end of cluster
    largestClustVec <- c(timeVecSec[valIdx[startIdx]], timeVecSec[valIdx[stopIdx]]) # save
    cat(paste0("Largest cluster from ", largestClustVec[1], " - ", largestClustVec[2], " sec.\n"))
    yMin <- min(tVec, na.rm = T)
    yMax <- max(tVec, na.rm = T)
    if (yMin < 0 & yMax < 0){yMax <- 0}
    if (yMin > 0 & yMax > 0){yMin <- 0}
    plot(timeVecSec, tVec, type = "l", ylim = c(yMin, yMax),
         xlab = "Time (in sec.)", ylab = "t-value",
         main = "Empirical t-value time course")
    abline(h = 1*tThresh, lty = 2)
    abline(h = -1*tThresh, lty = 2)
    
    # ------------------------------------------------------------------------ #
    ### Compute maximum cluster mass in actual data implementation:
    
    cat(paste0("Compute t-value for ", nIter, " permutations...\n"))
    tPermVec <- rep(NA, nIter)
    
    ### Loop over iterations, permute:
    for (iIter in 1:nIter){ # iIter <- 1
      
      ## Progress bar?
      nCheck <- 100 # check every 10%
      nStep <- nIter / nCheck
      if(iIter/nStep == round(iIter/nStep)){ProgressBar(iIter, nIter)}
      # if(iIter/nStep == round(iIter/nStep)){cat(paste0("Progress: ", iIter/nCheck, "%\n"))}
      
      ## Flip sign randomly for each subject:
      signVec <- sample(c(1, -1), nSub, replace = TRUE)
      permMat <- difMat * replicate(nSample, signVec)
      # plot(difMat[3, ])
      # plot(permMat[3, ])
      
      ## Compute mean/ SE/ t for difference across subjects:
      permMeanVec <- colMeans(permMat, na.rm = T) # compute mean across subjects
      permSeVec <- apply(permMat, 2, sd, na.rm = T) / sqrt(nSub) # compute standard error across subjects 
      permTVec <- permMeanVec / permSeVec # compute t-statistics
      # plot(permTVec)
      
      ## Compute maximum cluster mass in iteration:
      # tPermVec[iIter] <- computeClusterMass(permTVec, tThresh)
      # tPermVec[iIter] <- computeClusterMass2(permTVec, tThresh)
      # tPermVec[iIter] <- computeClusterMass3(permTVec, tThresh)
      
      # ---------------------------------------------------------------------- #
      ### Detect cluster mass via cum-sum:
      
      ## Exclude NAs:  
      tVecValid <- permTVec[!(is.na(permTVec))]
      # plot(tVecValid)
      
      ## Detect clusters above threshold:
      aboveThreshVec <- as.numeric(abs(tVecValid) > tThresh)
      # plot(aboveThreshVec)
      
      ## Cumulative sum over elements:
      cumMassVec <- cumsum(abs(tVecValid)*aboveThreshVec)
      # plot(cumMassVec)
      
      ## Detect end of clusters:
      stopVec <- diff(c(aboveThreshVec, 0)) == -1 # pad by 0 at end; starts and stops of clusters detectable as +1/-1
      
      ## Extract cluster mass at end of clusters:
      clustMassVec <- cumMassVec[stopVec] # extract cluster mass at jumps from 1 to 0 (i.e. -1) --> total sum of that cluster
      clustMassVec <- diff(c(0, clustMassVec)) # correct for cluster mass of previous cluster
      tPermVec[iIter] <- max(c(abs(clustMassVec), 0), na.rm = T)
      
    }
    
    # ---------------------------------------------------------------------- #
    ### Compute p-value and print to console:
    
    ## Stop time:
    end.time <- Sys.time(); beep()
    dif <- difftime(end.time,start.time); print(dif)
    
    ## Print to console:
    cat(paste0("Empirical t = ", tVal, "; maximal t in permutation distribution = ", max(abs(tPermVec)), "\n"))
    cat(paste0("Largest cluster from ", largestClustVec[1], " - ", largestClustVec[2], " sec.\n"))
    pVal <- mean(abs(tPermVec) > abs(tVal))
    cat(paste0("p = ", pVal, "\n"))
    
    ## Plot permutation distribution:  
    print(densityplot(tPermVec, 
                      xlim = c(0, max(c(tVal, tPermVec)) + 1),
                      xlab = "Permutation t-values",
                      main = paste0("Permutation t-distribution; empirical t-value in red;\np = ", pVal, "; ", nIter, " iterations"),
                      panel = function(x,...){
                        panel.densityplot(x,...)
                        panel.abline(v = tVal, col.line = "red", lwd = 3, lty = 2) 
                      }))
    
  }
  
  return(pVal)
  
}

# END OF FILE.
