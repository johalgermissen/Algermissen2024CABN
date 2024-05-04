#!/usr/bin/env Rscript
# ============================================================================ #
## 00_mgngunc_functions_preprocess_pupil.R
## MGNGUncertainty functions for pre-processing of pupillometry data.
# Johannes Algermissen, 2023.

## Implementing a pipeline used in Tobias Donner's lab, as specified in:
# Colizoli, O., De Gee, J. W., Urai, A. E., & Donner, T. H. (2018). Task-evoked pupil responses reflect internal belief states. Scientific reports, 8(1), 13702.
# de Gee, J. W., Knapen, T., & Donner, T. H. (2014). Decision-related pupil dilation reflects upcoming choice and individual bias. Proceedings of the National Academy of Sciences, 111(5), E618-E625.
# de Gee, J. W., Colizoli, O., Kloosterman, N. A., Knapen, T., Nieuwenhuis, S., & Donner, T. H. (2017). Dynamic modulation of decision biases by brainstem arousal systems. Elife, 6, e23232.
# Urai, A. E., Braun, A., & Donner, T. H. (2017). Pupil-linked arousal is driven by decision uncertainty and alters serial choice bias. Nature communications, 8, 14637.

# ============================================================================ #
#### 01) Read data into trials ####

readEyelink <- function(fileName){
  #' Just a wrapper for function \code{read.asc} from package \code{eyelinker},
  #' so everything is together.
  #' @param fileName string argument, file name ending on \code{.asc}.
  #' @return data object read in with \code{read.asc}, has fields:
  #'  \code{raw}, \code{msg}, \code{sacc}, \code{fix}, \code{blinks}, \code{info}.
  
  require(eyelinker)
  # https://cran.r-project.org/web/packages/eyelinker/vignettes/basics.html
  # https://cran.r-project.org/web/packages/eyelinker/eyelinker.pdf
  
  ## Load data:
  cat(paste0("Read file ", fileName, " of size ", file.size(fileName)/1000000, " MB\n"))
  data <- read.asc(fileName)
  cat(paste0("Output size: ", object.size(data)/1000000, " MB\n"))
  
  return(data)
}

# ============================================================================ #
#### 02) Segment pupil into trials: ####

readTrials <- function(rawData, segmentMessage, beforeEvent, afterEvent, 
                       varMessages = NULL, recode2Num = NULL, varNames = NULL, isVerbose = FALSE){
  #' Epoch data into trials based on occurrence of \code{segmentMessage} as time point zero, 
  #' with \code{beforeEvent} as positive integer of samples before \code{segmentMessage} 
  #' and \code{afterEvent} as positive integer of samples after \code{segmentMessage};
  #' can also read other variables if written in eyelink file
  #' returns epoched data with all trials features retrieved.
  #' @param rawData raw data read-in with read.asc from package \code{Eyelinker}.
  #' @param segmentMessage string, message to be used for segmenting data into trials, used to set time point 0.
  #' @param beforeEvent positive integer, number of samples to be included before \code{segmentMessage}.
  #' @param afterEvent positive integer, number of samples to be included after \code{segmentMessage}.
  #' @param varMessages vector of strings, contains messages for which values will be retrieved (assumes value one row below message).
  #' @param recode2Num vector of booleans, specify for each variable whether to recode or not.
  #' @param varNames vector of strings, new variable names for extracted variables in output object; default is NULL, which means \code{variablesMessage} will be used as variable names.
  #' @param isVerbose Boolean, print detailed messages about every trial to console or not (default: FALSE).
  #' @return data frame with epoched data.
  #' @examples 
  #' rawData <- read.asc("s_01.asc")
  #' data <- readTrials(rawData, segmentMessage = "StartMask", beforeEvent = 1000, afterEvent = 3000)
  #' (Consider letting messages for other data to-be-retrieved be specified in input??)
  
  
  cat("Epoch data into trials\n")
  
  # -------------------------------------------------------------------------- #
  ### Check whether varMessages, recode2Num, and varNames are of same length:
  if (!(is.null(varMessages))){
    
    if(length(varMessages) != length(recode2Num)){
      stop("Arguments varMessages and recode2Num are of different lengths")
    }
    
    if(is.null(varNames)){
      varNames <- varMessages
      cat("Argument varNames unspecified, will use varMessages as variable names\n")
    }
    
    if(length(varMessages) != length(varNames)){
      stop("Arguments varMessages and varNames are of different lengths")
    }
  } # end of if-clause if varMessages exist
  
  # -------------------------------------------------------------------------- #
  ### Extract raw data:
  
  raw <- rawData$raw
  
  # -------------------------------------------------------------------------- #
  ### Determine time of start of forward mask:
  
  maskRow <- which(rawData$msg$text == segmentMessage) # rows in msg where selected message ocurred
  maskTime <- rawData$msg$time[maskRow] # timing of those messages
  nTrial <- length(maskTime)
  startTime <- maskTime - beforeEvent # start beforeEvent ms before segmentMessage
  stopTime <- maskTime + afterEvent # start beforeEvent ms before segmentMessage
  
  # -------------------------------------------------------------------------- #
  ### Extract selected variables for each trial:
  
  if (!(is.null(varMessages))){
    
    varList <- list() # initialize
    nVar <- length(varMessages) # count
    
    ## Loop through variables:
    for (iVar in 1:nVar){ # iVar <- 1
      
      # Read variable values (appears one row below variable name):
      varVec <- rawData$msg$text[which(rawData$msg$text == varMessages[iVar])+1]
      
      if (recode2Num[iVar] == T){ # if numeric: recode into numeric
        
        varVec <- as.numeric(varVec)
        
      } else { # else: assume factor, only remove unused levels
        
        varVec <- droplevels(factor(varVec))
        
      } # end if-clause for recoding to numeric
      
      varList[[iVar]] <- varVec # save variable vector into list
      if (isVerbose & sum(is.na(varVec)) > 0){cat(paste0("Variable ", varMessages[iVar], " contains NAs!!!!!!!!!!!!!!!!!\n"))}
      
    } # end for-loop iVar
  } # end if-clause varMessages are not empty
  
  # -------------------------------------------------------------------------- #
  ### Loop over trials, extract pupil size during trial:
  
  trialList <- list() # initialize
  
  for (iTrial in 1:nTrial){ # iTrial = 1
    
    
    if(isVerbose){cat(paste0("Read trial ", iTrial, "\n"))}
    
    ## Determine rows where to start and stop extracting data:
    startRow <- which(raw$time == startTime[iTrial])
    stopRow <- which(raw$time == stopTime[iTrial])
    
    # Check if rows exist:
    if (length(stopRow) == 0){
      if(isVerbose){cat(paste0("Trial ", iTrial, ": startRow undefined, consider shortening beforeEvent"))}
    }
    if (length(stopRow) == 0){
      if(isVerbose){cat(paste0("Trial ", iTrial, ": stopRow undefined, consider shortening afterEvent"))}
    }
    
    ## Extract pupil as vector:
    pupil <- raw$ps[startRow:stopRow]
    
    ## Initialize vectors with trial number and timing for this trial:
    trialDur <- length(pupil) # number of pupil samples determines trial duration
    trialnr <- rep(iTrial, trialDur) # trial number for this trial
    absTime <- raw$time[startRow:stopRow] # extract total time for this trial
    trialTime <- seq(from = -1*beforeEvent, to = (trialDur - beforeEvent - 1), by = 1) # time relative to segmentMessage 
    
    ## Trial features:
    if (!(is.null(varMessages))){
      
      varMat <- as.data.frame(matrix(NA, trialDur, nVar)) # initialize empty data frame
      
      ## Loop through variables, extract values for this trial for selected variables:
      for (iVar in 1:nVar){
        varVec <- varList[[iVar]] # extract this variable
        varVal <- varVec[iTrial] # extract value of variable for this trial 
        varMat[, iVar] <- rep(varVal, trialDur) # repeat trial feature as often as samples within trial:
      } # end for-loop iVar
      
      colnames(varMat) <- varNames # set variable names
      
      # Concatenate into data frame with varMat, save n list:
      trialList[[iTrial]] <- as.data.frame(cbind(trialnr, absTime, trialTime, pupil, varMat))
      
    } else {
      
      # Concatenate into data frame without varMat:
      trialList[[iTrial]] <- as.data.frame(cbind(trialnr, absTime, trialTime, pupil))
      
    } # end if-clause varMessages are not empty
  } # end for-loop iTrial
  
  # -------------------------------------------------------------------------- #
  ### Append data frames of all trials:
  
  data <- do.call(rbind, trialList)
  
  ## Warn about warning messages when converting factors to numeric:
  if (!(is.null(recode2Num)) & sum(recode2Num) > 0){
    if(isVerbose){cat(paste0("Expect warning about 'NAs introduced by coercion' for every variable converted to numeric, i.e. ", sum(recode2Num), " warnings\n"))}
  } # end if recode2Num > 0
  
  return(data)
} # end function

# ============================================================================ #
#### 03) Response-lock data: ####

respLock <- function(data, rtVar = "RT", timeVar = "trialTime"){
  #' Plots epoched data on trial-by-trial basis.
  #' @param data data-frame, trial-epoched data (with variables \code{trialnr}, \code{trialTime}, and \code{pupil})
  #' @param selVar string, variable, with RTs, default "RT".
  #' @return no output, just plotting.
  
  ## Fixed settings:
  # data <- trialData
  trialVar <- "trialnr"
  timeVar <- "trialTime"
  
  respVar <- "response"
  valenceVar <- "valence"
  
  rtVar <- "RT";
  interpRtVar <- "RT_interp"
  
  cat("Response-lock pupil data\n")
  
  # ------------------------------------------------------------------------- #
  ### Interpolate RTs for NoGo trials:
  
  ## Retain one row per trial:
  perTrialData <- data[which(data[, timeVar] == 0), ]
  
  ## Mean RTs for Go2Win:
  G2WtrlIdx <- which(perTrialData[, respVar] == 1 & perTrialData[, valenceVar] == 1)
  G2WRT <- mean(perTrialData[G2WtrlIdx, rtVar], na.rm = T) 
  cat(paste0("Mean RTs for Go2Win is ", G2WRT, " sec.\n"))
  
  ## Mean RTs for Go2Avoid:
  G2AtrlIdx <- which(perTrialData[, respVar] == 1 & perTrialData[, valenceVar] == 0)
  G2ART <- mean(perTrialData[G2AtrlIdx, rtVar], na.rm = T) 
  cat(paste0("Mean RTs for Go2Avoid is ", G2ART, " sec.\n"))
  
  ## Identify NoGo trials to be interpreted:  
  NG2WtrlIdx <- which(perTrialData[, respVar] == 0 & perTrialData[, valenceVar] == 1)
  stopifnot(all(is.na(perTrialData[NG2WtrlIdx, rtVar])))
  cat(paste0("Found ", length(NG2WtrlIdx), " trials with NoGo responses to Win cues\n"))
  
  NG2AtrlIdx <- which(perTrialData[, respVar] == 0 & perTrialData[, valenceVar] == 0)
  stopifnot(all(is.na(perTrialData[NG2AtrlIdx, rtVar])))
  cat(paste0("Found ", length(NG2AtrlIdx), " trials with NoGo responses to Avoid cues\n"))
  
  ## Interpolate RTs for NoGos:
  cat(paste0("Interpolate RTs of ", round(G2WRT*1000), " for NG2W trials and ", round(G2ART*1000), " for NG2A trials\n"))
  data[, interpRtVar] <- data[, rtVar] # copy over, contains RTs for G2W and G2A
  data[which(data[, trialVar] %in% NG2WtrlIdx), interpRtVar] <- G2WRT
  data[which(data[, trialVar] %in% NG2AtrlIdx), interpRtVar] <- G2ART
  stopifnot(all(!is.na(data[, interpRtVar])))
  
  # ------------------------------------------------------------------------- #
  ### Loop over trials, subtract RT from time:
  # if RT is 500 ms, subtract 500 ms from all time points, so former 500 ms becomes new 0 ms
  
  nTrial <- max(data[, trialVar])
  for(iTrial in 1:nTrial){ # iTrial <- 1
    rowIdx <- which(data[, trialVar] == iTrial) # rows for this trial
    trlRT <- round(data[rowIdx[1], interpRtVar] * 1000) # RT for this trial, in ms, round to integers
    data[rowIdx, timeVar] <- data[rowIdx, timeVar] - trlRT
  }
  
  return(data)  
}

# ============================================================================ #
#### 04) Plot epoched data: ####

plotTrials <- function(data, selVar = "pupil"){
  #' Plots epoched data on trial-by-trial basis.
  #' @param data data-frame, trial-epoched data (with variables \code{trialnr}, \code{trialTime}, and \code{pupil})
  #' @param selVar selVar to-be-plotted.
  #' @return no output, just plotting.
  
  nTrial <- length(unique(data$trialnr)) # count trials
  varIdx <- which(colnames(data) == selVar) # index of selected variable
  
  cat("Plot trials\n")
  
  # -------------------------------------------------------------------------- #
  ## Loop over trials and plot:
  for (iTrial in 1:nTrial){ # iTrial <- 1
    
    cat(paste0("Plot trial ", iTrial, "\n"))
    
    rowIdx <- which(data$trialnr == iTrial)
    plot(data$trialTime[rowIdx], data[rowIdx, varIdx], type = "l", 
         xlab = "Relative trial time (in ms)", ylab = selVar, main = paste0("Trial ", iTrial))
    readline(prompt = "Press [enter] to continue")
  }
  
  cat("Finished! :-)\n")
  
}  

# ============================================================================ #
#### 05) Delete blinks: ####

deleteBlkSac <- function(rawData, trialData, selVar = "pupil",
                         removeBlinks = TRUE, removeSaccades = TRUE,
                         zeropadBlinks = 150, zeropadSaccades = 0, percentReject, isVerbose = FALSE){
  #' Extract timing of blinks and saccades automatically detected by Eyelink from \code{dat},
  #' zero-pad by \code{zeropad} ms before and after blink/saccade,
  #' find blinks and saccades back in \code{data} and set samples to \code{NA}.
  #' Automatically reject trials if more than \code{percentReject} % deleted.  
  #' @param rawData raw data read in via read.asc from eyelinker (needed to extract blinks and saccades).
  #' @param trialData data-frame, trial-epoched data.
  #' @param zeropad positive integer, also delete samples \code{zeropad} ms before and after blink/saccade (150 ms recommended by SR Research).
  #' @param percentReject positive integer between 0 and 100, set trial to invalid if more than \code{percentReject} % of samples set to NA.
  #' @param isVerbose Boolean, print detailed messages about every trial to console or not (default: FALSE).
  #' @return trialData data frame with \code{selVar} set to NA for blinks;
  #' variable \code{valid} set to 0 for trials with more than \code{percentReject} set to NA;
  #' variable \code{isInterpreted} set to 1 for rows with samples of \code{selVar} deleted .
  
  require(tidyr) # for pipe
  require(intervals)
  
  cat(paste0("Delete samples identified as blinks (padded by ", zeropadBlinks, " ms) and saccades (padded by ", zeropadSaccades, " ms) in variable \'", selVar, "\'\n"))
  
  # -------------------------------------------------------------------------- #
  ### Initialize book-keeping variables for interpolation (valid, isInterp)
  
  ## Set all trials to valid as default:
  if (!("valid" %in% colnames(trialData))){
    trialData$valid <- 1
  }
  
  ## Set all trials to not-interpolated:
  if (!("isInterp" %in% colnames(trialData))){
    trialData$isInterp <- 0
  }
  
  # -------------------------------------------------------------------------- #
  ### Retrieve any detected blinks or saccades
  
  if(removeBlinks == TRUE){
    
    ## a) Extract blinks from eyelink data as intervals: 
    Blk <- cbind(rawData$blinks$stime, rawData$blinks$etime) #Define a set of intervals
    ## Remove extremely long blinks:
    blkDur <- Blk[, 2] - Blk[, 1]
    # sum(blkDur > 1000)
    Blk <- Blk[blkDur < 1000, ]
    ## Zero-pad blinks:
    BlinkSuspect <- Intervals(Blk) %>% intervals::expand(zeropadBlinks, "absolute")
    
  }
  
  # plot(trialData$pupil[trialData$absTime > 1235450 & trialData$absTime < 1642742])
  
  if(removeSaccades == TRUE){
    
    ## b) Extract saccades from eyelink data as intervals: 
    Sac <- cbind(rawData$sac$stime,rawData$sac$etime) #Define a set of intervals
    ## Remove extremely long blinks:
    sacDur <- Sac[, 2] - Sac[, 1]
    # sum(sacDur > 1000)
    Sac <- Sac[sacDur < 1000, ]
    ## Zero-pad saccades:
    SacSuspect <- Intervals(Sac) %>% intervals::expand(zeropadSaccades, "absolute")
    
  }
  
  # -------------------------------------------------------------------------- #
  ### Loop over trials, delete blinks and saccades:
  
  nTrial <- length(unique(trialData$trialnr)) # count trials
  nTrialInvalid <- 0 # initialize
  
  for (iTrial in 1:nTrial){ # iTrial <- 50
    
    if(isVerbose){cat(paste0("Trial ", iTrial, ": Locate blinks\n"))}
    # plot(trialData$pupil[trialData$trialnr == iTrial])
    
    # iTrial <- 50
    # Retrieve start and end of trial, put into interval:
    rowIdx <- which(trialData$trialnr == iTrial)
    startTrial <- min(trialData$absTime[rowIdx])
    endTrial <- max(trialData$absTime[rowIdx])
    trialDur <- endTrial - startTrial
    trialInterval <- Intervals(c(startTrial, endTrial)) # concatenate into Intervals object
    
    # ------------------------------------------------------------------------ #
    ## a) Compare trial with blinks:
    if(removeBlinks == TRUE){
      
      allBlinks <- interval_overlap(trialInterval, BlinkSuspect, check_valid = TRUE)[[1]]
      nBlink <- length(allBlinks) # count blinks
      
      ## If any blinks during trial detected:
      if (nBlink > 0){
        
        if(isVerbose){cat(paste0("Trial ", iTrial, ": Detected ", nBlink, " blinks\n"))}
        blinkDur <- rep(NA, nBlink) # store duration of removed blinks
        
        for (iBlink in 1:nBlink){ # iBlink <- 1
          
          ## Determine start and end of blink within trial:
          startBlink <- max(BlinkSuspect[allBlinks[iBlink], 1], startTrial) # keep later of both events
          endBlink <- min(BlinkSuspect[allBlinks[iBlink], 2], endTrial) # keep earlier of both events
          blinkDur[iBlink] <- endBlink - startBlink # compute duration
          
          ## Determine rows to be set to NA:
          startRow <- which(trialData$absTime == startBlink)
          endRow <- which(trialData$absTime == endBlink)
          trialData[startRow:endRow, selVar] <- NA # set to NA
          trialData$isInterp[startRow:endRow] <- 1 # set to 1
        } # end for-loop iBlink
        
        ## Report number of detected blinks and saccades:
        blinkPerc <- sum(blinkDur)/trialDur*100
        if(isVerbose){cat(paste0("Trial ", iTrial, ": Removed ", round(blinkPerc, 2), "% of trial due to blinks\n"))}
        
      } # end if clause for blinks
    } # end if clause for removeBlinks
    
    # ------------------------------------------------------------------------ #
    ### b) Compare trial with saccades:
    if(removeSaccades == TRUE){
      
      allSacs <- interval_overlap(trialInterval, SacSuspect, check_valid = TRUE)[[1]]
      nSac <- length(allSacs) # count saccades
      
      ## If any saccades during trial detected:
      if (nSac > 0){
        
        if(isVerbose){cat(paste0("Trial ", iTrial, ": Detected ", nSac, " saccades\n"))}
        sacDur <- rep(NA, nSac) # store duration of removed saccades
        
        ## Loop through saccades:
        for (iSac in 1:nSac){ # iSac <- 1
          
          # Determine start and end of saccade within trial:
          startSac <- max(SacSuspect[allSacs[iSac], 1], startTrial) # keep later of both events
          endSac <- min(SacSuspect[allSacs[iSac] ,2], endTrial) # keep earlier of both events
          sacDur[iSac] <- endSac - startSac # compute duration
          
          # Determine rows to be set to NA:
          startRow <- which(trialData$absTime == startSac)
          endRow <- which(trialData$absTime == endSac)
          trialData[startRow:endRow, selVar] <- NA # set to NA
          trialData$isInterp[startRow:endRow] <- 1 # set to 1
        } # end for-loop iSac
        
        ## Report number of detected saccades:
        sacPerc <- sum(sacDur)/trialDur*100
        if(isVerbose){cat(paste0("Trial ", iTrial, ": Removed ", round(sacPerc, 2), "% of trial due to saccades\n"))}
        
      } # end if clause for Sac
    } # end if clause removeSaccades
    
    # -------------------------------------------------------------------------- #
    ### Determine total amount of deleted samples:
    
    rowIdx <- which(trialData$trialnr == iTrial) # samples belonging to this trial
    
    # a) If complete trial is NA:
    if (sum(is.na(trialData$pupil[rowIdx])) == length(rowIdx)){
      if(isVerbose){cat(paste0("Trial ", iTrial,": only NA, set to invalid!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"))}
      trialData$valid[rowIdx] <- 0 # set to 0
      trialData$pupil[rowIdx] <- 0 # set to 0 (otherwise hard to pre-process at later stages)
      nTrialInvalid <- nTrialInvalid + 1
    }
    
    # b) Determine proportion interpolated trials:
    rejectPerc <- mean(trialData$isInterp[rowIdx])*100
    if (rejectPerc > percentReject){
      if(isVerbose){cat(paste0("Trial ", iTrial, ": Set to invalid!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"))}
      trialData$valid[rowIdx] <- 0 # set to zero
      nTrialInvalid <- nTrialInvalid + 1 # increment
    } # end if-loop
    
    # Sys.sleep(3)
    
  } # end for-loop iTrial
  if(isVerbose){cat("Finished!\n")}
  
  ## Print number of invalid trials:
  if(isVerbose){cat(paste0("Number of invalid trials: ", nTrialInvalid, "\n"))}
  
  return(trialData)
  
} # end function

# ============================================================================ #
#### 06) Compute derivative, delete strong changes: ####

deleteHighDerivative <- function(data, selVar = "pupil", 
                                 derivThresh = 2, clusterThresh = 10, isInteractive = FALSE, isVerbose = FALSE){
  #' Compute first derivative to detect points of fast (i.e. biologically implausible) velocity
  #' Delete samples with absolute z-standardized derivative above \code{derivThresh},
  #' also delete clusters of samples of maximal length \code{clusterThresh} between samples with high velocity.
  #' @param data data-frame, trial-epoched data featuring variables  \code{trial_nr} and \code{selVar}.
  #' @param derivThresh positive integer, cut-off for absolute z-standardized derivative above which sample are deleted.
  #' @param clusterThresh positive integer, maximal size of clusters of samples between samples of high velocity which also get deleted.
  #' @param isInteractive boolean, plot data before (black) and after (green) deletion of high-velocity samples and clusters of samples between them.
  #' @param isVerbose Boolean, print detailed messages about every trial to console or not (default: FALSE).
  #' @return data-frame with variable \code{selVar} updated,
  #' also variables \code{deriv} and \code{deriv_z} added.
  # https://stackoverflow.com/questions/14082525/how-to-calculate-first-derivative-of-time-series
  
  cat(paste0("Delete samples with absolute values of standardized derivative > ", derivThresh, " in variable \'", selVar, "\'\n"))
  
  ### Count trials:
  nTrial <- length(unique(data$trialnr))
  
  # -------------------------------------------------------------------------- #
  ### Initialize derivative variable:
  
  if (!("deriv" %in% colnames(data))){
    data$deriv <- NA
  }
  if (!("deriv_z" %in% colnames(data))){
    data$deriv_z <- NA
  }
  
  # -------------------------------------------------------------------------- #
  ### Compute derivative for each trial:
  for (iTrial in 1:nTrial){ # iTrial <- 1
    
    if(isVerbose){cat(paste0("Trial ", iTrial, ": Compute derivative\n"))}
    rowIdx <- which(data$trialnr == iTrial) 
    
    ## Takes difference to previous (lag = 1) value:
    deriv <- c(0, diff(data[rowIdx, selVar]) / diff(data$trialTime[rowIdx])) # difference in pupil divided by difference in time
    
    ## z-standardize:
    deriv_z <- as.numeric(scale(deriv))
    
    # Store in data frame:
    data$deriv[rowIdx] <- deriv
    data$deriv_z[rowIdx] <- deriv_z
    
    ## Detect values for which abs(deriv_z) > cutoff:
    markIdx <- which(abs(data$deriv_z[rowIdx]) > derivThresh)
    
    ## Detect gaps > cluster_cutoff between marked samples, also delete those:
    derivMarkIdx <- c(clusterThresh*2, diff(markIdx)) # difference between consecutive marked samples (add 1 sample at beginning; delete later again)
    
    # -------------------------------------------------------------------------- #
    ### Loop over marked samples, check if only clusterThresh away from next sample:
    
    for (iIdx in 1:length(derivMarkIdx)){
      if (derivMarkIdx[iIdx] < clusterThresh){ # if difference between consecutive marked samples lower than clusterTresh
        markIdx <- c(markIdx, markIdx[iIdx-1]:markIdx[iIdx]) # append all samples in-between to existing vector
      }
    }
    markIdx <- sort(unique(markIdx)) # clean doublings and order
    
    if(isVerbose){cat(paste0("Trial ", iTrial,": Remove ", length(markIdx), " samples with high velocity\n"))}
    
    # ------------------------------------------------------------------------ #
    ### Plot variable before cleaning:
    
    if (isInteractive == TRUE){
      plot(data$trialTime[rowIdx], data[rowIdx, selVar], 
           type = "l", main = paste0("Trial ", iTrial, " - before velocity delection"))
      readline(prompt = "Press [enter] to continue")
    }
    
    # ------------------------------------------------------------------------ #
    ### Delete marked values:
    
    delIdx <- rowIdx[markIdx]
    data[delIdx, selVar] <- NA
    
    # ------------------------------------------------------------------------ #
    ### Plot variable after cleaning:
    
    if (isInteractive == TRUE){
      plot(data$trialTime[rowIdx], data[rowIdx, selVar], 
           type = "l", col = "green", main = paste0("Trial ", iTrial, " - after velocity delection"))
      readline(prompt = "Press [enter] to continue")
    } # end if for plotting
  } # end for-loop for iTrial
  if(isVerbose){cat("Finished!\n")}
  
  return(data)
  
}

# ============================================================================ #
#### 07) Interpolate all NAs: #####

applyInterpolation <- function(data, selVar = "pupil", padEdges = FALSE, isInteractive = FALSE, isVerbose = FALSE){
  #' Performs linear interpolation using na.approx from package zoo.
  #' @param padEdges boolean, fill leading and trailing NAs of trials with first/ valid non-NA.
  #' @param isInteractive boolean, plot each trial before/ after interpolation (after padding edges) or not.
  #' @param isVerbose Boolean, print detailed messages about every trial to console or not (default: FALSE).
  #' @return data frame, variable \code{selVar} manipulated.
  require(zoo) # for na.approx
  
  cat(paste0("Interpolate missing data in variable \'", selVar, "\'\n"))
  
  nTrial <- length(unique(data$trialnr)) # number of trials
  
  # ------------------------------------------------------------------------- #
  ### Loop over trials, interpolate:
  for (iTrial in 1:nTrial){ # iTrial <- 1
    
    if(isVerbose){cat(paste0("Trial ", iTrial, ": Interpolate NAs\n"))}
    rowIdx <- which(data$trialnr == iTrial) # determine samples belonging to this trial
    
    ## If interactive, then plot before interpolation, proceed upon button press:
    if (isInteractive == TRUE){
      plot(data$trialTime[rowIdx], data[rowIdx, selVar], 
           type = "l", col = "black", main = paste0("Trial ", iTrial, " - before interpolation"))
      readline(prompt = "Press [enter] to continue")
    }
    
    ## Apply interpolation:
    data[rowIdx, selVar] <- na.approx(data[rowIdx, selVar], na.rm = FALSE) # replace NAs with interpolated value
    
    ## If interactive, then plot after interpolation, proceed upon button press:
    if (isInteractive == TRUE){
      plot(data$trialTime[rowIdx], data[rowIdx, selVar], 
           type = "l", col = "green", main = paste0("Trial ", iTrial, " - after interpolation"))
      readline(prompt = "Press [enter] to continue")
    }
    
    ## Pad edges to remove artifacts after interpolation:
    if (padEdges == TRUE){
      data <- applyPadEdges(data, rowIdx, iTrial, isInteractive)
    } # end if for padding edges
  } # end for-loop iTrial
  if(isVerbose){cat("Finished!\n")}
  
  return(data)
  
}

# ============================================================================ #
#### 08) Pad edges: ####

applyPadEdges <- function(data, selVar = "pupil", rowIdx, iTrial, isInteractive = FALSE, isVerbose = FALSE){
  #' Replaces leading NAs by first non-NA value,
  #' Replaces trailing NAs by last non-NA value.
  #' @param data data frame, trial-epoched data with variable \code{ selVar}.
  #' @param selVar string argument, variable to be padded.
  #' @param rowIdx vector of integers, vector of rows of trial to be padded.
  #' @param iTrial integer, trial number (for title in plot).
  #' @param isInteractive boolean, plot trial time series after padding edges (in blue) or not.
  #' @param isVerbose Boolean, print detailed messages about every trial to console or not (default: FALSE).
  #' @return no return, just plotting
  #' Should only be used if enough empty time before baseline/ after trial 
  
  require(zoo) # for na.approx
  
  if(isVerbose){cat(paste0("Trial ", iTrial, ": Pad edges of variable \'", selVar, "\'\n"))}
  
  ## First interpolate:
  data[rowIdx, selVar] <- na.approx(data[rowIdx, selVar], na.rm = FALSE) # replace NAs with interpolated value
  
  ## Determine first and last sample:
  firstRowIdx <- head(rowIdx, 1)
  lastRowIdx <- tail(rowIdx, 1)
  
  ## Determine non-NA samples:
  filledIdx <- which(!(is.na(data[rowIdx, selVar]))) # non-NAs
  firstFilled <- rowIdx[filledIdx[1]] # first non-NA
  lastFilled <- rowIdx[filledIdx[length(filledIdx)]] # last non-NA
  
  ## Replace leading NAs with first non-NA value:
  if(firstFilled != firstRowIdx){
    data[firstRowIdx:(firstFilled - 1), selVar] <- data[firstFilled, selVar]
  }
  
  ## Replace trailing NAs with last non-NA value:
  if(lastFilled != lastRowIdx){
    data[(lastFilled + 1):lastRowIdx, selVar] <- data[lastFilled, selVar]
  }
  
  ## Plot after removal of NAs:
  if (isInteractive == TRUE){
    plot(data$trialTime[rowIdx], data[rowIdx, selVar], 
         type = "l", col = "blue", main = paste0("Trial ", iTrial, " - after padding edges"))
    readline(prompt = "Press [enter] to continue")
  }
  
  return(data)
  
}

# ============================================================================ #
#### 09) Filter data: ####

applyFilter <- function(data, selVar = "pupil", filtType = "Butterworth", order, filtPos = "low", cutoff, sf = 1000, 
                        deleteEdges = F, lengthEdges = 0, padEdges = F, isInteractive = FALSE, isVerbose = FALSE){
  #' Use butterworth or FIR filter from package \code{signal} to filter signal;
  #' @param data data-frame, trial-epoched data, needs to contain variables \code{trialnr}, \code{ selVar}.
  #' @param filtType string argument specifying filter type, either \code{Butterworth} or \code{FIR}.
  #' @param order integer, order of filter.
  #' @param filtPos string argument, either "low", "high", "stop", "pass".
  #' @param cutoff specifying upper/lower bands: if type = "low" or "high: positive integer, if type = "pass" or "stop": two-element vector.
  #' @param deleteEdges Boolean, delete edges after filter or not (recommended for Butterworth filter).
  #' @param lengthEdges integer, length edges to-be-deleted after filtering (depends on length of epoch, roughly 5% at each edge looks like a good idea).
  #' @param padEdges Boolean, interpolate NAs at start and end of trial by setting to first/ last valid value (default: FALSE).
  #' @param isVerbose Boolean, print detailed messages about every trial to console or not (default: FALSE).
  
  ## Define butterworth filter:
  # buttord(0.01,30,1,1)
  # https://stackoverflow.com/questions/7105962/how-do-i-run-a-high-pass-or-low-pass-filter-on-data-points-in-r
  # https://octave.sourceforge.io/signal/function/butter.html
  
  require(signal)
  
  cat(paste0("Filter variable \'", selVar, "\' with ", order, "th order ", filtPos, "-pass ", filtType, " filter at ", cutoff, " Hz\n"))
  
  ## Count trials:
  nTrial <- length(unique(data$trialnr))
  
  ## Divide cut-offs by Nyqvist frequency
  Nyqvist <- sf/2 # filter from package signal work in fractions of the Nyqvist frequency
  
  ## Initialize variable in case of removing edges:
  if(padEdges == FALSE){
    data$keepSegment <- 1 
  }
  
  # -------------------------------------------------------------------------- #
  ### Initialize filter function:
  
  if (filtType == "Butterworth"){
    myFilt <- butter(order, cutoff/Nyqvist, type = filtPos, plane = "z") # Define low-pass filter like Allen 2016
  } else if (filtType == "FIR"){
    myFilt <- fir1(order, cutoff/Nyqvist, type = filtPos, window = hamming(order+1), scale = TRUE)
  }  else {
    stop("Invalid filter filtType specified")
  }
  
  # -------------------------------------------------------------------------- #
  ### Loop over trials, filter each trial:
  
  for (iTrial in 1:nTrial){ # filter separately for each trial: iTrial <- 5
    
    if(isVerbose){cat(paste0("Trial ", iTrial, ": Apply filter\n"))}
    rowIdx <- which(data$trialnr == iTrial)
    
    ## If interactive: plot before filtering:
    if (isInteractive == TRUE){
      yMin <- min(data[rowIdx, selVar], na.rm = T)
      yMax <- max(data[rowIdx, selVar], na.rm = T)
      yLim <- c(yMin, yMax)
      plot(data$trialTime[rowIdx], data[rowIdx, selVar], ylim = yLim,
           type = "l", main = paste0("Trial ", iTrial, " - before filtering"))
      readline(prompt = "Press [enter] to continue")
    }
    
    ## Extract selected variable and filter using filtfilt from package signal:
    x <- data[rowIdx, selVar] # Extract selected variable fir thosi trial
    fillIdx <- which(!is.na(x)) # indices of non-NA signals relative to selected trial
    fillRowIdx <- rowIdx[fillIdx] # indices of non-NA signals relative to entire data set of subject
    x <- x[fillIdx] # extract only non-NA signal
    y <- signal::filtfilt(myFilt, x) # use filtfilt to prevent phase-shift
    
    ## Delete edges:
    effLengthEdges <- min(lengthEdges, length(y)) # to not make y longer than it is
    if (deleteEdges == T){
      y[1:(effLengthEdges)] <- NA # delete leading edge
      y[(length(y)-effLengthEdges):length(y)] <- NA # delete trailing edge
    }
    
    ## Save filtered signal back into data frame under selVar:
    data[fillRowIdx, selVar] <- y # put into data frame
    
    ## If interactive: plot after filtering:
    if (isInteractive == TRUE){
      plot(data$trialTime[rowIdx],data[rowIdx, selVar], ylim = yLim, 
           type = "l", col = "green", main = paste0("Trial ", iTrial, " - after filtering"))
      readline(prompt="Press [enter] to continue")
    }
    
    ## Option A: Interpolate also from start and until end of window: Pad edges (repeat first/last valid sample):
    if (padEdges == TRUE){ # interpolate
      data <- applyPadEdges(data, rowIdx, iTrial, isInteractive) # interpolate
    } # end if for padding edges
    
    ## Option B: Crop edges (only NAs) to shorten this trial:
    if(padEdges == FALSE){
      newMinTime <- min(data$trialTime[rowIdx]) + lengthEdges # new start of trial after removing edges
      newMaxTime <- max(data$trialTime[rowIdx]) - lengthEdges # new end of trial after removing edges
      if(isVerbose){cat(paste0("Crop data at edges to eliminate NAs; keep data from ", newMinTime, " - ", newMaxTime," \n"))}
      data$keepSegment[rowIdx] <- data$trialTime[rowIdx] %in% c(newMinTime:newMaxTime)
    }
    
    
  } # end for-loop iTrial
  
  if(padEdges == FALSE){
    data <- data[which(data$keepSegment == 1), ]
  }
  
  if(isVerbose){cat("Finished!\n")}
  
  return(data)
}

# ============================================================================ #
#### 10) Define blocks based on trial number: ####

defineBlocks <- function(data, nBlock, isVerbose = FALSE){
  #' Define block number based on trial number
  #' Assumes that blocks have the same number of trials
  #' No definition performed if variable \code{block} already exists
  #' @param data data-frame, trial-epoched data featuring variable \code{trialnr}
  #' @param isVerbose Boolean, print detailed messages about every block to console or not (default: FALSE).
  #' @return data-frame of trial-epoched data, variable \code{block} added
  
  cat("Define blocks\n")
  
  if ("block" %in% names(data)){
    cat("Variable \'block\' already exists, overwrite\n")
  }
  
  nTrial <- length(unique(data$trialnr)) # count total number of trials 
  nTrialPerBlock <- nTrial/nBlock # compute number of trials per block
  data$block <- NA # initialize
  
  ### Loop over blocks:
  for (iBlock in 1:nBlock){
    
    if(isVerbose){cat(paste0("Define block ", iBlock, "\n"))}
    trlIdx <- (nTrialPerBlock*(iBlock - 1) + 1):(nTrialPerBlock*iBlock) # indices of trials in this block
    blockRowIdx <- which(data$trialnr %in% trlIdx) # indices of rows in this block
    data$block[blockRowIdx] <- iBlock # set to block ID
    
  } # end for-loop
  if(isVerbose){cat("Finished!\n")}
  
  return(data) # end if variable already exists or not
  
}

# ============================================================================ #
#### 11) Standardize data set per run: ####

applyRunStandardization <- function(data, selVar = "pupil", stdType, isVerbose = FALSE){
  #' Standardize data per run
  #' either using z-standardization, maximum-normalization, or percent signal change.
  #' @param data epoched data, must contain variables \code{block}, \code{ selVar}, and \code{valid}.
  #' @param stdType string specifying transformation type, options: \code{zTransform}, \code{maxNormalization}, \code{percentSignalChange}.
  #' @param isVerbose Boolean, print detailed messages about every trial to console or not (default: FALSE).
  #' @return epoched data, added variables \code{varValid} and \code{varStd}.
  
  cat(paste0("Standardize variable \'", selVar, "\' using ", stdType, "\n"))
  
  # -------------------------------------------------------------------------- #
  ### Initialize new variables:
  
  ## Initialize book-keeping variables:
  varValid <- paste0(selVar, "_valid") 
  varStd <- paste0(selVar, "_std")
  
  ## Initialize variable with pupil samples of only non-rejected trials:
  if (!(varValid %in% colnames(data))){
    data[, varValid] <- data[, selVar] # copy over
    invalidIdx <- which(data$valid == 0) # rows with invalid samples
    data[invalidIdx, varValid] <- NA # delete invalid samples
  }
  
  ## Initialize variable with standardized pupil samples:
  if (!(varStd %in% colnames(data))){
    data[, varStd] <- NA # initialize
  }
  
  ## Initialize block variable:
  if (!("block" %in% names(data))){
    targetNBlock <- 4
    cat(paste0("Variable \'block\' not found--create blocks assuming trials are split into ", targetNBlock, "blocks\n"))
    data <- defineBlocks(data, nBlock = targetNBlock)
  }
  
  nBlock <- length(unique(data$block)) # count blocks:
  validIdx <- which(data$valid == 1) # initialize valid samples:
  
  # -------------------------------------------------------------------------- #
  ### Loop over blocks:
  for (iBlock in 1:nBlock){ # filter separately for each trial: iBlock <- 1
    
    rowIdx <- which(data$block == iBlock)
    if(isVerbose){cat(paste0("Block ", iBlock, ": ", mean(data$valid[rowIdx]*100), "% of trials will be included\n"))}
    
    ## Standardize only based on samples for valid trials:
    if (stdType == "zTransform"){
      
      if(isVerbose){cat(paste0("Block ", iBlock, ": Z-standardize per run\n"))}
      data[rowIdx, varStd] <- scale(data[rowIdx, varValid])
      
    } else if (stdType == "maxNormalization"){
      
      if(isVerbose){cat(paste0("Block ", iBlock, ": Maximum-normalization per run\n"))}
      grandMin <- min(data[rowIdx, varValid], na.rm = T)
      grandMax <- max(data[rowIdx, varValid], na.rm = T)
      data[rowIdx, varStd] <- (data[rowIdx, varValid] - grandMin) / grandMax # standardize into range [0, 1]
      
    } else if (stdType == "percentSignalChange"){
      
      if(isVerbose){cat(paste0("Block ", iBlock, ": Compute percent signal change per run\n"))}
      grandMean <- mean(data[rowIdx,varValid],na.rm = T)
      data[rowIdx, varStd] <- data[rowIdx, varValid]/grandMean*100 - 100
      
    } else {
      
      stop("Unknown standardization method in stdType")
      
    } # end if standardization type 
    
  } # end for-loop over blocks
  if(isVerbose){cat("Finished!\n")}
  
  return(data)
}

# ============================================================================ #
#### 12) Deconfound for certain variable: ####

applyDeconfound <- function(data, DV = "pupil", IV = NULL, regType = "linear", isVerbose = FALSE){
  #' Apply linear or logistic regression to remove variance associated with certain variable.
  #' @param data data frame with trial-epoched data, needs to contain \code{DV} and \code{IV}.
  #' @param DV string argument, dependent variable, default is \code{pupil}.
  #' @param IV string argument, single independent variable, no default.
  #' @param regType string argument, either \code{linear} or \code{logistic}, specifies type of regression.
  #' @param isVerbose Boolean, print detailed messages about every trial to console or not (default: FALSE).
  #' @return data-frame with trial-epoched data, DV only contains residuals after removing influence if \code{IV}.
  
  cat(paste0("Defound ", DV, " from ", IV, " by performing ", regType, " regression and saving residuals\n"))
  
  ## Check presence of IV:
  if(is.null(IV)){
    stop("Forgot to specify IV")
  }
  
  nTrial <- length(unique(data$trialnr)) # count trials
  
  # ------------------------------------------------------------------------- #
  ### Loop over trials:
  
  for (iTrial in 1:nTrial){ # deconfound separately per trial; iTrial = 1
    
    if(isVerbose){cat(paste0("Trial ", iTrial, ": Deconfound for variable", IV, "\n"))}
    rowIdx <- which(data$trialnr == iTrial) # Detect samples belonging to this trial
    y <- data[rowIdx, DV] # extract DV
    x <- data[rowIdx, IV] # extract IV
    
    ## Distinguish regression type:
    if (regType == "linear"){
      
      mod <- lm(y ~ x)
      
    } else if (regType == "logistic"){
      
      mod <- glm(y ~ x, family = "binomial")
      
    } else {
      
      stop("Invalid regType argument")
      
    } 
    
    ## Store residuals back in DV:
    data[rowIdx, DV] <- as.numeric(resid(mod)) # insert back in data frame
    
  } # end for-loop iTrial
  if(isVerbose){cat("Finished!\n")}
  
  return(data)
  
}

# ============================================================================ #
#### 13) Downsample using signal::decimate ####

applyDownsampling <- function(data, varNames, factor, filtType = "iir", deleteEdges = F, lengthEdges = 20, isVerbose = FALSE){
  #' Downsample variables \code{varNames} by factor \code{factor}.
  #' Either applies \code{fir} or \code{iir} filter.
  #' Not possible with internal NAs in varNames.
  #' For method type=\code{fir}, \code{factor} must be a multiple of 4.
  #' @param data data-frame, trial-epoched data with variables \code{trialnr}, \{absTime}, \{trialTime}, and variables \code{varNames}) to-be-downsampled
  #' @param varNames vector of strings, names of variables to-be-downsampled.
  #' @param factor positive integer, factor to downsample by.
  #' @param filtType string argument, either "iir" or "fir".
  #' @param deleteEdges boolean, delete edges after using \code{iir} filter or not.
  #' @param lengthEdges positive integer, length of edges to-be-deleted in downsampled samples, default is 20.
  #' @param isVerbose Boolean, print detailed messages about every trial to console or not (default: FALSE).
  #' @return downsampled, trial-epoched data.
  
  cat(paste0("Downsample data by factor ", factor, "\n"))
  
  require(signal)
  
  nTrial <- length(unique(data$trialnr)) # count trials
  trialList <- list() # initialize list for saving
  
  # ------------------------------------------------------------------------- #
  ### Loop over trials, downsample:
  
  for (iTrial in 1:nTrial){ # downsample separately for each trial: iTrial <- 1
    
    ## Determine start and end of trial:
    rowIdx <- which(data$trialnr == iTrial)
    rowIdxStart <- min(rowIdx)
    rowIdxStop <- max(rowIdx)
    
    ## Determine trial-time related variables: 
    trialnr <- data$trialnr[rowIdxStart] # extract trial number for this trial
    absTime <- seq(from = data$absTime[rowIdxStart], to = data$absTime[rowIdxStop], by = factor)
    trialTime <- seq(from = data$trialTime[rowIdxStart], to = data$trialTime[rowIdxStop], by = factor)
    
    ## Prepare data-frame for variables to-be-downsampled: 
    nDownVar <- length(varNames) # count number variables to downsample
    downVar <- as.data.frame(matrix(NA, length(trialTime), nDownVar)) # initialize data frame
    colnames(downVar) <- varNames # rename
    
    # ----------------------------------------------------------------------- #
    ## Loop over variables, extract variables, downsample:
    
    for (iVar in 1:nDownVar){
      
      ## Pad edges:
      data <- applyPadEdges(data, selVar = varNames[iVar], rowIdx = rowIdx, iTrial = iTrial, isInteractive = FALSE)
      
      ## Extract variable of interest:
      if(isVerbose){cat(paste0("Trial ", iTrial, ": Downsample variable ", varNames[iVar], " by factor ", factor, "\n"))}
      xLong <- data[rowIdx, varNames[iVar]] # extract old version of vector
      xShort <- decimate(x = xLong, q = factor, ftype = filtType) # downsample
      
      ## Delete edges:
      if (deleteEdges == T){
        if(isVerbose){cat(paste0("Trial ", iTrial, ": Cut edges of variable ", varNames[iVar], " by ", lengthEdges, " samples\n"))}
        xShort[1:(1 + lengthEdges)] <- NA # delete leading edge
        xShort[(length(xShort) - lengthEdges):length(xShort)] <- NA # delete trailing edge
      }
      
      # plot(xShort)
      downVar[, iVar] <- xShort # return into data-frame
      
    } # end for-loop iVar
    
    # ------------------------------------------------------------------------- #
    ## Extract values for other variables:  
    
    restVarNames <- colnames(data) # copy over names
    varIdx <- which(!(restVarNames %in% c("trialnr", "absTime", "trialTime", varNames))) # which ones NOT among time variable and varNames
    restVars <- data[rowIdxStart, varIdx] # just extract value of first sample
    trialList[[iTrial]] <- as.data.frame(cbind(trialnr, absTime, trialTime, downVar, restVars)) # concatenate into data frame
    
  } # end for-loop iTrial
  if(isVerbose){cat("Finished!\n")}
  
  
  ## Append all trials below each other:
  downData <- do.call(rbind, trialList) # concatenate all trials
  return(downData) # return
  
}

# ============================================================================ #
#### 14) Compute baseline and dilation (aggregate per trial) using plyr ####

applyAggregation <- function(data, selVar = "pupil_std", baselineLim = c(-500, 0), trialLim = c(0, 2966), iSub, 
                             extraVarNames = NULL){
  #' Aggregate data per trial (using ddply from plyr)
  #' Get baseline as mean of window within baselineWindow (positive number) till 0
  #' Get dilation as maximum within trialWindow (positive number)
  #' @param data Epoched-data, needs to have all variables to-be-aggregated (so also task factors)
  #' @param selVar string argument, variable to use for aggregation
  #' @param baselineLim vector of two integers, start and end of baseline window in ms to average over (default: -500 before segMessage till 0 i.e. at segMessage).
  #' @param trialLim vector of two integers, start and end of trial window in ms within which to find maximum dilation (default: 0 i.e. at segMessage till 2,966 ms after segMessage).
  #' @param iSub Subject number to-be-added to data frame
  #' @param extraVarNames specify vector of variables names for which first entry per trial will be retained
  #' @return Aggregated data frame with one row per trial
  
  require(plyr)
  
  cat("Aggregate data on trial-by-trial level, compute baseline pupil and pupil dilation\n")
  
  ## Initialize vectors where to extract baseline and TEPR:
  baselineWindow <- seq(baselineLim[1], baselineLim[2], 1) # vector from baselineStart till 0
  trialWindow <- seq(trialLim[1], trialLim[2], 1) # vector from 0 till trialEnd
  
  aggrData <- ddply(data, .(trialnr), function(x){
    
    ## Copy over trial-identifying variables:
    subject <- iSub
    trialnr <- x$trialnr[1]
    
    ### Extract first entry for optional variables:
    
    if (!(is.null(extraVarNames))){
      extraVars <- rep(NA, length(extraVarNames)) # initialize empty output vector
      for (iVar in extraVarNames){
        extraVars[iVar] <- x[1, extraVarNames[iVar]]
      }
    }
    
    ## Extract baseline:
    baseline <- mean(x[which(x$trialTime %in% baselineWindow), selVar], na.rm = T)
    
    ## Extract maximal pupil dilation:
    pupilmax <- max(x[which(x$trialTime %in% trialWindow), selVar], na.rm = T)
    
    ## Compute dilation corrected for baseline:
    dilation <- (pupilmax - baseline) 
    
    if (!(is.null(extraVarNames))){
      
      output <- data.frame(subject, trialnr, extraVars, baseline, dilation)
      names(output) <- c("subject", "trialnr", extraVarNames, "baseline", "dilation")
      
    } else {
      
      output <- data.frame(subject, trialnr, baseline, dilation)
      
    }
    return(output)
    dev.off()})
  
  ## Set infinite values to NA:
  aggrData$baseline[!is.finite(aggrData$baseline)] <- NA
  aggrData$dilation[!is.finite(aggrData$dilation)] <- NA
  
  ## Report number of trials with NAs:
  cat(paste0("Number of trials with baseline == NA: ", sum(is.na(aggrData$baseline)), "\n"))
  cat(paste0("Number of trials with dilation == NA: ", sum(is.na(aggrData$dilation)), "\n"))
  
  return(aggrData)
}

# END OF FILE.
