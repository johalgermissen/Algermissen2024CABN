#!/usr/bin/env Rscript
# ============================================================================ #
## 00_mgngunc_functions_preprocess_gaze.R
## MGNGUncertainty functions for pre-processing of gaze data.
# Johannes Algermissen, 2023.

# ============================================================================ #
#### Read asc file and return size ####

readEyelink <- function(fileName){
  #' Just a wrapper for function \code{read.asc} from package \code{eyelinker},
  #' so everything is together
  #' @param fileName string argument, file name ending on \code{.asc}
  #' @return object read in with \code{read.asc}, has fields
  #'  \code{raw}, \code{msg}, \code{sacc}, \code{fix}, \code{blinks}, \code{info}
  
  require(eyelinker)
  
  cat(paste0("Read file ",fileName, " of size ",file.size(fileName)/1000000, " MB\n"))
  dat <- read.asc(fileName, parse_all = T)
  cat(paste0("Output size: ", object.size(dat)/1000000, " MB\n"))
  
  return(dat)
}

# ============================================================================ #
#### Segment data into trials, read additional trial-by-trial variables ####

readTrials <- function(dat, segmentMessage, beforeEvent, afterEvent, 
                       variablesMessages = NULL, recode2Num = NULL, variablesNames = NULL){
  #' Epoch data into trials based on occurrence of \code{segmentMessage} as time point zero, 
  #' with \code{beforeEvent} as positive integer of samples before \code{segmentMessage} 
  #' and \code{afterEvent} as positive integer of samples after \code{segmentMessage};
  #' can also read other variables if written in eyelink file
  #' returns epoched data with all trials features retrieved
  #' @param dat raw data read-in with read.asc from package \code{Eyelinker}
  #' @param segmentMessage string, message to be used for segmenting data into trials, used to set time point 0
  #' @param beforeEvent positive integer, number of samples to be included before \code{segmentMessage}
  #' @param afterEvent positive integer, number of samples to be included after \code{segmentMessage}
  #' @param variablesMessages vector of strings, contains messages for which values will be retrieved (assumes value one row below message)
  #' @param recode2Num vector of booleans, specify for each variable whether to recode or not
  #' @param variablesNames vector of strings, new variable names for extracted variables in output object; default is NULL, which means \code{variablesMessage} will be used as variable names
  #' @return data frame with epoched data
  #' @examples 
  #' dat <- read.asc("s_01.asc")
  #' data <- readTrials(dat, segmentMessage="StartMask", beforeEvent=1000, afterEvent=3000)
  #' (Consider letting messages for other data to-be-retrieved be specified in input??)
  
  ## Check whether variablesMessages, recode2Num, and variablesNames are of same length:
  if (!(is.null(variablesMessages))){
    if(length(variablesMessages)!=length(recode2Num)){
      stop("Arguments variablesMesssage and recode2Num are of different lengths")
    }
    if(is.null(variablesNames)){
      variablesNames <- variablesMessages
      cat("Argument variablesNames unspecified, will use variablesMessages as variable names\n")
    }
    if(length(variablesMessages)!=length(variablesNames)){
      stop("Arguments variablesMesssage and variablesNames are of different lengths")
    }
  } # end of if-clause if variablesMessages exist
  
  ## Extract raw data:
  raw <- dat$raw
  
  ## Determine time of start:
  msgIdx <- which(dat$msg$text == segmentMessage) # indices where message for segmentation presented
  msgTime <- dat$msg$time[msgIdx] # timings of all trials at those indices 
  nTrials <- length(msgTime) # number trials
  startTime <- msgTime - beforeEvent # start segmentation xx ms before message
  stopTime <- msgTime + afterEvent # stop segmentation xx ms after message
  
  ## Determine trial details:
  if (!(is.null(variablesMessages))){
    
    variablesList <- list() # initialize empty list
    nVariables <- length(variablesMessages) # length to loop through
    
    for (iVariable in 1:nVariables){
      
      # Read variable values (one row below variable names appears):
      thisVariable <- dat$msg$text[which(dat$msg$text == variablesMessages[iVariable])+1]
      
      if (recode2Num[iVariable] == T){ # if input says convert to numeric: recode into numeric
        thisVariable <- as.numeric(levels(thisVariable))[thisVariable]
        # https://stackoverflow.com/questions/3418128/how-to-convert-a-factor-to-integer-numeric-without-loss-of-information
        
      } else { # else: assume factor, only remove unused levels
        thisVariable <- droplevels(factor(thisVariable))
        
      } # end if-clause for recoding to numeric
      
      variablesList[[iVariable]] <- thisVariable # add to list
      
    } # end for-loop iVariable
  } # end if-clause variablesMessages are not empty
  
  ## Initialize data frame:
  trialList <- list()
  
  ## Loop over trials:
  for (iTrial in 1:nTrials){ # iTrial = 1
    cat(paste0("Read trial ", iTrial, "\n"))
    
    # Retrieve start and stop time for this trial:
    startRow <- which(raw$time == startTime[iTrial])
    stopRow <- which(raw$time == stopTime[iTrial])
    
    # Check whether computed rows exist, otherwise warn:
    if (length(startRow) == 0){
      cat(paste0("Trial ", iTrial, ": startRow undefined, consider shortening beforeEvent\n"))
    }
    if (length(stopRow) == 0){
      cat(paste0("Trial ", iTrial, ": stopRow undefined, consider shortening afterEvent\n"))
    }
    
    ## Extract eye-tracking measure for this trial:
    xp <- raw$xp[startRow:stopRow]
    yp <- raw$yp[startRow:stopRow]
    pupil <- raw$ps[startRow:stopRow]
    
    ## Trial number and timing:
    trialDuration <- length(pupil) # number of pupil samples determines trial duration
    trialnr <- rep(iTrial, trialDuration)
    absTime <- raw$time[startRow:stopRow] # absolute timing relative to start of eye-tracker
    trialTime <- seq(from = -1*beforeEvent, to = (trialDuration - beforeEvent - 1), by = 1) # relative time within trial
    
    ## Extract trial features from messages:
    if (!(is.null(variablesMessages))){
      
      allVariables <- as.data.frame(matrix(NA, trialDuration, nVariables)) # initialize empty data frame
      
      for (iVariable in 1:nVariables){
        variableAllValues <- variablesList[[iVariable]] # extract this variable
        variableValue <- variableAllValues[iTrial] # extract value of variable for this trial 
        allVariables[, iVariable] <- rep(variableValue, trialDuration) # repeat trial feature as often as samples within trial:
      } # end for-loop iVariable
      
      colnames(allVariables) <- variablesNames # add variable names given as inputs
      
      # Concatenate eye-tracking variables and trial features (allVariables) into data frame:
      trialList[[iTrial]] <- as.data.frame(cbind(trialnr, absTime, trialTime, xp, yp, pupil, allVariables))
      
    } else {
      
      # Concatenate into data frame without allVariables:
      trialList[[iTrial]] <- as.data.frame(cbind(trialnr, absTime, trialTime, xp, yp, pupil))
      
    } # end if-clause variablesMessages are not empty
  } # end for-loop iTrial
  
  ## Append data frames of all trials:
  data <- do.call(rbind, trialList)
  
  ## Warn about warning messages that can occur when converting factors to numeric:
  if (!(is.null(recode2Num)) & sum(recode2Num) > 0){ # if recode2Num defined at at least one True
    cat(paste0("Expect warning about 'NAs introduced by coercion' for every variable converted to numeric, i.e. ", sum(recode2Num), " warnings\n"))
  } # end if recode2Num > 0
  
  cat("Done :-)\n")
  return(data)
} # end function

# ============================================================================ #
#### Recode factors to numeric ####

recode_behavioral_variables <- function(data){
  #' Recode factors into numeric
  #' @param data data frame, trial-epoched data with variable \code{variable}
  #' @return data data frame with recoded variables
  
  ## Define simple function:
  as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}

  ## Factor to numeric:
  cat("Convert factors to numeric\n")
  data$Block <- as.numeric.factor(data$Block)
  data$TrialNr <- as.numeric.factor(data$TrialNr)
  data$Trialnr <- data$TrialNr
  data$Condition <- as.numeric.factor(data$Condition)
  data$Valence <- as.numeric.factor(data$Valence)
  data$`Required Action` <- as.numeric.factor(data$`Required Action`)
  data$Manipulation <- as.numeric.factor(data$Manipulation)
  data$GoValidity <- as.numeric.factor(data$GoValidity)
  data$NoGoValidity <- as.numeric.factor(data$NoGoValidity)
  data$ISI <- as.numeric.factor(data$ISI)
  data$ITI <- as.numeric.factor(data$ITI)
  
  ## Rename required action:
  cat("Rename required action variable\n")
  data$ReqAction <- data$`Required Action`
  
  ## Combine validities:
  cat("Combine Go and NoGo validities\n")
  data$Validity <- ifelse(data$Response == 1, data$GoValidity, data$NoGoValidity)
  
  cat("Finished :-)\n")
  
  return(data)
}

# ============================================================================ #
#### Determine amount of missingness ####

# data1 <- subBehavData; data2 <- subEyeData;

compare_data_sets <- function(data1, data2, varVec = NULL){
  
  #' Determine whether 2 data sets are identical with regard to selected variables
  #' @param data1 data frame, first data set, use only rows where data1 is not NA
  #' @param data2 data frame, second data set
  #' @param varVec vector of strings, variables to compare (default: all of data1)
  #' @return TRUE (all variables the same) or FALSE (at least one variable not the same)
  
  if(is.null(varVec)){
    cat("No selection of variables provided, use all variables of first data set\n")
    varVec <- names(data1)
  }
  
  cat("Variables data set 1:\n", names(data1), "\n")
  cat("Variables data set 2:\n", names(data2), "\n")
  
  allSame <- T
  
  for (iVar in varVec){ # iVar <- 1
    valIdx <- !(is.na(data1[, iVar]))
    mSame <- mean(data1[valIdx, iVar] == data2[valIdx, iVar])
    if (is.na(mSame)){
      cat("Variable", iVar, ": NAs in data2 detected\n")
      allSame <- F
    } else if (mSame == 1){
      cat("Variable", iVar, ": data sets correspond\n")
    } else {
      misMatchIdx <- which(data1[, iVar] != data2[, iVar])
      cat("Variable", iVar, ": data sets mismatch in rows", misMatchIdx, "\n")
      allSame <- F
    }
  }
  return(allSame)
}


# ============================================================================ #
#### Determine amount of missingness: ####

count_missingness <- function(data, variable){
  #' Determine how much (in absolute and relative terms) of variable \code{variable} is NA.
  #' @param data data frame, trial-epoched data with variable \code{variable}
  #' @param variable string argument, variable to be interpolated
  #' just print how much is missing, no object returned
  
  nNA <- sum(is.na(trialData$xp))
  nRow <- nrow(data)
  cat(paste0("Variable ", variable, ": ", nNA, " out of ", nRow, " samples (", round(nNA/nRow, 2), "%) missing\n"))
}

# ============================================================================ #
#### 2) Delete blinks ####

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
#### 3) Compute derivative, delete strong changes: ####

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
    
    ## Detect gaps < cluster_cutoff between marked samples, also delete those:
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
#### 4) Linearly interpolate all NAs: #####

# data <- trialData; variable <- "xp"; maxGap = 100; maxKeepTrial <- NULL; padEdges = TRUE; interactive = TRUE;
applyInterpolation <- function(data, variable = "xp", maxGap = NULL, maxKeepTrial = NULL, padEdges = FALSE, interactive = FALSE){
  #' Performs linear interpolation using na.approx from package zoo
  #' @param data data frame, trial-epoched data with variable \code{variable}
  #' @param variable string argument, variable to be interpolated
  #' @param maxGap numeric, maximum length of NA sequences to be interpolated (in samples)
  #' @param maxKeepTrial numeric, maximum number of NAs to still interpolate trial, otherwise set entire trial to NA (in samples or as percentage)
  #' @param padEdges boolean, fill leading and trailing NAs of trials with first/ valid non-NA
  #' @param interactive boolean, plot each trial before/ after interpolation (after padding edges) or not
  #' @return data frame, variable \code{variable} manipulated
  require(zoo) # for na.approx
  
  # Defaults:
  timeVar <- "trialTime"
  trialVar <- "trialnr"
  nTrials <- length(unique(data[, trialVar]))
  
  NATrials <- 0
  
  ## Loop over trials:
  for (iTrial in 1:nTrials){ # iTrial <- 253
    cat(paste0("Trial ", iTrial, ": Interpolate NAs\n"))
    rowIdx <- which(data[, trialVar] == iTrial)
    
    ## Trial length and number of missing samples (in samples):    
    nSamples <- length(data[rowIdx, variable]) # length of trial in samples
    nMissing <- sum(is.na(data[rowIdx, variable])) # number of missing samples in this trial
    
    ## Set maxKeepTrial and maxGap:
    if (is.null(maxKeepTrial)){ # if maxKeepTrial not set: set to total length of trial
      maxKeepTrial <- nSamples
    }
    if (maxKeepTrial < 1){ # if given as percentage:
      maxKeepTrial <- maxKeepTrial*nSamples
    }
    
    if (is.null(maxGap)){ # if maxGap not set: set to total length of trial
      maxGap <- nSamples
    }
    if (maxGap < 1){ # if given as percentage:
      maxGap <- maxGap*nSamples
    }
    
    ## Check if any data non-NA:
    if(nMissing >= maxKeepTrial){
      
      data[rowIdx, variable] <- NA # set entire trial to NA ()
      cat(paste0("Trial ", iTrial, ": ", nMissing, " out of ", nSamples, " samples are NA, no interpolation\n"))
      NATrials <- NATrials + 1
      
    } else { # if any non-NA data points:
      
      ## Plot before interpolation:
      if (interactive == TRUE){
        plot(data[rowIdx, timeVar], data[rowIdx, variable], type = "l", 
             ylim = c(scrx/2 - tolFix*plotFactor, scrx/2 + tolFix*plotFactor),
             main = paste0("Trial ", iTrial, " - before interpolation"))
        readline(prompt = "Press [enter] to continue")
      }
      
      # isNA <- is.na(data[rowIdx, variable])
      # which(isNA)
      # diff(which(isNA))
      
      ## Apply interpolation:
      data[rowIdx, variable] <- na.approx(data[rowIdx, variable], maxgap = maxGap, na.rm = FALSE) # replace NAs with interpolated value
      
      ## Plot after interpolation:
      if (interactive == TRUE){
        plot(data[rowIdx, timeVar], data[rowIdx, variable], type = "l", col = "green",
             ylim = c(scrx/2-tolFix*plotFactor, scrx/2+tolFix*plotFactor), 
             main = paste0("Trial ", iTrial, " - after interpolation"))
        readline(prompt="Press [enter] to continue")
      }
      
      ## Interpolate also from start and until end of window: Pad edges:
      if (padEdges == TRUE){
        data <- applyPadEdges(data, variable, rowIdx, iTrial, interactive)
      } # end if for padding edges
      
    } # end if loop
    
  } # end for-loop iTrial
  
  cat("Found", nATrials, "trials without any data\n")
  
  return(data)
}

# ============================================================================ #
#### Interpolate leading and trailing NAs at beginning & end of trial ####

applyPadEdges <- function(data, variable = "xp", rowIdx, iTrial, interactive = F){
  #' Replaces leading NAs by first non-NA value,
  #' Replaces trailing NAs by last non-NA value.
  #' @param data data frame, trial-epoched data with variable \code{variable}
  #' @param variable string argument, variable to be padded
  #' @param rowIdx vector of integers, vector of rows of trial to be padded
  #' @param iTrial integer, trial number (for title in plot)
  #' @param interactive boolean, plot trial time series after padding edges (in blue) or not
  #' @return no return, just plotting
  #' Should only be used if enough empty time before baseline/ after trial 
  require(zoo)
  cat(paste0("Trial ", iTrial, ": Pad edges of variable ", variable, "\n"))
  
  timeVar <- "trialTime"
  trialVar <- "trialnr"
  
  # First interpolate:
  data[rowIdx, variable] <- na.approx(data[rowIdx, variable], na.rm = FALSE) # replace NAs with interpoloated value
  
  # Determine first and last sample:
  firstRowIdx <- head(rowIdx, 1)
  lastRowIdx <- tail(rowIdx, 1)
  
  # Determine non-NA samples:
  filledIdx <- which(!(is.na(data[rowIdx, variable]))) # non-NAs
  firstfilled <- rowIdx[head(filledIdx, 1)] # first non-NA
  lastfilled <- rowIdx[tail(filledIdx, 1)] # last non-NA
  
  # Replace leading NAs:
  if(firstfilled!=firstRowIdx){
    data[firstRowIdx:(firstfilled-1), variable] <- data[firstfilled, variable]
  }
  
  # Replace trailing NAs:
  if(lastfilled!=lastRowIdx){
    data[(lastfilled+1):lastRowIdx, variable] <- data[lastfilled, variable]
  }
  
  # Plot after removal of NAs:
  if (interactive == TRUE){
    plot(data[rowIdx, timeVar], data[rowIdx, variable], type = "l", col = "blue", 
         ylim = c(scrx/2-tolFix*plotFactor, scrx/2+tolFix*plotFactor),
         main = paste0("Trial ", iTrial, " - after padding edges"))
    readline(prompt="Press [enter] to continue")
  }
  return(data)
}

# ============================================================================ #
#### Plot one coordinate over time with ROIs: ####

plot_coordinate_ROI <- function(data, xVar = "trialTime", yVar = "xp",
                                leftROI, rightROI, radius){
  #' Line plot of variable (over time) per trial, highlight ROIs.
  #' @param data data frame, trial-epoched data with variable \code{variable} and 'trialnr'
  #' @param xVar variable on x-axis, should be e.g. time within trial
  #' @param yVar variable on y-axis, should be variable of interest
  #' @param leftROI relevant coordinate of center of left ROI
  #' @param rightROI relevant coordinate of center of right ROI
  #' @param radius radius of both ROIs
  #' No output, just plotting.
  
  # Determines edges of ROIs:
  leftROImin <- leftROI - radius
  leftROImax <- leftROI + radius
  rightROImin <- rightROI - radius
  rightROImax <- rightROI + radius
  
  # Defaults:
  timeVar <- xVar
  trialVar <- "trialnr"
  nTrial <- length(unique(data[, trialVar]))
  RtVar <- "RT"
  
  # Determine average RT:
  avgRT <- mean(data[which(data[, timeVar] == 1), rtVar], na.rm = T) # retrieve only first sample, mean with NAs
  
  # Loop over trials:
  for (iTrial in 1:nTrial){ # iTrial <- 1
    rowIdx <- which(data[, trialVar] == iTrial) # row indices of this trial
    
    # Time variable:
    timeLine <- data[rowIdx, xVar] # extract time variable for this trial
    timeMin <- timeLine[1] # first entry
    timeMax <- timeLine[length(timeLine)] # last entry
    
    # Check whether trial is Go trial:
    isGo <- data[rowIdx[1], "response"]
    
    # Add RT:
    if (isGo == 1){ # If trial is Go trial: retrieve RT
      RT <- data[rowIdx[1], rtVar]
      RTcol <- "green" # green for Go
      # RTs recorded relevant to StartCueOutcome (when outcomes appear), so same timing as recoded trial time
    } else { # else if NoGo trial: take average RT
      RT <- avgRT
      RTcol <- "red" # red for NoGo
    }# end if
    
    # Make line plot:
    plot(timeLine, data[rowIdx, yVar], type = "l", lwd = 2,
         cex.main = 2, cex.lab = 1.5, cex.axis = 1.25, 
         ylim = c(leftROImin, rightROImax),
         xlab = "Time (in ms)", ylab = "X coordinate (in pixels)", main = paste0("Trial ", iTrial))
    
    # Vertical line for RT;
    abline(v = RT*1000,col = RTcol, lwd = 2, lty = 2)
    
    # Make polygons for ROIs:
    xx <- c(timeLine[1], timeLine[length(timeLine)], timeLine[length(timeLine)], timeLine[1], timeLine[1])
    yy <- c(leftROImin, leftROImin, leftROImax, leftROImax, leftROImin) # left polygon
    polygon(xx, yy, density = 0, col = "#FF8000") # left ROI
    yy <- c(rightROImin, rightROImin, rightROImax, rightROImax, rightROImin) # right polygon
    polygon(xx, yy, density = 0, col = "#0000FF") # right ROI
    
    # Press to continue:
    readline(prompt = "Press [enter] to continue")
  }
}

# -------------------------------------------------------------------------------------------------- #
#### Identify which samples are before response: ####

mark_before_RT <- function(data, tolerance = 0){
  #' Mark for trials with Go responses whether sample occurrred \code{tolerance}
  #' seconds before response (1) or not (0). NA for trials with NoGo responses.
  #' @param data data frame, trial-epoched data 
  #' @param tolerance numeric, how many milliseconds before RT to stop (default: 0)
  #' @return data with new variable "isBeforeRT" indicating whether sample is before RT or not
  
  # data <- trialData
  
  require(stringr)
  
  # Defaults:
  timeVar <- "trialTime"
  trialVar <- "trialnr"
  RtVar <- "RT"
  nTrial <- length(unique(data[, trialVar]))
  
  newVarName <- paste0("isBeforeRT_tol", str_pad(tolerance, 3, pad = "0"))
  
  # Initialize variable:
  data[, newVarName] <- NA
  
  # Determine average RT:
  avgRT <- mean(data[which(data[, timeVar] == 1), rtVar], na.rm = T) # retrieve only first sample, mean with NAs
  
  # Loop over trials:
  for (iTrial in 1:nTrial){ # iTrial <- 1
    
    cat(paste0("Trial ", iTrial, ": Mark samples before (RT - ", tolerance, " ms)\n"))
    
    # Retrieve row indices of trial:
    rowIdx <- which(data[, trialVar] == iTrial) # row indices of this trial
    # length(rowIdx)
    
    # --------------------------------------------------------------- #
    ### Retrieve RT:
    
    # Check whether trial is Go trial:
    isGo <- data[rowIdx[1], "response"]
    
    if (isGo == 1){ # If trial is Go trial: retrieve RT
      RT <- data[rowIdx[1], rtVar]
      # RTs recorded relevant to StartCueOutcome (when outcomes appear), so same timing as recoded trial time
    } else { # else if NoGo trial: take average RT
      RT <- avgRT
    }# end if
    
    # Detect samples before RT:
    idx <- which(data[rowIdx, timeVar]*.001 <= RT - tolerance*.001) # indices relative to trial onset
    # length(idx)
    beforeRTidx <- rowIdx[idx] # retrieve indices absolute to entire data frame
    # length(beforeRTidx)
    data[beforeRTidx, newVarName] <- 1
    
    # Detect samples after RT:
    idx <- which(!(rowIdx %in% beforeRTidx))
    afterRTidx <- rowIdx[idx] # retrieve indices absolute to entire data frame
    # length(afterRTidx)
    data[afterRTidx, newVarName] <- 0
    
    # View(data[rowIdx,c("trialnr", "trialTime", "response", "RT", newVarName)])
    
  } # end iTrial
  # View(data[,c("trialnr", "trialTime", "response", "RT", newVarName)])
  
  # Statistics:
  # cat("After classifying samples into before/after RT:\n")
  # cat(table(data[, newVarName]))
  
  return(data)
}

# -------------------------------------------------------------------------------------------------- #
#### Plot average time course per condition: ####

plot_avg_cond <- function(data1, data2, yVar="xp", timeVar="trialTime", 
                          yLabel = "X coordinate (in pixels)", mainLabel = "condition",
                          colorVec=c("green", "red"), legendVec=c("Cond 1", "Cond 2"),
                          isSave=FALSE, splitType="condition", suffix=""){
  #' Aggregates time course within a trial separately for two conditions, plots them
  #' @param data1 data frame, trial-epoched data, first condition 
  #' @param data2 data frame, trial-epoched data, second condition 
  #' @param yVar string, dependent variable to plot (default: "xp")
  #' @param timeVar string, variable containing time within trial (default: "trialTime")
  #' @param yLabel string, y-axis label (default: "X coordinate (in pixels)")
  #' @param mainLabel string, input for caption specifying conditions (default: "condition")
  #' @param colorVec vector of 2 strings, containing colors to use (default: green and red)
  #' @param legendVec vector of 2 strings, condition names for legend (default: cond1 and cond2) 
  #' @return data with new variable "isBeforeRT" indicating whether sample is before RT or not
  
  # Assume that plotDir yLim and and subID are set globally
  
  # Defaults:
  RtVar <- "RT"
  
  # Time variable:
  xVar <- unique(data1[, timeVar])
  
  # Aggregate relevant variable across trials:
  yVar1 <- tapply(data1[, yVar], data1[, timeVar], mean, na.rm=T) # average per sample within trial
  yVar2 <- tapply(data2[, yVar], data2[, timeVar], mean, na.rm=T) # average per sample within trial
  
  # Start saving:
  if(isSave == T){
    png(paste0(plotDir, "MGNGFreeView_xp_trial_", splitType, suffix, "_sub", str_pad(subID, 2, "left", "0"), ".png"),width = 480, height = 480)
  }
  
  # Plot average time courses:
  plot(xVar, yVar1, "l",col=colorVec[1], lwd = 2,
       cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.25, 
       xlab = "Time (in ms)", ylab = yLabel, ylim = yLim,
       main = paste0("Sub ", subID, ": Average per ", mainLabel))
  lines(xVar, yVar2, "l",col=colorVec[2], lwd = 2)
  
  # Plot middle of screen:
  abline(h=scrx/2, lwd=1.5, lty=2)
  
  # Add legend:
  legend("topleft", legend=legendVec,
         col=colorVec, lty=1, cex=1)
  
  # End saving:
  if(isSave == T){
    dev.off()
  }
  
}

# -------------------------------------------------------------------------------------------------- #
#### Identify which samples are before response: ####

reepoch_resplocked <- function(data, variables, startTime = -500, stopTime = 500){
  #' Mark for trials with Go responses whether sample occurrred \code{tolerance}
  #' seconds before response (1) or not (0). NA for trials with NoGo responses.
  #' @param data data frame, trial-epoched data 
  #' @param variables vector of string, variables to change (other variables will keep value at timeVar=)
  #' @param startTime numeric, how many ms before RT to begin reepoched trials (fill in NA if necessary) 
  #' @param endTime numeric, how many ms after RT to end reepoched trials  (fill in NA if necessary)
  #' @return newData re-epoched to RT
  
  # data <- trialData; 
  # variables <- c("xp", "yp", "pupil", "dist_left", "fix_left", "dist_right", "fix_right", "isBeforeRT_tol000")
  
  
  # Defaults:
  timeVar <- "trialTime"
  maxDur <- max(data[, timeVar])
  trialVar <- "trialnr"
  nTrial <- length(unique(data[, trialVar]))
  RtVar <- "RT"
  nVar <- length(variables)
  trialDuration <- stopTime - startTime
  
  # Determine other(variables where to just keep value at first sample):
  allVariables <- names(data)
  otherVariables <- allVariables[!(allVariables %in% variables)]
  
  # Determine average RT:
  avgRT <- mean(data[which(data[, timeVar] == 1), rtVar], na.rm = T) # retrieve only first sample, mean with NAs
  
  # New trial time:
  trialTime <- startTime:(stopTime-1) # account for crossing 0
  
  ## Initialize new data frame:
  trialList <- list()
  
  # ----------------------------------------------------------------- #
  # Loop over trials:
  for (iTrial in 1:nTrial){ # iTrial <- 209
    cat(paste0("Trial ", iTrial, ": Relock to RT\n"))
    
    # Retrieve row indices of trial:
    rowIdx <- which(data[, trialVar] == iTrial) # row indices of this trial
    # length(rowIdx)
    
    # --------------------------------------------------------------- #
    ### Retrieve RT:
    
    # Check whether trial is Go trial:
    isGo <- data[rowIdx[1], "response"]
    
    if (isGo == 1){ # If trial is Go trial: retrieve RT
      RT <- data[rowIdx[1], rtVar]
      # RTs recorded relevant to StartCueOutcome (when outcomes appear), so same timing as recoded trial time
    } else { # else if NoGo trial: take average RT
      RT <- avgRT
    }# end if
    
    RT <- round(RT,3) # round to ms precision
    
    # --------------------------------------------------------------- #
    ### Determine when trial starts/ stops in old and new data:
    
    ## Samples to extract from old data:
    startSample <- RT*1000+startTime # where to start this re-epoched trial
    stopSample <- RT*1000+stopTime-1 # where to stop this re-epoched trial
    # make stopSample one sample shorter become 0 included
    
    ## Samples to fill in in new data:
    newStartSample <- 1
    newStopSample <- trialDuration
    
    ## Correct if necessary:
    if (startSample < 1){
      newStartSample <- newStartSample - (startSample-1)  # set new start sample further forward
      startSample <- 1 # reset old start sample to minimally zero
      cat(paste0("Old start sample of ", startSample, " too early, reset to 1 \n"))
    }
    if (stopSample > maxDur){
      cat("Old start sample too late, reset \n")
      newStopSample <- newStopSample - (stopSample-maxDur) # set new start sample further back
      stopSample <- maxDur # reset old start sample to minimally zero
      cat(paste0("Old stop sample of ", stopSample, " too late, reset to", maxDur, "\n"))
    }
    
    # length(startSample:stopSample)
    # length(newStartSample:newStopSample)
    # --------------------------------------------------------------- #
    # Extract selected variables in selected time frame:
    varFrame <- data.frame(matrix(nrow = trialDuration, ncol = nVar)) # initialize empty data frame
    
    for(iVar in 1:nVar){ # iVar <- 1
      varName <- variables[iVar] # retrieve variable name
      varFrame[newStartSample:newStopSample, iVar]  <- data[rowIdx[startSample:stopSample], varName]
    }
    
    # --------------------------------------------------------------- #
    # Extract only first sample for other variables:
    
    otherValues <- data[rowIdx[1],otherVariables]
    
    # --------------------------------------------------------------- #
    # Put into one data frame per trial in list:
    
    trialList[[iTrial]] <- as.data.frame(cbind(iTrial, trialTime, varFrame,otherValues))
    
  }
  
  ## Append data frames of all trials:
  newData <- do.call(rbind, trialList)
  
  ## Add names back in:
  names(newData) <- c(trialVar, timeVar, variables,otherVariables)
  
  # head(newData)
  # str(newData)
  
  return(newData)
}

# -------------------------------------------------------------------------------------------------- #
#### Identify adjacent fixations within ROI as clusters, give index ####

identify_fixations <- function(data, variables, interactive = F){
  
  #' Loop through trials, identify clusters of adjacent fixations within an ROI, give indices to clusters.
  #' @param data data frame, trial-epoched data with variables \code{variables}
  #' @param variables vector of variable names for which to identify clusters (each variable needs to have values 1 and 0)
  #' @param interactive whether to print specifics on fixation to console or not
  #' @return data with new variable "_fixIdx" for each input variable indicating index of each cluster
  
  # Defaults:
  nVariables <- length(variables) # length to loop through
  trialVar <- "trialnr"
  nTrials <- length(unique(data[, trialVar]))
  
  for (iVariable in 1:nVariables){ # for each variable
    
    thisVariable <- variables[iVariable] # retrieve variable
    cat(paste0("Start for variable ", thisVariable), "\n")
    
    newVariable <- paste0(thisVariable, "_fixIdx") # create name for fixation index
    data[, newVariable] <- NA # initialize new variable
    
    for (iTrial in 1:nTrials){ # iTrial <- 1
      cat(paste0("Trial ", iTrial, ": Identify fixations\n"))
      rowIdx <- which(data[, trialVar] == iTrial)
      nSamples <- length(rowIdx)
      
      ## Initialize:
      sampleCount <- 0
      fixCount <- 0
      memory <- 0
      
      ## Loop through samples:
      for (iSample in 1:nSamples){
        sampleIdx <- rowIdx[iSample]
        
        ## Retrieve whether fixation or not:
        sampleIsFix <- data[sampleIdx, thisVariable] # retrieve whether fixation or not
        
        ## Increment fixation count at start of new fixation:
        if(sampleIsFix == 1 & memory == 0){
          fixCount <- fixCount + 1
        }
        
        ## Add fixation index to newVariable:
        data[sampleIdx, newVariable] <- ifelse(sampleIsFix == 1,fixCount, nA)
        
        ## Print end of fixation:
        if(sampleIsFix == 0 & memory == 1 & interactive == TRUE){
          cat(paste0("Trial ", iTrial, ": Finished fixation no. ", fixCount, " with ", length(which(data[rowIdx, newVariable] == fixCount)), " samples\n"))
        }
        
        ## Update memory:
        memory <- sampleIsFix
        
      } # end iSample 
      if (interactive == TRUE){
        cat(paste0("Trial ", iTrial, ": Detected ", fixCount, " fixations\n"))
      }
      
    } # end iTrial
    
  } # end iVariable
  
  return(data)
}

# -------------------------------------------------------------------------------------------------- #
#### Discard too short fixations: ####

remove_short_fixations <- function(data, isFixVar, numFixVar, threshold){
  #' Loop through trials, remove identified clusters shorter than threshold.
  #' @param data data frame, trial-epoched data with variables \code{isFixVar} and \code{numFixVar}
  #' @param isFixVar contains indicator whether fixation is within ROI (1) or not (0)
  #' @param numFixVar contains indices attached to each cluster of adjacent fixations of ROI
  #' @param threshold threshold duration below which fixations of clusters should be deleted in \code{isFixVar} 
  #' @return data where \code{ixFixVar} is cleaned for too short fixations
  
  # If name for numFixVar not defined:
  if(is.null(numFixVar)){
    numFixVar <- paste0(numFixVar, "_fixIdx")
    cat(paste0("Argument numFixVar unspecified, now called ", numFixVar, "\n"))
  }
  
  # Defaults:
  trialVar <- "trialnr"
  nTrials <- length(unique(data[, trialVar])) # number trials
  
  ## Loop through trials: 
  for (iTrial in 1:nTrials){ # iTrial <- 1
    
    cat(paste0("Start trial ", iTrial, "\n"))
    
    rowIdx <- which(data[, trialVar] == iTrial) # row indices of trial 
    nFix <- max(data[rowIdx, numFixVar], na.rm = T) # number fixations of trial
    if(nFix == -Inf){nFix <- 0}
    # cat(paste0("Number fixations is ", nFix, "\n"))
    
    ## Loop through fixations:
    if (nFix > 0){
      
      for (iFix in 1:nFix){ # loop through fixations
        
        fixIdx <- which(data[, trialVar] == iTrial & data[, numFixVar] == iFix) # row indices of fixation
        durFix <- length(fixIdx) # fixation duration
        # cat(paste0("Fixation ", iFix, ": Duration is ", durFix, "\n"))
        
        if (durFix < threshold){ # if shorter than tolerance
          data[fixIdx, isFixVar] <- 0 # set status of fixation back to zero
          data[fixIdx, numFixVar] <- NA # delete cluster index
          cat(paste0("Delete fixation ", iFix, ": Only ", durFix, " samples\n"))
        }
      }
    } # end iFix
  } # end iTrial
  
  return(data)
}

# ============================================================================ #
#### Compute gaze deviation relative to pre-cue baseline: ####

# subjectID = 1;
# segMessage = "StartCue"; beforeSamples = 866; afterSamples = 1300
# baselineTime = c(-866, -367)

# aggrMethod = "distanceDeriv"

aggregate_gaze <- function(subjectID, aggrMethod, segMessage, beforeSamples, afterSamples, baselineTime){
  
  require(stringr) # for str_pad
  require(zoo) # for na.approx
  
  ## Fixed settings:
  eyeMeasure <- "gaze"
  # aggrMethod <- "distanceBaselineTrial"
  scrx <- 516*2 # screen size x dimension
  scry <- 392*2 # screen size y dimension
  
  ## Subject ID as string:  
  cat("* ------------------------------------------------------------ *\n")
  subStr <- str_pad(subjectID, width = 2, side = "left", pad = "0")
  cat(paste0(">>> Subject ", subStr, ": Start processing\n"))
  
  # --------------------------------------------------------------------------------------------- #
  ### Determine input directory:
  
  beforeSamplesStr <- str_pad(beforeSamples, width = 4, side = "left", pad = "0")
  afterSamplesStr <- str_pad(afterSamples, width = 4, side = "left", pad = "0")
  newDir <- paste0("timecourse_", segMessage, "_", beforeSamplesStr, "_", afterSamplesStr, "ms/")
  cat(paste0(">>> Load time course data per trial from ", newDir, "\n"))
  xpDir <- paste0(dirs$processedDataDir, eyeMeasure, "/", "xp_", newDir)
  ypDir <- paste0(dirs$processedDataDir, eyeMeasure, "/", "yp_", newDir)

  # --------------------------------------------------------------------------------------------- #
  ### Determine output directory and file name:
  
  ## Output directory:
  outputDir <- paste0(dirs$processedDataDir, eyeMeasure, "/", aggrMethod, "_timecourse_", 
                      segMessage, "_", beforeSamplesStr, "_", afterSamplesStr, "ms/")
  dir.create(outputDir, showWarnings = FALSE, recursive = TRUE)
  
  ## Output file:
  outFileName <- paste0("MGNGUNC_Sub_", subStr, "_", eyeMeasure) 
  fullFileName <- paste0(outputDir, outFileName, ".csv")
  
  # --------------------------------------------------------------------------------------------- #
  ### Load raw data:
  
  ## Load time course data for this subject:
  xpMat <- as.matrix(read.csv(paste0(xpDir, "MGNGUNC_Sub_", subStr, "_gaze.csv")))
  ypMat <- as.matrix(read.csv(paste0(ypDir, "MGNGUNC_Sub_", subStr, "_gaze.csv")))
  
  # --------------------------------------------------------------------------------------------- #
  ### Select method and perform aggregation:
  
  if (aggrMethod == "distanceGroupMean"){
    
    ## Compute Euclidean distance per trial:
    cat(paste0(">>> Compute distance from center of screen at ", scrx/2, ", ", scry/2, "\n"))
    outMat <- sqrt((xpMat - scrx/2)^2 + (ypMat - scry/2)^2)  
    
    
    # ------------------------------------------------------------------------------------------- #
  } else if (aggrMethod == "distanceBaselineTrial"){
    
    ## Determine baseline window and segment indices:
    cat(paste0(">>> Baseline trials from ", baselineTime[1], " to ", baselineTime[2], "\n"))
    timeVec <- seq(-1*beforeSamples, afterSamples, 1) # create time vector
    if(length(timeVec) != ncol(xpMat)){stop("Time vector based on beforeSamples/ afterSamples and columns of xpMat do not match")}
    baselineVec <- seq(1*baselineTime[1], baselineTime[2], 1) # in units of ms
    baselineIdx <- which(baselineVec %in% timeVec) # in index units
    
    ## Initialize:
    nTrial <- nrow(xpMat)
    nTime <- ncol(xpMat)
    xpDistMat <- matrix(NA, nrow = nTrial, ncol = nTime)
    ypDistMat <- matrix(NA, nrow = nTrial, ncol = nTime)
    
    ## Loop over trials, subtract trial-by-trial baseline:
    
    for (iTrial in 1:nTrial){ # iTrial <- 1
      
      ## xp:
      baseVec <- xpMat[iTrial, baselineIdx] # extract baseline
      baseVec_z <- as.numeric(scale(baseVec)) # z-standardize
      validIdx <- which(abs(baseVec_z) <= 2) # reject outliers at 2*z in computation of baseline
      nReject <- length(baseVec) - length(validIdx)
      cat(paste0(">>> Trial ", iTrial, ", xp: Reject ", str_pad(nReject, width = 3, side = "left", pad = "0"), 
                 " out of ", length(baseVec), " samples\n"))
      baselineMeanTrial <- mean(baseVec[validIdx], na.rm = T) # compute mean in baseline window
      if (is.na(baselineMeanTrial)){cat(paste0("Trial ", iTrial, ": Baseline is NA"))} # interpolate of no mean at all/ if NA: whole trial NA
      xpDistMat[iTrial, ] <- xpMat[iTrial, ] - baselineMeanTrial # subtract mean
      
      ## yp:
      baseVec <- ypMat[iTrial, baselineIdx] # extract baseline
      baseVec_z <- as.numeric(scale(baseVec)) # z-standardize
      validIdx <- which(abs(baseVec_z) <= 2) # reject outliers at 2*z in computation of baseline
      baselineMeanTrial <- mean(baseVec[validIdx], na.rm = T) # compute mean in baseline window
      ypDistMat[iTrial, ] <- ypMat[iTrial, ] - baselineMeanTrial
    } # end iTrial

    ## Compute Euclidean distance per trial:
    outMat <- sqrt(xpDistMat^2 + ypDistMat^2)  
    
    ## Inspect:
    # iTrial <- 1
    # selTime <- 1:20
    # plot(xpDistMat[iTrial, selTime])
    # plot(ypDistMat[iTrial, selTime])
    # plot(outMat[iTrial, selTime])
    
    # ------------------------------------------------------------------------------------------- #
  } else if (aggrMethod == "distanceAbsDeriv"){
  
    ## Compute Euclidean distance per trial:
    cat(paste0(">>> Compute distance from center of screen at ", scrx/2, ", ", scry/2, "\n"))
    euclDistMat <- sqrt((xpMat - scrx/2)^2 + (ypMat - scry/2)^2)  
    
    ## Loop over trials, subtract trial-by-trial baseline:
    cat(">>> Compute derivative on Euclidean distance from center\n")
    
    nTrial <- nrow(xpMat)
    nTime <- ncol(xpMat)
    derivMat <- matrix(NA, nrow = nTrial, ncol = nTime) # initialize
    for (iTrial in 1:nTrial){ # iTrial <- 1
      derivMat[iTrial, 2:nTime] <- diff(euclDistMat[iTrial, ])
    }
    
    ## Inspect:
    # iTrial <- 1
    # selTime <- 1:20
    # plot(euclDistMat[iTrial, selTime])
    # plot(derivMat[iTrial, selTime])
    
    outMat <- abs(derivMat) # save absolute derivative
    
    # ------------------------------------------------------------------------------------------- #
  } else if (aggrMethod == "distanceAbsCumSum"){
    
    ## Compute Euclidean distance per trial:
    cat(paste0(">>> Compute distance from center of screen at ", scrx/2, ", ", scry/2, "\n"))
    euclDistMat <- sqrt((xpMat - scrx/2)^2 + (ypMat - scry/2)^2)  
    
    ## Loop over trials, subtract trial-by-trial baseline:
    cat(">>> Compute derivative on Euclidean distance from center\n")
    cat(">>> Interpolate NAs\n")
    cat(">>> Estimate SD of noise per trial, add random noise with estimated SD to interpolated data\n")
    cat(">>> Compute cumulative sum of (interpolated) derivative\n")
    
    nTrial <- nrow(xpMat)
    nTime <- ncol(xpMat)
    derivMat <- matrix(NA, nrow = nTrial, ncol = nTime) # initialize
    cumSumMat <- matrix(NA, nrow = nTrial, ncol = nTime) # initialize
    
    for (iTrial in 1:nTrial){ # iTrial <- 10
      
      ## Compute derivative:
      derivMat[iTrial, 2:nTime] <- diff(euclDistMat[iTrial, ])

      ## Extract:
      # selTime <- 1:2000
      # tmp1 <- derivMat[iTrial, selTime]
      tmp1 <- derivMat[iTrial, ]
      # plot(tmp1)
      
      ## Interpolate:
      nReject <- sum(is.na(tmp1))
      cat(paste0(">>> Trial ", iTrial, ": Interpolate ", str_pad(nReject, width = 3, side = "left", pad = "0"), 
                 " out of ", length(tmp1), " samples\n"))
      tmp2 <- na.approx(tmp1, na.rm = FALSE) 
      # plot(tmp2)
      
      ## Add noise to interpolated data:
      naIdx <- which(is.na(tmp1)) # detect which samples were interpolated
      sdTrial <- sd(tmp1, na.rm = T) # estimate SD based on valid samples
      tmp2[naIdx] <- rnorm(length(naIdx), 0, sdTrial) # add noise
      cat(paste0(">>> Trial ", iTrial, ": Add noise with SD = ", round(sdTrial, 3), ", to ", nReject, " interpolated samples\n"))
      # plot(tmp2)
      
      ## Detect residual (leading/ trailing) NAs:
      # sum(is.na(tmp2))
      # which(is.na(tmp2))
      validIdx <- which(!(is.na(tmp2)))
      
      ## Cumulative sum of absolute derivative:
      cat(paste0(">>> Trial ", iTrial, ": Compute cumulative sum of absolute derivative for ", length(validIdx), " samples available\n"))
      tmp3 <- cumsum(abs(tmp2[validIdx]))
      # plot(tmp3)
      
      ## Store:
      cumSumMat[iTrial, validIdx] <- tmp3
    }
    
    outMat <- cumSumMat
    
  } else {
    stop(paste0("Unknown aggrMethod ", aggrMethod))
  }
  
  # --------------------------------------------------------------------------------------------- #
  ### Save:
  
  ## Load time course data for this subject:
  cat(paste0(">>> Subject ", subStr, ": Save data\n"))
  write.csv(outMat, fullFileName, row.names = F)
  cat(paste0(">>> Subject ", subStr, ": Finished :-)\n"))
  
}  

# END OF FILE.
