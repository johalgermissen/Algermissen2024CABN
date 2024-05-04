# README General

## Content ##
This collection contains the following files:
- **analyses**: scripts for all behavioral, gaze, pupil, and questionnaire data required to reproduce the reported results.
- **data**: raw and pre-processed behavioral, eye-tracking (pupil and gaze), and questionnaire data.
- **task**: scripts and files for running the task.

## Root directory ##
Note that the **root directory** of this project needs to be adjusted inside analyses/regression/helpers/set_dirs.R. 

## Analyses folder ##
This folder contains all analyses scripts required to reproduce the results reported in the main results and the supplementary materials.

### Analysis software used ###
All analyses were performed in _R_. See the file /analyses/regression/helpers/package_manager.R_ for all respective package version numbers and a print out of both _sessionInfo()_ and _loadedNamespaces()_.

### helpers ###
Helper functions are located under analyses/regression/helpers.
- *package_manager.R*: Loads packages. Also contains version numbers and print-outs of _sessionInfo()_ and _loadedNamespaces()_), set _scipen_ and _factor codings_.
- *set_dirs.R*: Set directories for code and data. Mind adjusting the root directory to your file structure.


### functions ###
Custom functions are located under analyses/regression/functions.
Contains functions common to preprocessing and analyzing behavioral, gaze, pupil, and questionnaire data:
- *00_mgngunc_functions_preprocess_gaze.R*: Contains several functions for pre-processing gaze data (called by *01_mgngunc_preprocess_gaze.R*).
- *00_mgngunc_functions_preprocess_pupil.R*: Contains several functions for pre-processing pupil data (called by *01_mgngunc_preprocess_pupil.R*).
- *00_functions_regression.R*: Contains several functions for pre-processing behavioral data, plotting, and fitting linear and additive mixed-effects models (called by *02_mgngunc_plot.R*, *02_mgngunc_regression.R*, *05_mgngunc_gamm.R*).
- *00_mgngunc_functions_timecourse.R*: Contains several functions for loading, aggregating, plotting, and analyzing ms-by-ms time course data (gaze and pupil; called by *03_mgngunc_gaze.R* and *04_mgngunc_timecourse.R*).

### main scripts ###
Main scripts are located under analyses/regression/.
They comprise scripts to pre-process and analyze eye-tracking data for each respective sample. Run files in numerical order:

- *01_mgngunc_preprocess_gaze.R*: Pre-process gaze data.
- *01_mgngunc_preprocess_pupil.R*: Pre-process pupil data.
- *02_mgngunc_plot.R*: Create raw data plots of responses, RTs, and pupil data.
- *02_mgngunc_regression.R*: Fit mixed-effects linear/ logistic regression models to responses, RTs, pupil baselines, pupil dilations.
- *03_mgngunc_gaze.R*: Create plots and test differences in ms-by-ms gaze time course data.
- *04_mgngunc_timecourse.R*: Create plots and test differences in ms-by-ms pupil time course data.
- *05_mgngunc_gamm.R*: Fit generalized additive mixed effects models (GAMMs) to pupil baseline and pupil dilation data.
- *06_mgngunc_questionnaires.R*: Correlate regression coefficients from linear mixed-effects models with questionnaire scores.
- *07_mgngunc_power.R*: Estimate post-hoc power based on the ICC estimated from intercept-only models on response, RT, and dilation data.
- *08_mgngunc_plots_final.R*: Create final plots used in manuscript.

## Data folder ##
This folder contains raw (and pre-processed) behavioral, gaze, pupil, and questionnaire data.

#### rawData ####
- *behavior*: behavioral raw data (separate file per subject) with the following variables:
    - *Subject*: integer, subject number.
    - *Age*: integer, subject age.
    - *Gender*: string (male, female, other), subject gender.
    - *Hand*: string (leftHand, rightHand), subject dominant hand.
    - *Eye*: string (leftEye, rightEye), subject dominant eye.
    - *Block*: integer between 1 and 4, block number.
    - *Trialnr*: integer between 1 and 256, trial number.
    - *Stimulus*: string, cue identifier within block, letter (A-D) indicates block, number (1-4) cue identity (see task folder).
    - *Condition*: integer between 1 and 4, cue conditions: 1 = Go-to-Win, 2 = Go-to-Avoid, 3 = NoGo-to-Win, 4 = NoGo-to-Avoid.
    - *ReqAction*: 1 or 0, required response, either Go (1) or NoGo (0).
    - *Valence*: 1 or 0, cue valence, either Win (1) or Avoid (0).
    - *Manipulation*: 1 or 0, arousal manipulation, 1 = high/ angry face, 0 = low/ neutral face.
    - *Validity*: 1 or 0, feedback valdidity, either valid (1) or invalid (0).
    - *Response*: 1 or 0, either Go response (1) or NoGo response (0).
    - *ACC*: 1 or 0, either correct (1) or incorrect (0) response.
    - *RT*: float, RT of any Go response (in sec.), NA for NoGo responses.
    - *Outcome*: integer between -1 and 1, outcome obtained, 1 = reward, 0 = neutral, -1 = punishment.
    - For further information, see the *wrapper_preprocessing()* function in the *00_mgngunc_functions_regression.R* script.
- *pupil*: eye-tracking raw data (separate file per subject; both as .edf and converted to .asc) collected with an Eyelink 1000. Data can be conveniently read in with the *eyelinker* package in R. Data contains four columns (interrupted by messages); columns are:
    - sample identifier (in ms relative to recording start).
    - x-coordinate of gaze.
    - y-coordinate of gaze.
    - pupil diameter.
    - For more information, see the respective file header and the Eyelink 1000 documentation.
- *questionnaires*: questionnaire data read out from Limesurvey:
    - limesurvey_survey_676118.lss: readout from Limesurvey.
    - limesurvey_survey_676118.txt: readout from Limesurvey.
    - QuestionnaireRaw.sav: questionnaire raw data in SPSS file (see *load_preprocess_questionnaires()* function in the *00_mgngunc_functions_regression.R* script for more information).
    - survey_576118_en.xml: file with original questions read out from Limesurvey (see *load_preprocess_questionnaires()* function in the *00_mgngunc_functions_regression.R* script for more information).
    - survey_archive_576118.lsa: readout from Limesurvey.

#### processedData ####
This folder contains data sets of pre-processed and aggregated gaze and pupil data.
- folder *gaze*: 
    - *distanceBaselineTrial_timecourse_StartCue_0866_1300ms*: file per subject with Euclidean distance of gaze position from pre-trial baseline from 866 ms before until 1300 ms after cue onset; rows are trials, columns are time samples (ms resolution).
    - *xp_timecourse_StartCue_0866_1300ms*: file per subject with x-position of gaze from 866 ms before until 1300 ms after cue onset; rows are trials, columns are time samples (ms resolution).
    - *yp_timecourse_StartCue_0866_1300ms*: file per subject with y-position of gaze from 866 ms before until 1300 ms after cue onset; rows are trials, columns are time samples (ms resolution).
- folder *pupil*:
    - *pupil_timecourse_StartMask_1000_2966ms*: file per subject with pupil width from 1000 ms before until 2966 ms after mask onset; rows are trials, columns are time samples (ms resolution).
    - *pupil_timecourse_StartOutcome_1000_2000s*: file per subject with pupil width from 1000 ms before until 2000 ms after outcome onset; rows are trials, columns are time samples (ms resolution).
    - *pupil_trial_StartMask_1000_2966ms*: file per subject with subject ID, trial number, pupil baseline, and pupil dilation (columns) per trial (rows) in time window from 0 ms till 2966 ms after mask onset.
    - *pupil_trial_StartOutcome_1000_2000ms*: file per subject with subject ID, trial number, pupil baseline, and pupil dilation (columns) per trial (rows) in time window from 0 ms till 2000 ms after outcome onset.

## Task folder ##
This folder contains all the code and files necessary to run the task in PsychoPy.

- *MGNGUNC_Task.py*: script to execute task.
- *DataOutput*: target directory to store behavioral output data.
- *FaceStimuli*: prime stimuli (angry and neutral face), foreward mask (pixels scrambled), backward mask (another neutral face).
- *Instructions*: instructions stored as .png files to be loaded in task; created in file *Instructions_MGNGjohalg.pptx*.
- *TaskStimuli*: contains task stimuli (cues displayed during practice trials and actual task; all in grayscale).
- *TrialHandler*: files with input variables per trial:
    - *MGNGUnCreateTrialSequence.R*: script used to create individual sequence with constraints per subject ID.
    - *MGNGUnCreatePracticeTrials.R*: script used to create sequences for four practice trial blocks (fixed across subject IDs).
    - *stimuluslist_pract_PartX.csv*: trial sequences for four practice blocks.
    - *stimuluslist_test_stimNum_256_sub_XXX_PartX.csv*: individual trial sequences per block per subject ID.

END of file.

