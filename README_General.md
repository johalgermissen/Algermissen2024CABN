# README General

## Content ##
This collection contains the following files:
- **analyses**: scripts for all behavioral, gaze, pupil, and questionnaire data required to reproduce the reported results.
- **task**: scripts and files for running the task.

Code for the entire paper will be maintained under https://github.com/johalgermissen/Algermissen2024CABN, with a permanent copy of the code at the time of publication under https://github.com/denoudenlab/Algermissen2024CABN.

**(Raw data and processed data** for this project are available under https://data.ru.nl/collections/di/dcc/DSC_2019.00030_944 under a CC-BY-4.0 licence.

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
All data is available under https://data.ru.nl/collections/di/dcc/DSC_2019.00030_944 under a CC-BY-4.0 licence. 

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

