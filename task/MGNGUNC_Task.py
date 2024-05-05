#!/usr/bin/env python2
# -*- coding: utf-8 -*-
__prj__ = 'Motivational Go/NoGo Task with pupil'
__license__ = 'See attached licence'
__author__ = 'Johannes Algermissen'
__email__ = 'j.algermissen@donders.ru.nl'
__date__ = '2019/04/15 '

###########################################################################
# Load modules
###########################################################################

# eye-tracking modules
from pylink import *
import time
import gc
import sys
import os

# my own modules
import os.path # for absolute path
from psychopy import core, visual, event, data, gui
from Tkinter import *
# import pygame
import random

# track logging
import logging
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)

global res

useEyeTracker = False 
###########################################################################
# Request participant number, gender, and age via GUI
###########################################################################

def function_gui_info():
    guiDict = {"1. Subject":999,"2. Age":999,"3. Gender":["male","female","other"],"4. Hand":["leftHand","rightHand"],
    "5. Eye":["leftEye","rightEye"],"6. Use eyetracker":["1: Yes","2: No"], "7. Version":["1: Full","2: Demo"], "8. nRep":16} # define dictionary with infos to be asked for
    dialog = gui.DlgFromDict(dictionary = guiDict, title = "Go/NoGo Task", fixed = ("options")) # ask for infos
    if dialog.OK: # if info received: define variables
        subject = guiDict["1. Subject"]
        age = int(guiDict["2. Age"])
        gender = guiDict["3. Gender"]
        hand = guiDict["4. Hand"]
        eye = guiDict["5. Eye"]
        if guiDict["6. Use eyetracker"] == "1: Yes":
            useEyeTracker = True
        else: 
            useEyeTracker = False
        if guiDict["7. Version"] == "1: Full":
            isFull = True
        else: 
            isFull = False
        nTrials = guiDict["8. nRep"] * 16;
        exists = os.path.isfile("TrialHandler/stimuluslist_test_stimNum_{}_sub_{:0>3d}_Part1.csv".format(nTrials,subject))        
        if not exists:
            print u"Stop -- trialHandler file for subject {} with {} trials does not exist".format(subject,nTrials)
            core.quit()
        return (subject, age, gender, hand, eye, useEyeTracker, isFull, nTrials)
    else: # otherwise quit experiment and display reason
        print "Stop -- missing participant information"
        core.quit()
        

###########################################################################
# Waitkeys in pygame
###########################################################################

def waitforbutton(waitTime):
    core.wait(waitTime) # make sure that participants see instructions for a certain amount of time, and do not just click through it        
    keyPress = event.waitKeys(keyList = ['space','escape'], maxWait = 3000)
    if keyPress[0][0] == 'escape': # If abort:
        if useEyeTracker & (block == 'Test'):
            getEYELINK().sendMessage("EXPERIMENT ABORTED")
            EYELINK.stopRecording();
            EYELINK.setOfflineMode();
            msecDelay(500); # in ms
            #Close the file and transfer it to Display PC
            EYELINK.closeDataFile()
            EYELINK.receiveDataFile(edfFileName, edfFileName)
        core.quit()

###########################################################################
# Overall function displaying slide of instructions
###########################################################################

def function_display_instructions(file, waitTime = 0.3): # default 1.5 sec
    instr_image.setImage(file) # retrieve new image
    instr_image.draw()
    win.flip()
    waitforbutton(waitTime = waitTime)
    win.flip() # clean screen 
    core.wait(0.3) # wait for 300 ms until next screen follows (avoid rapid transitions)
    
###########################################################################
# Display countdown at beginning of each block
###########################################################################

def function_countdown():
    for i in ["3","2","1","Start!"]:
        if 'escape' in event.getKeys(): # exit with Escape
            core.quit()
        center_fix.setText(i)
        center_fix.draw()
        win.flip() 
        core.wait(1) # make count down slower
    
###########################################################################
# Run trials -- either for practice or for test
###########################################################################

def function_run_trials(trials, block, part):
    # Display countdown
    function_countdown()
    # Initial fixation cross to calm down pupil:
    center_fix.setText('+')
    center_fix.draw()
    win.flip()
    # Intialize trial
    trialcount = 0
    countACC = 0
    # Continue recording
    if useEyeTracker & (block == 'Test'):
        error = EYELINK.startRecording(1,1,0,0) # 
        if error: core.quit() #return error;    
        pylink.beginRealTimeMode(100)
        core.wait(2.0) # wait to get baseline for first trial
    # Run through trials
    for trial in trials:
        ###########################################################
        # Prepare trial
        trialcount += 1 # update trial
        print(u'Trialcount is {} ####################################################################'. format(str(trialcount)))
        if event.getKeys(['escape']): # exit with Escape
            core.quit()
        # Retrieve picture and answers
        trialnr         = trial['trialnr'] # either retrieve or count yourself
        stimulus        = trial['stimulus']
        condition       = trial['condition']
        valence         = trial['valence']
        reqAction       = trial['reqAction']
        manipulation    = trial['manipulation']
        goValidity      = trial['goValidity']
        nogoValidity    = trial['nogoValidity']
        ISI             = trial['ISI']
        ITI             = trial['ITI']
        print(u'condition is {}'. format(str(condition)))
        print(u'valence is {}'. format(str(valence)))
        print(u'reqAction is {}'. format(str(reqAction)))
        print(u'manipulation is {}'. format(str(manipulation)))
        print(u'goValidity is {}'. format(str(goValidity)))
        print(u'nogoValidity is {}'. format(str(nogoValidity)))
        print(u'ISI is {}'. format(str(ISI)))
        print(u'ITI is {}'. format(str(ITI)))
        ###########################################################
        # Write to eye-tracker at the beginning of trial
        if useEyeTracker & (block == 'Test'):
            getEYELINK().sendMessage("StartTrial")
            getEYELINK().sendMessage("Block")
            getEYELINK().sendMessage("{}".format(part))
            getEYELINK().sendMessage("TrialNr")
            getEYELINK().sendMessage("{}".format(trialnr))
            getEYELINK().sendMessage("Stimulus")
            getEYELINK().sendMessage("{}".format(stimulus))
            getEYELINK().sendMessage("Condition")
            getEYELINK().sendMessage("{}".format(condition))
            getEYELINK().sendMessage("Valence")
            getEYELINK().sendMessage("{}".format(valence))
            getEYELINK().sendMessage("Required Action")
            getEYELINK().sendMessage("{}".format(reqAction))
            getEYELINK().sendMessage("Manipulation")
            getEYELINK().sendMessage("{}".format(manipulation))
            getEYELINK().sendMessage("GoValidity")
            getEYELINK().sendMessage("{}".format(goValidity))
            getEYELINK().sendMessage("NoGoValidity")
            getEYELINK().sendMessage("{}".format(nogoValidity))
            getEYELINK().sendMessage("ISI")
            getEYELINK().sendMessage("{}".format(ISI))
            getEYELINK().sendMessage("ITI")
            getEYELINK().sendMessage("{}".format(ITI))
            getEYELINK().sendMessage("StartMask")
        ###########################################################
        # Set manipulation:
        # a) Forward mask: Phase-scrambled
        prime_image.setImage('FaceStimuli/Final_forward_mask.png')
        prime_image.draw()
        win.flip()
        core.wait(0.250)
        # b) Prime:
        if block == 'Test':
            if manipulation == 1:
                prime_image.setImage('FaceStimuli/Final_prime_a.png')
            else:
                prime_image.setImage('FaceStimuli/Final_prime_n.png')
            prime_image.draw()
            win.flip()
            core.wait(0.016) # refresh-rate 120 Hz, so use 8.333 ms, 16.6666 ms, 25 ms
            # Micah Allen: 0.016 (2 frames); 0.032 (4 frames); too long: 0.050
        # c) Stimulus Onset asynchrony
        # core.wait(0.016) 
        # d) Backward mask: neutral
        prime_image.setImage('FaceStimuli/Final_backward_mask.png')
        prime_image.draw()
        win.flip()
        core.wait(0.100)
        ###########################################################
        # Write to eye-tracker that stimulus presentation starts (i.e. fixations end, which matters more; afterwards only dilation important)
        if useEyeTracker & (block == 'Test'):
            getEYELINK().sendMessage("EndMask")
            getEYELINK().sendMessage("StartCue")
        ###########################################################
        # Present cue
        cue_image.setImage('TaskStimuli/Final_{}_grayscale.png'.format(stimulus))
        cue_image.draw()
        win.flip()
        # Register response
        rtClock = core.Clock()
        keyPress = event.waitKeys(keyList = ['space','escape'], maxWait = 1.3, timeStamped = rtClock)
        time_dif = rtClock.getTime()
        if keyPress:
            if keyPress[0][0] == 'escape': # If abort:
                if useEyeTracker & (block == 'Test'):
                    getEYELINK().sendMessage("EXPERIMENT ABORTED")
                core.quit()
            elif keyPress[0][0] == 'space': # if Go
                response = 1
                RT = time_dif
                print(u'Go!!!!!!!!!')
                print(u'RT is {}'.format(RT))
                core.wait(1.3 - RT) # wait till end of 1.3 sec
            else: 
                print(u'Error in response coding')
                core.quit()
        else: # if NoGo
            response = 0
            RT = 'NA'
        win.flip()
        if useEyeTracker & (block == 'Test'):
            getEYELINK().sendMessage("StopCue")
        ###########################################################
        # Check accuracy
        if reqAction == response:
            ACC = 1
            countACC += 1 # one more correct response
        else:
            ACC = 0
        print(u'response is {}.'.format(str(response)))
        print(u'ACC is {}.'.format(str(ACC)))
        ###########################################################
        # Inter-stimulus interval (ISI) between cue and outcome: 1300 - 1700 ms
        center_fix.setText('+')
        center_fix.draw()
        win.flip()
        core.wait(ISI)
        win.flip()
        ###########################################################
        # Determine validity:
        if response == 1:
            validity = goValidity
        else:
            validity = nogoValidity
        # Select outcome based on cue validity
        if valence == 1 and ACC == 1 and validity == 1:
            outcome = 1;
        elif valence == 1 and ACC == 1 and validity == 0:
            outcome = 0;
        elif valence == 1 and ACC == 0 and validity == 1:
            outcome = 0;
        elif valence == 1 and ACC == 0 and validity == 0:
            outcome = 1;
        elif valence == 0 and ACC == 1 and validity == 1:
            outcome = 0;
        elif valence == 0 and ACC == 1 and validity == 0:
            outcome = -1;
        elif valence == 0 and ACC == 0 and validity == 1:
            outcome = -1;
        elif valence == 0 and ACC == 0 and validity == 0:
            outcome = 0;
        else: 
            outputtext = (u'Error: Cannot determine outcome given valence = {}, ACC = {}, validity = {}.').format(str(valence),str(ACC),str(validity))
            print(outputtext)
            core.quit()
        ###########################################################
        # Present outcome 
        if outcome == 1:
            cue_image.setImage('TaskStimuli/Final_reward_grayscale.png')
        elif outcome == -1:
            cue_image.setImage('TaskStimuli/Final_punishment_grayscale.png')
        else:
            cue_image.setImage('TaskStimuli/Final_neutral_grayscale.png')
        if useEyeTracker & (block == 'Test'):
            getEYELINK().sendMessage("StartOutcome")
        cue_image.draw()
        win.flip()
        core.wait(0.7) # Present for 700 ms
        win.flip()
        if useEyeTracker & (block == 'Test'):
            getEYELINK().sendMessage("StopOutcome")
        ###########################################################
        # Write to eye-tracker at the end of trial
        if useEyeTracker & (block == 'Test'):
            getEYELINK().sendMessage("EndTrial")
            getEYELINK().sendMessage("Response")
            getEYELINK().sendMessage("{}".format(response))
            getEYELINK().sendMessage("ACC")
            getEYELINK().sendMessage("{}".format(ACC))
            getEYELINK().sendMessage("RT")
            getEYELINK().sendMessage("{}".format(RT))
            getEYELINK().sendMessage("Outcome")
            getEYELINK().sendMessage("{}".format(outcome))
        # write results (but only of actual trials, not practice trials)
        if block == 'Test':
            file_trials = open('DataOutput/MGNGUN_{:0>2d}_{:0>2d}.csv'.format(subject,outputAddon),'a') # open output file
            file_trials.write("{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}\n".format(subject,age,gender,hand,eye,part,trialnr,stimulus,condition,reqAction,valence,manipulation,validity,response,ACC,RT,outcome))
            file_trials.close()
        ###########################################################
        # Intertrial interval (ITI): 1800-2300 ms
        center_fix.setText('+')
        center_fix.draw()
        win.flip()
        core.wait(ITI) # wait for ITI (in seconds) so pupil returns to baseline levels
        win.flip()
        ###########################################################
        # End of trials
    # Pause eye-tracker
    if useEyeTracker & (block == 'Test'):
#        pylink.endRealTimeMode();    
        EYELINK.stopRecording();
    # Calculate global performance statistics
    print(u'countACC is {}.'.format(str(countACC)))
    percentACC = round(float(countACC)/trialcount*100,0)
    return(percentACC)

###########################################################################
# Run block including instructions and practice stimuli
###########################################################################

def function_block(block = None, part = None):
    # Load stimuli with trialhandler
    if block == 'Pract':
        stimulusList = data.importConditions("TrialHandler/stimuluslist_pract_Part{}.csv".format(part)) # Load test stimuli   
    elif block == 'Test':
        stimulusList = data.importConditions("TrialHandler/stimuluslist_test_stimNum_{}_sub_{:0>3d}_Part{}.csv".format(nTrials,subject, part)) # Load test stimuli   
    else:
        outputtext = (u'Error: Invalid settings for coming block.')
        print(outputtext)
        core.quit()        
    # Prepare trials:    
    trials = data.TrialHandler(stimulusList, nReps = 1, method = "sequential")
    # Assign cue pictures to cue numbers (1-8):
    # Do task:
    (percentACC) = function_run_trials(trials = trials, block = block, part = part)
    return(percentACC)
###################################################################################
# Main function (allows for running the program if executed from the file itself)
###################################################################################

def function_main():
    # Display instructions at beginning
    if isFull == True:
        # General instructions:
        function_display_instructions('Instructions/MGNGUNC_Instructions_General_01.jpg')
        function_display_instructions('Instructions/MGNGUNC_Instructions_General_02.jpg')
        function_display_instructions('Instructions/MGNGUNC_Instructions_General_03.jpg')
        function_display_instructions('Instructions/MGNGUNC_Instructions_General_04.jpg')
        function_display_instructions('Instructions/MGNGUNC_Instructions_General_05.jpg')
        function_display_instructions('Instructions/MGNGUNC_Instructions_General_06.jpg')
        function_display_instructions('Instructions/MGNGUNC_Instructions_General_07.jpg')
        function_display_instructions('Instructions/MGNGUNC_Instructions_General_08.jpg')
        # Practice rounds:
        function_display_instructions('Instructions/MGNGUNC_Instructions_Pract_G2W01.jpg')
        function_block(block = 'Pract', part = 1) # 5 Go2Win trials
        function_display_instructions('Instructions/MGNGUNC_Instructions_Pract_G2W02.jpg')
        function_display_instructions('Instructions/MGNGUNC_Instructions_Pract_NG2W01.jpg')
        function_block(block = 'Pract', part = 2) # 5 NoGo2Win trials
        function_display_instructions('Instructions/MGNGUNC_Instructions_Pract_NG2W02.jpg')
        function_display_instructions('Instructions/MGNGUNC_Instructions_Pract_G2A01.jpg')
        function_block(block = 'Pract', part = 3) # Go2Avoid trials
        function_display_instructions('Instructions/MGNGUNC_Instructions_Pract_G2A02.jpg')
        function_display_instructions('Instructions/MGNGUNC_Instructions_Pract_NG2A01.jpg')
        function_block(block = 'Pract', part = 4) # NoGo2Avoid trials
        function_display_instructions('Instructions/MGNGUNC_Instructions_Pract_NG2A02.jpg')
#   # Test Rounds:
    allACC = []; # empty list that will contain ACC of each block
    function_display_instructions('Instructions/MGNGUNC_Instructions_BeforeTest_01.jpg')
    function_display_instructions('Instructions/MGNGUNC_Instructions_BeforeTest_02.jpg')
    # Block 1:
    function_display_instructions('Instructions/MGNGUNC_Instructions_BeforeTest_Block01.jpg')
    function_display_instructions('Instructions/MGNGUNC_Instructions_Remember.jpg')
    (percentACC) = function_block(block = 'Test', part = 1)
    allACC.append(percentACC)
    function_display_instructions('Instructions/MGNGUNC_Instructions_Break_Block01.jpg')
    if isFull == True:
        # Block 2:
        function_display_instructions('Instructions/MGNGUNC_Instructions_BeforeTest_Block02.jpg')
        function_display_instructions('Instructions/MGNGUNC_Instructions_Remember.jpg')
        (percentACC) = function_block(block = 'Test', part = 2)
        allACC.append(percentACC)
        function_display_instructions('Instructions/MGNGUNC_Instructions_Break_Block02.jpg')
        # Block 3:
        function_display_instructions('Instructions/MGNGUNC_Instructions_BeforeTest_Block03.jpg')
        function_display_instructions('Instructions/MGNGUNC_Instructions_Remember.jpg')
        (percentACC) = function_block(block = 'Test', part = 3)
        allACC.append(percentACC)
        function_display_instructions('Instructions/MGNGUNC_Instructions_Break_Block03.jpg')
        # Block 4:
        function_display_instructions('Instructions/MGNGUNC_Instructions_BeforeTest_Block04.jpg')
        function_display_instructions('Instructions/MGNGUNC_Instructions_Remember.jpg')
        (percentACC) = function_block(block = 'Test', part = 4)
    allACC.append(percentACC)
    # Compute overall accuracy:
    meanACC = sum(allACC)/len(allACC);
    totalACC = round(float(meanACC),0)
    # Displays scores
    # First print (for security in case participants advance task too fast)
    outputtext = (u'The task is over now. You have responded correctly in {} percent of the trials.\nPlease contact the experimenter now.').format(str(totalACC))
    print(outputtext)
    # Then on screen
    center_text.setPos((0,.2)) # first line of text
    center_text.setText(u'The task is over now. You responded correctly in {} percent of the trials.'.format(totalACC))
    center_text.draw()
    center_text.setPos((0,0)) # second line of text
    if percentACC >= 75:
        center_text.setText(u'This gives you an extra reward!')
    else:
        center_text.setText(u'Unfortunately, that is not enough for an extra reward.')
    center_text.draw()
    center_text.setPos((0,-.2)) # third line of text
    center_text.setText(u'Please contact the experimenter now.')
    center_text.draw()
    win.flip()
    waitforbutton(waitTime = 2) # prevent participants from advancing task prematurely 
    win.flip() # clean screen 
    function_display_instructions('Instructions/MGNGUNC_Instructions_End01.jpg')
###################################################################################################################################################################################
###################################################################################################################################################################################
###################################################################################################################################################################################

###########################################################################
# Get initial participant info
(subject, age, gender, hand, eye, useEyeTracker, isFull, nTrials) = function_gui_info()
# Create output file addon
outputAddon = random.randrange(0,100)

###########################################################################
# connect to Eyelink 
if useEyeTracker:
#    edfFileName = "DataOutput/s" + str(subject) + "_" + str(outputAddon) + ".EDF"
    edfFileName = "s_{:0>2d}_{:0>2d}.EDF".format(subject,outputAddon) # write headers of file: 

    
    EYELINK = EyeLink('100.1.1.1') #EyeLinkListener();
    # If calibration worked:
#    print "do tracker setup"
#    EYELINK.doTrackerSetup()
#    print "TrackerSetup Done"
    
    EYELINK.openDataFile(edfFileName)
    
    # Actually start recording:
    error = EYELINK.startRecording(1,1,0,0) # 
    
###########################################################################
# Create global variables

# Eye-tracker
# Window
# Mind: color in PsychoPy is between -1 and 1, so use rgb/355*2-1 for conversion
win=visual.Window(fullscr=True, color=[.30,.30,.30], winType='pygame', units = 'norm') # check gray background color with PowerPoint; is 166/255=0.6509804 in original stimuli
win.setMouseVisible(False)
# Central text (for fixation cross)
center_fix  = visual.TextStim(win=win,text="", pos = (0,0), color=[0.03,0.03,0.03], wrapWidth = 0.8, alignHoriz = 'center', height = .2) # Test out color: .80? .30
center_text = visual.TextStim(win=win,text="", pos = (0,0), color=[1,1,1], wrapWidth = 0.8, alignHoriz = 'center', height = .1) # Test out color: white enough? used for final feedback
#center_text = visual.TextBox(window=win,text="", pos = (0,0), font_size=32, font_color=[.8,.8,.8]) # Test out color: white enough? used for fixation cross and final feedback
# Image
instr_image = visual.ImageStim(win=win,image='Instructions/MGNGUNC_Instructions_General_01.jpg', units = 'norm')
prime_image = visual.ImageStim(win=win,image='FaceStimuli/Final_forward_mask.png', pos = (0,0), size = (281, 381), units = 'pix') 
cue_image = visual.ImageStim(win=win,image='TaskStimuli/Final_A1_grayscale.png', pos = (0,0), size = (300, 300), units = 'pix') # original 720, 720 pixels
 
# Define trial output file
file_trials = open('DataOutput/MGNGUN_{:0>2d}_{:0>2d}.csv'.format(subject,outputAddon),'w') # write headers of file: 
file_trials.write("{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}\n".format('Subject','Age','Gender','Hand','Eye','Block','Trialnr','Stimulus','Condition','ReqAction','Valence','Manipulation','Validity','Response','ACC','RT','Outcome'))
file_trials.close()

###########################################################################
# Run experiment
###########################################################################

function_main()

###########################################################################
# End experiment
###########################################################################

# ---------------------------------------------
#---- stop recording and disconnect from iViewX
# ---------------------------------------------

# Stop recording and save idf file
if useEyeTracker:
    EYELINK.setOfflineMode();
    msecDelay(500); # in ms
    #Close the file and transfer it to Display PC
    EYELINK.closeDataFile()
    EYELINK.receiveDataFile(edfFileName, edfFileName)
# Quit
core.quit()
# END of program