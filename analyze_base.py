# /usr/bin/python

#Author: Ben Tannenwald
#Date: April 30, 2019
#Purpose: Worker code to process/analyze testbeam data (loosely based on Andrea Benaglia's analyze_FNAL.C code)

import os, sys, argparse, getpass
import shutil

# *** 0. setup parser for command line
parser = argparse.ArgumentParser()
parser.add_argument("--bar", help="which bar to process, e.g. box1, box2, box3, single")
parser.add_argument("--firstRun", help="first run to analyze")
parser.add_argument("--lastRun", help="last run to analyze")
parser.add_argument("--timeAlgo", help="time fitting algorithm", default="LP2_30")
parser.add_argument("--biasVoltage", help="bias voltage [V]")
parser.add_argument("--isCondor", help="flag to turn off display of plots for condor submission", default=False, action='store_true')
args = parser.parse_args()

if len(sys.argv)==1: # no options, easy way to get help
    os.system('python analyze_base.py -h')
    quit()

# ** A. Test bar option and exit if unexpected
if(args.bar is None):
    print "#### Need to bar to analyze using --bar <box1/box2/box3/single> ####\nEXITING"
    quit()
else:
    if (args.bar == "box1" or args.bar == "box2" or args.bar == "box3" or args.bar == "single")==False:
        print "#### Need to bar to analyze using --bar <box1/box2/box3/single> ####\nEXITING"
        quit()
    else:
        print '-- Setting  bar = {0}'.format(args.bar)

# ** B. first run
if(args.firstRun is None):
    print "#### Need to specify first run to analyze --firstRun <run number> ###\nEXITING"
    quit()
else:
    if (args.firstRun).isdigit == False:
        print "#### Need to specify first run to analyze --firstRun <run number> ###\nEXITING"
        quit()
    else:
        print '-- Setting firstRun = {0}'.format(args.firstRun)

# ** C. last run
if(args.lastRun is None):
    print "#### Need to specify last run to analyze --lastRun <run number> ###\nEXITING"
    quit()
else:
    if (args.lastRun).isdigit == False:
        print "#### Need to specify last run to analyze --lastRun <run number> ###\nEXITING"
        quit()
    else:
        print '-- Setting lastRun = {0}'.format(args.lastRun)

# ** D. time Algo
if(args.timeAlgo == "IL_50"):
    print '-- Default timeAlgo = {0}'.format(args.timeAlgo)
else:
    print '-- Setting lastRun = {0}'.format(args.timeAlgo)

# ** E. bias voltage
if(args.biasVoltage is None):
    print "#### Need to specify bias voltage to analyze --biasVoltage <bias voltage> ###\nEXITING"
    quit()
else:
    if (args.lastRun).isdigit == False:
        print "#### Need to specify bias voltage to analyze --biasVoltage <bias voltage> ###\nEXITING"
        quit()
    else:
        print '-- Setting biasVoltage = {0}'.format(args.biasVoltage)

# ** F. Condor submission --> turn on batch mode
if args.isCondor==True:
    sys.argv.append( '-b' )


## *** 1. Parse together output directory and create if does not exist
## ** A. Top level directory of which runs
outputDir = 'Runs_'+args.firstRun+'to'+args.lastRun
if ( not os.path.exists(outputDir) ):
    print "Specified run directory {0} DNE.\nCREATING NOW".format(outputDir)
    os.system("mkdir {0}".format(outputDir))
shutil.copy('lib/index.php',outputDir+'/index.php')
# ** B. Sub-level directory of what bias
outputDir = outputDir + '/' + args.biasVoltage + 'V'
if ( not os.path.exists(outputDir) ):
    print "Specified bias voltage sub-directory {0} DNE.\nCREATING NOW".format(outputDir)
    os.system("mkdir {0}".format(outputDir))
shutil.copy('lib/index.php',outputDir+'/index.php')
# ** C. Sub-level directory of which bar
outputDir = outputDir + '/' + args.bar
if ( not os.path.exists(outputDir) ):
    print "Specified bar sub-directory {0} DNE.\nCREATING NOW".format(outputDir)
    os.system("mkdir {0}".format(outputDir))
shutil.copy('lib/index.php',outputDir+'/index.php')
# ** D. Sub-level directory of which time algorithm
outputDir = outputDir + '/' + args.timeAlgo + '/'
if ( not os.path.exists(outputDir) ):
    print "Specified timeAlgo sub-directory {0} DNE.\nCREATING NOW".format(outputDir)
    os.system("mkdir {0}".format(outputDir))
shutil.copy('lib/index.php',outputDir+'/index.php')
print '-- Setting outputDir = {0}'.format(outputDir)


###############################################################################################

###############################################################################################


import numpy as np
import ROOT as rt
import root_numpy as rtnp
from lib.channelHandler import *
from lib.fileHandler import *


#outputFitInfo = open( outputDir+'/timeAndAmpFitInfo.txt', 'w' )
dataFolder = '/eos/uscms/store/group/cmstestbeam/2019_04_April_CMSTiming/VME/RecoData/RecoWithTracks/v3/' # BBT, 4-30-19, LPC EOS

if args.isCondor==True:
    dataFolder = 'root://cmseos.fnal.gov//store/group/cmstestbeam/2019_04_April_CMSTiming/VME/RecoData/RecoWithTracks/v3/' # BBT, 4-30-19, LPC EOS

ch = setChannelData( args.bar, int(args.firstRun), int(args.lastRun))
if ch is None:
   exit()

pulse = returnChain( dataFolder, args.firstRun, args.lastRun)

# set booleans for user-specific studies
doPositionResiduals = False
if getpass.getuser() == 'benjtann':
    doPositionResiduals = True

# run analysis code
timePeak, timeSigma, mipPeak, mipPeak_err, h_ampCut = calculatePeakPositionMIP( pulse, ch, args.timeAlgo, outputDir)
ampWalkCorr = calculateAmpWalkCorrection( pulse, ch, args.timeAlgo, outputDir, timePeak, timeSigma, mipPeak )
ampCorrectedMeasurements = applyAmpWalkCorrection( pulse, ch, args.timeAlgo, outputDir, timePeak, timeSigma, mipPeak, h_ampCut, ampWalkCorr )
applyPositionCorrection( pulse, ch, args.timeAlgo, outputDir, timePeak, timeSigma, mipPeak, mipPeak_err, h_ampCut, ampWalkCorr, ampCorrectedMeasurements )
if doPositionResiduals:
    calculatePositionResiduals(pulse, ch, args.timeAlgo, outputDir, timePeak, timeSigma, mipPeak, mipPeak_err, h_ampCut, ampWalkCorr, ampCorrectedMeasurements)

#outputFitInfo.close()
