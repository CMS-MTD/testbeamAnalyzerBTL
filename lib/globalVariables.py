# /usr/bin/python

#Author: Ben Tannenwald
#Date: April 30, 2019
#Purpose: Store global variables for channel analysis


from ROOT import TLatex

NCH = 3
ampch_id   = [0 for i in range(0,NCH)]
timech_id  = [0 for i in range(0,NCH)]
namech     = ['' for i in range(0,NCH)]
latexch    = [TLatex() for i in range(0,NCH)]
ampmin_cut = [0 for i in range(0,NCH)]
ampmax_cut = [0 for i in range(0,NCH)]

rel_amp_cut_low = 0.85 #  low amp cut in fraction of MIP peak
rel_amp_cut_hig = 4.   # high amp cut in fraction of MIP peak

lowerTimeCut = 20. #  low time cut in ns
upperTimeCut = 60. # high time cut in ns
nSigmaTimeCut = 2. # n of sigma on time cut

minX = -99. # range of X in mm
maxX = -99. # range of X in mm
centerX = -99. # hodoscope X coordinate of crystal center in mm
BSX = -99.     # half-size of beam spot selection around the center in mm
minY = -99.    # range of Y in mm
maxY = -99.    # range of Y in mm
centerY = -99. # hodoscope Y coordinate of crystal center in mm
BSY = -99.     # half-size of beam spot selection around the center in mm

NPOSCUTSX = 8
lowerPosCutX = -99
upperPosCutX = -99
NPOSCUTSY = 6
lowerPosCutY = -99
upperPosCutY = -99

NBSCUTS = 6 # number of beam spot bins
BScut = [10, 7, 6, 3, 2, 1] #in mm around the center

sigma_ref = 0.015 # MCP time resolution in ps

nAmpBins = 250
ampMin = 0.
ampMax = 1000.
nTimeBins = 500
minTime = 0.
maxTime = 200.
nDeltatBins = 500
minDeltat = -20.
maxDeltat = +20.
    
