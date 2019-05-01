# /usr/bin/python

#Author: Ben Tannenwald
#Date: April 30, 2019
#Purpose: Store functions+variables for channel analysis

from lib.globalVariables import *
import yaml
from ROOT import TLatex

def setChannelData( _bar, _firstRun):
    with open('config/config_v1.yml', 'r') as ymlfile:
        cfg = yaml.load(ymlfile)    

    print('channelSet: ', type(nSigmaTimeCut))
    print('channelSet: ', type(maxX), type(minY), type(maxY))
    #print( cfg[_firstRun][_bar] )
    ampch_id[0]   = int(cfg[_firstRun][_bar]['photek_amp_digiIndex'])
    timech_id[0]  = int(cfg[_firstRun][_bar]['photek_time_digiIndex']) 
    namech[0]     = cfg[_firstRun][_bar]['photek_name']
    ampmin_cut[0] = int(cfg[_firstRun][_bar]['photek_lowAmpCut'])
    ampmax_cut[0] = int(cfg[_firstRun][_bar]['photek_highAmpCut'])
    ampch_id[1]   = int(cfg[_firstRun][_bar]['left_amp_digiIndex'])
    timech_id[1]  = int(cfg[_firstRun][_bar]['left_time_digiIndex']) 
    namech[1]     = cfg[_firstRun][_bar]['left_name']
    ampmin_cut[1] = int(cfg[_firstRun][_bar]['left_lowAmpCut'])
    ampmax_cut[1] = int(cfg[_firstRun][_bar]['left_highAmpCut'])
    ampch_id[2]   = int(cfg[_firstRun][_bar]['right_amp_digiIndex'])
    timech_id[2]  = int(cfg[_firstRun][_bar]['right_time_digiIndex']) 
    namech[2]     = cfg[_firstRun][_bar]['right_name']
    ampmin_cut[2] = int(cfg[_firstRun][_bar]['right_lowAmpCut'])
    ampmax_cut[2] = int(cfg[_firstRun][_bar]['right_highAmpCut'])
    minX          = float(cfg[_firstRun][_bar]['minX'])
    maxX          = float(cfg[_firstRun][_bar]['maxX'])
    centerX       = float(cfg[_firstRun][_bar]['bar_center_x'])
    BSX           = float(cfg[_firstRun][_bar]['bar_halfwidth_x'])
    minY          = float(cfg[_firstRun][_bar]['y_min'])
    maxY          = float(cfg[_firstRun][_bar]['y_max'])
    centerY       = float(cfg[_firstRun][_bar]['bar_center_y'])
    BSY           = float(cfg[_firstRun][_bar]['bar_halfwidth_y'])

    print('channelSet: ', minX, maxX, minY, maxY)
    print('channelSet: ', type(minX), type(maxX), type(minY), type(maxY))
    lowerPosCutX = centerX-BSX
    upperPosCutX = centerX+BSX  
    lowerPosCutY = centerY-BSY
    upperPosCutY = centerY+BSY
    
    for _iCh in range(0, NCH):
        latexch[_iCh] = TLatex(0.16,0.96, namech[_iCh])
        latexch[_iCh].SetNDC()
        latexch[_iCh].SetTextFont(42)
        latexch[_iCh].SetTextSize(0.03)

    
#def getChannelData( _bar, _biasVoltage, _firstRun, _lastRun):
#    if _bar == "box1" or _bar=="box2" or _bar=="box3":
        

