# /usr/bin/python

#Author: Ben Tannenwald
#Date: April 30, 2019
#Purpose: Store functions+variables for channel analysis

#from globalVariables import *
import lib.globalVariables as ch
#from lib.globalVariables import minX, maxX, minY, maxY
import yaml, os
from ROOT import TLatex

def findFirstRunNumberOfRange(_firstRun):
    _firstRunNumberOfAllRanges = [7682, 7850]
    _firstRunOfRange = 0

    for iRunNumber in _firstRunNumberOfAllRanges:
        if _firstRun <= iRunNumber:
            _firstRunOfRange = iRunNumber
    
    return _firstRunOfRange


def findLastRunNumberOfRange(_lasstRun):
    _lastRunNumberOfAllRanges = [7800, 7877]
    _lastRunOfRange = 0

    for iRunNumber in _lastRunNumberOfAllRanges:
        if _lastRun <= iRunNumber:
            _lastRunOfRange = iRunNumber
    
    return _lastRunOfRange


def setChannelData( _bar, _firstRun, _firstRun):
    
    _firstRunOfRange = findFirstRunNumberOfRange( _firstRun )
    _lastRunOfRange = findFirstRunNumberOfRange( _lastRun )
    _configName = 'config/config_{0}to{1}.yml'.format(_firstRunOfRange, _lastRunOfRange)

    if ( not os.path.exists(_configName) ):
        print("!!!! NO CONFIG FILE FOUND FOR RUNS [{0}, {1}]. EXITING".format(_firstRun, _lastRun))
        return None
    else:
        print("!!!! Using config file: {0}".format( _configName) )
            
    with open(_configName, 'r') as ymlfile:
        cfg = yaml.load(ymlfile)    

    #print( cfg[_firstRun][_bar] )
    ch.ampch_id[0]   = int(cfg[_firstRun][_bar]['photek_amp_digiIndex'])
    ch.timech_id[0]  = int(cfg[_firstRun][_bar]['photek_time_digiIndex']) 
    ch.namech[0]     = cfg[_firstRun][_bar]['photek_name']
    ch.ampmin_cut[0] = int(cfg[_firstRun][_bar]['photek_lowAmpCut'])
    ch.ampmax_cut[0] = int(cfg[_firstRun][_bar]['photek_highAmpCut'])
    ch.ampch_id[1]   = int(cfg[_firstRun][_bar]['left_amp_digiIndex'])
    ch.timech_id[1]  = int(cfg[_firstRun][_bar]['left_time_digiIndex']) 
    ch.namech[1]     = cfg[_firstRun][_bar]['left_name']
    ch.ampmin_cut[1] = int(cfg[_firstRun][_bar]['left_lowAmpCut'])
    ch.ampmax_cut[1] = int(cfg[_firstRun][_bar]['left_highAmpCut'])
    ch.ampch_id[2]   = int(cfg[_firstRun][_bar]['right_amp_digiIndex'])
    ch.timech_id[2]  = int(cfg[_firstRun][_bar]['right_time_digiIndex']) 
    ch.namech[2]     = cfg[_firstRun][_bar]['right_name']
    ch.ampmin_cut[2] = int(cfg[_firstRun][_bar]['right_lowAmpCut'])
    ch.ampmax_cut[2] = int(cfg[_firstRun][_bar]['right_highAmpCut'])
    ch.minX          = float(cfg[_firstRun][_bar]['x_min'])
    ch.maxX          = float(cfg[_firstRun][_bar]['x_max'])
    ch.centerX       = float(cfg[_firstRun][_bar]['bar_center_x'])
    ch.BSX           = float(cfg[_firstRun][_bar]['bar_halfwidth_x'])
    ch.minY          = float(cfg[_firstRun][_bar]['y_min'])
    ch.maxY          = float(cfg[_firstRun][_bar]['y_max'])
    ch.centerY       = float(cfg[_firstRun][_bar]['bar_center_y'])
    ch.BSY           = float(cfg[_firstRun][_bar]['bar_halfwidth_y'])


    ch.lowerPosCutX = ch.centerX - ch.BSX
    ch.upperPosCutX = ch.centerX + ch.BSX  
    ch.lowerPosCutY = ch.centerY - ch.BSY
    ch.upperPosCutY = ch.centerY + ch.BSY
    
    for _iCh in range(0, ch.NCH):
        ch.latexch[_iCh] = TLatex(0.16,0.96, ch.namech[_iCh])
        ch.latexch[_iCh].SetNDC()
        ch.latexch[_iCh].SetTextFont(42)
        ch.latexch[_iCh].SetTextSize(0.03)

    print('--- channel configuration loaded')

    return ch
    
#def getChannelData( _bar, _biasVoltage, _firstRun, _lastRun):
#    if _bar == "box1" or _bar=="box2" or _bar=="box3":
        

