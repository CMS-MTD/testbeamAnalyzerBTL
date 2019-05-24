# /usr/bin/python

#Author: Ben Tannenwald
#Date: April 30, 2019
#Purpose: Store functions for opening files

from ROOT import TFile, TChain, TTree, TProfile, TH1F, TH2F, TProfile2D, TCanvas, TLine, TBranch, TF1, TLegend, TGraphErrors, TLatex
from ROOT import kRed, kBlack, kBlue, kOrange, kMagenta, kYellow, kViolet
from ROOT import gPad, gStyle
import lib.dataQualityHandler as dq
from time import sleep
import numpy as np

gStyle.SetPadTopMargin(0.05);
gStyle.SetPadBottomMargin(0.13);
gStyle.SetPadLeftMargin(0.13);
gStyle.SetPadRightMargin(0.17);
gStyle.SetTitleOffset(1.20, "X");
gStyle.SetTitleOffset(1.80, "Y");
gStyle.SetTitleOffset(1.60, "Z");

   
def printCanvas( _canvas, _outputDir, _timeAlgo):

    _canvas.Draw();
    _algoLabel = TLatex(0.74, 0.96, _timeAlgo)
    _algoLabel.SetNDC(True)
    _algoLabel.SetTextFont(42)
    _algoLabel.SetTextSize(0.03)
    _algoLabel.Draw("same");
    _canvas.Print( '{0}/{1}_{2}.png'.format( _outputDir, _canvas.GetTitle(), _timeAlgo ) )

    return
               

def returnChain( _dataFolder, _firstRun, _lastRun):

    myTree = TChain("pulse", "pulse")
    f0 = TFile()
    trackingIsSynced = False
    for _iRun in range( int(_firstRun), int(_lastRun)+1):
        #f0 = TFile::Open( '{0}/RawDataSaver0CMSVMETiming_Run{1}_0_Raw.root'.format(_dataFolder, _iRun) ) # necessary format for condor... may need to revist later
        f0 = TFile.Open( '{0}/RawDataSaver0CMSVMETiming_Run{1}_0_Raw.root'.format(_dataFolder, _iRun) )
        if  f0 == None :
            print ( '!! FILE NOT FOUND {0}/RawDataSaver0CMSVMETiming_Run{1}_0_Raw.root'.format( _dataFolder, _iRun ))
            continue
    

        trackingIsSynced = dq.checkTrackingSynchronization( '{0}/RawDataSaver0CMSVMETiming_Run{1}_0_Raw.root'.format(_dataFolder, _iRun), _iRun)

        if trackingIsSynced == True: 
            myTree.Add( '{0}/RawDataSaver0CMSVMETiming_Run{1}_0_Raw.root'.format(_dataFolder, _iRun) ) # BBT, 4-30-19 
            print( 'adding run {0} {1}/RawDataSaver0CMSVMETiming_Run{0}_0_Raw.root'.format(_iRun, _dataFolder))
    
        else:
            print( '!! TRACKING DE-SYNCED, skiping run {0} {1}/RawDataSaver0CMSVMETiming_Run{0}_0_Raw.root'.format(_iRun, _dataFolder) )
    
  
    nEntries = myTree.GetEntries()
    print ( '>>> Events read: {0}'.format(nEntries) )

    return myTree


def getTimeAlgorithmBranch( _chain, _timeAlgo ):
    """ return branch of time algorithm"""  
    time = []
    
    #this is a very inefficient way to do this, but adaptive branch naming is hard - BBT 05/02/19
    if (_timeAlgo == 'LP2_20'):
        time = _chain.LP2_20
    elif (_timeAlgo == 'LP2_30'):
        time = _chain.LP2_30            
    elif (_timeAlgo == 'LP2_50'):
        time = _chain.LP2_50
    elif (_timeAlgo == 'IL_20'):
        time = _chain.IL_20
    elif (_timeAlgo == 'IL_30'):
        time = _chain.IL_30
    elif (_timeAlgo == 'IL_50'):
        time = _chain.IL_50
    
    return time

def eventHasGoodHit( _analysisStage, _timeAlgo, _myX, _myY, _chain=[], _channelInfo=[], _timePeak=[], _timeSigma=[], _mipPeak=[]):
    """ check quality of hits. depends on stage of analysis"""

    _isGoodBar = True

    ## 0. Cut on hit tracking
    if _myX == -999 or _myY == -999:
        _isGoodBar = False
    if _analysisStage == 0:
        return _isGoodBar
        
    ## 1. Cut on BS
    if( abs(_myX-_channelInfo.centerX) > _channelInfo.BSX or abs(_myY-_channelInfo.centerY) > _channelInfo.BSY ):
        _isGoodBar = False
    if _analysisStage == 1:
        return _isGoodBar
    
    ## 2. MCP
    # a. Cut on MCP amp
    _time_ref = _chain.gaus_mean[_channelInfo.timech_id[0]]
    if  _chain.amp[_channelInfo.ampch_id[0]] < _channelInfo.ampmin_cut[0] or _chain.amp[_channelInfo.ampch_id[0]] > _channelInfo.ampmax_cut[0] :
        _isGoodBar = False

    # b. Cut on MCP time
    if ( _time_ref < max(_timePeak[0] - _timeSigma[0]*_channelInfo.nSigmaTimeCut,_channelInfo.lowerTimeCut) or
         _time_ref > min(_timePeak[0] + _timeSigma[0]*_channelInfo.nSigmaTimeCut,_channelInfo.upperTimeCut) ):
        _isGoodBar = False
    if _analysisStage == 2:
        return _isGoodBar


    ## 3. Selected bar
    _time = getTimeAlgorithmBranch( _chain, _timeAlgo )
    if( _chain.amp[_channelInfo.ampch_id[1]] < max(_mipPeak[1]*_channelInfo.rel_amp_cut_low,_channelInfo.ampmin_cut[1]) ): _isGoodBar = False
    if( _chain.amp[_channelInfo.ampch_id[1]] > min(_mipPeak[1]*_channelInfo.rel_amp_cut_hig,_channelInfo.ampmax_cut[1]) ): _isGoodBar = False
    if( _chain.amp[_channelInfo.ampch_id[2]] < max(_mipPeak[2]*_channelInfo.rel_amp_cut_low,_channelInfo.ampmin_cut[2]) ): _isGoodBar = False
    if( _chain.amp[_channelInfo.ampch_id[2]] > min(_mipPeak[2]*_channelInfo.rel_amp_cut_hig,_channelInfo.ampmax_cut[2]) ): _isGoodBar = False
    if( _time[_channelInfo.timech_id[1]] < max(_timePeak[1]-_timeSigma[1]*_channelInfo.nSigmaTimeCut,_channelInfo.lowerTimeCut) ): _isGoodBar = False
    if( _time[_channelInfo.timech_id[1]] > min(_timePeak[1]+_timeSigma[1]*_channelInfo.nSigmaTimeCut,_channelInfo.upperTimeCut) ): _isGoodBar = False
    if( _time[_channelInfo.timech_id[2]] < max(_timePeak[2]-_timeSigma[2]*_channelInfo.nSigmaTimeCut,_channelInfo.lowerTimeCut) ): _isGoodBar = False
    if( _time[_channelInfo.timech_id[2]] > min(_timePeak[2]+_timeSigma[2]*_channelInfo.nSigmaTimeCut,_channelInfo.upperTimeCut) ): _isGoodBar = False

    #if _analysisStage == 3: #drop condition for final stage
    return _isGoodBar    


def calculatePeakPositionMIP( _chain, _ch, _timeAlgo, _outputDir ):
    """(1) first event loop - to calculate mip peak position"""

    # define histograms
    h_beamXY = TH2F("h_beamXY","h_beamXY", 100 ,_ch.minX, _ch.maxX, 100, _ch.minY, _ch.maxY)
    h_beamX  = TH1F("h_beamX", "h_beamX", 100, _ch.minX, _ch.maxX)
    h_beamY  = TH1F("h_beamY", "h_beamY", 100, _ch.minY, _ch.maxY)
  
    h_amp        = [TH1F() for i in range(0,_ch.NCH)]
    h_amp_cut    = [TH1F() for i in range(0,_ch.NCH)]
    p2_amp_vs_XY = [TProfile2D() for i in range(0,_ch.NCH)]
    
    h_time       = [TH1F() for i in range(0,_ch.NCH)]


    # make histograms
    for _iCh in range(0,_ch.NCH):
        h_amp[_iCh]        = TH1F('h_amp_'+_ch.namech[_iCh],'h_amp_'+_ch.namech[_iCh],_ch.nAmpBins,_ch.ampMin,_ch.ampMax);
        h_amp_cut[_iCh]    = TH1F('h_amp_cut_'+_ch.namech[_iCh],'h_amp_cut_'+_ch.namech[_iCh],_ch.nAmpBins,_ch.ampMin,_ch.ampMax);    
        p2_amp_vs_XY[_iCh] = TProfile2D('p2_amp_vs_XY_'+_ch.namech[_iCh],'p2_amp_vs_XY_'+_ch.namech[_iCh],100,_ch.minX,_ch.maxX,100,_ch.minY,_ch.maxY);
        h_time[_iCh]       = TH1F('h_time_'+_ch.namech[_iCh],'h_time_'+_ch.namech[_iCh],_ch.nTimeBins,_ch.minTime,_ch.maxTime);


    # event loop
    nSelectedEntries = 0;
    chainCounter = 0
    for entry in _chain:
        if chainCounter%5000 ==0:
            print( '>>> 1st loop: reading entry {0}/{1}'.format(chainCounter, _chain.GetEntries()) )
        chainCounter += 1

        myX = _chain.x_dut[0]
        myY = _chain.y_dut[0]
        
        # check if tracking good
        if not eventHasGoodHit( 0, _timeAlgo, myX, myY):
            continue
    
        h_beamX.Fill( myX )
        h_beamY.Fill( myY )
        h_beamXY.Fill( myX,myY )

        for _iCh in range(0, _ch.NCH):
            h_amp[_iCh].Fill( _chain.amp[_ch.ampch_id[_iCh]] );    
            p2_amp_vs_XY[_iCh].Fill( myX, myY, _chain.amp[_ch.ampch_id[_iCh]] );
            
            # check if tracking good + cut on BS
            if not eventHasGoodHit( 1, _timeAlgo, myX, myY, _chain, _ch):
                continue
      
            #get branch for user-specified time algorithm
            time = getTimeAlgorithmBranch( _chain, _timeAlgo )

            h_time[_iCh].Fill( time[ _ch.timech_id[_iCh]] )
            
            if (_chain.amp[_ch.ampch_id[_iCh]] > _ch.ampmin_cut[_iCh]) and (_chain.amp[_ch.ampch_id[_iCh]] < _ch.ampmax_cut[_iCh]) :
                h_amp_cut[_iCh].Fill( _chain.amp[_ch.ampch_id[_iCh]] );
                
        nSelectedEntries += 1
    
    print('>>> 1st loop: selected entries {0}'.format(nSelectedEntries))

    # ---------------------------- drawing beam histos ----------------------------

    c_beamXY = TCanvas("c_beamXY","c_beamXY",500,500)
    c_beamXY.cd()
    h_beamXY.SetStats(0)
    h_beamXY.SetTitle(";X [mm]; Y[mm];entries")
    h_beamXY.Draw("COLZ")

    printCanvas( c_beamXY, _outputDir, _timeAlgo)

    x_BS_min = TLine(_ch.centerX - _ch.BSX, _ch.minY, _ch.centerX - _ch.BSX, _ch.maxY)
    x_BS_max = TLine(_ch.centerX + _ch.BSX, _ch.minY, _ch.centerX + _ch.BSX, _ch.maxY)
    y_BS_min = TLine(_ch.minX, _ch.centerY - _ch.BSY, _ch.maxX, _ch.centerY-_ch.BSY)
    y_BS_max = TLine(_ch.minX, _ch.centerY + _ch.BSY, _ch.maxX, _ch.centerY+_ch.BSY)

    #print(_ch.centerX - _ch.BSX, _ch.minY, _ch.centerX - _ch.BSX, _ch.maxY)
    #print(_ch.minX, _ch.centerY + _ch.BSY, _ch.maxX, _ch.centerY + _ch.BSY)

    x_BS_min.SetLineColor(kRed)
    x_BS_min.SetLineWidth(2)
    x_BS_max.SetLineColor(kRed)
    x_BS_max.SetLineWidth(2)
    y_BS_min.SetLineColor(kRed)
    y_BS_min.SetLineWidth(2)
    y_BS_max.SetLineColor(kRed)
    y_BS_max.SetLineWidth(2)

    c_amp_vs_XY = TCanvas("c_amp_vs_XY","c_amp_vs_XY",1000,500*((_ch.NCH-1))/2)
    c_amp_vs_XY.Divide(2,(_ch.NCH-1)/2)
    for _iCh in range(1,_ch.NCH):
        c_amp_vs_XY.cd(_iCh)
        p2_amp_vs_XY[_iCh].Draw("COLZ")
        p2_amp_vs_XY[_iCh].SetStats(0)
        p2_amp_vs_XY[_iCh].SetTitle(";X [mm];Y [mm];amplitude [mV]")
        gPad.SetLogz()
        
        x_BS_min.Draw("same")
        x_BS_max.Draw("same")
        y_BS_min.Draw("same")
        y_BS_max.Draw("same")
        
        _ch.latexch[_iCh].Draw("same")        

    printCanvas( c_amp_vs_XY, _outputDir, _timeAlgo)
  
    # ---------------------------- fitting and drawing mip peak ----------------------------
    mip_peak     = [0 for i in range(0,_ch.NCH)]
    mip_peak_err = [0 for i in range(0,_ch.NCH)]
  
    c_amp = TCanvas("c_amp","c_amp",1000,500*(_ch.NCH+1)/2)
    c_amp.Divide(2,(_ch.NCH+1)/2)
    l_peakAmp = [TLatex() for i in range(0,_ch.NCH)]

    for _iCh in range(0,_ch.NCH):

        c_amp.cd(_iCh+1)
        gPad.SetLogy()

        h_amp[_iCh].SetStats(0)
        h_amp[_iCh].SetLineColor(kBlack)
        h_amp[_iCh].SetLineWidth(2)
        h_amp[_iCh].Draw()
        h_amp[_iCh].SetTitle(";max. amp [mV];entries");
        h_amp_cut[_iCh].SetFillColor(kOrange-9)
        h_amp_cut[_iCh].SetLineColor(kBlack)
        h_amp_cut[_iCh].SetLineWidth(2)
        h_amp_cut[_iCh].Draw("same")

        if _iCh > 0:
            fitMipPeak = TF1 ("fitMipPeak","landau",_ch.ampmin_cut[_iCh],_ch.ampmax_cut[_iCh])
            fitMipPeak.SetParameter(1,h_amp_cut[_iCh].GetBinCenter(int(h_amp_cut[_iCh].GetMaximum())))
            fitMipPeak.SetNpx(10000)
            fitMipPeak.SetLineColor(kRed)
            fitMipPeak.SetLineWidth(2)
            h_amp_cut[_iCh].Fit(fitMipPeak,"SQR")
            mip_peak[_iCh] = fitMipPeak.GetParameter(1)
            mip_peak_err[_iCh] = fitMipPeak.GetParError(1)
            #print( 'Fit Peak [{0}]: {1} +/- {2} mV'.format( _iCh, round(mip_peak[_iCh],2), round(mip_peak_err[_iCh],2)) )
            #txtOutputFitInfo << Form("Angle: %d, ", angleScan) << namech[_iCh].c_str() <<", amp peak[" << _iCh << "] = " << Form("%.2f #pm %.2f mV", mip_peak[_iCh], mip_peak_err[_iCh]) << std::endl
            print( 'Fit Peak [{0}]: {1} +/- {2} mV'.format( _ch.namech[_iCh], round(mip_peak[_iCh],2), round(mip_peak_err[_iCh],2)) )
            #outputFitInfo.Write( 'Fit Peak [{0}]: {1} +/- {2} mV'.format( _iCh, round(mip_peak[_iCh],2), round(mip_peak_err[_iCh],2)) )
            
            #TLine* lowcut = new TLine(std::max(rel_amp_cut_low*mip_peak[_iCh],ampmin_cut[_iCh]),0.,std::max(rel_amp_cut_low*mip_peak[_iCh],ampmin_cut[_iCh]),h_amp[_iCh].GetMaximum())
            #TLine* higcut = new TLine(std::min(rel_amp_cut_hig*mip_peak[_iCh],ampmax_cut[_iCh]),0.,std::min(rel_amp_cut_hig*mip_peak[_iCh],ampmax_cut[_iCh]),h_amp[_iCh].GetMaximum())
            #lowcut.Draw("same")
            #higcut.Draw("same")
            
            l_peakAmp[_iCh] = TLatex(0.45, 0.8, 'Fit Peak: {0} #pm {1} mV'.format( round(mip_peak[_iCh],2), round(mip_peak_err[_iCh],2)))
            l_peakAmp[_iCh].SetNDC()
            l_peakAmp[_iCh].SetTextFont(42)
            l_peakAmp[_iCh].SetTextSize(0.03)
            l_peakAmp[_iCh].Draw("same")
            
        _ch.latexch[_iCh].Draw("same")

        c_amp.Update()

    c_amp.Update()
    #c_amp.Print("c_amp.png")
    printCanvas( c_amp, _outputDir, _timeAlgo)

    # ---------------------------- fitting and drawing time peak ----------------------------
    fitTimePeak = TF1("fitTimePeak","gaus",_ch.lowerTimeCut,_ch.upperTimeCut)
    time_peak   = [0 for i in range(0,_ch.NCH)]
    time_sigma  = [0 for i in range(0,_ch.NCH)]
    
    c_time = TCanvas("c_time","c_time",1000,500*(_ch.NCH+1)/2)
    c_time.Divide(2,(_ch.NCH+1)/2)
  
    for _iCh in range(0,_ch.NCH):

        c_time.cd(_iCh+1)
        gPad.SetLogy()

        h_time[_iCh].SetStats(0)
        h_time[_iCh].SetFillColor(kOrange-9)
        h_time[_iCh].SetLineColor(kBlack)
        h_time[_iCh].Draw()
        h_time[_iCh].SetTitle(";time[ns];entries")
    
        fitTimePeak.SetParameters(1.,20.,5.)
        fitTimePeak.SetNpx(10000)
        fitTimePeak.SetLineColor(kRed)
        fitTimePeak.SetLineWidth(2)
        h_time[_iCh].Fit(fitTimePeak,"QRS")
        time_peak[_iCh] = fitTimePeak.GetParameter(1)
        time_sigma[_iCh] = fitTimePeak.GetParameter(2)
        print( 'Time Peak [{0}]: {1} +/- {2} ns'.format( _iCh, round(time_peak[_iCh],2), round(time_sigma[_iCh],2)) )
        #outputFitInfo.Write( 'Time Peak [{0}]: {1} +/- {2} ns'.format( _iCh, round(time_peak[_iCh],2), round(time_sigma[_iCh],2)) )
        #double quantileTimePeak = GetEffSigma( h_time[iCh] )
        #std::cout << "(Quantile) time_peak[" << iCh << "] = " << quantileTimePeak <<  std::endl            

        #lowcut = new TLine(std::max(time_peak[iCh]-time_sigma[iCh]*nSigmaTimeCut,lowerTimeCut),0,std::max(time_peak[iCh]-time_sigma[iCh]*nSigmaTimeCut,lowerTimeCut),h_time[iCh].GetMaximum())
        #higcut = new TLine(std::min(time_peak[iCh]+time_sigma[iCh]*nSigmaTimeCut,upperTimeCut),0,std::min(time_peak[iCh]+time_sigma[iCh]*nSigmaTimeCut,upperTimeCut),h_time[iCh].GetMaximum())
        #lowcut.Draw("same")
        #higcut.Draw("same")
    
        l_time = TLatex(0.45, 0.8, 'Fit Time: {0} #pm {1} ns'.format( round(time_peak[_iCh],2), round(time_sigma[_iCh],2)))
        l_time.SetNDC()
        l_time.SetTextFont(42)
        l_time.SetTextSize(0.03)
        _ch.latexch[_iCh].Draw("same")

    c_time.Update()
    #c_time.Print("c_time.png")
    printCanvas( c_time, _outputDir, _timeAlgo)
    
    return time_peak, time_sigma, mip_peak, mip_peak_err, h_amp_cut


def calculateAmpWalkCorrection( _chain, _ch, _timeAlgo, _outputDir, _timePeak, _timeSigma, _mipPeak ):
    """(2) second loop events -  to calculate amp-walk correction"""
    
    # define histograms
    h_deltat         = [TH1F() for i in range(0,_ch.NCH)]
    h2_deltat_vs_amp = [TH2F() for i in range(0,_ch.NCH)]
    p_deltat_vs_amp  = [TProfile() for i in range(0,_ch.NCH)]
    
    for _iCh in range(0,_ch.NCH):
        h_deltat[_iCh] = TH1F( 'h_deltat_{0}'.format(_ch.namech[_iCh]), 'h_deltat_{0}'.format(_ch.namech[_iCh]),_ch.nDeltatBins,_ch.minDeltat, _ch.maxDeltat)
        h2_deltat_vs_amp[_iCh] = TH2F( 'h2_deltat_vs_amp_{0}'.format(_ch.namech[_iCh]), 'h2_deltat_vs_amp_{0}'.format(_ch.namech[_iCh]), _ch.nAmpBins, _ch.ampMin, _ch.ampMax, _ch.nDeltatBins, _ch.minDeltat, _ch.maxDeltat);
        p_deltat_vs_amp[_iCh] = TProfile( 'p_deltat_vs_amp_{0}'.format(_ch.namech[_iCh]), 'p_deltat_vs_amp_{0}'.format(_ch.namech[_iCh]), _ch.nAmpBins, _ch.ampMin, _ch.ampMax)
        
  
    # event loop
    nSelectedEntries = 0;
    nCountedEntries = 0;
    chainCounter = 0
    for entry in _chain:
        if chainCounter%5000 ==0:
            print( '>>> 2nd loop: reading entry {0}/{1}'.format(chainCounter, _chain.GetEntries()) )
        chainCounter += 1
  
        myX = _chain.x_dut[0]
        myY = _chain.y_dut[0]
    
        # check if tracking good + cut on BS + cut on MCP
        if not eventHasGoodHit( 2, _timeAlgo, myX, myY, _chain, _ch, _timePeak, _timeSigma):
            continue

        # set time algorithm
        time = getTimeAlgorithmBranch( _chain, _timeAlgo )
        time_ref = _chain.gaus_mean[_ch.timech_id[0]]

        # loop over channels
        for _iCh in range(0,_ch.NCH):
            if( _chain.amp[_ch.ampch_id[_iCh]] > max(_mipPeak[_iCh]*_ch.rel_amp_cut_low,_ch.ampmin_cut[_iCh]) and
                _chain.amp[_ch.ampch_id[_iCh]] < min(_mipPeak[_iCh]*_ch.rel_amp_cut_hig,_ch.ampmax_cut[_iCh]) and
                time[_ch.timech_id[_iCh]] > max(_timePeak[_iCh]-_timeSigma[_iCh]*_ch.nSigmaTimeCut, _ch.lowerTimeCut) and
                time[_ch.timech_id[_iCh]] < min(_timePeak[_iCh]+_timeSigma[_iCh]*_ch.nSigmaTimeCut, _ch.upperTimeCut) ):
                
                h_deltat[_iCh].Fill( time[_ch.timech_id[_iCh]]-time_ref )
                h2_deltat_vs_amp[_iCh].Fill( _chain.amp[_ch.ampch_id[_iCh]],time[_ch.timech_id[_iCh]]-time_ref )
                p_deltat_vs_amp[_iCh].Fill( _chain.amp[_ch.ampch_id[_iCh]],time[_ch.timech_id[_iCh]]-time_ref )
                nCountedEntries += 1
        nSelectedEntries+= 1
  
    print( '>>> 2nd loop: selected entries {0}'.format(nSelectedEntries))
    print( '>>> 2nd loop: counted entries {0}'.format(nCountedEntries))

    # ---------------------------- fitting and drawing time walk correction ----------------------------
    fitAmpCorr = [TF1() for i in range(0,_ch.NCH)]

    c_time_vs_amp = TCanvas("c_time_vs_amp","c_time_vs_amp",1000,500*((_ch.NCH-1))/2)
    c_time_vs_amp.Divide(2,(_ch.NCH-1)/2)
  
    for _iCh in range(1,_ch.NCH):
        c_time_vs_amp.cd(_iCh)
    
        h2_deltat_vs_amp[_iCh].SetStats(0)
        h2_deltat_vs_amp[_iCh].GetYaxis().SetRangeUser(h_deltat[_iCh].GetMean()-5.*h_deltat[_iCh].GetRMS(),h_deltat[_iCh].GetMean()+5.*h_deltat[_iCh].GetRMS())
        h2_deltat_vs_amp[_iCh].SetTitle(";#Deltat [ns]; max. amplitude [mV]")
        h2_deltat_vs_amp[_iCh].Draw("COLZ")
        p_deltat_vs_amp[_iCh].SetMarkerStyle(20)
        p_deltat_vs_amp[_iCh].SetMarkerSize(0.7)
        p_deltat_vs_amp[_iCh].SetMarkerColor(kMagenta)
        p_deltat_vs_amp[_iCh].Draw("same")
    
        #fitAmpCorr[_iCh] = new TF1(Form("fitAmpCorr_%s", namech[_iCh].c_str()),"[0]*log([1]*x)+[2]",0.,1000.)
        #fitAmpCorr[_iCh].SetParameters(-0.3,6e-13,-10.)
        fitAmpCorr[_iCh] = TF1( 'fitAmpCorr_{0}'.format(_ch.namech[_iCh]),"pol4",0.,1000.)
        #fitAmpCorr[_iCh] = new TF1(Form("fitAmpCorr_%s", namech[_iCh].c_str()),"pol4",0.,200.)

        p_deltat_vs_amp[_iCh].Fit(fitAmpCorr[_iCh],"QNRS+")
        fitAmpCorr[_iCh].SetLineColor(kMagenta)
        fitAmpCorr[_iCh].Draw("same")
        
        _ch.latexch[_iCh].Draw("same")

    c_time_vs_amp.Update()
    printCanvas( c_time_vs_amp, _outputDir, _timeAlgo)
    
    return fitAmpCorr

def applyAmpWalkCorrection( _chain, _ch, _timeAlgo, _outputDir, _timePeak, _timeSigma, _mipPeak, _hAmpCut, _ampWalkCorr ):
    """(3) third loop events --> to apply amp-walk correction"""

    # define histograms
    h_deltat_avg = TH1F('h_deltat_avg','h_deltat_avg',6000, _ch.minDeltat, _ch.maxDeltat)
    h_deltat_avg_ampCorr = TH1F('h_deltat_avg_ampCorr','h_deltat_avg_ampCorr',6000, _ch.minDeltat, _ch.maxDeltat)
    
    h_deltat_avg_ampCorr_comb = TH1F('h_deltat_avg_ampCorr_comb','h_deltat_avg_ampCorr_comb',6000, _ch.minDeltat, _ch.maxDeltat)
    
    h_deltat_left  = TH1F('h_deltat_left','h_deltat_left',6000, _ch.minDeltat, _ch.maxDeltat)
    h_deltat_right = TH1F('h_deltat_right','h_deltat_right',6000, _ch.minDeltat, _ch.maxDeltat)
    h_deltat_diff  = TH1F('h_deltat_diff','h_deltat_diff',6000, _ch.minDeltat, _ch.maxDeltat)
    
    h_deltat_left_ampCorr  = TH1F('h_deltat_left_ampCorr','h_deltat_left_ampCorr',6000, _ch.minDeltat, _ch.maxDeltat)
    h_deltat_right_ampCorr = TH1F('h_deltat_right_ampCorr','h_deltat_right_ampCorr',6000, _ch.minDeltat, _ch.maxDeltat)
    h_deltat_diff_ampCorr  = TH1F('h_deltat_diff_ampCorr','h_deltat_diff_ampCorr',6000, _ch.minDeltat, _ch.maxDeltat)
    
    p_left_vs_X  = TProfile('p_left_vs_X','p_left_vs_X',400, _ch.minX, _ch.maxX)
    p_right_vs_X = TProfile('p_right_vs_X','p_right_vs_X',400, _ch.minX, _ch.maxX)  
    p_avg_vs_X   = TProfile('p_avg_vs_X','p_avg_vs_X',400, _ch.minX, _ch.maxX)
    p_diff_vs_X  = TProfile('p_diff_vs_X','p_diff_vs_X',400, _ch.minX, _ch.maxX)
    
    p_left_vs_diff  = TProfile('p_left_vs_diff','p_left_vs_diff',6000,-10,10)
    p_right_vs_diff = TProfile('p_right_vs_diff','p_right_vs_diff',6000,-10,10)  
    p_avg_vs_diff   = TProfile('p_avg_vs_diff','p_avg_vs_diff',6000,-10,10)
    
    h_deltat_avg_ampCorr_posCutX = [TH1F() for i in range(0,_ch.NPOSCUTSX)]
    h_deltat_left_ampCorr_posCutX = [TH1F() for i in range(0,_ch.NPOSCUTSX)]
    h_deltat_right_ampCorr_posCutX = [TH1F() for i in range(0,_ch.NPOSCUTSX)]
     
    h_deltat_avg_ampCorr_posCutY = [TH1F() for i in range(0,_ch.NPOSCUTSY)] 
    h_deltat_left_ampCorr_posCutY = [TH1F() for i in range(0,_ch.NPOSCUTSY)] 
    h_deltat_right_ampCorr_posCutY = [TH1F() for i in range(0,_ch.NPOSCUTSY)]
    
    h_deltat_avg_ampCorr_posCutXY = [ [TH1F() for i in range(0,_ch.NPOSCUTSY)] for i in range(0,_ch.NPOSCUTSX)]

    for _iPosCut in range(0,_ch.NPOSCUTSX):
        h_deltat_avg_ampCorr_posCutX[_iPosCut]   = TH1F('h_deltat_ampCorr_posCut_{0}'.format(_iPosCut), 'h_deltat_ampCorr_posCut_{0}'.format(_iPosCut), 3000, _ch.minDeltat, _ch.maxDeltat)
        h_deltat_left_ampCorr_posCutX[_iPosCut]  = TH1F('h_deltat_left_ampCorr_posCut_{0}'.format(_iPosCut), 'h_deltat_left_ampCorr_posCut_{0}'.format(_iPosCut), 3000, _ch.minDeltat, _ch.maxDeltat)
        h_deltat_right_ampCorr_posCutX[_iPosCut] = TH1F('h_deltat_right_ampCorr_posCut_{0}'.format(_iPosCut), 'h_deltat_right_ampCorr_posCut_{0}'.format(_iPosCut), 3000, _ch.minDeltat, _ch.maxDeltat)
        for _iPosCutY in range(0,_ch.NPOSCUTSY):
            h_deltat_avg_ampCorr_posCutXY[_iPosCut][_iPosCutY]  = TH1F('h_deltat_ampCorr_posCutXY_{0}_{1}'.format(_iPosCut,_iPosCutY), 'h_deltat_ampCorr_posCutXY_{0}_{1}'.format(_iPosCut,_iPosCutY), 3000, _ch.minDeltat, _ch.maxDeltat)

    for _iPosCut in range(0,_ch.NPOSCUTSY):
        h_deltat_avg_ampCorr_posCutY[_iPosCut]   = TH1F('h_deltat_ampCorr_posCutY_{0}'.format(_iPosCut), 'h_deltat_ampCorr_posCutY_{0}'.format(_iPosCut), 3000, _ch.minDeltat, _ch.maxDeltat)
        h_deltat_left_ampCorr_posCutY[_iPosCut]  = TH1F('h_deltat_left_ampCorr_posCutY_{0}'.format(_iPosCut), 'h_deltat_left_ampCorr_posCutY_{0}'.format(_iPosCut), 3000, _ch.minDeltat, _ch.maxDeltat)
        h_deltat_right_ampCorr_posCutY[_iPosCut] = TH1F('h_deltat_right_ampCorr_posCutY_{0}'.format(_iPosCut), 'h_deltat_right_ampCorr_posCutY_{0}'.format(_iPosCut), 3000, _ch.minDeltat, _ch.maxDeltat)

    h_deltat_avg_ampCorr_BSCut = [TH1F() for i in range(0,_ch.NBSCUTS)]
    h_deltat_left_ampCorr_BSCut = [TH1F() for i in range(0,_ch.NBSCUTS)]
    h_deltat_right_ampCorr_BSCut = [TH1F() for i in range(0,_ch.NBSCUTS)]
    for _iBSCut in range(0,_ch.NBSCUTS):
        h_deltat_avg_ampCorr_BSCut[_iBSCut]   = TH1F('h_deltat_avg_ampCorr_BSCut_{0}'.format(_iBSCut), 'h_deltat_avg_ampCorr_BSCut_{0}'.format(_iBSCut), _ch.nDeltatBins, _ch.minDeltat, _ch.maxDeltat);
        h_deltat_left_ampCorr_BSCut[_iBSCut]  = TH1F('h_deltat_left_ampCorr_BSCut_{0}'.format(_iBSCut), 'h_deltat_left_ampCorr_BSCut_{0}'.format(_iBSCut), _ch.nDeltatBins, _ch.minDeltat, _ch.maxDeltat);
        h_deltat_right_ampCorr_BSCut[_iBSCut] = TH1F('h_deltat_right_ampCorr_BSCut_{0}'.format(_iBSCut), 'h_deltat_right_ampCorr_BSCut_{0}'.format(_iBSCut), _ch.nDeltatBins, _ch.minDeltat, _ch.maxDeltat);


    # event loop
    nSelectedEntries = 0
    chainCounter = 0
    for entry in _chain:
        if chainCounter%5000 ==0:
            print( '>>> 3rd loop: reading entry {0}/{1}'.format(chainCounter, _chain.GetEntries()) )
        chainCounter += 1
     
        myX = _chain.x_dut[0]
        myY = _chain.y_dut[0]
    
        # check if tracking good + cut on BS + cut on MCP + bar selection
        if not eventHasGoodHit( 3, _timeAlgo, myX, myY, _chain, _ch, _timePeak, _timeSigma, _mipPeak):
            continue

        # set time algorithm
        time = getTimeAlgorithmBranch( _chain, _timeAlgo )
        time_ref = _chain.gaus_mean[_ch.timech_id[0]]

        amp1 = _chain.amp[_ch.ampch_id[1]]
        amp2 = _chain.amp[_ch.ampch_id[2]]
        time1 = time[_ch.timech_id[1]]
        time2 = time[_ch.timech_id[2]]
        time1_ampCorr = time1 - _ampWalkCorr[1].Eval(amp1) + _ampWalkCorr[1].Eval(_hAmpCut[1].GetMean())
        time2_ampCorr = time2 - _ampWalkCorr[2].Eval(amp2) + _ampWalkCorr[2].Eval(_hAmpCut[2].GetMean())
    
        deltat_avg = 0.5*(time1+time2) - time_ref
        deltat_avg_ampCorr = 0.5*(time1_ampCorr+time2_ampCorr) - time_ref
    
        h_deltat_left.Fill( time1 - time_ref )
        h_deltat_right.Fill( time2 - time_ref )
        h_deltat_diff.Fill( time2 - time1 )
        h_deltat_avg.Fill( deltat_avg )
        
        h_deltat_left_ampCorr.Fill( time1_ampCorr - time_ref )
        h_deltat_right_ampCorr.Fill( time2_ampCorr - time_ref )      
        h_deltat_diff_ampCorr.Fill( time2_ampCorr - time1_ampCorr )    
        h_deltat_avg_ampCorr.Fill( deltat_avg_ampCorr )
        
        # filling plots vs position
        p_left_vs_X.Fill( myX,time1_ampCorr-time_ref )
        p_right_vs_X.Fill( myX,time2_ampCorr-time_ref )
        p_avg_vs_X.Fill( myX,0.5*(time1_ampCorr+time2_ampCorr)-time_ref )
        p_diff_vs_X.Fill( myX,time2_ampCorr-time1_ampCorr )
        
        
        # filling plots vd t_diff
        p_left_vs_diff.Fill( time2_ampCorr-time1_ampCorr,time1_ampCorr-time_ref)
        p_right_vs_diff.Fill( time2_ampCorr-time1_ampCorr,time2_ampCorr-time_ref)
        p_avg_vs_diff.Fill( time2_ampCorr-time1_ampCorr,0.5*(time1_ampCorr+time2_ampCorr)-time_ref)
        
    
        # X dependency  
        for _iPosCutX in range(0,_ch.NPOSCUTSX):
            stepX = (_ch.upperPosCutX-_ch.lowerPosCutX) / _ch.NPOSCUTSX

            if( myX > _ch.lowerPosCutX + _iPosCutX*stepX and  myX < _ch.lowerPosCutX + (_iPosCutX+1)*stepX ):
                h_deltat_avg_ampCorr_posCutX[_iPosCutX].Fill( deltat_avg_ampCorr )
                h_deltat_left_ampCorr_posCutX[_iPosCutX].Fill( time1_ampCorr - time_ref )
                h_deltat_right_ampCorr_posCutX[_iPosCutX].Fill( time2_ampCorr - time_ref )
                
            for _iPosCutY in range(0,_ch.NPOSCUTSY):
                stepY = (_ch.upperPosCutY-_ch.lowerPosCutY) / _ch.NPOSCUTSY
                if( myX > _ch.lowerPosCutX + _iPosCutX*stepX and  myX < _ch.lowerPosCutX + (_iPosCutX+1)*stepX and myY >_ch.lowerPosCutY + _iPosCutY*stepY and myY < _ch.lowerPosCutY+(_iPosCut+1)*stepY):
                    h_deltat_avg_ampCorr_posCutXY[_iPosCutX][_iPosCutY].Fill( deltat_avg_ampCorr )
        
        # Y dependency
        for _iPosCutY in range(0,_ch.NPOSCUTSY):
            stepY = (_ch.upperPosCutY-_ch.lowerPosCutY) / _ch.NPOSCUTSY

            if( myY >_ch.lowerPosCutY + _iPosCutY*stepY and myY < _ch.lowerPosCutY+(_iPosCut+1)*stepY ):
                h_deltat_avg_ampCorr_posCutY[_iPosCutY].Fill( deltat_avg_ampCorr )
                h_deltat_left_ampCorr_posCutY[_iPosCutY].Fill( time1_ampCorr - time_ref )
                h_deltat_right_ampCorr_posCutY[_iPosCutY].Fill( time2_ampCorr - time_ref )
    
        
    
        for _iBSCut in range(0, _ch.NBSCUTS):
            if ( abs(myX-_ch.centerX) < _ch.BScut[_iBSCut] ):
                h_deltat_avg_ampCorr_BSCut[_iBSCut].Fill( 0.5*(time1_ampCorr+time2_ampCorr) - time_ref )
                h_deltat_left_ampCorr_BSCut[_iBSCut].Fill( time1_ampCorr - time_ref )
                h_deltat_right_ampCorr_BSCut[_iBSCut].Fill( time2_ampCorr - time_ref )
    
    
        nSelectedEntries += 1

    print('>>> 3rd loop: selected entries {0}'.format(nSelectedEntries))
   
    # ---------------------------- position plots ----------------------------

    c_time_vs_X = TCanvas('c_time_vs_X','c_time_vs_X',1000,500)
    c_time_vs_X.Divide(2,1)
    
    c_time_vs_X.cd(1)
    gPad.SetGridy()
    
    fitFunc_corrX = TF1("fitFunc_corrX","pol2",_ch.centerX-2.*_ch.BSX,_ch.centerX+2.*_ch.BSX)
    p_avg_vs_X.Fit(fitFunc_corrX,"QNR")
    
    p_avg_vs_X.GetXaxis().SetRangeUser(_ch.centerX-2.*_ch.BSX,_ch.centerX+2.*_ch.BSX)
    p_avg_vs_X.GetYaxis().SetRangeUser(h_deltat_avg.GetMean()-15.*h_deltat_avg.GetRMS(),
                                         h_deltat_avg.GetMean()+15.*h_deltat_avg.GetRMS())
    p_avg_vs_X.Draw()
    p_avg_vs_X.SetStats(0)
    p_avg_vs_X .SetTitle(";X [mm];#Deltat [ns]")
    fitFunc_corrX.Draw("same")
    p_avg_vs_X.SetMarkerStyle(21)
    p_left_vs_X.SetLineColor(kRed+1)
    p_left_vs_X.SetMarkerColor(kRed+1)
    p_left_vs_X.SetMarkerStyle(20)
    p_left_vs_X.Draw("same")
    p_right_vs_X.SetLineColor(kBlue+1)
    p_right_vs_X.SetMarkerColor(kBlue+1)
    p_right_vs_X.SetMarkerStyle(20)
    p_right_vs_X.Draw("same")
    p_diff_vs_X.SetLineColor(kYellow+2)
    p_diff_vs_X.SetMarkerColor(kYellow+2)
    p_diff_vs_X.SetMarkerStyle(22)
    p_diff_vs_X.SetLineStyle(7)
    p_diff_vs_X.Draw("same")
    
    leg = TLegend(0.71,0.73,0.86,0.93,'','brNDC')
    leg.SetBorderSize(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.03)
    leg.SetLineColor(1)
    leg.SetLineStyle(1)
    leg.SetLineWidth(1)
    leg.SetFillColor(0)  
    leg.AddEntry(p_left_vs_X,  "t_{left} - t_{MCP}", "lpe")     
    leg.AddEntry(p_right_vs_X, "t_{right} - t_{MCP}", "lpe")     
    leg.AddEntry(p_avg_vs_X,   "t_{avg} - t_{MCP}", "lpe")     
    leg.AddEntry(p_diff_vs_X,  "t_{left} - t_{right}", "lpe")     
    leg.Draw("same")
    
    c_time_vs_X.cd(2)
    gPad.SetGridy()
    
    fitFunc_corrDiff = TF1("fitFunc_corrDiff","pol2",h_deltat_diff.GetMean()-5.*h_deltat_diff.GetRMS(),h_deltat_diff.GetMean()+5.*h_deltat_diff.GetRMS())
    p_avg_vs_diff.Fit(fitFunc_corrDiff,"QNR")
    p_avg_vs_diff.GetXaxis().SetRangeUser(h_deltat_diff.GetMean()-5.*h_deltat_diff.GetRMS(), h_deltat_diff.GetMean()+5.*h_deltat_diff.GetRMS())
    p_avg_vs_diff.GetYaxis().SetRangeUser(h_deltat_avg.GetMean()-15.*h_deltat_avg.GetRMS(),  h_deltat_avg.GetMean()+15.*h_deltat_avg.GetRMS())
    p_avg_vs_diff.Draw()
    p_avg_vs_diff.SetStats(0)
    p_avg_vs_diff.SetTitle(";t_{left} - t_{right} [ns];#Deltat [ns]")
    fitFunc_corrDiff.Draw("same")
    p_avg_vs_diff.SetMarkerStyle(20)
    p_left_vs_diff.SetLineColor(kRed+1)
    p_left_vs_diff.SetMarkerColor(kRed+1)
    p_left_vs_diff.SetMarkerStyle(20)
    p_left_vs_diff.Draw("same")
    p_right_vs_diff.SetLineColor(kBlue+1)
    p_right_vs_diff.SetMarkerColor(kBlue+1)
    p_right_vs_diff.SetMarkerStyle(20)
    p_right_vs_diff.Draw("same")
    
    leg = TLegend(0.71,0.78,0.86,0.93,'','brNDC')
    leg.SetBorderSize(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.03)
    leg.SetLineColor(1)
    leg.SetLineStyle(1)
    leg.SetLineWidth(1)
    leg.SetFillColor(0)
    leg.AddEntry(p_left_vs_diff,  "t_{left} - t_{MCP}", "lpe")     
    leg.AddEntry(p_right_vs_diff, "t_{right} - t_{MCP}", "lpe")     
    leg.AddEntry(p_avg_vs_diff,   "t_{avg} - t_{MCP}", "lpe")     
    leg.Draw("same")
    printCanvas( c_time_vs_X, _outputDir, _timeAlgo)
    
    # ---------------------------- compare left, right, sum time resolution ----------------------------
    c_timeRes_comp = TCanvas("c_timeRes_comp","c_timeRes_comp",1000,500)
    c_timeRes_comp.Divide(2,1)
    
    c_timeRes_comp.cd(1)
    h_deltat_avg.SetStats(0)
    h_deltat_avg.SetTitle(";#Deltat (no amp-walk corr.) [ns];entries")
    h_deltat_avg.SetLineColor(kBlack)
    h_deltat_avg.SetLineWidth(2)
    h_deltat_avg.GetXaxis().SetRangeUser(h_deltat_avg.GetMean()-15.*h_deltat_avg.GetRMS(), h_deltat_avg.GetMean()+15.*h_deltat_avg.GetRMS())
    h_deltat_avg.Draw()
    h_deltat_left.Draw("same")
    h_deltat_left.SetLineColor(kRed+1)
    h_deltat_left.SetLineWidth(2)
    h_deltat_right.Draw("same")
    h_deltat_right.SetLineColor(kBlue+1)
    h_deltat_right.SetLineWidth(2)
    
    fitdeltat_left = TF1("fitdeltat_left", "gaus", h_deltat_left.GetMean()-h_deltat_left.GetRMS()*2, h_deltat_left.GetMean()+h_deltat_left.GetRMS()*2)
    fitdeltat_left.SetLineColor(kRed+1)
    h_deltat_left.Fit(fitdeltat_left, "QR")
    fitdeltat_right = TF1("fitdeltat_right", "gaus", h_deltat_right.GetMean()-h_deltat_right.GetRMS()*2, h_deltat_right.GetMean()+h_deltat_right.GetRMS()*2)
    fitdeltat_right.SetLineColor(kBlue+1)
    h_deltat_right.Fit(fitdeltat_right, "QR")
    fitdeltat_avg = TF1("fitdeltat_avg", "gaus", h_deltat_avg.GetMean()-h_deltat_avg.GetRMS()*2, h_deltat_avg.GetMean()+h_deltat_avg.GetRMS()*2)
    fitdeltat_avg.SetLineColor(kBlack)
    h_deltat_avg.Fit(fitdeltat_avg, "QR")
    
    sigmaLeft  = np.sqrt( pow(fitdeltat_left.GetParameter(2),2)  - pow(_ch.sigma_ref,2) )    
    sigmaRight = np.sqrt( pow(fitdeltat_right.GetParameter(2),2) - pow(_ch.sigma_ref,2) )    
    sigmaAvg   = np.sqrt( pow(fitdeltat_avg.GetParameter(2),2)   - pow(_ch.sigma_ref,2) )    
 
    leg = TLegend(0.65,0.78,0.80,0.93,'',"brNDC")
    leg.SetBorderSize(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.03)
    leg.SetLineColor(1)
    leg.SetLineStyle(1)
    leg.SetLineWidth(1)
    leg.SetFillColor(0)
    leg.SetFillStyle(1000)
    leg.AddEntry(h_deltat_left,  '#sigma^{{left}}_{{t}} = {0} ps'.format(round(sigmaLeft*1000),1), "l")
    leg.AddEntry(h_deltat_right, '#sigma^{{right}}_{{t}} = {0} ps'.format(round(sigmaRight*1000),1), "l")
    leg.AddEntry(h_deltat_avg,   '#sigma^{{avg}}_{{t}} = {0} ps'.format(round(sigmaAvg*1000),1), "l")
    leg.Draw("same")
    
    c_timeRes_comp.cd(2)
    h_deltat_avg_ampCorr.SetStats(0)
    h_deltat_avg_ampCorr.SetTitle(";#Deltat (amp-walk corr.) [ns];entries")
    h_deltat_avg_ampCorr.SetLineColor(kBlack)
    h_deltat_avg_ampCorr.SetLineWidth(2)
    h_deltat_avg_ampCorr.GetXaxis().SetRangeUser(h_deltat_avg.GetMean()-15.*h_deltat_avg.GetRMS(), h_deltat_avg.GetMean()+15.*h_deltat_avg.GetRMS())
    h_deltat_avg_ampCorr.Draw()
    h_deltat_left_ampCorr.SetLineColor(kRed+1)
    h_deltat_left_ampCorr.SetLineWidth(2)
    h_deltat_left_ampCorr.Draw("same")
    h_deltat_right_ampCorr.SetLineColor(kBlue+1)
    h_deltat_right_ampCorr.SetLineWidth(2)
    h_deltat_right_ampCorr.Draw("same")
    
    fitdeltat_left_ampCorr = TF1("fitdeltat_left_ampCorr", "gaus", h_deltat_left_ampCorr.GetMean()-h_deltat_left_ampCorr.GetRMS()*2, h_deltat_left_ampCorr.GetMean()+h_deltat_left_ampCorr.GetRMS()*2)
    fitdeltat_left_ampCorr.SetLineColor(kRed+1)
    h_deltat_left_ampCorr.Fit(fitdeltat_left_ampCorr, "QR")
    fitdeltat_right_ampCorr = TF1("fitdeltat_right_ampCorr", "gaus", h_deltat_right_ampCorr.GetMean()-h_deltat_right_ampCorr.GetRMS()*2, h_deltat_right_ampCorr.GetMean()+h_deltat_right_ampCorr.GetRMS()*2)
    fitdeltat_right_ampCorr.SetLineColor(kBlue+1)
    h_deltat_right_ampCorr.Fit(fitdeltat_right_ampCorr, "QR")
    fitdeltat_avg_ampCorr = TF1("fitdeltat_avg_ampCorr", "gaus", h_deltat_avg_ampCorr.GetMean()-h_deltat_avg_ampCorr.GetRMS()*2, h_deltat_avg_ampCorr.GetMean()+h_deltat_avg_ampCorr.GetRMS()*2)
    fitdeltat_avg_ampCorr.SetLineColor(kBlack)
    h_deltat_avg_ampCorr.Fit(fitdeltat_avg_ampCorr, "QR")
    
    sigmaLeftCorr  = np.sqrt(pow(fitdeltat_left_ampCorr.GetParameter(2),2)  - pow(_ch.sigma_ref,2) )    
    sigmaRightCorr = np.sqrt(pow(fitdeltat_right_ampCorr.GetParameter(2),2) - pow(_ch.sigma_ref,2) )    
    sigmaAvgCorr   = np.sqrt(pow(fitdeltat_avg_ampCorr.GetParameter(2),2)   - pow(_ch.sigma_ref,2) )    
    sigmaLeftCorr_err  = fitdeltat_left_ampCorr.GetParError(2)
    sigmaRightCorr_err = fitdeltat_right_ampCorr.GetParError(2)
    sigmaAvgCorr_err   = fitdeltat_avg_ampCorr.GetParError(2)

    printCanvas( c_timeRes_comp, _outputDir, _timeAlgo )
    
    # ---------------------------- BS cut plots ----------------------------
    fitdeltat_ampCorr_BSCut_L = TF1("fitdeltat_ampCorr_BSCut_L", "gaus", h_deltat_left_ampCorr.GetMean() -h_deltat_left_ampCorr.GetRMS()*2,  h_deltat_left_ampCorr.GetMean() +h_deltat_left_ampCorr.GetRMS()*2)
    fitdeltat_ampCorr_BSCut_R = TF1("fitdeltat_ampCorr_BSCut_R", "gaus", h_deltat_right_ampCorr.GetMean()-h_deltat_right_ampCorr.GetRMS()*2, h_deltat_right_ampCorr.GetMean()+h_deltat_right_ampCorr.GetRMS()*2)
    fitdeltat_ampCorr_BSCut   = TF1("fitdeltat_ampCorr_BSCut",   "gaus", h_deltat_avg_ampCorr.GetMean()  -h_deltat_avg_ampCorr.GetRMS()*2,   h_deltat_avg_ampCorr.GetMean()  +h_deltat_avg_ampCorr.GetRMS()*2)
    
    gdeltat_vs_BS_L = TGraphErrors()
    gdeltat_vs_BS_R = TGraphErrors()
    gdeltat_vs_BS   = TGraphErrors()
    
    myPoint = 0
    for _iBSCut in range(0,_ch.NBSCUTS):
        #/left
        h_deltat_left_ampCorr_BSCut[_iBSCut].Fit(fitdeltat_ampCorr_BSCut_L, "QNR")
        tempSigma_L =  np.sqrt(pow(fitdeltat_ampCorr_BSCut_L.GetParameter(2),2) - pow(_ch.sigma_ref, 2))  
        sigmaErr_L = fitdeltat_ampCorr_BSCut_L.GetParError(2)
          
        #right
        h_deltat_right_ampCorr_BSCut[_iBSCut].Fit(fitdeltat_ampCorr_BSCut_R, "QNR")      
        tempSigma_R =  np.sqrt(pow(fitdeltat_ampCorr_BSCut_R.GetParameter(2),2) - pow(_ch.sigma_ref, 2))  
        sigmaErr_R = fitdeltat_ampCorr_BSCut_R.GetParError(2)
    
        #avg
        h_deltat_avg_ampCorr_BSCut[_iBSCut].Fit(fitdeltat_ampCorr_BSCut, "QNR")
        tempSigma =  np.sqrt(pow(fitdeltat_ampCorr_BSCut.GetParameter(2),2) - pow(_ch.sigma_ref, 2))  
        sigmaErr = fitdeltat_ampCorr_BSCut.GetParError(2)
    
        if (tempSigma>0 and tempSigma<0.5 and  h_deltat_avg_ampCorr_BSCut[_iBSCut].GetEntries()>30):
            gdeltat_vs_BS_L.SetPoint(myPoint, _ch.BScut[_iBSCut]*2, tempSigma_L)            
            gdeltat_vs_BS_L.SetPointError(myPoint,0, sigmaErr_L)    
            
            gdeltat_vs_BS_R.SetPoint(myPoint, _ch.BScut[_iBSCut]*2, tempSigma_R)            
            gdeltat_vs_BS_R.SetPointError(myPoint,0, sigmaErr_R)    
            
            gdeltat_vs_BS.SetPoint(myPoint, _ch.BScut[_iBSCut]*2, tempSigma)            
            gdeltat_vs_BS.SetPointError(myPoint,0, sigmaErr)    
            
            myPoint += 1


    c_timeRes_vs_BS = TCanvas("c_timeRes_vs_BS","c_timeRes_vs_BS",500,500)
    c_timeRes_vs_BS.cd()
    gPad.SetGridy()
    gPad.SetLogx()
    #hPad = (TH1F*)( gPad.DrawFrame(BScut[NBSCUTS-1],0.02,3.*BScut[0],0.1) )
    hPad = ( gPad.DrawFrame(_ch.BScut[_ch.NBSCUTS-1],0.02,3.*_ch.BScut[0],0.1) )
    hPad.SetTitle(";beam spot width [mm];#sigma_{t} [ns]")
    hPad.Draw()
    gdeltat_vs_BS.Draw("PLE,same")
    gdeltat_vs_BS.SetMarkerStyle(20)
    gdeltat_vs_BS_L.SetLineColor(kRed+1)
    gdeltat_vs_BS_L.SetMarkerColor(kRed+1)
    gdeltat_vs_BS_L.SetMarkerStyle(21)
    gdeltat_vs_BS_L.Draw("same LPE")
    gdeltat_vs_BS_R.SetLineColor(kBlue+1)
    gdeltat_vs_BS_R.SetMarkerColor(kBlue+1)
    gdeltat_vs_BS_R.SetMarkerStyle(21)
    gdeltat_vs_BS_R.Draw("same LPE")
    
    leg = TLegend(0.71,0.78,0.86,0.93,'',"brNDC")
    leg.SetBorderSize(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.03)
    leg.SetLineColor(1)
    leg.SetLineStyle(1)
    leg.SetLineWidth(1)
    leg.SetFillColor(0)  
    leg.AddEntry(gdeltat_vs_BS_L, 'left only', "lpe")     
    leg.AddEntry(gdeltat_vs_BS_R, 'right only', "lpe")     
    leg.AddEntry(gdeltat_vs_BS,   'avgrage', "lpe")     
    leg.Draw("same") 
    
    printCanvas( c_timeRes_vs_BS, _outputDir, _timeAlgo )
    # ---------------------------- plots in position bins ----------------------------

    g_timeRes_left_vs_X = TGraphErrors()
    g_timeRes_right_vs_X = TGraphErrors()
    g_timeRes_avg_vs_X = TGraphErrors()
    
    g_timeRes_left_vs_Y = TGraphErrors()
    g_timeRes_right_vs_Y = TGraphErrors()
    g_timeRes_avg_vs_Y   = TGraphErrors()
    
    h2_timeRes_avg_vs_XY = TH2F("h2_timeRes_avg_vs_XY","h2_timeRes_avg_vs_XY",_ch.NPOSCUTSX,_ch.lowerPosCutX,_ch.upperPosCutX,_ch.NPOSCUTSY,_ch.lowerPosCutY,_ch.upperPosCutY)
    
    # X position
    # avgrage
    myPoint = 0
    for _iPosCutX in range(0,_ch.NPOSCUTSX):
        stepX = ( _ch.upperPosCutX - _ch.lowerPosCutX ) / _ch.NPOSCUTSX
        
        fitdeltat_ampCorr_posCut = TF1("fitdeltat_avg_ampCorr_posCut","gaus", h_deltat_avg_ampCorr_posCutX[_iPosCutX].GetMean()-h_deltat_avg_ampCorr_posCutX[_iPosCutX].GetRMS()*2.,
                                       h_deltat_avg_ampCorr_posCutX[_iPosCutX].GetMean()+h_deltat_avg_ampCorr_posCutX[_iPosCutX].GetRMS()*2.)
        h_deltat_avg_ampCorr_posCutX[_iPosCutX].Fit(fitdeltat_ampCorr_posCut, "QNR")
        selPos = _ch.lowerPosCutX+(stepX*_iPosCutX)
        tempSigma =  np.sqrt( pow(fitdeltat_ampCorr_posCut.GetParameter(2),2) - pow(_ch.sigma_ref,2) )  
        sigmaErr = fitdeltat_ampCorr_posCut.GetParError(2)      
        if( tempSigma > 0. and tempSigma < 0.5 and h_deltat_avg_ampCorr_posCutX[_iPosCutX].GetEntries() > 20 ):
            g_timeRes_avg_vs_X.SetPoint(myPoint,selPos,tempSigma)            
            g_timeRes_avg_vs_X.SetPointError(myPoint,stepX/2/np.sqrt(12),sigmaErr)    
            myPoint += 1

    # left only
    myPoint = 0
    for _iPosCutX in range(0,_ch.NPOSCUTSX):
        stepX = ( _ch.upperPosCutX - _ch.lowerPosCutX ) / _ch.NPOSCUTSX
        fitdeltat_ampCorr_posCutL = TF1("fitdeltat_left_ampCorr_posCut","gaus", h_deltat_left_ampCorr.GetMean()-h_deltat_left_ampCorr.GetRMS()*2.,
                                    h_deltat_left_ampCorr.GetMean()+h_deltat_left_ampCorr.GetRMS()*2.)
        h_deltat_left_ampCorr_posCutX[_iPosCutX].Fit(fitdeltat_ampCorr_posCutL, "QNR")
        selPos = _ch.lowerPosCutX+(stepX*_iPosCutX)
        tempSigma =  np.sqrt(pow(fitdeltat_ampCorr_posCutL.GetParameter(2),2) - pow(_ch.sigma_ref, 2))
        sigmaErr = fitdeltat_ampCorr_posCutL.GetParError(2)
        if( tempSigma > 0. and tempSigma < 0.5 and h_deltat_left_ampCorr_posCutX[_iPosCutX].GetEntries() > 20 ):
            g_timeRes_left_vs_X.SetPoint(myPoint, selPos, tempSigma)
            g_timeRes_left_vs_X.SetPointError(myPoint,stepX/2, sigmaErr)
            myPoint += 1

    # right only
    myPoint = 0
    for _iPosCutX in range(0,_ch.NPOSCUTSX):
        stepX = ( _ch.upperPosCutX - _ch.lowerPosCutX ) / _ch.NPOSCUTSX
        fitdeltat_ampCorr_posCutR = TF1("fitdeltat_right_ampCorr_posCut","gaus", h_deltat_right_ampCorr.GetMean()-h_deltat_right_ampCorr.GetRMS()*2.,
                                       h_deltat_right_ampCorr.GetMean()+h_deltat_right_ampCorr.GetRMS()*2.)
        h_deltat_right_ampCorr_posCutX[_iPosCutX].Fit(fitdeltat_ampCorr_posCutR, "QNR")
        selPos = _ch.lowerPosCutX+(stepX*_iPosCutX)
        tempSigma =  np.sqrt(pow(fitdeltat_ampCorr_posCutR.GetParameter(2),2) - pow(_ch.sigma_ref, 2))
        sigmaErr = fitdeltat_ampCorr_posCutR.GetParError(2)
        if (tempSigma > 0. and tempSigma < 0.5 and  h_deltat_right_ampCorr_posCutX[_iPosCutX].GetEntries() > 20 ):
            g_timeRes_right_vs_X.SetPoint(myPoint,selPos,tempSigma)
            g_timeRes_right_vs_X.SetPointError(myPoint,stepX/2,sigmaErr)
            myPoint += 1


    c_timeRes_vs_X_Y = TCanvas("c_timeRes_vs_X_Y","c_timeRes_vs_X_Y",1000,500)
    c_timeRes_vs_X_Y.Divide(2,1)
    c_timeRes_vs_X_Y.cd(1)
    gPad.SetGridy()
    g_timeRes_avg_vs_X.GetYaxis().SetRangeUser(0, 0.2)
    g_timeRes_avg_vs_X.SetTitle(";X [mm];#sigma_{t} [ns]")
    g_timeRes_avg_vs_X.SetMarkerStyle(20)
    g_timeRes_avg_vs_X.Draw("ALPE")
    g_timeRes_left_vs_X.SetLineColor(kRed+1)
    g_timeRes_left_vs_X.SetMarkerColor(kRed+1)
    g_timeRes_right_vs_X.SetLineColor(kBlue+1)
    g_timeRes_right_vs_X.SetMarkerColor(kBlue+1)
    g_timeRes_left_vs_X.Draw("same LPE")
    g_timeRes_right_vs_X.Draw("same LPE")
    
    leg = TLegend(0.70,0.78,0.85,0.93,'',"brNDC")
    leg.SetBorderSize(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.03)
    leg.SetLineColor(1)
    leg.SetLineStyle(1)
    leg.SetLineWidth(1)
    leg.SetFillColor(0)
    leg.SetFillStyle(1000)
    leg.AddEntry(g_timeRes_left_vs_X,  "#sigma^{left}_{t}", "el")
    leg.AddEntry(g_timeRes_right_vs_X, "#sigma^{right}_{t}", "el")
    leg.AddEntry(g_timeRes_avg_vs_X,   "#sigma^{avg}_{t}",  "pl")
    leg.Draw("same")
   
    # Y position
    # average
    myPoint = 0
    for _iPosCutY in range(0,_ch.NPOSCUTSY):
        stepY = ( _ch.upperPosCutY - _ch.lowerPosCutY ) / _ch.NPOSCUTSY

        fitdeltat_ampCorr_posCut_avgY = TF1("fitdeltat_avg_ampCorr_posCut","gaus", h_deltat_avg_ampCorr_posCutY[_iPosCutY].GetMean()-h_deltat_avg_ampCorr_posCutY[_iPosCutY].GetRMS()*2.,
                                       h_deltat_avg_ampCorr_posCutY[_iPosCutY].GetMean()+h_deltat_avg_ampCorr_posCutY[_iPosCutY].GetRMS()*2.)
        h_deltat_avg_ampCorr_posCutY[_iPosCutY].Fit(fitdeltat_ampCorr_posCut_avgY, "QNR")
        selPos = _ch.lowerPosCutY+(stepY*_iPosCutY)
        tempSigma =  np.sqrt( pow(fitdeltat_ampCorr_posCut_avgY.GetParameter(2),2) - pow(_ch.sigma_ref,2) )  
        sigmaErr = fitdeltat_ampCorr_posCut_avgY.GetParError(2)      
        if( tempSigma > 0. and tempSigma < 0.5 and h_deltat_avg_ampCorr_posCutY[_iPosCutY].GetEntries() > 20 ):
            g_timeRes_avg_vs_Y.SetPoint(myPoint,selPos,tempSigma)            
            g_timeRes_avg_vs_Y.SetPointError(myPoint,stepY/2/np.sqrt(12),sigmaErr)    
            myPoint += 1

    # left only
    myPoint = 0
    for _iPosCutY in range(0,_ch.NPOSCUTSY):
        stepY = ( _ch.upperPosCutY - _ch.lowerPosCutY ) / _ch.NPOSCUTSY
        fitdeltat_ampCorr_posCut_leftY = TF1("fitdeltat_left_ampCorr_posCut","gaus", h_deltat_left_ampCorr.GetMean()-h_deltat_left_ampCorr.GetRMS()*2.,
                                       h_deltat_left_ampCorr.GetMean()+h_deltat_left_ampCorr.GetRMS()*2.)
        h_deltat_left_ampCorr_posCutY[_iPosCutY].Fit(fitdeltat_ampCorr_posCut_leftY, "QNR")
        selPos = _ch.lowerPosCutY+(stepY*_iPosCutY)
        tempSigma =  np.sqrt(pow(fitdeltat_ampCorr_posCut_leftY.GetParameter(2),2) - pow(_ch.sigma_ref, 2))
        sigmaErr = fitdeltat_ampCorr_posCut_leftY.GetParError(2)
        if( tempSigma > 0. and tempSigma < 0.5 and h_deltat_left_ampCorr_posCutY[_iPosCutY].GetEntries() > 20 ):
            g_timeRes_left_vs_Y.SetPoint(myPoint, selPos, tempSigma)
            g_timeRes_left_vs_Y.SetPointError(myPoint,stepY/2, sigmaErr)
            myPoint += 1

    # right only
    myPoint = 0
    for _iPosCutY in range(0,_ch.NPOSCUTSY):
        stepY = ( _ch.upperPosCutY - _ch.lowerPosCutY ) / _ch.NPOSCUTSY
        fitdeltat_ampCorr_posCut_rightY = TF1("fitdeltat_right_ampCorr_posCut","gaus", h_deltat_right_ampCorr.GetMean()-h_deltat_right_ampCorr.GetRMS()*2.,
                                       h_deltat_right_ampCorr.GetMean()+h_deltat_right_ampCorr.GetRMS()*2.)
        h_deltat_right_ampCorr_posCutY[_iPosCutY].Fit(fitdeltat_ampCorr_posCut_rightY, "QNR")
        selPos = _ch.lowerPosCutY+(stepY*_iPosCutY)
        tempSigma =  np.sqrt(pow(fitdeltat_ampCorr_posCut_rightY.GetParameter(2),2) - pow(_ch.sigma_ref, 2))
        sigmaErr = fitdeltat_ampCorr_posCut_rightY.GetParError(2)
        if (tempSigma > 0. and tempSigma < 0.5 and  h_deltat_right_ampCorr_posCutY[_iPosCutY].GetEntries() > 20 ):
            g_timeRes_right_vs_Y.SetPoint(myPoint,selPos,tempSigma)
            g_timeRes_right_vs_Y.SetPointError(myPoint,stepY/2,sigmaErr)
            myPoint += 1

  
    c_timeRes_vs_X_Y.cd(2)
    gPad.SetGridy()
    g_timeRes_avg_vs_Y.GetYaxis().SetRangeUser(0, 0.2)
    g_timeRes_avg_vs_Y.SetTitle(";Y [mm];#sigma_{t} [ns]")
    g_timeRes_avg_vs_Y.SetMarkerStyle(20)
    g_timeRes_avg_vs_Y.Draw("ALPE")
    g_timeRes_left_vs_Y.SetLineColor(kRed+1)
    g_timeRes_left_vs_Y.SetMarkerColor(kRed+1)
    g_timeRes_right_vs_Y.SetLineColor(kBlue+1)
    g_timeRes_right_vs_Y.SetMarkerColor(kBlue+1)
    g_timeRes_left_vs_Y.Draw("same LPE")
    g_timeRes_right_vs_Y.Draw("same LPE")
    
    leg = TLegend(0.70,0.78,0.85,0.93,'',"brNDC")
    leg.SetBorderSize(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.03)
    leg.SetLineColor(1)
    leg.SetLineStyle(1)
    leg.SetLineWidth(1)
    leg.SetFillColor(0)
    leg.SetFillStyle(1000)
    leg.AddEntry(g_timeRes_left_vs_Y,  "#sigma^{left}_{t}", "el")
    leg.AddEntry(g_timeRes_right_vs_Y, "#sigma^{right}_{t}","el")
    leg.AddEntry(g_timeRes_avg_vs_Y,   "#sigma^{avg}_{t}",  "pl")
    leg.Draw("same")
  
    printCanvas( c_timeRes_vs_X_Y, _outputDir, _timeAlgo)
  
    # XY position
    #right only
    for _iPosCutX in range(0,_ch.NPOSCUTSX):
        stepX = (_ch.upperPosCutX - _ch.lowerPosCutX) / _ch.NPOSCUTSX
        for _iPosCutY in range(0,_ch.NPOSCUTSY):
            stepY = (_ch.upperPosCutY - _ch.lowerPosCutY) / _ch.NPOSCUTSY
    
            if( h_deltat_avg_ampCorr_posCutXY[_iPosCutX][_iPosCutY].GetEntries() < 20 ): continue

            fitdeltat_ampCorr_posCut_XY = TF1("fitdeltat_right_ampCorr_posCut","gaus", h_deltat_avg_ampCorr_posCutXY[_iPosCutX][_iPosCutY].GetMean()-h_deltat_avg_ampCorr_posCutXY[_iPosCutX][_iPosCutY].GetRMS()*2.,
                                           h_deltat_avg_ampCorr_posCutXY[_iPosCutX][_iPosCutY].GetMean()+h_deltat_avg_ampCorr_posCutXY[_iPosCutX][_iPosCutY].GetRMS()*2.)
            h_deltat_avg_ampCorr_posCutXY[_iPosCutX][_iPosCutY].Fit("fitdeltat_right_ampCorr_posCut","QNR")
            selPosX = _ch.lowerPosCutX+(stepX*_iPosCutX)
            selPosY = _ch.lowerPosCutY+(stepY*_iPosCutY)
            tempSigma = np.sqrt( pow(fitdeltat_ampCorr_posCut_XY.GetParameter(2),2) - pow(_ch.sigma_ref,2) )
            sigmaErr = fitdeltat_ampCorr_posCut_XY.GetParError(2)
      
            if (tempSigma > 0. and tempSigma < 0.5 and h_deltat_avg_ampCorr_posCutXY[_iPosCutX][_iPosCutY].GetEntries() > 20 ):
                h2_timeRes_avg_vs_XY.Fill(selPosX,selPosY,tempSigma)

    c_timeRes_vs_XY = TCanvas ("c_timeRes_vs_XY","c_timeRes_vs_XY",500,500)
    c_timeRes_vs_XY.cd()
    gStyle.SetPaintTextFormat(".3f")
    h2_timeRes_avg_vs_XY.SetStats(0)
    h2_timeRes_avg_vs_XY.Draw("COLZ,text")
    h2_timeRes_avg_vs_XY.GetZaxis().SetRangeUser(0, 0.08)
    h2_timeRes_avg_vs_XY.SetTitle(";x [mm];y [mm];#sigma_{t} [ns]")

    printCanvas( c_timeRes_vs_XY, _outputDir, _timeAlgo)

    ampCorrectedMeasurements  = [sigmaLeft, sigmaRight, sigmaAvg, fitFunc_corrX, fitFunc_corrDiff, h_deltat_diff, sigmaLeftCorr, sigmaLeftCorr_err, sigmaRightCorr, sigmaRightCorr_err, sigmaAvgCorr, sigmaAvgCorr_err, h_deltat_left_ampCorr, h_deltat_right_ampCorr, h_deltat_avg_ampCorr, c_timeRes_comp]

    return ampCorrectedMeasurements


def applyPositionCorrection( _chain, _ch, _timeAlgo, _outputDir, _timePeak, _timeSigma, _mipPeak, _mipPeak_err, _hAmpCut, _ampWalkCorr, _ampCorrectedMeas ):
    """(4) fourth loop events --> to apply pos correction"""

    # define histograms
    h_deltat_wei_ampCorr = TH1F("h_deltat_wei_ampCorr", "h_deltat_wei_ampCorr", 6000, _ch.minDeltat, _ch.maxDeltat)
    h_deltat_avg_ampCorr_diffCorr = TH1F("h_deltat_avg_ampCorr_diffCorr", "h_deltat_avg_ampCorr_diffCorr", 6000, _ch.minDeltat, _ch.maxDeltat)
    h_deltat_avg_ampCorr_posCorr = TH1F("h_deltat_avg_ampCorr_posCorr", "h_deltat_avg_ampCorr_posCorr", 6000, _ch.minDeltat, _ch.maxDeltat)

    # event loop
    nSelectedEntries = 0
    chainCounter = 0
    for entry in _chain:
        if chainCounter%5000 ==0:
            print( '>>> 4th loop: reading entry {0}/{1}'.format(chainCounter, _chain.GetEntries()) )
        chainCounter += 1
     
        myX = _chain.x_dut[0]
        myY = _chain.y_dut[0]
    
        # check if tracking good + cut on BS + cut on MCP + bar selection
        if not eventHasGoodHit( 3, _timeAlgo, myX, myY, _chain, _ch, _timePeak, _timeSigma, _mipPeak):
            continue

        # set time algorithm
        time = getTimeAlgorithmBranch( _chain, _timeAlgo )
        time_ref = _chain.gaus_mean[_ch.timech_id[0]]

        amp1 = _chain.amp[_ch.ampch_id[1]]
        amp2 = _chain.amp[_ch.ampch_id[2]]
        time1 = time[_ch.timech_id[1]]
        time2 = time[_ch.timech_id[2]]
        time1_ampCorr = time1 - _ampWalkCorr[1].Eval(amp1) + _ampWalkCorr[1].Eval(_hAmpCut[1].GetMean())
        time2_ampCorr = time2 - _ampWalkCorr[2].Eval(amp2) + _ampWalkCorr[2].Eval(_hAmpCut[2].GetMean())
    
        sigmaLeft        = _ampCorrectedMeas[0]
        sigmaRight       = _ampCorrectedMeas[1]
        sigmaAvg         = _ampCorrectedMeas[2]
        fitFunc_corrX    = _ampCorrectedMeas[3]
        fitFunc_corrDiff = _ampCorrectedMeas[4]
        h_deltat_diff    = _ampCorrectedMeas[5]

        deltat_avg = 0.5*(time1+time2) - time_ref
        deltat_avg_ampCorr = 0.5*(time1_ampCorr+time2_ampCorr) - time_ref
        deltat_wei_ampCorr = ( time1_ampCorr/pow(sigmaLeft,2) + time2_ampCorr/pow(sigmaRight,2) ) / ( 1/pow(sigmaLeft,2) + 1/pow(sigmaRight,2) ) - time_ref

        posCorr = -1.*fitFunc_corrX.Eval(myX) + fitFunc_corrX.Eval(_ch.centerX)
        diffCorr = -1.*fitFunc_corrDiff.Eval(time2_ampCorr-time1_ampCorr) + fitFunc_corrDiff.Eval(h_deltat_diff.GetMean())

        h_deltat_wei_ampCorr.Fill( deltat_wei_ampCorr )
        h_deltat_avg_ampCorr_posCorr.Fill( deltat_avg_ampCorr + posCorr )
        h_deltat_avg_ampCorr_diffCorr.Fill( deltat_avg_ampCorr + diffCorr )
    
        nSelectedEntries += 1

    print ('>>> 4th loop: selected entries {0}'.format(nSelectedEntries))
    
    c_timeRes_comp = _ampCorrectedMeas[-1]
    c_timeRes_comp.cd(2)

    h_deltat_wei_ampCorr.SetLineColor(kOrange+1)
    h_deltat_wei_ampCorr.SetLineWidth(2)
    h_deltat_wei_ampCorr.Draw("same")
    h_deltat_avg_ampCorr_posCorr.SetLineColor(kViolet+1)
    h_deltat_avg_ampCorr_posCorr.SetLineWidth(2)
    h_deltat_avg_ampCorr_posCorr.Draw("same")  
    
    fitdeltat_wei_ampCorr = TF1("fitdeltat_wei_ampCorr", "gaus", h_deltat_wei_ampCorr.GetMean()-h_deltat_wei_ampCorr.GetRMS()*2, h_deltat_wei_ampCorr.GetMean()+h_deltat_wei_ampCorr.GetRMS()*2)
    fitdeltat_wei_ampCorr.SetLineColor(kOrange+1)
    h_deltat_wei_ampCorr.Fit(fitdeltat_wei_ampCorr, "QR")
    
    fitdeltat_avg_ampCorr_posCorr = TF1("fitdeltat_avg_ampCorr_posCorr", "gaus", h_deltat_avg_ampCorr_posCorr.GetMean()-h_deltat_avg_ampCorr_posCorr.GetRMS()*2, h_deltat_avg_ampCorr_posCorr.GetMean()+h_deltat_avg_ampCorr_posCorr.GetRMS()*2)
    fitdeltat_avg_ampCorr_posCorr.SetLineColor(kViolet+1)
    h_deltat_avg_ampCorr_posCorr.Fit(fitdeltat_avg_ampCorr_posCorr, "QR")
    
    sigmaLeftCorr          = _ampCorrectedMeas[6]
    sigmaLeftCorr_err      = _ampCorrectedMeas[7]
    sigmaRightCorr         = _ampCorrectedMeas[8]
    sigmaRightCorr_err     = _ampCorrectedMeas[9]
    sigmaAvgCorr           = _ampCorrectedMeas[10]
    sigmaAvgCorr_err       = _ampCorrectedMeas[11]
    h_deltat_left_ampCorr  = _ampCorrectedMeas[12]
    h_deltat_right_ampCorr = _ampCorrectedMeas[13]
    h_deltat_avg_ampCorr   = _ampCorrectedMeas[14]
    sigmaWeiCorr        = np.sqrt(pow(fitdeltat_wei_ampCorr.GetParameter(2),2)         - pow(_ch.sigma_ref,2) )
    sigmaAvgCorrPos     = np.sqrt(pow(fitdeltat_avg_ampCorr_posCorr.GetParameter(2),2) - pow(_ch.sigma_ref,2) )
    sigmaWeiCorr_err    = fitdeltat_wei_ampCorr.GetParError(2)
    sigmaAvgCorrPos_err = fitdeltat_avg_ampCorr_posCorr.GetParError(2)
    
    leg = TLegend(0.5,0.68,0.80,0.93,'',"brNDC")
    leg.SetBorderSize(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.03)
    leg.SetLineColor(1)
    leg.SetLineStyle(1)
    leg.SetLineWidth(1)
    leg.SetFillColor(0)  
    leg.AddEntry(h_deltat_left_ampCorr,       '#sigma^{{left}}_{{t}} = {0} #pm {1} ps'.format(round(sigmaLeftCorr*1000,1), round(sigmaLeftCorr_err*1000,1)), "l")
    leg.AddEntry(h_deltat_right_ampCorr,      '#sigma^{{right}}_{{t}} = {0} #pm {1} ps'.format(round(sigmaRightCorr*1000,1), round(sigmaLeftCorr_err*1000,1)), "l")
    leg.AddEntry(h_deltat_avg_ampCorr,        '#sigma^{{avg}}_{{t}} = {0} #pm {1} ps'.format(round(sigmaAvgCorr*1000,1), round(sigmaAvgCorr_err*1000,1)), "l")
    leg.AddEntry(h_deltat_wei_ampCorr,        '#sigma^{{wei}}_{{t}} = {0} #pm {1} ps'.format(round(sigmaWeiCorr*1000,1), round(sigmaWeiCorr_err*1000,1)), "l")
    leg.AddEntry(h_deltat_avg_ampCorr_posCorr,'#sigma^{{avg+pos}}_{{t}} = {0} #pm {1} ps'.format(round(sigmaAvgCorrPos*1000,1), round(sigmaAvgCorrPos_err*1000,1)), "l")
    leg.Draw("same")
    
    print( 'TimeReso : {0} +/- {1} ns'.format( round(sigmaWeiCorr*1000,1), round(sigmaWeiCorr_err*1000,1)) )

    printCanvas( c_timeRes_comp, _outputDir, _timeAlgo )

    """
    txtOutputFitInfo << Form("Angle: %d, #sigma^{left}_{t} = %.1f #pm %.1f ps", angleScan, sigmaLeftCorr*1000, sigmaLeftCorr_err*1000) << " (" << timeAlgo.c_str() << ")"<< std::endl;
    txtOutputFitInfo << Form("Angle: %d, #sigma^{right}_{t} = %.1f #pm %.1f ps", angleScan,  sigmaRightCorr*1000, sigmaRightCorr_err*1000) << " (" << timeAlgo.c_str() << ")"<< std::endl;
    txtOutputFitInfo << Form("Angle: %d, #sigma^{avg}_{t} = %.1f #pm %.1f ps", angleScan,  sigmaAvgCorr*1000, sigmaAvgCorr_err*1000) << " (" << timeAlgo.c_str() << ")"<< std::endl;
    txtOutputFitInfo << Form("Angle: %d, #sigma^{avg+pos}_{t} = %.1f #pm %.1f ps", angleScan,  sigmaAvgCorrPos*1000, sigmaAvgCorrPos_err*1000) << " (" << timeAlgo.c_str() << ")"<< std::endl;
    txtOutputFitInfo << Form("Angle: %d, #sigma^{wei}_{t} = %.1f #pm %.1f ps", angleScan,  sigmaWeiCorr*1000, sigmaWeiCorr_err*1000) << " (" << timeAlgo.c_str() << ")"<< std::endl;
    
    double quantileWei    = GetEffSigma( h_deltat_avg_ampCorr );
    double quantileAvg    = GetEffSigma( h_deltat_wei_ampCorr );
    double quantileAvgPos = GetEffSigma( h_deltat_avg_ampCorr_posCorr );
    
    cout << Form("Angle: %d, #sigma^{avg}_{t} = %.1f #pm %.1f ps", angleScan,  sigmaAvgCorr*1000, sigmaAvgCorr_err*1000) << " (" << timeAlgo.c_str() << ")"<< std::endl;
    cout << Form("Angle: %d, #sigma^{avg+pos}_{t} = %.1f #pm %.1f ps", angleScan,  sigmaAvgCorrPos*1000, sigmaAvgCorrPos_err*1000) << " (" << timeAlgo.c_str() << ")"<< std::endl;
    cout << Form("Angle: %d, #sigma^{wei}_{t} = %.1f #pm %.1f ps", angleScan,  sigmaWeiCorr*1000, sigmaWeiCorr_err*1000) << " (" << timeAlgo.c_str() << ")"<< std::endl;
    
    std::cout << Form("[QUANTILE] Angle: %d, #sigma^{avg}_{t} = %.1f ps", angleScan,  quantileAvg*1000) <<  std::endl;            
    std::cout << Form("[QUANTILE] Angle: %d, #sigma^{wei}_{t} = %.1f ps", angleScan,  quantileWei*1000) <<  std::endl;            
    std::cout << Form("[QUANTILE] Angle: %d, #sigma^{avg+pos}_{t} = %.1f ps", angleScan,  quantileAvgPos*1000) <<  std::endl;            
    """
    outputFitInfo = open( _outputDir+'/timeAndAmpFitInfo.txt', 'w' )
    outputFitInfo.write('OutPut {0} with {1}\n'.format(_outputDir, _timeAlgo))
    for _iCh in range(0,_ch.NCH):
        if _iCh > 0:
		#outputFitInfo.write('Fit Peak [{0} ch {1}]: {2} +/- 0.0 mV\n'.format( _ch.namech[_iCh], _ch.ampch_id[_iCh], round(_mipPeak[_iCh],2)) )
		outputFitInfo.write('Fit Peak [{0} ch {1}]: {2} +/- {3} mV\n'.format( _ch.namech[_iCh], _ch.ampch_id[_iCh], round(_mipPeak[_iCh],2), round(_mipPeak_err[_iCh],2)) )
    outputFitInfo.write('#sigma^{{left}}_{{t}} = {0} #pm {1} ps\n'.format(round(sigmaLeftCorr*1000,1), round(sigmaLeftCorr_err*1000,1)))
    outputFitInfo.write('#sigma^{{right}}_{{t}} = {0} #pm {1} ps\n'.format(round(sigmaRightCorr*1000,1), round(sigmaRightCorr_err*1000,1)))
    outputFitInfo.write('#sigma^{{avg}}_{{t}} = {0} #pm {1} ps\n'.format(round(sigmaAvgCorr*1000,1), round(sigmaAvgCorr_err*1000,1)))
    outputFitInfo.write('#sigma^{{wei}}_{{t}} = {0} #pm {1} ps\n'.format(round(sigmaWeiCorr*1000,1), round(sigmaWeiCorr_err*1000,1)))
    outputFitInfo.write('#sigma^{{avg+pos}}_{{t}} = {0} #pm {1} ps\n'.format(round(sigmaAvgCorrPos*1000,1), round(sigmaAvgCorrPos_err*1000,1)))
    outputFitInfo.close()

    return


def calculatePositionResiduals( _chain, _ch, _timeAlgo, _outputDir, _timePeak, _timeSigma, _mipPeak, _mipPeak_err, _hAmpCut, _ampWalkCorr):
    """(5) fifth loop events --> to calculate position residuals """

    # define histograms
    h_deltat_leftVright         = TH1F("h_deltat_leftVright", "h_deltat_leftVright", 120, -600, 600)
    h_deltat_leftVright_ampCorr = TH1F("h_deltat_leftVright_ampCorr", "h_deltat_leftVright_ampCorr", 120, -600, 600)
    h_deltax_leftVright         = TH1F("h_deltax_leftVright", "h_deltax_leftVright", 40, -100, 100)
    h_deltax_leftVright_ampCorr = TH1F("h_deltax_leftVright_ampCorr", "h_deltax_leftVright_ampCorr", 40, -100, 100)
    h_x1_ampCorr                = TH1F("h_x1_ampCorr", "h_x1_ampCorr", 50, -1000, -750)
    h_x2_ampCorr                = TH1F("h_x2_ampCorr", "h_x2_ampCorr", 50, -1000, -750)
    h2_x1x2                     = TH2F("h2_x1x2", "h2_x1x2", 50, -950, -800, 50, -950, -800)
    h2_x1x2_ampCorr             = TH2F("h2_x1x2_ampCorr", "h2_x1x2_ampCorr", 50, -950, -800, 50, -950, -800)

    # event loop
    nSelectedEntries = 0
    chainCounter = 0
    for entry in _chain:
        if chainCounter%5000 ==0:
            print( '>>> 5th loop: reading entry {0}/{1}'.format(chainCounter, _chain.GetEntries()) )
        chainCounter += 1
     
        myX = _chain.x_dut[0]
        myY = _chain.y_dut[0]
    
        # check if tracking good + cut on BS + cut on MCP + bar selection
        if not eventHasGoodHit( 3, _timeAlgo, myX, myY, _chain, _ch, _timePeak, _timeSigma, _mipPeak):
            continue

        # set time algorithm
        time = getTimeAlgorithmBranch( _chain, _timeAlgo )
        time_ref = _chain.gaus_mean[_ch.timech_id[0]]

        # get left/right SiPM data
        amp1 = _chain.amp[_ch.ampch_id[1]]
        amp2 = _chain.amp[_ch.ampch_id[2]]
        time1 = time[_ch.timech_id[1]]
        time2 = time[_ch.timech_id[2]]
        time1_ampCorr = time1 - _ampWalkCorr[1].Eval(amp1) + _ampWalkCorr[1].Eval(_hAmpCut[1].GetMean())
        time2_ampCorr = time2 - _ampWalkCorr[2].Eval(amp2) + _ampWalkCorr[2].Eval(_hAmpCut[2].GetMean())

        #fill histograms
        deltat         = (time2 - time1)*1000 # [ps]
        deltat_ampCorr = (time2_ampCorr - time1_ampCorr) * 1e3 # [ps]
        t1_ampCorr = (time1_ampCorr - time_ref) * 1e3 # [ps]
        t2_ampCorr = (time2_ampCorr - time_ref) * 1e3 # [ps]
        t1 = (time1 - time_ref) * 1e3 # [ps]
        t2 = (time2 - time_ref) * 1e3 # [ps]
        
        c_lyso = 3e8/1.8 * 1e-12 * 1e3 # [mm/ps] ,  c/n_lyso * s/ps * mm/m
        #print("[AMP CORR] Event {0}, t1 = {1}, t2 = {2}, x1 = {3}, x2 = {4}".format(chainCounter, round(t1_ampCorr,2), round(t2_ampCorr,2), round(t1_ampCorr*c_lyso,2), round(t2_ampCorr*c_lyso,2) ))
        #print("Event {0}, t1 = {1}, t2 = {2}, x1 = {3}, x2 = {4}".format(chainCounter, round(t1,2), round(t2,2), round(t1*c_lyso,2), round(t2*c_lyso,2) ))
        
        h_deltat_leftVright.Fill( deltat )
        h_deltat_leftVright_ampCorr.Fill( deltat_ampCorr )
        h_deltax_leftVright.Fill( deltat*c_lyso )
        h_deltax_leftVright_ampCorr.Fill( deltat_ampCorr*c_lyso )
        h_x1_ampCorr.Fill( t1_ampCorr*c_lyso )
        h_x2_ampCorr.Fill( t2_ampCorr*c_lyso )
        h2_x1x2_ampCorr.Fill( t1_ampCorr*c_lyso, t2_ampCorr*c_lyso )
        h2_x1x2.Fill( t1*c_lyso, t2*c_lyso )

        nSelectedEntries += 1

    print ('>>> 5th loop: selected entries {0}'.format(nSelectedEntries))

    # make and fill canvas
    c_leftVright_deltaT_deltaX = TCanvas ("c_leftVright_deltat_deltaX","c_leftVright_deltat_deltaX",1500,500)
    c_leftVright_deltaT_deltaX.Divide(3,1)
    c_leftVright_deltaT_deltaX.cd(1)
    h_deltat_leftVright.SetLineColor(kBlack)
    h_deltat_leftVright_ampCorr.SetLineColor(kRed)
    h_deltat_leftVright.SetTitle(";t_{left} - t_{right} [ps];Entries / 10 ps")
    h_deltat_leftVright.Draw()
    h_deltat_leftVright_ampCorr.Draw("same")

    c_leftVright_deltaT_deltaX.cd(2)
    h_deltax_leftVright.SetLineColor(kBlack)
    h_deltax_leftVright_ampCorr.SetLineColor(kRed)
    h_deltax_leftVright.SetTitle(";x_{left} - x_{right} [mm];Entries / 5 mm")
    h_deltax_leftVright.Draw()
    h_deltax_leftVright_ampCorr.Draw("same")

    c_leftVright_deltaT_deltaX.cd(3)
    h_x1_ampCorr.SetLineColor(kBlack)
    h_x2_ampCorr.SetLineColor(kRed)
    h_x1_ampCorr.SetTitle(";t_{hit}*c_{LYSO} [mm];Entries / 5 mm")
    h_x1_ampCorr.Draw()
    h_x2_ampCorr.Draw("same")

    c_leftVright_x1Vx2 = TCanvas ("c_leftVright_x1Vx2","c_leftVright_x1Vx2",1000,500)
    c_leftVright_x1Vx2.Divide(2,1)
    c_leftVright_x1Vx2.cd(1)
    h2_x1x2.SetTitle(";x_{left} [mm];x_{right} [mm]")
    h2_x1x2.Draw("colz")
    label = TLatex(0.45, 0.8, 'Uncorrected')
    label.SetNDC()
    label.SetTextFont(42)
    label.SetTextSize(0.03)
    label.Draw("same")

    c_leftVright_x1Vx2.cd(2)
    h2_x1x2_ampCorr.SetTitle(";x_{left} [mm];x_{right} [mm]")
    h2_x1x2_ampCorr.Draw("colz")
    label = TLatex(0.45, 0.8, 'Amp-Walk Corrected')
    label.Draw("same")

    print("[NO Correction] Correlation between x1 and x2 = {0}".format( round(h2_x1x2.GetCorrelationFactor(),3) ))
    print("[Amp Corrected] Correlation between x1_ampCorr and x2_ampCorr = {0}".format( round(h2_x1x2_ampCorr.GetCorrelationFactor(),3) ))

    printCanvas( c_leftVright_deltaT_deltaX, _outputDir, _timeAlgo)
    printCanvas( c_leftVright_x1Vx2, _outputDir, _timeAlgo)

    return
