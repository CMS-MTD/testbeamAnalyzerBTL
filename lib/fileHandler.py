# /usr/bin/python

#Author: Ben Tannenwald
#Date: April 30, 2019
#Purpose: Store functions for opening files

from ROOT import TFile, TChain, TTree, TProfile, TH1F, TH2F, TProfile2D, TCanvas, TLine, TBranch, TF1
from ROOT import kRed, kBlack, kBlue, kOrange
from ROOT import gPad, gStyle
from lib.globalVariables import *
import lib.dataQualityHandler as dq
from time import sleep


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
    _c0.Print( '{0}/{1}_{2}.png'.format( _outputDir, c0.GetTitle, _timeAlgo ) )

    return
               

def returnChain( _dataFolder, _firstRun, _lastRun):

    myTree = TChain("pulse", "pulse")
    f0 = TFile()
    trackingIsSynced = False
    for _iRun in range( int(_firstRun), int(_lastRun)):
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
        
        if myX == -999 or myY == -999:
            continue
    
        h_beamX.Fill( myX )
        h_beamY.Fill( myY )
        h_beamXY.Fill( myX,myY )

        for _iCh in range(0, _ch.NCH):
            h_amp[_iCh].Fill( _chain.amp[_ch.ampch_id[_iCh]] );    
            p2_amp_vs_XY[_iCh].Fill( myX, myY, _chain.amp[_ch.ampch_id[_iCh]] );
            
            #cut on BS
            if myX == -999 or myY == -999:
                continue
            if (abs(myX-_ch.centerX) > _ch.BSX) or (abs(myY-_ch.centerY) > _ch.BSY) :
                continue
      
            #this is a very inefficient way to do this, but adaptive branch naming is hard - BBT 05/02/19
            if (_timeAlgo == 'LP2_20'):
                h_time[_iCh].Fill( _chain.LP2_20[ _ch.timech_id[_iCh]] )
            elif (_timeAlgo == 'LP2_30'):
                h_time[_iCh].Fill( _chain.LP2_30[ _ch.timech_id[_iCh]] )
            elif (_timeAlgo == 'LP2_50'):
                h_time[_iCh].Fill( _chain.LP2_50[ _ch.timech_id[_iCh]] )
            elif (_timeAlgo == 'IL_20'):
                h_time[_iCh].Fill( _chain.IL_20[ _ch.timech_id[_iCh]] )
            elif (_timeAlgo == 'IL_30'):
                h_time[_iCh].Fill( _chain.IL_30[ _ch.timech_id[_iCh]] )
            elif (_timeAlgo == 'IL_50'):
                h_time[_iCh].Fill( _chain.IL_50[ _ch.timech_id[_iCh]] )
            
            if (_chain.amp[_ch.ampch_id[_iCh]] > _ch.ampmin_cut[_iCh]) and (_chain.amp[_ch.ampch_id[_iCh]] < _ch.ampmax_cut[_iCh]) :
                h_amp_cut[_iCh].Fill( _chain.amp[_ch.ampch_id[_iCh]] );
                
        nSelectedEntries += 1
    
    print('\n>>> 1st loop: selected entries {0}'.format(nSelectedEntries))

    # ---------------------------- drawing beam histos ----------------------------

    c_beamXY = TCanvas("c_beamXY","c_beamXY",500,500)
    c_beamXY.cd()
    h_beamXY.SetStats(0)
    h_beamXY.SetTitle(";X [mm]; Y[mm];entries")
    h_beamXY.Draw("COLZ")

    x_BS_min = TLine(_ch.centerX - _ch.BSX, _ch.minY, _ch.centerX - _ch.BSX, _ch.maxY)
    x_BS_max = TLine(_ch.centerX + _ch.BSX, _ch.minY, _ch.centerX + _ch.BSX, _ch.maxY)
    y_BS_min = TLine(_ch.minX, _ch.centerY - _ch.BSY, _ch.maxX, _ch.centerY-_ch.BSY)
    y_BS_max = TLine(_ch.minX, _ch.centerY + _ch.BSY, _ch.maxX, _ch.centerY+_ch.BSY)

    print(_ch.centerX - _ch.BSX, _ch.minY, _ch.centerX - _ch.BSX, _ch.maxY)
    print(_ch.minX, _ch.centerY + _ch.BSY, _ch.maxX, _ch.centerY + _ch.BSY)

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
  
    # ---------------------------- fitting and drawing mip peak ----------------------------
    mip_peak     = [0 for i in range(0,_ch.NCH)]
    mip_peak_err = [0 for i in range(0,_ch.NCH)]
  
    c_amp = TCanvas("c_amp","c_amp",1000,500*(_ch.NCH+1)/2)
    c_amp.Divide(2,(_ch.NCH+1)/2)
    l_peakAmp =  TLatex()

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
            print( 'Fit Peak [{0}]: {1} +/- {2} mV'.format( _iCh, round(mip_peak[_iCh],2), round(mip_peak_err[_iCh],2)) )
            #txtOutputFitInfo << Form("Angle: %d, ", angleScan) << namech[_iCh].c_str() <<", amp peak[" << _iCh << "] = " << Form("%.2f #pm %.2f mV", mip_peak[_iCh], mip_peak_err[_iCh]) << std::endl
            
            #TLine* lowcut = new TLine(std::max(rel_amp_cut_low*mip_peak[_iCh],ampmin_cut[_iCh]),0.,std::max(rel_amp_cut_low*mip_peak[_iCh],ampmin_cut[_iCh]),h_amp[_iCh].GetMaximum())
            #TLine* higcut = new TLine(std::min(rel_amp_cut_hig*mip_peak[_iCh],ampmax_cut[_iCh]),0.,std::min(rel_amp_cut_hig*mip_peak[_iCh],ampmax_cut[_iCh]),h_amp[_iCh].GetMaximum())
            #lowcut.Draw("same")
            #higcut.Draw("same")
            
            l_peakAmp = TLatex(0.45, 0.8, 'Fit Peak: {0} #pm {1} mV'.format( round(mip_peak[_iCh],2), round(mip_peak_err[_iCh],2)))
            l_peakAmp.SetNDC()
            l_peakAmp.SetTextFont(42)
            l_peakAmp.SetTextSize(0.03)
            l_peakAmp.Draw("same")
            
        _ch.latexch[_iCh].Draw("same")

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

    sleep(15)
