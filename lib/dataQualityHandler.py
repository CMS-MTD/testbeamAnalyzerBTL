# /usr/bin/python

#Author: Ben Tannenwald
#Date: May 3, 2019
#Purpose: Store functions for data quality monitoring

from ROOT import TFile, TChain, TTree, TProfile, TH1F, TH2F, TProfile2D, TCanvas, TLine, TBranch, TF1

def returnNumberOfBinsAboveAmpThreshold( _profile, _threshold = 0.35 ):

    _binsAboveThreshold = 0 

    for _iBin in range(1, _profile.GetNbinsX()):
        if _profile.GetBinContent(_iBin) > _threshold and _profile.GetBinContent(_iBin-1) > _threshold: 
            _binsAboveThreshold += 1
  
    #print (_binsAboveThreshold)
    return _binsAboveThreshold


def checkTrackingSynchronization( _filename, _iRun):
    _fileToCheck = TFile.Open( _filename );
    if _fileToCheck.GetListOfKeys().Contains("pulse")==False:
        return False

    _pulse = _fileToCheck.Get("pulse");
  
    h_box1_L = TProfile("h_box1_L", "h_box1_L", 30, 10, 35);
    h_box1_R = TProfile("h_box1_R", "h_box1_R", 30, 10, 35);
    h_box2_L = TProfile("h_box2_L", "h_box2_L", 30, 10, 35);
    h_box2_R = TProfile("h_box2_R", "h_box2_R", 30, 10, 35);
    h_box3_L = TProfile("h_box3_L", "h_box3_L", 30, 10, 35);
    h_box3_R = TProfile("h_box3_R", "h_box3_R", 30, 10, 35);
    h_single_L = TProfile("h_single_L", "h_single_L", 30, 10, 35);
    h_single_R = TProfile("h_single_R", "h_single_R", 30, 10, 35);
    
    _pulse.Draw("amp[1] > 400 :  y_dut[0]>>+h_box1_L", "ntracks==1");
    _pulse.Draw("amp[2] > 400 :  y_dut[0]>>+h_box2_L", "ntracks==1");
    _pulse.Draw("amp[3] > 400 :  y_dut[0]>>+h_box3_L", "ntracks==1");
    _pulse.Draw("amp[4] > 400 :  y_dut[0]>>+h_box1_R", "ntracks==1");
    _pulse.Draw("amp[5] > 400 :  y_dut[0]>>+h_box2_R", "ntracks==1");
    _pulse.Draw("amp[6] > 400 :  y_dut[0]>>+h_box3_R", "ntracks==1");
    _pulse.Draw("amp[19] > 200 : y_dut[0]>>+h_single_L", "ntracks==1");
    _pulse.Draw("amp[20] > 200 : y_dut[0]>>+h_single_R", "ntracks==1");

    _nBinsForSyncedFlag = 2;
    found_box1_L   = True if returnNumberOfBinsAboveAmpThreshold(h_box1_L) >= _nBinsForSyncedFlag else False
    found_box2_L   = True if returnNumberOfBinsAboveAmpThreshold(h_box2_L) >= _nBinsForSyncedFlag else False
    found_box3_L   = True if returnNumberOfBinsAboveAmpThreshold(h_box3_L) >= _nBinsForSyncedFlag else False
    found_box1_R   = True if returnNumberOfBinsAboveAmpThreshold(h_box1_R) >= _nBinsForSyncedFlag else False
    found_box2_R   = True if returnNumberOfBinsAboveAmpThreshold(h_box2_R) >= _nBinsForSyncedFlag else False
    found_box3_R   = True if returnNumberOfBinsAboveAmpThreshold(h_box3_R) >= _nBinsForSyncedFlag else False
    found_single_L = True if returnNumberOfBinsAboveAmpThreshold(h_single_L) >= _nBinsForSyncedFlag else False
    found_single_R = True if returnNumberOfBinsAboveAmpThreshold(h_single_R) >= _nBinsForSyncedFlag else False

    # ** A. This is the tightest "goodness" condition possible. looser is potentially an option
    #if found_box1_L == True and found_box2_L == True and found_box3_L == True and found_box1_R == True and found_box2_R == True and found_box3_R == True and found_single_L == True and found_single_R == True:
    #** B. This checks tracking using box2+box3+single but skipping box1 (top bar in box)
    #if found_box2_L == True and found_box3_L == True and found_box2_R == True and found_box3_R == True and found_single_L == True and found_single_R == True:
    #** C. This checks tracking using box1+box3 but skipping box2 (middle bar and single bar)
    if found_box1_L == True and found_box3_L == True and found_box1_R == True and found_box3_R == True and found_single_L == True and found_single_R == True:
        return True
  
    return False
