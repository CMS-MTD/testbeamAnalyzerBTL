# testbeamAnalyzerBTL
Python code for analyzing BTL testbeam data

## Setup

>$ cmsrel CMSSW_9_4_0
>$ cd CMSSW_9_4_0/src
>$ cmsenv
>$ git clone git@github.com:CMS-MTD/testbeamAnalyzerBTL.git


# Commands
## Local running
python analyze_base.py [-h] [--bar BAR] [--firstRun FIRSTRUN] [--lastRun LASTRUN] [--timeAlgo TIMEALGO] [--biasVoltage BIASVOLTAGE]

## Condor submission
python submitRunsToCondor.py [-h] [--bar BAR] [--firstRun FIRSTRUN] [--lastRun LASTRUN] [--timeAlgo TIMEALGO] [--biasVoltage BIASVOLTAGE]
