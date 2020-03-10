# AnalysisCommon


## Getting started
```shell
cmsrel CMSSW_10_2_15   (if not already done)
cd CMSSW_10_2_15/src
cmsenv
git cms-init
```


## Get locally full HiggsAnalysis-CombinedLimit tools, to easily export libraries
link from https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit
```shell
git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
cd HiggsAnalysis/CombinedLimit
source env_standalone.sh 
make -j 8; make # second make fixes compilation error of first
```


## Get the scripts to run the fitting chain
```shell
cd CMSSW_10_2_15/src
git@github.com:CMSBParking/AnalysisCommon.git
cd AnalysisCommon
```


## Foreseen sequences
* run on MC resonant and non-resonant to compute AxE and save signal model in a workspace
Need to adapt format for input files and selections to other's ntuples (currently based on Otto's)
```shell
root AxE_MC_withFits_root.C
```
* run on MC K* to get the model for the partially reconstructed background and save in a workspace
Need to adapt format for input files to other's ntuples
```shell
python makePlot_fitPeak_unbinned.py -p -i /afs/cern.ch/user/k/klau/myWorkspace/public/ForArabella/03Mar2020_pf/RootTree_2020Jan16_BdToKstarJpsi_ToKPiee_BToKEEAnalyzer_2020Feb18_fullq2_pf_isoPFMVADphiptImb_weighted_pauc02_mva.root -o part_workspace.root
```
Model WP dependent, q2 bin dependent and for both resonant and non-resonant
```shell
root modelFrom_partialRecoBkg.C
```



* run the simultaneous fit on the resonant (data) and non-resonant (toy MC) to extract expected S/sqrt(S+B) with error
  - import models for signal (MC resonant) and partially reconstructed background. Build resonant_PDF (signal + partialRecoBkg + expo for combinatorial bkg)
  - fit data in the resonant q2 bin and get expected nSignalEvents
  - fit data in the low q2 bin blinded (sidebands only) with PDF (partialRecoBkg + expo) to get expected background
  - rescale nSignalEvents_resonant to expected nSignalEvents_lowq2 by (BRfraction * AxEfraction) 
  - generate toy MC with background from lowq2 sidebands and expected nSignalEvents_lowq2
  - perform a simultaneous fit to toy MC (lowq2) and data (resonant) to have a better estimate of nEvents with errors
Currently set for Otto's ntuples format and signal selection. Need to be updated to other's ntuples and other WPs.
```shell
root makeSimultaneousFits.C
```