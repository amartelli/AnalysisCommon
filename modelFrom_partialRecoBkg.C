//root modelFrom_partialRecoBkg.C'(1, 0)'


//roofit
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooChebychev.h"
#include "RooCBShape.h"
#include "RooArgusBG.h"
#include "RooFFTConvPdf.h"
#include "RooAddPdf.h"
#include "RooWorkspace.h"
#include "RooFitResult.h"
#include "RooFormulaVar.h"
#include "TText.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TH1F.h"
#include "TH1D.h"
#include "RooPlot.h"

 #include <iostream>
 #include <string>
 #include <vector>
 #include <utility>


using namespace RooFit;
void modelFrom_partialRecoBkg(int isEE=1, int isResonant=0, float mvaCut = 12.68){

  gSystem->Load("../HiggsAnalysis/CombinedLimit/lib/libHiggsAnalysisCombinedLimit.so"); 


  int nBins = 3;
  std::vector<std::pair<float, float>> mll_bins_Range;
  mll_bins_Range.push_back(std::pair<float, float>(1.1, 2.4));
  mll_bins_Range.push_back(std::pair<float, float>(2.4, 4.));
  mll_bins_Range.push_back(std::pair<float, float>(1.1, 4.));

 
  //current files from Otto, need to adapt to other's outputs
  std::string inputFileList = "/afs/cern.ch/user/k/klau/myWorkspace/public/ForArabella/03Mar2020_pf/";
  if(!isResonant && isEE)
    inputFileList += "RootTree_2020Jan16_BdToKstaree_BToKEEAnalyzer_2020Feb18_fullq2_pf_isoPFMVADphiptImb_weighted_pauc02_mva.root";
  else if(isResonant && isEE)
    inputFileList += "RootTree_2020Jan16_BdToKstarJpsi_ToKPiee_BToKEEAnalyzer_2020Feb18_fullq2_pf_isoPFMVADphiptImb_weighted_pauc02_mva.root";
  //need to add for muons


  std::cout << " isEE = " << isEE << " isResonant = " << isResonant << std::endl;
  std::cout << "inputFileList = " << inputFileList << std::endl;


  TChain* tree = new TChain("tree");
  tree->Add(inputFileList.c_str());

  int nEntriesTot = tree->GetEntries();
  std::cout << " loaded nEvents = " << nEntriesTot  << std::endl;

  float B_mll_fullfit;
  float B_fit_mass;
  float B_bdt_weight;

  tree->SetBranchStatus("*", 0);
  if(isEE){
    tree->SetBranchStatus("BToKEE_mll_fullfit", 1);    tree->SetBranchAddress("BToKEE_mll_fullfit", &B_mll_fullfit);
    tree->SetBranchStatus("BToKEE_fit_mass", 1);    tree->SetBranchAddress("BToKEE_fit_mass", &B_fit_mass);
    tree->SetBranchStatus("BToKEE_xgb", 1);    tree->SetBranchAddress("BToKEE_xgb", &B_bdt_weight);
  }
  else{
    tree->SetBranchStatus("BToKMuMu_mll_fullfit", 1);    tree->SetBranchAddress("BToKMuMu_mll_fullfit", &B_mll_fullfit);
    tree->SetBranchStatus("BToKMuMu_fit_mass", 1);    tree->SetBranchAddress("BToKMuMu_fit_mass", &B_fit_mass);
    tree->SetBranchStatus("BToKMuMu_xgb", 1);    tree->SetBranchAddress("BToKMuMu_xgb", &B_bdt_weight);
  }

  //loop over tree load  RooDataSet
  RooWorkspace w("partialModel");
  RooRealVar x("x", "", 4.5, 6.);
  RooDataSet* data[3];

  int ijC = 0;
  for(auto ij: mll_bins_Range){
    data[ijC] = new RooDataSet(Form("data_mll_%.1f-%.1f", ij.first, ij.second),
                               Form("data_mll_%.1f-%.1f", ij.first, ij.second), RooArgSet(x));
    ++ijC;
  }
  for(int ij=0; ij<nEntriesTot; ++ij){
    tree->GetEntry(ij);

    if(B_bdt_weight <= mvaCut) continue;
    int ijC = 0;
    for(auto ijB: mll_bins_Range){
      if(B_mll_fullfit >= ijB.first && B_mll_fullfit < ijB.second) {
        x = B_fit_mass;
        data[ijC]->add(RooArgSet(x));
      }
      ++ijC;
    }
  }
  w.import(x);


  for(auto ij=0; ij<nBins; ++ij){

    w.import(*data[ij]);
    std::string suffixName = Form("_mll_%.1f-%.1f", mll_bins_Range[ij].first, mll_bins_Range[ij].second);

    std::string pippo = "KeysPdf::partial"+suffixName+"(x, data"+suffixName+", MirrorLeft, 2.0)";
    std::cout << pippo.c_str() << std::endl;

    w.factory(("KeysPdf::partial"+suffixName+"(x, data"+suffixName+", MirrorLeft, 2.0)").c_str());
    //w.factory(("KeysPdf::partial"+suffixName+"(x, *data[ij], MirrorLeft, 2.0)").c_str());
    RooAbsPdf* partialModel =w.pdf(("partial"+suffixName).c_str());
    
    RooPlot* plot = w.var("x")->frame();
    if(isEE){
      if(isResonant) plot->SetXTitle("K*(JPsi)ee mass (GeV)");
      else plot->SetXTitle("Kee mass (GeV)");
    }
    plot->SetTitle("");
    data[ij]->plotOn(plot, Binning(50));
    partialModel->plotOn(plot, LineColor(kRed)); //MoveToBack());

    float chi2 = plot->chiSquare(6);


    TCanvas * cc = new TCanvas();
    cc->SetLogy(0);
    plot->Draw();
    
    TLatex tL;
    tL.SetNDC();
    tL.SetTextSize(0.04);
    tL.SetTextFont(42);
    tL.DrawLatex(0.65,0.8, Form("chi2 = %.1f ",chi2));
    TLatex tL1;
    tL1.SetNDC();
    tL1.SetTextSize(0.05);
    tL1.SetTextFont(42);
    tL1.DrawLatex(0.65,0.7, Form("mll in %.1f-%.1f", mll_bins_Range[ij].first, mll_bins_Range[ij].second));
    std::string outName = "plots/PartiallyReco_MC/";
    if(isEE) outName += "Kee";
    else outName += "Kmm";
    outName += Form("_mll_%.1f-%.1f", mll_bins_Range[ij].first, mll_bins_Range[ij].second);
    cc->Print((outName+".png").c_str(), "png");
    cc->Print((outName+".pdf").c_str(), "pdf");
    cc->Print((outName+".root").c_str(), "root");
  }

  w.Print();
  TFile outWS(Form("partiallyRecoBkgModel_MC_isResonant%d_isEE%d_workspace.root", isResonant, isEE), "recreate");
  outWS.cd();
  w.Write();
  outWS.Close();
}
