//root AxE_MC_withFits_root.C'(1, 1, 1)'

//look at MC to compute AxE = nEvents surviving a selection in a q2 bin / original nEvents in MC
//compute AxE by counting (events in q2 bin) and by fitting
//results are AxE and signal model
//save results in workspace to be loaded in the final fits

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
void AxE_MC_withFits_root(int isMC=1, int isEE=1, int isResonant=0, float mvaCut = 12.68){

  gSystem->Load("../HiggsAnalysis/CombinedLimit/lib/libHiggsAnalysisCombinedLimit.so"); 

  // float LumiNonReso = 1.53E+04; // fb-1
  // float LumiReso = 45.16919828; // fb-1

  float NOriginalNonReso = 7587169;
  float NOriginalReso = 2122456;


  int nBins = 3;
  std::vector<std::pair<float, float>> mll_bins_Range;
  mll_bins_Range.push_back(std::pair<float, float>(1.1, 2.4));
  mll_bins_Range.push_back(std::pair<float, float>(2.4, 4.));
  mll_bins_Range.push_back(std::pair<float, float>(1.1, 4.));


  //current files from Otto, need to adapt to other's outputs
  std::string inputFileList = "/afs/cern.ch/user/k/klau/myWorkspace/public/ForArabella/03Mar2020_pf/";
  if(!isResonant && isEE) 
    inputFileList += "RootTree_2020Jan16_BuToKee_BToKEEAnalyzer_2020Feb18_fullq2_pf_isoPFMVADphiptImb_weighted_pauc02_mva.root";
  else if(isResonant && isEE) 
    inputFileList += "RootTree_2020Jan16_BuToKJpsi_Toee_BToKEEAnalyzer_2020Feb18_fullq2_pf_isoPFMVADphiptImb_weighted_pauc02_mva.root";
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

  //loop over tree load hostos and RooDataSet
  TH1F* hHisto[3];
  RooWorkspace w("signalModel");
  RooRealVar x("x", "", 4.5, 6.);
  RooDataSet* data[3];

  int ijC = 0;
  for(auto ij: mll_bins_Range){
    hHisto[ijC] = new TH1F(Form("hHisto_mll_%.1f-%.1f", ij.first, ij.second), "", 500, 4.5, 6.);

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
	hHisto[ijC]->Fill(B_fit_mass);
      }
      ++ijC;
    }
  }
  w.import(x);

  //not really needed. Just to take a look at histos
  TFile out(Form("AxE_isEE%d_isResonant%d.root", isEE, isResonant), "recreate");
  out.cd();
  for(auto ij=0; ij<nBins; ++ij){
    //std::cout << " histo name = " << hHisto[ij]->GetName() << std::endl;
    hHisto[ij]->Write(hHisto[ij]->GetName());
  }
  out.Close();


  //now fitting part

  float nEv_postFit[3];
  float nEvError_postFit[3];
  float meanVal_postFit[3];
  float sigmaVal_postFit[3];

  for(auto ij=0; ij<nBins; ++ij){

    std::string suffixName = Form("_mll_%.1f-%.1f", mll_bins_Range[ij].first, mll_bins_Range[ij].second);

    w.factory(("nSig"+suffixName+"[1.e5, 0.0, 1.e6]").c_str());
    //ele double CB
    if(isEE)
      w.factory(("RooDoubleCBFast::smodel"+suffixName+"(x, mu"+suffixName+"[5.3,4.5,6], sigma"+suffixName+"[0.04,0.001,0.08], aL"+suffixName+"[1.5, 0, 5.], nL"+suffixName+"[2, 0, 10], aR"+suffixName+"[1.5, 0, 5.], nR"+suffixName+"[10, 0, 15])").c_str());
    //muon double CB
    else
      w.factory(("RooDoubleCBFast::smodel"+suffixName+"(x, mu"+suffixName+"[5.3,4.5,6], sigma"+suffixName+"[0.02,0.001,0.05], aL"+suffixName+"[1.5, 0, 3.], nL"+suffixName+"[2, 0, 5], aR"+suffixName+"[1.5, 0, 5.], nR"+suffixName+"[10, 5, 15])").c_str());

    w.factory(("SUM::model"+suffixName+"(nSig"+suffixName+" * smodel"+suffixName+")").c_str());
    RooAbsPdf * model = w.pdf(("model"+suffixName).c_str());


    //fit Bmass                                                                       
    RooFitResult * r = model->fitTo(*(data[ij]), Minimizer("Minuit2"),Save(true));
    
    RooPlot * plot = w.var("x")->frame();
    if(isEE){
      if(isResonant) plot->SetXTitle("K(JPsi)ee mass (GeV)");
      else plot->SetXTitle("Kee mass (GeV)");
    }
    plot->SetTitle("");
    data[ij]->plotOn(plot, Binning(200));
    model->plotOn(plot);
    model->paramOn(plot, RooFit::Layout(0.2,0.4,0.9),RooFit::Format("NEA",AutoPrecision(1)));
    plot->getAttLine()->SetLineColorAlpha(kWhite, 0.2);
    plot->getAttText()->SetTextSize(0.03);
    plot->getAttText()->SetTextFont(42);

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
    std::string outName = "plots/Bmass_MC/";
    if(isEE) outName += "Kee";
    else outName += "Kmm";
    outName += Form("_mll_%.1f-%.1f", mll_bins_Range[ij].first, mll_bins_Range[ij].second);
    cc->Print((outName+".png").c_str(), "png");
    cc->Print((outName+".pdf").c_str(), "pdf");
    cc->Print((outName+".root").c_str(), "root");

    RooRealVar* parS = (RooRealVar*) r->floatParsFinal().find(("nSig"+suffixName).c_str());
    nEv_postFit[ij] = parS->getValV();
    nEvError_postFit[ij] = parS->getError();
    RooRealVar* parMean = (RooRealVar*) r->floatParsFinal().find(("mu"+suffixName).c_str());
    RooRealVar* parSigma = (RooRealVar*) r->floatParsFinal().find(("sigma"+suffixName).c_str());
    meanVal_postFit[ij] = parMean->getValV();
    sigmaVal_postFit[ij] = parSigma->getValV();

  }

  // print AxE results
  for(auto ij=0; ij<nBins; ++ij){

    std::cout << "\n mll bin = " << mll_bins_Range[ij].first << " - " << mll_bins_Range[ij].second << std::endl;

    std::cout << " by fitting: " << std::endl; 
    std::cout << "    mean = " << meanVal_postFit[ij] << " sigma = " << sigmaVal_postFit[ij] << std::endl;
    std::cout << "    n SigEvents = " << nEv_postFit[ij] << " +/- " << nEvError_postFit[ij] << std::endl;

    RooRealVar num("num", "", nEv_postFit[ij]);
    RooRealVar den("den", "", (isResonant ? NOriginalReso : NOriginalNonReso)); 
    RooFormulaVar axe(Form("axe_mll_%.1f-%.1f", mll_bins_Range[ij].first, mll_bins_Range[ij].second), "", " 1. * @0 / @1", RooArgList(num, den));
    RooRealVar axeVal(Form("axeVal_mll_%.1f-%.1f", mll_bins_Range[ij].first, mll_bins_Range[ij].second), "", axe.getVal());

    std::cout << "    /nOriginal  = " << 1.*nEv_postFit[ij] / (isResonant ? NOriginalReso : NOriginalNonReso)
	      << " +/- " << 1.*nEvError_postFit[ij] / (isResonant ? NOriginalReso : NOriginalNonReso) << std::endl;
  
    
    std::cout << " by counting: " << std::endl; 
    std::cout << "    n SigEvents = " << hHisto[ij]->Integral() << std::endl;
    std::cout << "    /nOriginal  = " << 1.*hHisto[ij]->Integral() / (isResonant ? NOriginalReso : NOriginalNonReso) << std::endl;

    w.import(axeVal); 
  }

  w.Print();
  TFile outWS(Form("signalModel_MC_isResonant%d_isEE%d_workspace.root", isResonant, isEE), "recreate");
  outWS.cd();
  w.Write();
  outWS.Close();
}
