//#include "RooArgusBG.h"
#include "RooAddPdf.h"
#include "RooCBShape.h"
#include "RooChebychev.h"
#include "RooChi2Var.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooFFTConvPdf.h"
#include "RooFitResult.h"
#include "RooFormulaVar.h"
#include "RooGaussian.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "RooSimultaneous.h"
#include "RooWorkspace.h"

#include "TFile.h"
#include "TText.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TPaveText.h"


using namespace RooFit;


void makeSimultaneousFits(int isEE=1, float mvaCut = 12.68){

  gSystem->Load("../HiggsAnalysis/CombinedLimit/lib/libHiggsAnalysisCombinedLimit.so");
  gROOT->Macro("/afs/cern.ch/user/a/amartell/public/setStyle.C");

  RooWorkspace w("w");
  RooRealVar x("x", "", 4.5, 6.);


  //read workspace for signal model
  TFile* inF1 = TFile::Open("signalModel_MC_isResonant1_isEE1_workspace.root");
  RooWorkspace* wSignal = (RooWorkspace*)inF1->Get("signalModel"); //->Clone("wSignal");

  //non resonant for AxE
  TFile* inF1nnR = TFile::Open("signalModel_MC_isResonant0_isEE1_workspace.root");
  RooWorkspace* wSignalnnR = (RooWorkspace*)inF1nnR->Get("signalModel"); //->Clone("wSignal"); 

  //read workspace for partially reco bkg
  TFile* inF2 = TFile::Open("part_workspace.root");
  RooWorkspace* wPartial = (RooWorkspace*)inF2->Get("myPartialWorkSpace"); //->Clone("wPartial");

  wSignal->Print();
  wPartial->Print();


  float BmassFit = wSignal->var("mu_mll_2.4-4.0")->getVal();
  float BsigmaFit = wSignal->var("sigma_mll_2.4-4.0")->getVal();

  float Blow_sideband = BmassFit - 3.*BsigmaFit;
  float Bhigh_sideband = BmassFit + 3.*BsigmaFit;

  float Blow_blindmin = x.getBinning().lowBound();
  float Bhigh_blindmax = x.getBinning().highBound();

  std::cout << " Blow_blindmin = " << Blow_blindmin << " Bhigh_blindmax = " << Bhigh_blindmax << std::endl;

  //load data
  std::string inputFileList = "/afs/cern.ch/user/k/klau/myWorkspace/public/ForArabella/03Mar2020_pf/";
  if(isEE)
    inputFileList += "RootTree_2020Jan16_Run2018ABCD_BToKEEAnalyzer_2020Feb18_fullq2_pf_isoPFMVADphiptImb_weighted_pauc02_mvaCut02.root";

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
    tree->SetBranchStatus("BToKEE_fit_mass", 1);       tree->SetBranchAddress("BToKEE_fit_mass", &B_fit_mass);
    tree->SetBranchStatus("BToKEE_xgb", 1);            tree->SetBranchAddress("BToKEE_xgb", &B_bdt_weight);
  }
  else{
    tree->SetBranchStatus("BToKMuMu_mll_fullfit", 1);  tree->SetBranchAddress("BToKMuMu_mll_fullfit", &B_mll_fullfit);
    tree->SetBranchStatus("BToKMuMu_fit_mass", 1);     tree->SetBranchAddress("BToKMuMu_fit_mass", &B_fit_mass);
    tree->SetBranchStatus("BToKMuMu_xgb", 1);          tree->SetBranchAddress("BToKMuMu_xgb", &B_bdt_weight);
  }

  //for the moment JPsi and lowq2
  std::vector<std::pair<float, float>> mll_bins_Range;
  mll_bins_Range.push_back(std::pair<float, float>(1.1, 2.4));
  mll_bins_Range.push_back(std::pair<float, float>(2.4, 4.));

  RooDataSet* data[2];
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
	//blind lowq2
	if(ijC == 0 && 
	   B_fit_mass >= Blow_sideband && B_fit_mass <= Bhigh_sideband) continue;
        x = B_fit_mass;
        data[ijC]->add(RooArgSet(x));
      }
      ++ijC;
    }
  }
  w.import(x);


  //Resonant category: background model + signal model
  w.factory("nsigR[5.e3, 0., 1.e4]");
  w.factory("nbkgR[1.e4, 0., 1.e5]");

  w.factory("Exponential::combBkgR(x, tauR[-3.0, -100.0, -1.0e-5])");
  w.import("part_workspace.root:myPartialWorkSpace:partial", RenameAllNodes("Bkg"));

  RooRealVar smodelR_mu("smodelR_mu", "", (wSignal->var("mu_mll_2.4-4.0"))->getVal());
  RooRealVar smodelR_sigma("smodelR_sigma", "", (wSignal->var("sigma_mll_2.4-4.0"))->getVal());
  RooRealVar smodelR_aL("smodelR_aL", "", (wSignal->var("aL_mll_2.4-4.0"))->getVal());
  RooRealVar smodelR_nL("smodelR_nL", "", (wSignal->var("nL_mll_2.4-4.0"))->getVal());
  RooRealVar smodelR_aR("smodelR_aR", "", (wSignal->var("aR_mll_2.4-4.0"))->getVal());
  RooRealVar smodelR_nR("smodelR_nR", "", (wSignal->var("nR_mll_2.4-4.0"))->getVal()); 
  w.import(smodelR_mu);
  w.import(smodelR_sigma);
  w.import(smodelR_aL);
  w.import(smodelR_nL);
  w.import(smodelR_aR);
  w.import(smodelR_nR);
  w.factory("RooDoubleCBFast::smodelR(x, smodelR_mu, smodelR_sigma, smodelR_aL, smodelR_nL, smodelR_aR, smodelR_nR)");

  w.factory("SUM::modelBkg(f1[0.5, 0., 10.] * partial_Bkg, combBkgR)");
  w.factory("SUM::modelResonant(nsigR * smodelR, nbkgR * modelBkg)");

  RooAbsPdf * modelResonant = w.pdf("modelResonant");


  //fit Bmass in resonant category
  RooFitResult* rR = modelResonant->fitTo(*(data[1]), Extended(true), Minimizer("Minuit2"),Save(true));

  w.var("x")->setRange("signalRegioncut", Blow_sideband, Bhigh_sideband);
  RooAbsReal* backgroundFraction = w.pdf("modelBkg")->createIntegral(*(w.var("x")), RooFit::NormSet(*(w.var("x"))), Range("signalRegioncut"));
  RooAbsReal* signalFraction = w.pdf("smodelR")->createIntegral(*(w.var("x")), RooFit::NormSet(*(w.var("x"))), Range("signalRegioncut"));

  float nBkgInsignalRegion = backgroundFraction->getVal() * w.var("nbkgR")->getVal();
  float nSignalInsignalRegion = signalFraction->getVal() * w.var("nsigR")->getVal();

  std::cout << " resonant q2: nSignal = " << nSignalInsignalRegion << " nBkg = " << nBkgInsignalRegion 
	    << " S/sqrt(S+B) = " << nSignalInsignalRegion / sqrt(nSignalInsignalRegion+nBkgInsignalRegion) << std::endl;
  std::cout << "\n\n " << std::endl;


  RooPlot* plotR = w.var("x")->frame();
  plotR->SetXTitle("K(JPsi)ee mass (GeV)");
  plotR->SetTitle("");
  data[1]->plotOn(plotR, Binning(50));
  modelResonant->plotOn(plotR);
  modelResonant->plotOn(plotR, Name("fitBkg"), Components("modelBkg"),LineStyle(kDashed), LineColor(kBlue));
  modelResonant->plotOn(plotR, Name("combBkg"), Components("combBkgR"), FillColor(42), LineColor(42), DrawOption("F"), MoveToBack());
  modelResonant->plotOn(plotR, Name("fitPartialBkg"), Components("partial_Bkg"), FillColor(40), LineColor(40), DrawOption("F"), AddTo("combBkg"));
  modelResonant->plotOn(plotR, Name("combBkg"), Components("combBkgR"), FillColor(42), LineColor(42), DrawOption("F"), MoveToBack());
  modelResonant->plotOn(plotR, Name("fitSig"), Components("smodelR"),LineColor(kRed+1));
  data[1]->plotOn(plotR, Binning(50));

  TCanvas* ccR = new TCanvas();
  ccR->SetLogy(0);
  plotR->Draw();
  TPaveText pt(0.7,0.70,0.9,0.90,"nbNDC");
  pt.SetFillColor(0);
  pt.SetLineWidth(0);
  pt.SetBorderSize(1);
  pt.SetTextFont(42);
  pt.SetTextSize(0.04);
  pt.SetTextAlign(12);
  pt.AddText(Form("MVA cut = %.2f", mvaCut));
  pt.AddText(Form("S = %.0f #pm %.f", nSignalInsignalRegion, w.var("nsigR")->getError()));
  pt.AddText(Form("B = %.0f #pm %.f", nBkgInsignalRegion, w.var("nbkgR")->getError()));
  pt.AddText(Form("S/#sqrt{S+B} = %.1f", nSignalInsignalRegion / sqrt(nSignalInsignalRegion+nBkgInsignalRegion)));
  pt.Draw("same");
  ccR->Update();

  TLegend tl(0.70,0.5,0.90,0.70);
  tl.AddEntry(plotR->findObject("fitPartialBkg"), "partial","l");
  tl.AddEntry(plotR->findObject("fitBkg"), "all bkg","l");
  tl.AddEntry(plotR->findObject("fitSig"), "signal","l");
  tl.SetTextFont(42);
  tl.SetTextSize(0.04);
  tl.Draw("same");

  ccR->Print("modelResonant.png");
  ccR->Print("modelResonant.pdf");

  //Lowq2 
  w.factory("nbkgLowq2[500, 0., 1.e5]");

  w.factory("Exponential::combBkgL(x, tauL[-3.0, -100.0, -1.0e-5])");
  w.import("part_workspace.root:myPartialWorkSpace:partial", RenameAllNodes("BkgL"));

  w.factory("SUM::modelBkgL(f2[0.5,0.,1.] * partial_BkgL, combBkgL)");
  w.factory("SUM::modelBkgLowq2(nbkgLowq2 * modelBkgL)");

  RooAbsPdf * modelBkgLowq2 = w.pdf("modelBkgLowq2");


  //fit Bmass in lowq2 category
  w.var("x")->setRange("Slow", Blow_blindmin, Blow_sideband);
  w.var("x")->setRange("Shigh", Bhigh_sideband, Bhigh_blindmax);

  RooFitResult* rL = modelBkgLowq2->fitTo(*(data[0]), Extended(true), Minimizer("Minuit2"), Save(true), Range("Slow,Shigh"));

  RooAbsReal* backgroundFractionL = w.pdf("modelBkgLowq2")->createIntegral(*(w.var("x")), RooFit::NormSet(*(w.var("x"))), Range("signalRegioncut"));
  float nBkgInsignalRegionL = backgroundFractionL->getVal() * w.var("nbkgLowq2")->getVal();

  RooPlot* plotL = w.var("x")->frame();
  plotL->SetXTitle("Kee mass (GeV)");
  plotL->SetTitle("");
  data[0]->plotOn(plotL, Binning(50));
  modelBkgLowq2->plotOn(plotL);
  modelBkgLowq2->plotOn(plotL, Name("fitBkg"), Components("modelBkgLowq2"),LineStyle(kDashed), LineColor(kBlue));
  modelBkgLowq2->plotOn(plotL, Name("fitPartialBkgL"), Components("partial_BkgL"), FillColor(40), LineColor(40), DrawOption("F"));
  data[0]->plotOn(plotL, Binning(50));


  TCanvas* ccL = new TCanvas();
  ccL->SetLogy(0);
  plotL->Draw();
  TPaveText ptL(0.7,0.70,0.9,0.90,"NBNDC");
  ptL.SetFillColor(0);
  ptL.SetLineWidth(0);
  ptL.SetBorderSize(1);
  ptL.SetTextFont(42);
  ptL.SetTextSize(0.04);
  ptL.SetTextAlign(12);
  ptL.AddText(Form("MVA cut = %.2f", mvaCut));
  ptL.AddText(Form("B = %.0f #pm %.f", nBkgInsignalRegionL, w.var("nbkgLowq2")->getError()));
  ptL.Draw("same");
  ccL->Update();

  TLegend tlL(0.70,0.5,0.90,0.70);
  tlL.AddEntry(plotL->findObject("fitPartialBkgL"), "partial","l");
  tlL.AddEntry(plotL->findObject("fitBkg"), "all bkg","l");
  tlL.SetTextFont(42);
  tlL.SetTextSize(0.04);
  tlL.Draw("same");

  ccL->Print("modelBkgLowq2.png");
  ccL->Print("modelBkgLowq2.pdf");


  //now build the full lowq2 model with toys: dataset from combined signal(reascaled from resonant) + bkg
  RooRealVar nResonantEvents ("nResonantEvents", "", nSignalInsignalRegion);
  RooRealVar BR_Kll ("BR_Kll", "", 4.51e-7);
  RooRealVar BR_KJPsill ("BR_KJPsill", "", 1.01e-3 * 0.0597);

  RooRealVar AxE_KJPsill_mll_2p4_4 ("AxE_KJPsill_mll_2p4_4", "", wSignal->var("axeVal_mll_2.4-4.0")->getValV());
  RooRealVar AxE_KJPsill_mll_1p1_4 ("AxE_KJPsill_mll_1p1_4", "", wSignal->var("axeVal_mll_1.1-4.0")->getValV());
  RooRealVar AxE_Kll_mll_1p1_2p4 ("AxE_Kll_mll_1p1_2p4", "", wSignalnnR->var("axeVal_mll_1.1-2.4")->getValV());
  RooRealVar AxE_Kll_mll_1p1_4 ("AxE_Kll_mll_1p1_4", "", wSignalnnR->var("axeVal_mll_1.1-4.0")->getValV());

  //  RooFormulaVar nsigLowq2("nsigLowq2", "", "@0 * @1/@2 * @3/@4 / (@5/@6) ", RooArgList(*parNSignalResonant, BR_Kll, BR_KJPsill, AxE_Kll_mll_1p1_2p4, AxE_Kll_mll_1p1_4, AxE_KJPsill_mll_2p4_4, AxE_KJPsill_mll_1p1_4));

  RooFormulaVar nsigLowq2("nsigLowq2", "", "@0 * @1/@2 * @3/@4 ", RooArgList(nResonantEvents, BR_Kll, BR_KJPsill, AxE_Kll_mll_1p1_2p4, AxE_KJPsill_mll_2p4_4));

  std::cout << " nResonant = " << nResonantEvents.getVal() << " AxE = " << AxE_Kll_mll_1p1_2p4.getVal()/AxE_KJPsill_mll_2p4_4.getVal()
	    << " BR ratio = " << BR_Kll.getVal() / BR_KJPsill.getVal() << " resonant to nn resonant = " << nsigLowq2.getVal() << std::endl;
  std::cout << "\n " << std::endl;
  w.import(nsigLowq2);


  w.factory("RooDoubleCBFast::smodelL(x, smodelR_mu, smodelR_sigma, smodelR_aL, smodelR_nL, smodelR_aR, smodelR_nR)");
  w.factory("SUM::modelToy(nsigLowq2 * smodelL, nbkgLowq2 * modelBkgLowq2)");
  RooAbsPdf* modelToy = w.pdf("modelToy");

  //generate toy roodataset 
  RooDataSet* toyDataRaw = modelToy->generate(RooArgSet(*(w.var("x"))), 100 * w.var("nbkgLowq2")->getVal());
  //weight for events(want same number of bkg events but reducing the stat fluctuations)
  RooFormulaVar evtWeight("evtWeight", "", " 0.01", *(w.var("x")));
  RooRealVar* evtWeightVal = (RooRealVar*) toyDataRaw->addColumn(evtWeight) ;
  RooDataSet toyData(toyDataRaw->GetName(), toyDataRaw->GetTitle(), toyDataRaw, *toyDataRaw->get(), 0, evtWeightVal->GetName());

  RooPlot* plotT = w.var("x")->frame();
  plotT->SetXTitle("Kee mass (GeV) - toy");
  plotT->SetTitle("");
  toyData.plotOn(plotT, Binning(50));

  modelToy->fitTo(toyData, Extended(true), Minimizer("Minuit2"), Save(true), SumW2Error(kTRUE));
  modelToy->plotOn(plotT);


  TCanvas* ccT = new TCanvas();
  ccT->SetLogy(0);
  plotT->Draw();
  ccT->Print("toyData.png");
  ccT->Print("toyData.pdf");


  std::cout << " fit toy nbkgLowq2 = " << w.var("nbkgLowq2")->getVal() << " +/- " <<  w.var("nbkgLowq2")->getError() << std::endl;
  std::cout << "\n\n " << std::endl;

  //fit simultaneously resonant (data) and non resonant (toy) to have correct estimate of errors
  //Define categories
  w.factory("nsigL[5.e3, 0., 1.e4]");
  w.factory("nbkgL[1.e4, 0., 1.e5]");
  w.factory("SUM::modelLowq2(nsigL * smodelL, nbkgL * modelBkgL)");
  RooAbsPdf* modelLowq2 = w.pdf("modelLowq2");


  RooCategory sample("sample", "sample");
  sample.defineType("lowq2");
  sample.defineType("resonant");

  //Construct combined dataset in (x, sample) with weights
  RooRealVar combW("combW", "", 0, 1);  

  RooArgSet varsPlusCat(*(w.var("x"))); varsPlusCat.add(sample);
  RooArgSet varsPlusWeight(varsPlusCat); varsPlusWeight.add(combW);
  RooDataSet combData("combData", "", varsPlusWeight, WeightVar(combW));
 
  sample.setLabel("lowq2");
  for(int i=0, n=toyData.numEntries(); i<n; ++i){
    varsPlusCat = *toyData.get(i);
    combData.add(varsPlusCat, toyData.weight());
  }
  sample.setLabel("resonant");
  for(int i=0, n=data[1]->numEntries(); i<n; ++i){
    varsPlusCat = *data[1]->get(i);
    combData.add(varsPlusCat, data[1]->weight());
  }

  //combined dataset without weights
  //RooDataSet combData("combData", "", *(w.var("x")), RooFit::Index(sample), RooFit::Import("lowq2", toyDataW), RooFit::Import("resonant", resonantData));

  //Construct a simultaneous pdf in (x, sample)
  //with category sample as index
  RooSimultaneous simPdf("simPdf", "", sample);

  //Associate each model to the corresponding sample
  simPdf.addPdf(*(w.pdf("modelLowq2")), "lowq2");
  simPdf.addPdf(*(w.pdf("modelResonant")), "resonant");

  RooFitResult * r = simPdf.fitTo(combData, Extended(true), Minimizer("Minuit2"),Save(true), SumW2Error(kTRUE));

  RooAbsReal* comb_backgroundFractionR = w.pdf("modelBkg")->createIntegral(*(w.var("x")), RooFit::NormSet(*(w.var("x"))), Range("signalRegioncut"));
  RooAbsReal* comb_signalFractionR = w.pdf("smodelR")->createIntegral(*(w.var("x")), RooFit::NormSet(*(w.var("x"))), Range("signalRegioncut"));
  RooAbsReal* comb_backgroundFractionL = w.pdf("modelBkgL")->createIntegral(*(w.var("x")), RooFit::NormSet(*(w.var("x"))), Range("signalRegioncut"));
  RooAbsReal* comb_signalFractionL = w.pdf("smodelL")->createIntegral(*(w.var("x")), RooFit::NormSet(*(w.var("x"))), Range("signalRegioncut"));

  float comb_nBkgInsignalRegionR = comb_backgroundFractionR->getVal() * w.var("nbkgR")->getVal();
  float comb_nSignalInsignalRegionR = comb_signalFractionR->getVal() * w.var("nsigR")->getVal();
  float comb_nBkgInsignalRegionL = comb_backgroundFractionR->getVal() * w.var("nbkgL")->getVal();
  float comb_nSignalInsignalRegionL = comb_signalFractionR->getVal() * w.var("nsigL")->getVal();

  std::cout << " resonant q2: nSignal = " << comb_nSignalInsignalRegionR << " nBkg = " << comb_nBkgInsignalRegionR 
	    << " S/sqrt(S+B) = " << comb_nSignalInsignalRegionR / sqrt(comb_nSignalInsignalRegionR+comb_nBkgInsignalRegionR) << std::endl;
  std::cout << " low q2: nSignal = " << comb_nSignalInsignalRegionL << " nBkg = " << comb_nBkgInsignalRegionL 
	    << " S/sqrt(S+B) = " << comb_nSignalInsignalRegionL / sqrt(comb_nSignalInsignalRegionL+comb_nBkgInsignalRegionL) << std::endl;
  std::cout << "\n\n " << std::endl;



  RooPlot* frame1 = w.var("x")->frame();
  frame1->SetXTitle("Kee mass (GeV) - toy");
  frame1->SetTitle("");
  RooPlot* frame2 = w.var("x")->frame();
  frame2->SetXTitle("K(JPsi)ee mass (GeV)");
  frame2->SetTitle("");

  combData.plotOn(frame1, Binning(50), RooFit::Cut("sample==sample::lowq2"));
  combData.plotOn(frame2, Binning(50), RooFit::Cut("sample==sample::resonant"));

  simPdf.plotOn(frame1, RooFit::Slice(sample, "lowq2"), RooFit::ProjWData(RooArgSet(sample), combData));  
  simPdf.plotOn(frame1, RooFit::Slice(sample, "lowq2"), RooFit::ProjWData(RooArgSet(sample), combData), Name("fitBkg"), Components("modelBkgL"), LineStyle(kDashed), LineColor(kBlue));
  simPdf.plotOn(frame1, RooFit::Slice(sample, "lowq2"), RooFit::ProjWData(RooArgSet(sample), combData), Name("combBkg"), Components("combBkgL"), FillColor(42), LineColor(42), MoveToBack(), DrawOption("F"));
  simPdf.plotOn(frame1, RooFit::Slice(sample, "lowq2"), RooFit::ProjWData(RooArgSet(sample), combData), Name("fitPartialBkg"), Components("partial_BkgL"), FillColor(40), LineColor(40), DrawOption("F"), AddTo("combBkg"));
  simPdf.plotOn(frame1, RooFit::Slice(sample, "lowq2"), RooFit::ProjWData(RooArgSet(sample), combData), Name("combBkg"), Components("combBkgL"), FillColor(42), LineColor(42), MoveToBack(), DrawOption("F"));
  simPdf.plotOn(frame1, RooFit::Slice(sample, "lowq2"), RooFit::ProjWData(RooArgSet(sample), combData), Name("fitSig"), Components("smodelL"), LineColor(kRed+1));
  combData.plotOn(frame1, Binning(50), RooFit::Cut("sample==sample::lowq2"));


  simPdf.plotOn(frame2, RooFit::Slice(sample, "resonant"), RooFit::ProjWData(RooArgSet(sample), combData));
  simPdf.plotOn(frame2, RooFit::Slice(sample, "resonant"), RooFit::ProjWData(RooArgSet(sample), combData), Name("fitBkg"), Components("modelBkg"), LineStyle(kDashed), LineColor(kBlue));
  simPdf.plotOn(frame2, RooFit::Slice(sample, "resonant"), RooFit::ProjWData(RooArgSet(sample), combData), Name("combBkg"), Components("combBkgR"), FillColor(42), LineColor(42), MoveToBack(), DrawOption("F"));
  simPdf.plotOn(frame2, RooFit::Slice(sample, "resonant"), RooFit::ProjWData(RooArgSet(sample), combData), Name("fitPartialBkg"), Components("partial_Bkg"), FillColor(40), LineColor(40), DrawOption("F"), AddTo("combBkg"));
  simPdf.plotOn(frame2, RooFit::Slice(sample, "resonant"), RooFit::ProjWData(RooArgSet(sample), combData), Name("combBkg"), Components("combBkgR"), FillColor(42), LineColor(42), MoveToBack(), DrawOption("F"));
  simPdf.plotOn(frame2, RooFit::Slice(sample, "resonant"), RooFit::ProjWData(RooArgSet(sample), combData), Name("fitSig"), Components("smodelR"), LineColor(kRed+1));
  combData.plotOn(frame2, Binning(50), RooFit::Cut("sample==sample::resonant"));

  

  TCanvas* tc = new TCanvas("tc", "", 1300, 600);
  tc->Divide(2);
  tc->cd(1);
  //gPad->SetLeftMargin(0.15);
  //  frame1->GetYaxis()->SetTitleOffset(1.4);
  frame1->Draw();

  TPaveText ptCombL(0.65,0.70,0.85,0.90,"nbNDC");
  ptCombL.SetFillColor(0);
  ptCombL.SetLineWidth(0);
  ptCombL.SetBorderSize(1);
  ptCombL.SetTextFont(42);
  ptCombL.SetTextSize(0.04);
  ptCombL.SetTextAlign(12);
  ptCombL.AddText(Form("MVA cut = %.2f", mvaCut));
  ptCombL.AddText(Form("S = %.0f #pm %.f", comb_nSignalInsignalRegionL, w.var("nsigL")->getError()));
  ptCombL.AddText(Form("B = %.0f #pm %.f", comb_nBkgInsignalRegionL, w.var("nbkgL")->getError()));
  ptCombL.AddText(Form("S/#sqrt{S+B} = %.1f", comb_nSignalInsignalRegionL / sqrt(comb_nSignalInsignalRegionL+comb_nBkgInsignalRegionL)));
  ptCombL.Draw("same");
  tc->Update();
  TLegend tlCombL(0.70,0.5,0.90,0.70);
  tlCombL.AddEntry(frame1->findObject("combBkg"), "combinatorial","l");
  tlCombL.AddEntry(frame1->findObject("fitPartialBkg"), "partial","l");
  tlCombL.AddEntry(frame1->findObject("fitBkg"), "all bkg","l");
  tlCombL.AddEntry(frame1->findObject("fitSig"), "signal","l");
  tlCombL.SetTextFont(42);
  tlCombL.SetTextSize(0.04);
  tlCombL.Draw("same");


  tc->cd(2);
  //gPad->SetLeftMargin(0.15);
  //  frame2->GetYaxis()->SetTitleOffset(1.4);
  frame2->Draw();

  TPaveText ptCombR(0.65,0.70,0.85,0.90,"nbNDC");
  ptCombR.SetFillColor(0);
  ptCombR.SetLineWidth(0);
  ptCombR.SetBorderSize(1);
  ptCombR.SetTextFont(42);
  ptCombR.SetTextSize(0.04);
  ptCombR.SetTextAlign(12);
  ptCombR.AddText(Form("MVA cut = %.2f", mvaCut));
  ptCombR.AddText(Form("S = %.0f #pm %.f", comb_nSignalInsignalRegionR, w.var("nsigR")->getError()));
  ptCombR.AddText(Form("B = %.0f #pm %.f", comb_nBkgInsignalRegionR, w.var("nbkgR")->getError()));
  ptCombR.AddText(Form("S/#sqrt{S+B} = %.1f", comb_nSignalInsignalRegionR / sqrt(comb_nSignalInsignalRegionR+comb_nBkgInsignalRegionR)));
  ptCombR.Draw("same");
  tc->Update();
  TLegend tlCombR(0.70,0.5,0.90,0.70);
  tlCombR.AddEntry(frame2->findObject("combBkg"), "combinatorial","l");
  tlCombR.AddEntry(frame2->findObject("fitPartialBkg"), "partial","l");
  tlCombR.AddEntry(frame2->findObject("fitBkg"), "all bkg","l");
  tlCombR.AddEntry(frame2->findObject("fitSig"), "signal","l");
  tlCombR.SetTextFont(42);
  tlCombR.SetTextSize(0.04);
  tlCombR.Draw("same");

  tc->Print("combined.png");
  tc->Print("combined.pdf");

  
}
