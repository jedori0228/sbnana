// Make a plot with cuts
#include "sbnana/CAFAna/Core/Binning.h"
#include "sbnana/CAFAna/Core/Spectrum.h"
#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Systs/SBNWeightSysts.h"
#include "sbnana/CAFAna/Core/EnsembleRatio.h"
#include "sbnana/CAFAna/Core/EnsembleSpectrum.h"

using namespace ana;

#include "sbnana/CAFAna/StandardRecord/Proxy/SRProxy.h"

#include "TCanvas.h"
#include "TH1.h"
#include "TFile.h"
#include "TString.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TSystem.h"

#include "canvas_margin.h"

#include "myMuonSelection.h"
#include "myProtonSelection.h"
#include "myEstimator.h"
#include "myTruth.h"
#include "myEventSelection.h"
#include "myTool.h"

#include "myFilelist.h"
#include "myFilelistShort.h"
#include "myFilelistReCAF.txt"

double TargetPOT = 6.6e20;
TString str_TargetPOT = "6.6e20 POT";
void runByParam(int xxx=0);

void myAnalyzer(){

  //runByParam(0);
  //runByParam(1);
  //runByParam(2);
  //runByParam(3);
  runByParam(4);

}

void runByParam(int xxx){

  setTDRStyle();

  std::string inputName = "";

  vector<string> goodInputs = removeZombie(inputFilesReCAF);
  SpectrumLoader loader(goodInputs);

  //SpectrumLoader loader(inputFiles);
  //SpectrumLoader loader(inputFilesShort);
  //SpectrumLoader loader(inputName);
  //SpectrumLoader loader("/pnfs/icarus/persistent/users/icaruspro/SBNworkshopApril2021/CAF/corsika_nue_BNB/0/0/*.root");
  //SpectrumLoader loader("/pnfs/icarus/persistent/users/icaruspro/SBNworkshopApril2021/CAF/corsika_nue_BNB/0/0/caf-01608068-6435-4a36-93b5-29ead574d963.root");

  Cut cutMyCut = kIsNuMuCC && cutHasMuonTrack && cutIsMuonTrackLong && cutHasProtonTrack && cutmyFMScore; // TODO cutmyFMScore test

  //==== cutTruthNoNeutron
  //==== cutTruthNoChargedPion
  //==== cutTruthNoPiZero
  TString outputDir = "./output/myAnalyzer/BindingE40MeV/";
  if(xxx==0){
    cutMyCut = cutMyCut;
    outputDir = outputDir+"/";
  }
  else if(xxx==1){
    cutMyCut = cutMyCut && cutTruthNoChargedPion;
    outputDir = outputDir+"/NoChargedPion/";
  }
  else if(xxx==2){
    cutMyCut = cutMyCut && cutTruthNoNeutron;
    outputDir = outputDir+"/NoNeutron/";
  }
  else if(xxx==3){
    cutMyCut = cutMyCut && cutTruthNoPiZero;
    outputDir = outputDir+"/NoPiZero/";
  }
  else if(xxx==4){
    cutMyCut = cutMyCut && cutTruthNoChargedPion && cutTruthNoNeutron && cutTruthNoPiZero;
    outputDir = outputDir+"/ProtonSel__ProtonChi2LT60/NoNeutron_NoChargedPion_NoPiZero/ReCAF/";
  }
  //==== For the efficiency
  else if(xxx==5){
    cutMyCut = kNoCut;
    outputDir = outputDir+"/NoNeutron_NoChargedPion_NoPiZero/";
  }
  //outputDir = "./output/Null/";
  gSystem->mkdir(outputDir, kTRUE);
  cout << "[runByParam] outputDir = " << outputDir << endl;

  TString outputName = "output.root";
  TFile *outputfile = new TFile(outputDir+outputName,"RECREATE");

  //==== some binnings
  double xEnergyMin = 0.; // GeV
  double xEnergyMax = 10.;
  double dxEnergy = 0.1;
  double dxEnergy2 = 0.05;
  const Binning binsEnergy = Binning::Simple( int( (xEnergyMax-xEnergyMin)/dxEnergy ), xEnergyMin, xEnergyMax );
  const Binning binsEnergy2 = Binning::Simple( int( (xEnergyMax-xEnergyMin)/dxEnergy2 ), xEnergyMin, xEnergyMax );
  const Binning binsEResidual = Binning::Simple(600, -3., 3.);
  const Binning binsEResidualFraction = Binning::Simple(600, -3., 3.);
  const Binning binsPDG = Binning::Simple(6000, -3000, 3000);

  //==== HistAxis to be converted to Spectrum
  //==== Muon
  //====   Reco
  const HistAxis axMuonTrackRangeP("MuonTrackRangeP", binsEnergy, varMuonTrackRangeP);
  const HistAxis axMuonTrackMCSP("MuonTrackMCSP", binsEnergy, varMuonTrackMCSP);
  const HistAxis axMuonTrackCombinedP("MuonTrackCombinedP", binsEnergy, varMuonTrackCombinedP);
  const HistAxis axMuonTrackCaloP("MuonTrackCaloP", binsEnergy, varMuonTrackCaloP);
  const HistAxis axMuonTrackPlane0CaloP("MuonTrackPlane0CaloP", binsEnergy, varMuonTrackPlane0CaloP);
  const HistAxis axMuonTrackPlane1CaloP("MuonTrackPlane1CaloP", binsEnergy, varMuonTrackPlane1CaloP);
  const HistAxis axMuonTrackPlane2CaloP("MuonTrackPlane2CaloP", binsEnergy, varMuonTrackPlane2CaloP);
  //====   Truth
  const HistAxis axMuonTrackMatchedTruthP("MuonTrackMatchedTruthP", binsEnergy, varMuonTrackMatchedTruthP);
  const HistAxis axMuonTrackMatchedTruthPDG("MuonTrackMatchedTruthPDG", binsPDG, varMuonTrackMatchedTruthPDG);
  //====   Residual
  const HistAxis axMuonTrackRangePResidual("MuonTrackRangePResidual", binsEResidual, varMuonTrackRangePResidual);
  const HistAxis axMuonTrackMCSPResidual("MuonTrackMCSPResidual", binsEResidual, varMuonTrackMCSPResidual);
  const HistAxis axMuonTrackCombinedPResidual("MuonTrackCombinedPResidual", binsEResidual, varMuonTrackCombinedPResidual);
  const HistAxis axMuonTrackCaloPResidual("MuonTrackCaloPResidual", binsEResidual, varMuonTrackCaloPResidual);
  const HistAxis axMuonTrackPlane0CaloPResidual("MuonTrackPlane0CaloPResidual", binsEResidual, varMuonTrackPlane0CaloPResidual);
  const HistAxis axMuonTrackPlane1CaloPResidual("MuonTrackPlane1CaloPResidual", binsEResidual, varMuonTrackPlane1CaloPResidual);
  const HistAxis axMuonTrackPlane2CaloPResidual("MuonTrackPlane2CaloPResidual", binsEResidual, varMuonTrackPlane2CaloPResidual);
  //====   ResidualFraction
  const HistAxis axMuonTrackRangePResidualFraction("MuonTrackRangePResidualFraction", binsEResidualFraction, varMuonTrackRangePResidualFraction);
  const HistAxis axMuonTrackMCSPResidualFraction("MuonTrackMCSPResidualFraction", binsEResidualFraction, varMuonTrackMCSPResidualFraction);
  const HistAxis axMuonTrackCombinedPResidualFraction("MuonTrackCombinedPResidualFraction", binsEResidualFraction, varMuonTrackCombinedPResidualFraction);
  const HistAxis axMuonTrackCaloPResidualFraction("MuonTrackCaloPResidualFraction", binsEResidualFraction, varMuonTrackCaloPResidualFraction);
  const HistAxis axMuonTrackPlane0CaloPResidualFraction("MuonTrackPlane0CaloPResidualFraction", binsEResidualFraction, varMuonTrackPlane0CaloPResidualFraction);
  const HistAxis axMuonTrackPlane1CaloPResidualFraction("MuonTrackPlane1CaloPResidualFraction", binsEResidualFraction, varMuonTrackPlane1CaloPResidualFraction);
  const HistAxis axMuonTrackPlane2CaloPResidualFraction("MuonTrackPlane2CaloPResidualFraction", binsEResidualFraction, varMuonTrackPlane2CaloPResidualFraction);
  //==== Proton
  //====   Reco
  const HistAxis axProtonTrackRangeP("ProtonTrackRangeP", binsEnergy, varProtonTrackRangeP);
  const HistAxis axProtonTrackMCSP("ProtonTrackMCSP", binsEnergy, varProtonTrackMCSP);
  const HistAxis axProtonTrackCombinedP("ProtonTrackCombinedP", binsEnergy, varProtonTrackCombinedP);
  const HistAxis axProtonTrackCaloP("ProtonTrackCaloP", binsEnergy, varProtonTrackCaloP);
  const HistAxis axProtonTrackPlane0CaloP("ProtonTrackPlane0CaloP", binsEnergy, varProtonTrackPlane0CaloP);
  const HistAxis axProtonTrackPlane1CaloP("ProtonTrackPlane1CaloP", binsEnergy, varProtonTrackPlane1CaloP);
  const HistAxis axProtonTrackPlane2CaloP("ProtonTrackPlane2CaloP", binsEnergy, varProtonTrackPlane2CaloP);
  //====   Truth
  const HistAxis axProtonTrackMatchedTruthP("ProtonTrackMatchedTruthP", binsEnergy, varProtonTrackMatchedTruthP);
  const HistAxis axProtonTrackMatchedTruthPDG("ProtonTrackMatchedTruthPDG", binsPDG, varProtonTrackMatchedTruthPDG);
  //====   Residual
  const HistAxis axProtonTrackRangePResidual("ProtonTrackRangePResidual", binsEResidual, varProtonTrackRangePResidual);
  const HistAxis axProtonTrackMCSPResidual("ProtonTrackMCSPResidual", binsEResidual, varProtonTrackMCSPResidual);
  const HistAxis axProtonTrackCombinedPResidual("ProtonTrackCombinedPResidual", binsEResidual, varProtonTrackCombinedPResidual);
  const HistAxis axProtonTrackCaloPResidual("ProtonTrackCaloPResidual", binsEResidual, varProtonTrackCaloPResidual);
  const HistAxis axProtonTrackPlane0CaloPResidual("ProtonTrackPlane0CaloPResidual", binsEResidual, varProtonTrackPlane0CaloPResidual);
  const HistAxis axProtonTrackPlane1CaloPResidual("ProtonTrackPlane1CaloPResidual", binsEResidual, varProtonTrackPlane1CaloPResidual);
  const HistAxis axProtonTrackPlane2CaloPResidual("ProtonTrackPlane2CaloPResidual", binsEResidual, varProtonTrackPlane2CaloPResidual);
  //====   ResidualFraction
  const HistAxis axProtonTrackRangePResidualFraction("ProtonTrackRangePResidualFraction", binsEResidualFraction, varProtonTrackRangePResidualFraction);
  const HistAxis axProtonTrackMCSPResidualFraction("ProtonTrackMCSPResidualFraction", binsEResidualFraction, varProtonTrackMCSPResidualFraction);
  const HistAxis axProtonTrackCombinedPResidualFraction("ProtonTrackCombinedPResidualFraction", binsEResidualFraction, varProtonTrackCombinedPResidualFraction);
  const HistAxis axProtonTrackCaloPResidualFraction("ProtonTrackCaloPResidualFraction", binsEResidualFraction, varProtonTrackCaloPResidualFraction);
  const HistAxis axProtonTrackPlane0CaloPResidualFraction("ProtonTrackPlane0CaloPResidualFraction", binsEResidualFraction, varProtonTrackPlane0CaloPResidualFraction);
  const HistAxis axProtonTrackPlane1CaloPResidualFraction("ProtonTrackPlane1CaloPResidualFraction", binsEResidualFraction, varProtonTrackPlane1CaloPResidualFraction);
  const HistAxis axProtonTrackPlane2CaloPResidualFraction("ProtonTrackPlane2CaloPResidualFraction", binsEResidualFraction, varProtonTrackPlane2CaloPResidualFraction);
  //==== Neutrino
  //====   Reco
  const HistAxis axNeutrinoCaloEnergy("NeutrinoCaloEnergy", binsEnergy, varNeutrinoCaloEnergy);
  const HistAxis axNeutrinoCombinedEnergy("NeutrinoCombinedEnergy", binsEnergy, varNeutrinoCombinedEnergy);
  const HistAxis axNeutrinoQE("NeutrinoQE", binsEnergy, varNeutrinoQE);
  const HistAxis axNeutrinoFakeRecoEnergy("NeutrinoFakeRecoEnergy", binsEnergy, varNeutrinoFakeRecoEnergy);
  //====   Truth
  const HistAxis axNeutrinoTruthE("NeutrinoTruthE", binsEnergy, varNeutrinoTruthE);
  //====   Residual
  const HistAxis axNeutrinoCaloEnergyResidual("NeutrinoCaloEnergyResidual", binsEResidual, varNeutrinoCaloEnergyResidual);
  const HistAxis axNeutrinoCombinedEnergyResidual("NeutrinoCombinedEnergyResidual", binsEResidual, varNeutrinoCombinedEnergyResidual);
  const HistAxis axNeutrinoQEResidual("NeutrinoQEResidual", binsEResidual, varNeutrinoQEResidual);
  //====   ResidualFraction
  const HistAxis axNeutrinoCaloEnergyResidualFraction("NeutrinoCaloEnergyResidualFraction", binsEResidualFraction, varNeutrinoCaloEnergyResidualFraction);
  const HistAxis axNeutrinoCombinedEnergyResidualFraction("NeutrinoCombinedEnergyResidualFraction", binsEResidualFraction, varNeutrinoCombinedEnergyResidualFraction);
  const HistAxis axNeutrinoQEResidualFraction("NeutrinoQEResidualFraction", binsEResidualFraction, varNeutrinoQEResidualFraction);

  //==== Spectrum
  //==== Muon
  //====   Reco
  Spectrum sMuonTrackRangeP(loader, axMuonTrackRangeP, kNoSpillCut, cutMyCut);
  Spectrum sMuonTrackMCSP(loader, axMuonTrackMCSP, kNoSpillCut, cutMyCut);
  Spectrum sMuonTrackCombinedP(loader, axMuonTrackCombinedP, kNoSpillCut, cutMyCut);
  Spectrum sMuonTrackCaloP(loader, axMuonTrackCaloP, kNoSpillCut, cutMyCut);
  Spectrum sMuonTrackPlane0CaloP(loader, axMuonTrackPlane0CaloP, kNoSpillCut, cutMyCut);
  Spectrum sMuonTrackPlane1CaloP(loader, axMuonTrackPlane1CaloP, kNoSpillCut, cutMyCut);
  Spectrum sMuonTrackPlane2CaloP(loader, axMuonTrackPlane2CaloP, kNoSpillCut, cutMyCut);
  //====   Truth
  Spectrum sMuonTrackMatchedTruthP(loader, axMuonTrackMatchedTruthP, kNoSpillCut, cutMyCut);
  Spectrum sMuonTrackMatchedTruthPDG(loader, axMuonTrackMatchedTruthPDG, kNoSpillCut, cutMyCut);
  //====   Truth vs Reco
  Spectrum sMuonTrack_MatchedTruthP_vs_CombinedP(loader, axMuonTrackMatchedTruthP, axMuonTrackCombinedP, kNoSpillCut, cutMyCut);
  //====   Residual
  Spectrum sMuonTrackRangePResidual(loader, axMuonTrackRangePResidual, kNoSpillCut, cutMyCut);
  Spectrum sMuonTrackMCSPResidual(loader, axMuonTrackMCSPResidual, kNoSpillCut, cutMyCut);
  Spectrum sMuonTrackCombinedPResidual(loader, axMuonTrackCombinedPResidual, kNoSpillCut, cutMyCut);
  Spectrum sMuonTrackCaloPResidual(loader, axMuonTrackCaloPResidual, kNoSpillCut, cutMyCut);
  Spectrum sMuonTrackPlane0CaloPResidual(loader, axMuonTrackPlane0CaloPResidual, kNoSpillCut, cutMyCut);
  Spectrum sMuonTrackPlane1CaloPResidual(loader, axMuonTrackPlane1CaloPResidual, kNoSpillCut, cutMyCut);
  Spectrum sMuonTrackPlane2CaloPResidual(loader, axMuonTrackPlane2CaloPResidual, kNoSpillCut, cutMyCut);
  //====   ResidualFraction
  Spectrum sMuonTrackRangePResidualFraction(loader, axMuonTrackRangePResidualFraction, kNoSpillCut, cutMyCut);
  Spectrum sMuonTrackMCSPResidualFraction(loader, axMuonTrackMCSPResidualFraction, kNoSpillCut, cutMyCut);
  Spectrum sMuonTrackCombinedPResidualFraction(loader, axMuonTrackCombinedPResidualFraction, kNoSpillCut, cutMyCut);
  Spectrum sMuonTrackCaloPResidualFraction(loader, axMuonTrackCaloPResidualFraction, kNoSpillCut, cutMyCut);
  Spectrum sMuonTrackPlane0CaloPResidualFraction(loader, axMuonTrackPlane0CaloPResidualFraction, kNoSpillCut, cutMyCut);
  Spectrum sMuonTrackPlane1CaloPResidualFraction(loader, axMuonTrackPlane1CaloPResidualFraction, kNoSpillCut, cutMyCut);
  Spectrum sMuonTrackPlane2CaloPResidualFraction(loader, axMuonTrackPlane2CaloPResidualFraction, kNoSpillCut, cutMyCut);
  //==== Proton
  //====   Reco
  Spectrum sProtonTrackRangeP(loader, axProtonTrackRangeP, kNoSpillCut, cutMyCut);
  Spectrum sProtonTrackMCSP(loader, axProtonTrackMCSP, kNoSpillCut, cutMyCut);
  Spectrum sProtonTrackCombinedP(loader, axProtonTrackCombinedP, kNoSpillCut, cutMyCut);
  Spectrum sProtonTrackCaloP(loader, axProtonTrackCaloP, kNoSpillCut, cutMyCut);
  Spectrum sProtonTrackPlane0CaloP(loader, axProtonTrackPlane0CaloP, kNoSpillCut, cutMyCut);
  Spectrum sProtonTrackPlane1CaloP(loader, axProtonTrackPlane1CaloP, kNoSpillCut, cutMyCut);
  Spectrum sProtonTrackPlane2CaloP(loader, axProtonTrackPlane2CaloP, kNoSpillCut, cutMyCut);
  //====   Truth
  Spectrum sProtonTrackMatchedTruthP(loader, axProtonTrackMatchedTruthP, kNoSpillCut, cutMyCut);
  Spectrum sProtonTrackMatchedTruthPDG(loader, axProtonTrackMatchedTruthPDG, kNoSpillCut, cutMyCut);
  //====   Truth vs Reco
  Spectrum sProtonTrack_MatchedTruthP_vs_CaloP(loader, axProtonTrackMatchedTruthP, axProtonTrackCaloP, kNoSpillCut, cutMyCut);
  //====   Residual
  Spectrum sProtonTrackRangePResidual(loader, axProtonTrackRangePResidual, kNoSpillCut, cutMyCut);
  Spectrum sProtonTrackMCSPResidual(loader, axProtonTrackMCSPResidual, kNoSpillCut, cutMyCut);
  Spectrum sProtonTrackCombinedPResidual(loader, axProtonTrackCombinedPResidual, kNoSpillCut, cutMyCut);
  Spectrum sProtonTrackCaloPResidual(loader, axProtonTrackCaloPResidual, kNoSpillCut, cutMyCut);
  Spectrum sProtonTrackPlane0CaloPResidual(loader, axProtonTrackPlane0CaloPResidual, kNoSpillCut, cutMyCut);
  Spectrum sProtonTrackPlane1CaloPResidual(loader, axProtonTrackPlane1CaloPResidual, kNoSpillCut, cutMyCut);
  Spectrum sProtonTrackPlane2CaloPResidual(loader, axProtonTrackPlane2CaloPResidual, kNoSpillCut, cutMyCut);
  //====   ResidualFraction
  Spectrum sProtonTrackRangePResidualFraction(loader, axProtonTrackRangePResidualFraction, kNoSpillCut, cutMyCut);
  Spectrum sProtonTrackMCSPResidualFraction(loader, axProtonTrackMCSPResidualFraction, kNoSpillCut, cutMyCut);
  Spectrum sProtonTrackCombinedPResidualFraction(loader, axProtonTrackCombinedPResidualFraction, kNoSpillCut, cutMyCut);
  Spectrum sProtonTrackCaloPResidualFraction(loader, axProtonTrackCaloPResidualFraction, kNoSpillCut, cutMyCut);
  Spectrum sProtonTrackPlane0CaloPResidualFraction(loader, axProtonTrackPlane0CaloPResidualFraction, kNoSpillCut, cutMyCut);
  Spectrum sProtonTrackPlane1CaloPResidualFraction(loader, axProtonTrackPlane1CaloPResidualFraction, kNoSpillCut, cutMyCut);
  Spectrum sProtonTrackPlane2CaloPResidualFraction(loader, axProtonTrackPlane2CaloPResidualFraction, kNoSpillCut, cutMyCut);
  //==== Neutrino
  //====   Reco
  Spectrum sNeutrinoCaloEnergy(loader, axNeutrinoCaloEnergy, kNoSpillCut, cutMyCut);
  Spectrum sNeutrinoCombinedEnergy(loader, axNeutrinoCombinedEnergy, kNoSpillCut, cutMyCut);
  Spectrum sNeutrinoQE(loader, axNeutrinoQE, kNoSpillCut, cutMyCut);
  Spectrum sNeutrinoFakeRecoEnergy(loader, axNeutrinoFakeRecoEnergy, kNoSpillCut, cutMyCut);
  //====   Truth
  Spectrum sNeutrinoTruthE(loader, axNeutrinoTruthE, kNoSpillCut, cutMyCut);
  //====   Truth vs Reco
  Spectrum sNeutrino_TruthE_vs_CombinedEnergy(loader, axNeutrinoTruthE, axNeutrinoCombinedEnergy, kNoSpillCut, cutMyCut);
  Spectrum sNeutrino_TruthE_vs_QE(loader, axNeutrinoTruthE, axNeutrinoQE, kNoSpillCut, cutMyCut);
  //====   Residual
  Spectrum sNeutrinoCaloEnergyResidual(loader, axNeutrinoCaloEnergyResidual, kNoSpillCut, cutMyCut);
  Spectrum sNeutrinoCombinedEnergyResidual(loader, axNeutrinoCombinedEnergyResidual, kNoSpillCut, cutMyCut);
  Spectrum sNeutrinoQEResidual(loader, axNeutrinoQEResidual, kNoSpillCut, cutMyCut);
  //====   ResidualFraction
  Spectrum sNeutrinoCaloEnergyResidualFraction(loader, axNeutrinoCaloEnergyResidualFraction, kNoSpillCut, cutMyCut);
  Spectrum sNeutrinoCombinedEnergyResidualFraction(loader, axNeutrinoCombinedEnergyResidualFraction, kNoSpillCut, cutMyCut);
  Spectrum sNeutrinoQEResidualFraction(loader, axNeutrinoQEResidualFraction, kNoSpillCut, cutMyCut);

/*
  //==== Syst
  std::vector<Var> weis;
  const std::vector<const ISyst*>& systs = GetSBNGenieWeightSysts();
  weis.reserve(systs.size());
  for(long unsigned int i = 0; i<systs.size(); ++i){
    weis.push_back(GetUniverseWeight(systs, i));
  }
  EnsembleSpectrum esTruthEnergy(loader, axEnergy, kNoSpillCut, cutMyCut, weis);
*/

  //==== This is the call that actually fills in the spectrum
  loader.Go();

  //==== Draw

  double y_max = -999.;
  double y_min = 999999999;

  TLatex latex_POT;
  latex_POT.SetNDC();
  latex_POT.SetTextSize(0.035);

  //==== Muon

  //====   Muon momentum (P)

  y_max = -999.;
  TCanvas *c_MuonP = new TCanvas("c_MuonP", "", 800, 800);
  canvas_margin(c_MuonP);
  c_MuonP->cd();

  TH1D *hist_dummy_MuonP = new TH1D("hist_dummy_MuonP", "", int( (xEnergyMax-xEnergyMin)/dxEnergy ), xEnergyMin, xEnergyMax);
  hist_axis(hist_dummy_MuonP);
  hist_dummy_MuonP->GetXaxis()->SetTitle("Muon momentum (GeV)");
  hist_dummy_MuonP->GetXaxis()->SetRangeUser(0.,5.0);
  hist_dummy_MuonP->Draw("hist");

  TH1 *hMuonTrackRangeP = sMuonTrackRangeP.ToTH1(TargetPOT);
  hMuonTrackRangeP->SetName("hMuonTrackRangeP");
  hMuonTrackRangeP->SetLineColor(kGreen);
  hMuonTrackRangeP->SetLineWidth(2);
  y_max = max( y_max, hMuonTrackRangeP->GetMaximum() );

  TH1 *hMuonTrackMCSP = sMuonTrackMCSP.ToTH1(TargetPOT);
  hMuonTrackMCSP->SetName("hMuonTrackMCSP");
  hMuonTrackMCSP->SetLineColor(kBlue);
  hMuonTrackMCSP->SetLineWidth(2);
  y_max = max( y_max, hMuonTrackMCSP->GetMaximum() );

  TH1 *hMuonTrackCombinedP = sMuonTrackCombinedP.ToTH1(TargetPOT);
  hMuonTrackCombinedP->SetName("hMuonTrackCombinedP");
  hMuonTrackCombinedP->SetLineColor(kBlack);
  hMuonTrackCombinedP->SetLineWidth(2);
  y_max = max( y_max, hMuonTrackCombinedP->GetMaximum() );

  TH1 *hMuonTrackCaloP = sMuonTrackCaloP.ToTH1(TargetPOT);
  hMuonTrackCaloP->SetName("hMuonTrackCaloP");
  hMuonTrackCaloP->SetLineColor(kViolet);
  hMuonTrackCaloP->SetLineWidth(2);
  y_max = max( y_max, hMuonTrackCaloP->GetMaximum() );

  TH1 *hMuonTrackPlane0CaloP = sMuonTrackPlane0CaloP.ToTH1(TargetPOT);
  hMuonTrackPlane0CaloP->SetName("hMuonTrackPlane0CaloP");
  hMuonTrackPlane0CaloP->SetLineColor(kViolet);
  hMuonTrackPlane0CaloP->SetLineWidth(2);
  hMuonTrackPlane0CaloP->SetLineStyle(3);
  y_max = max( y_max, hMuonTrackPlane0CaloP->GetMaximum() );

  TH1 *hMuonTrackPlane1CaloP = sMuonTrackPlane1CaloP.ToTH1(TargetPOT);
  hMuonTrackPlane1CaloP->SetName("hMuonTrackPlane1CaloP");
  hMuonTrackPlane1CaloP->SetLineColor(kViolet);
  hMuonTrackPlane1CaloP->SetLineWidth(2);
  hMuonTrackPlane1CaloP->SetLineStyle(5);
  y_max = max( y_max, hMuonTrackPlane1CaloP->GetMaximum() );

  TH1 *hMuonTrackPlane2CaloP = sMuonTrackPlane2CaloP.ToTH1(TargetPOT);
  hMuonTrackPlane2CaloP->SetName("hMuonTrackPlane2CaloP");
  hMuonTrackPlane2CaloP->SetLineColor(kViolet);
  hMuonTrackPlane2CaloP->SetLineWidth(2);
  hMuonTrackPlane2CaloP->SetLineStyle(7);
  y_max = max( y_max, hMuonTrackPlane2CaloP->GetMaximum() );

  TH1 *hMuonTrackMatchedTruthP = sMuonTrackMatchedTruthP.ToTH1(TargetPOT);
  hMuonTrackMatchedTruthP->SetName("hMuonTrackMatchedTruthP");
  hMuonTrackMatchedTruthP->SetLineColor(kRed);
  hMuonTrackMatchedTruthP->SetLineWidth(2);
  hMuonTrackMatchedTruthP->SetLineStyle(3);
  y_max = max( y_max, hMuonTrackMatchedTruthP->GetMaximum() );

  hMuonTrackRangeP->Draw("histsame");
  hMuonTrackMCSP->Draw("histsame");
  hMuonTrackCombinedP->Draw("histsame");
  hMuonTrackCaloP->Draw("histsame");
  hMuonTrackPlane0CaloP->Draw("histsame");
  hMuonTrackPlane1CaloP->Draw("histsame");
  hMuonTrackPlane2CaloP->Draw("histsame");
  hMuonTrackMatchedTruthP->Draw("histsame");

  TLegend *lg_MuonP = new TLegend(0.35, 0.70, 0.90, 0.90);
  lg_MuonP->SetBorderSize(0);
  lg_MuonP->SetFillStyle(0);
  lg_MuonP->AddEntry(hMuonTrackRangeP, "Range", "l");
  lg_MuonP->AddEntry(hMuonTrackMCSP, "MCS", "l");
  lg_MuonP->AddEntry(hMuonTrackCombinedP, "Range (Contained) + MCS (Exiting)", "l");
  lg_MuonP->AddEntry(hMuonTrackCaloP, "Bestplane calo", "l");
  lg_MuonP->AddEntry(hMuonTrackPlane0CaloP, "Ind. 1 calo", "l");
  lg_MuonP->AddEntry(hMuonTrackPlane1CaloP, "Ind. 2 calo", "l");
  lg_MuonP->AddEntry(hMuonTrackPlane2CaloP, "Coll. calo", "l");
  lg_MuonP->AddEntry(hMuonTrackMatchedTruthP, "Truth", "l");

  lg_MuonP->Draw();

  latex_POT.DrawLatex(0.20, 0.96, str_TargetPOT);

  hist_dummy_MuonP->GetYaxis()->SetRangeUser(0., y_max*1.1);
  c_MuonP->SaveAs(outputDir+"MuonP.pdf");

  outputfile->cd();
  c_MuonP->Write();

  c_MuonP->Close();

  //====   Truth vs Reco 2D : only save TH2, not canvas
  TH2 *hMuonTrack_MatchedTruthP_vs_CombinedP = sMuonTrack_MatchedTruthP_vs_CombinedP.ToTH2(TargetPOT);
  hMuonTrack_MatchedTruthP_vs_CombinedP->SetName("hMuonTrack_MatchedTruthP_vs_CombinedP");
  outputfile->cd();
  hMuonTrack_MatchedTruthP_vs_CombinedP->Write();

  //====    PDG
  TH1 *hMuonTrackMatchedTruthPDG = sMuonTrackMatchedTruthPDG.ToTH1(TargetPOT);
  hMuonTrackMatchedTruthPDG->SetName("hMuonTrackMatchedTruthPDG");
  outputfile->cd();
  hMuonTrackMatchedTruthPDG->Write();

  //====   Residual

  y_max = -999.;
  TCanvas *c_MuonPResidual = new TCanvas("c_MuonPResidual", "", 800, 800);
  canvas_margin(c_MuonPResidual);
  c_MuonPResidual->cd();

  TH1D *hist_dummy_MuonPResidual = new TH1D("hist_dummy_MuonPResidual", "", 600, -3., 3.);
  hist_axis(hist_dummy_MuonPResidual);
  hist_dummy_MuonPResidual->GetXaxis()->SetTitle("Muon momentum (Reco-Truth) (GeV)");
  hist_dummy_MuonPResidual->Draw("hist");

  TH1 *hMuonTrackRangePResidual = sMuonTrackRangePResidual.ToTH1(TargetPOT);
  hMuonTrackRangePResidual->SetName("hMuonTrackRangePResidual");
  hMuonTrackRangePResidual->SetLineColor(kGreen);
  hMuonTrackRangePResidual->SetLineWidth(2);
  y_max = max( y_max, hMuonTrackRangePResidual->GetMaximum() );

  TH1 *hMuonTrackMCSPResidual = sMuonTrackMCSPResidual.ToTH1(TargetPOT);
  hMuonTrackMCSPResidual->SetName("hMuonTrackMCSPResidual");
  hMuonTrackMCSPResidual->SetLineColor(kBlue);
  hMuonTrackMCSPResidual->SetLineWidth(2);
  y_max = max( y_max, hMuonTrackMCSPResidual->GetMaximum() );

  TH1 *hMuonTrackCombinedPResidual = sMuonTrackCombinedPResidual.ToTH1(TargetPOT);
  hMuonTrackCombinedPResidual->SetName("hMuonTrackCombinedPResidual");
  hMuonTrackCombinedPResidual->SetLineColor(kBlack);
  hMuonTrackCombinedPResidual->SetLineWidth(2);
  y_max = max( y_max, hMuonTrackCombinedPResidual->GetMaximum() );

  TH1 *hMuonTrackCaloPResidual = sMuonTrackCaloPResidual.ToTH1(TargetPOT);
  hMuonTrackCaloPResidual->SetName("hMuonTrackCaloPResidual");
  hMuonTrackCaloPResidual->SetLineColor(kViolet);
  hMuonTrackCaloPResidual->SetLineWidth(2);
  y_max = max( y_max, hMuonTrackCaloPResidual->GetMaximum() );

  TH1 *hMuonTrackPlane0CaloPResidual = sMuonTrackPlane0CaloPResidual.ToTH1(TargetPOT);
  hMuonTrackPlane0CaloPResidual->SetName("hMuonTrackPlane0CaloPResidual");
  hMuonTrackPlane0CaloPResidual->SetLineColor(kViolet);
  hMuonTrackPlane0CaloPResidual->SetLineWidth(2);
  hMuonTrackPlane0CaloPResidual->SetLineStyle(3);
  y_max = max( y_max, hMuonTrackPlane0CaloPResidual->GetMaximum() );

  TH1 *hMuonTrackPlane1CaloPResidual = sMuonTrackPlane1CaloPResidual.ToTH1(TargetPOT);
  hMuonTrackPlane1CaloPResidual->SetName("hMuonTrackPlane1CaloPResidual");
  hMuonTrackPlane1CaloPResidual->SetLineColor(kViolet);
  hMuonTrackPlane1CaloPResidual->SetLineWidth(2);
  hMuonTrackPlane1CaloPResidual->SetLineStyle(5);
  y_max = max( y_max, hMuonTrackPlane1CaloPResidual->GetMaximum() );

  TH1 *hMuonTrackPlane2CaloPResidual = sMuonTrackPlane2CaloPResidual.ToTH1(TargetPOT);
  hMuonTrackPlane2CaloPResidual->SetName("hMuonTrackPlane2CaloPResidual");
  hMuonTrackPlane2CaloPResidual->SetLineColor(kViolet);
  hMuonTrackPlane2CaloPResidual->SetLineWidth(2);
  hMuonTrackPlane2CaloPResidual->SetLineStyle(7);
  y_max = max( y_max, hMuonTrackPlane2CaloPResidual->GetMaximum() );

  hMuonTrackRangePResidual->Draw("histsame");
  hMuonTrackMCSPResidual->Draw("histsame");
  hMuonTrackCombinedPResidual->Draw("histsame");
  hMuonTrackCaloPResidual->Draw("histsame");
  hMuonTrackPlane0CaloPResidual->Draw("histsame");
  hMuonTrackPlane1CaloPResidual->Draw("histsame");
  hMuonTrackPlane2CaloPResidual->Draw("histsame");

  TLegend *lg_MuonPResidual = new TLegend(0.35, 0.70, 0.90, 0.90);
  lg_MuonPResidual->SetBorderSize(0);
  lg_MuonPResidual->SetFillStyle(0);
  lg_MuonPResidual->AddEntry(hMuonTrackRangePResidual, "Range", "l");
  lg_MuonPResidual->AddEntry(hMuonTrackMCSPResidual, "MCS", "l");
  lg_MuonPResidual->AddEntry(hMuonTrackCombinedPResidual, "Range (Contained) + MCS (Exiting)", "l");
  lg_MuonPResidual->AddEntry(hMuonTrackCaloPResidual, "Bestplane calo", "l");
  lg_MuonPResidual->AddEntry(hMuonTrackPlane0CaloPResidual, "Ind. 1 calo", "l");
  lg_MuonPResidual->AddEntry(hMuonTrackPlane1CaloPResidual, "Ind. 2 calo", "l");
  lg_MuonPResidual->AddEntry(hMuonTrackPlane2CaloPResidual, "Coll. calo", "l");

  lg_MuonPResidual->Draw();

  latex_POT.DrawLatex(0.20, 0.96, str_TargetPOT);

  hist_dummy_MuonPResidual->GetYaxis()->SetRangeUser(0., y_max*1.1);
  c_MuonPResidual->SaveAs(outputDir+"MuonPResidual.pdf");

  outputfile->cd();
  c_MuonPResidual->Write();

  c_MuonPResidual->Close();

  //====   ResidualFraction

  y_max = -999.;
  TCanvas *c_MuonPResidualFraction = new TCanvas("c_MuonPResidualFraction", "", 800, 800);
  canvas_margin(c_MuonPResidualFraction);
  c_MuonPResidualFraction->cd();

  TH1D *hist_dummy_MuonPResidualFraction = new TH1D("hist_dummy_MuonPResidualFraction", "", 600, -3., 3.);
  hist_axis(hist_dummy_MuonPResidualFraction);
  hist_dummy_MuonPResidualFraction->GetXaxis()->SetTitle("Muon momentum (Reco-Truth)/Truth");
  hist_dummy_MuonPResidualFraction->Draw("hist");

  TH1 *hMuonTrackRangePResidualFraction = sMuonTrackRangePResidualFraction.ToTH1(TargetPOT);
  hMuonTrackRangePResidualFraction->SetName("hMuonTrackRangePResidualFraction");
  hMuonTrackRangePResidualFraction->SetLineColor(kGreen);
  hMuonTrackRangePResidualFraction->SetLineWidth(2);
  y_max = max( y_max, hMuonTrackRangePResidualFraction->GetMaximum() );

  TH1 *hMuonTrackMCSPResidualFraction = sMuonTrackMCSPResidualFraction.ToTH1(TargetPOT);
  hMuonTrackMCSPResidualFraction->SetName("hMuonTrackMCSPResidualFraction");
  hMuonTrackMCSPResidualFraction->SetLineColor(kBlue);
  hMuonTrackMCSPResidualFraction->SetLineWidth(2);
  y_max = max( y_max, hMuonTrackMCSPResidualFraction->GetMaximum() );

  TH1 *hMuonTrackCombinedPResidualFraction = sMuonTrackCombinedPResidualFraction.ToTH1(TargetPOT);
  hMuonTrackCombinedPResidualFraction->SetName("hMuonTrackCombinedPResidualFraction");
  hMuonTrackCombinedPResidualFraction->SetLineColor(kBlack);
  hMuonTrackCombinedPResidualFraction->SetLineWidth(2);
  y_max = max( y_max, hMuonTrackCombinedPResidualFraction->GetMaximum() );

  TH1 *hMuonTrackCaloPResidualFraction = sMuonTrackCaloPResidualFraction.ToTH1(TargetPOT);
  hMuonTrackCaloPResidualFraction->SetName("hMuonTrackCaloPResidualFraction");
  hMuonTrackCaloPResidualFraction->SetLineColor(kViolet);
  hMuonTrackCaloPResidualFraction->SetLineWidth(2);
  y_max = max( y_max, hMuonTrackCaloPResidualFraction->GetMaximum() );

  TH1 *hMuonTrackPlane0CaloPResidualFraction = sMuonTrackPlane0CaloPResidualFraction.ToTH1(TargetPOT);
  hMuonTrackPlane0CaloPResidualFraction->SetName("hMuonTrackPlane0CaloPResidualFraction");
  hMuonTrackPlane0CaloPResidualFraction->SetLineColor(kViolet);
  hMuonTrackPlane0CaloPResidualFraction->SetLineWidth(2);
  hMuonTrackPlane0CaloPResidualFraction->SetLineStyle(3);
  y_max = max( y_max, hMuonTrackPlane0CaloPResidualFraction->GetMaximum() );

  TH1 *hMuonTrackPlane1CaloPResidualFraction = sMuonTrackPlane1CaloPResidualFraction.ToTH1(TargetPOT);
  hMuonTrackPlane1CaloPResidualFraction->SetName("hMuonTrackPlane1CaloPResidualFraction");
  hMuonTrackPlane1CaloPResidualFraction->SetLineColor(kViolet);
  hMuonTrackPlane1CaloPResidualFraction->SetLineWidth(2);
  hMuonTrackPlane1CaloPResidualFraction->SetLineStyle(5);
  y_max = max( y_max, hMuonTrackPlane1CaloPResidualFraction->GetMaximum() );

  TH1 *hMuonTrackPlane2CaloPResidualFraction = sMuonTrackPlane2CaloPResidualFraction.ToTH1(TargetPOT);
  hMuonTrackPlane2CaloPResidualFraction->SetName("hMuonTrackPlane2CaloPResidualFraction");
  hMuonTrackPlane2CaloPResidualFraction->SetLineColor(kViolet);
  hMuonTrackPlane2CaloPResidualFraction->SetLineWidth(2);
  hMuonTrackPlane2CaloPResidualFraction->SetLineStyle(7);
  y_max = max( y_max, hMuonTrackPlane2CaloPResidualFraction->GetMaximum() );

  hMuonTrackRangePResidualFraction->Draw("histsame");
  hMuonTrackMCSPResidualFraction->Draw("histsame");
  hMuonTrackCombinedPResidualFraction->Draw("histsame");
  hMuonTrackCaloPResidualFraction->Draw("histsame");
  hMuonTrackPlane0CaloPResidualFraction->Draw("histsame");
  hMuonTrackPlane1CaloPResidualFraction->Draw("histsame");
  hMuonTrackPlane2CaloPResidualFraction->Draw("histsame");

  TLegend *lg_MuonPResidualFraction = new TLegend(0.35, 0.70, 0.90, 0.90);
  lg_MuonPResidualFraction->SetBorderSize(0);
  lg_MuonPResidualFraction->SetFillStyle(0);
  lg_MuonPResidualFraction->AddEntry(hMuonTrackRangePResidualFraction, "Range", "l");
  lg_MuonPResidualFraction->AddEntry(hMuonTrackMCSPResidualFraction, "MCS", "l");
  lg_MuonPResidualFraction->AddEntry(hMuonTrackCombinedPResidualFraction, "Range (Contained) + MCS (Exiting)", "l");
  lg_MuonPResidualFraction->AddEntry(hMuonTrackCaloPResidualFraction, "Bestplane calo", "l");
  lg_MuonPResidualFraction->AddEntry(hMuonTrackPlane0CaloPResidualFraction, "Ind. 1 calo", "l");
  lg_MuonPResidualFraction->AddEntry(hMuonTrackPlane1CaloPResidualFraction, "Ind. 2 calo", "l");
  lg_MuonPResidualFraction->AddEntry(hMuonTrackPlane2CaloPResidualFraction, "Coll. calo", "l");


  lg_MuonPResidualFraction->Draw();

  latex_POT.DrawLatex(0.20, 0.96, str_TargetPOT);

  hist_dummy_MuonPResidualFraction->GetYaxis()->SetRangeUser(0., y_max*1.1);
  c_MuonPResidualFraction->SaveAs(outputDir+"MuonPResidualFraction.pdf");

  outputfile->cd();
  c_MuonPResidualFraction->Write();

  c_MuonPResidualFraction->Close();

  //==== Proton

  //====   Proton momentum (P)

  y_max = -999.;
  TCanvas *c_ProtonP = new TCanvas("c_ProtonP", "", 800, 800);
  canvas_margin(c_ProtonP);
  c_ProtonP->cd();

  TH1D *hist_dummy_ProtonP = new TH1D("hist_dummy_ProtonP", "", int( (xEnergyMax-xEnergyMin)/dxEnergy ), xEnergyMin, xEnergyMax);
  hist_axis(hist_dummy_ProtonP);
  hist_dummy_ProtonP->GetXaxis()->SetTitle("Proton momentum (GeV)");
  hist_dummy_ProtonP->GetXaxis()->SetRangeUser(0.,5.0);
  hist_dummy_ProtonP->Draw("hist");

  TH1 *hProtonTrackRangeP = sProtonTrackRangeP.ToTH1(TargetPOT);
  hProtonTrackRangeP->SetName("hProtonTrackRangeP");
  hProtonTrackRangeP->SetLineColor(kGreen);
  hProtonTrackRangeP->SetLineWidth(2);
  y_max = max( y_max, hProtonTrackRangeP->GetMaximum() );

  TH1 *hProtonTrackMCSP = sProtonTrackMCSP.ToTH1(TargetPOT);
  hProtonTrackMCSP->SetName("hProtonTrackMCSP");
  hProtonTrackMCSP->SetLineColor(kBlue);
  hProtonTrackMCSP->SetLineWidth(2);
  y_max = max( y_max, hProtonTrackMCSP->GetMaximum() );

  TH1 *hProtonTrackCombinedP = sProtonTrackCombinedP.ToTH1(TargetPOT);
  hProtonTrackCombinedP->SetName("hProtonTrackCombinedP");
  hProtonTrackCombinedP->SetLineColor(kBlack);
  hProtonTrackCombinedP->SetLineWidth(2);
  y_max = max( y_max, hProtonTrackCombinedP->GetMaximum() );

  TH1 *hProtonTrackCaloP = sProtonTrackCaloP.ToTH1(TargetPOT);
  hProtonTrackCaloP->SetName("hProtonTrackCaloP");
  hProtonTrackCaloP->SetLineColor(kViolet);
  hProtonTrackCaloP->SetLineWidth(2);
  y_max = max( y_max, hProtonTrackCaloP->GetMaximum() );

  TH1 *hProtonTrackPlane0CaloP = sProtonTrackPlane0CaloP.ToTH1(TargetPOT);
  hProtonTrackPlane0CaloP->SetName("hProtonTrackPlane0CaloP");
  hProtonTrackPlane0CaloP->SetLineColor(kViolet);
  hProtonTrackPlane0CaloP->SetLineWidth(2);
  hProtonTrackPlane0CaloP->SetLineStyle(3);
  y_max = max( y_max, hProtonTrackPlane0CaloP->GetMaximum() );

  TH1 *hProtonTrackPlane1CaloP = sProtonTrackPlane1CaloP.ToTH1(TargetPOT);
  hProtonTrackPlane1CaloP->SetName("hProtonTrackPlane1CaloP");
  hProtonTrackPlane1CaloP->SetLineColor(kViolet);
  hProtonTrackPlane1CaloP->SetLineWidth(2);
  hProtonTrackPlane1CaloP->SetLineStyle(5);
  y_max = max( y_max, hProtonTrackPlane1CaloP->GetMaximum() );

  TH1 *hProtonTrackPlane2CaloP = sProtonTrackPlane2CaloP.ToTH1(TargetPOT);
  hProtonTrackPlane2CaloP->SetName("hProtonTrackPlane2CaloP");
  hProtonTrackPlane2CaloP->SetLineColor(kViolet);
  hProtonTrackPlane2CaloP->SetLineWidth(2);
  hProtonTrackPlane2CaloP->SetLineStyle(7);
  y_max = max( y_max, hProtonTrackPlane2CaloP->GetMaximum() );

  TH1 *hProtonTrackMatchedTruthP = sProtonTrackMatchedTruthP.ToTH1(TargetPOT);
  hProtonTrackMatchedTruthP->SetName("hProtonTrackMatchedTruthP");
  hProtonTrackMatchedTruthP->SetLineColor(kRed);
  hProtonTrackMatchedTruthP->SetLineWidth(2);
  hProtonTrackMatchedTruthP->SetLineStyle(3);
  y_max = max( y_max, hProtonTrackMatchedTruthP->GetMaximum() );

  hProtonTrackRangeP->Draw("histsame");
  hProtonTrackMCSP->Draw("histsame");
  hProtonTrackCombinedP->Draw("histsame");
  hProtonTrackCaloP->Draw("histsame");
  hProtonTrackPlane0CaloP->Draw("histsame");
  hProtonTrackPlane1CaloP->Draw("histsame");
  hProtonTrackPlane2CaloP->Draw("histsame");
  hProtonTrackMatchedTruthP->Draw("histsame");

  TLegend *lg_ProtonP = new TLegend(0.35, 0.70, 0.90, 0.90);
  lg_ProtonP->SetBorderSize(0);
  lg_ProtonP->SetFillStyle(0);
  lg_ProtonP->AddEntry(hProtonTrackRangeP, "Range", "l");
  lg_ProtonP->AddEntry(hProtonTrackMCSP, "MCS", "l");
  lg_ProtonP->AddEntry(hProtonTrackCombinedP, "Range (Contained) + MCS (Exiting)", "l");
  lg_ProtonP->AddEntry(hProtonTrackCaloP, "Bestplane calo", "l");
  lg_ProtonP->AddEntry(hProtonTrackPlane0CaloP, "Ind. 1 calo", "l");
  lg_ProtonP->AddEntry(hProtonTrackPlane1CaloP, "Ind. 2 calo", "l");
  lg_ProtonP->AddEntry(hProtonTrackPlane2CaloP, "Coll. calo", "l");
  lg_ProtonP->AddEntry(hProtonTrackMatchedTruthP, "Truth", "l");

  lg_ProtonP->Draw();

  latex_POT.DrawLatex(0.20, 0.96, str_TargetPOT);

  hist_dummy_ProtonP->GetYaxis()->SetRangeUser(0., y_max*1.1);
  c_ProtonP->SaveAs(outputDir+"ProtonP.pdf");

  outputfile->cd();
  c_ProtonP->Write();

  c_ProtonP->Close();

  //====   Truth vs Reco 2D : only save TH2, not canvas
  TH2 *hProtonTrack_MatchedTruthP_vs_CaloP = sProtonTrack_MatchedTruthP_vs_CaloP.ToTH2(TargetPOT);
  hProtonTrack_MatchedTruthP_vs_CaloP->SetName("hProtonTrack_MatchedTruthP_vs_CaloP");
  outputfile->cd();
  hProtonTrack_MatchedTruthP_vs_CaloP->Write();

  //====    PDG
  TH1 *hProtonTrackMatchedTruthPDG = sProtonTrackMatchedTruthPDG.ToTH1(TargetPOT);
  hProtonTrackMatchedTruthPDG->SetName("hProtonTrackMatchedTruthPDG");
  outputfile->cd();
  hProtonTrackMatchedTruthPDG->Write();

  //====   Residual

  y_max = -999.;
  TCanvas *c_ProtonPResidual = new TCanvas("c_ProtonPResidual", "", 800, 800);
  canvas_margin(c_ProtonPResidual);
  c_ProtonPResidual->cd();

  TH1D *hist_dummy_ProtonPResidual = new TH1D("hist_dummy_ProtonPResidual", "", 600, -3., 3.);
  hist_axis(hist_dummy_ProtonPResidual);
  hist_dummy_ProtonPResidual->GetXaxis()->SetTitle("Proton momentum (Reco-Truth) (GeV)");
  hist_dummy_ProtonPResidual->Draw("hist");

  TH1 *hProtonTrackRangePResidual = sProtonTrackRangePResidual.ToTH1(TargetPOT);
  hProtonTrackRangePResidual->SetName("hProtonTrackRangePResidual");
  hProtonTrackRangePResidual->SetLineColor(kGreen);
  hProtonTrackRangePResidual->SetLineWidth(2);
  y_max = max( y_max, hProtonTrackRangePResidual->GetMaximum() );

  TH1 *hProtonTrackMCSPResidual = sProtonTrackMCSPResidual.ToTH1(TargetPOT);
  hProtonTrackMCSPResidual->SetName("hProtonTrackMCSPResidual");
  hProtonTrackMCSPResidual->SetLineColor(kBlue);
  hProtonTrackMCSPResidual->SetLineWidth(2);
  y_max = max( y_max, hProtonTrackMCSPResidual->GetMaximum() );

  TH1 *hProtonTrackCombinedPResidual = sProtonTrackCombinedPResidual.ToTH1(TargetPOT);
  hProtonTrackCombinedPResidual->SetName("hProtonTrackCombinedPResidual");
  hProtonTrackCombinedPResidual->SetLineColor(kBlack);
  hProtonTrackCombinedPResidual->SetLineWidth(2);
  y_max = max( y_max, hProtonTrackCombinedPResidual->GetMaximum() );

  TH1 *hProtonTrackCaloPResidual = sProtonTrackCaloPResidual.ToTH1(TargetPOT);
  hProtonTrackCaloPResidual->SetName("hProtonTrackCaloPResidual");
  hProtonTrackCaloPResidual->SetLineColor(kViolet);
  hProtonTrackCaloPResidual->SetLineWidth(2);
  y_max = max( y_max, hProtonTrackCaloPResidual->GetMaximum() );

  TH1 *hProtonTrackPlane0CaloPResidual = sProtonTrackPlane0CaloPResidual.ToTH1(TargetPOT);
  hProtonTrackPlane0CaloPResidual->SetName("hProtonTrackPlane0CaloPResidual");
  hProtonTrackPlane0CaloPResidual->SetLineColor(kViolet);
  hProtonTrackPlane0CaloPResidual->SetLineWidth(2);
  hProtonTrackPlane0CaloPResidual->SetLineStyle(3);
  y_max = max( y_max, hProtonTrackPlane0CaloPResidual->GetMaximum() );

  TH1 *hProtonTrackPlane1CaloPResidual = sProtonTrackPlane1CaloPResidual.ToTH1(TargetPOT);
  hProtonTrackPlane1CaloPResidual->SetName("hProtonTrackPlane1CaloPResidual");
  hProtonTrackPlane1CaloPResidual->SetLineColor(kViolet);
  hProtonTrackPlane1CaloPResidual->SetLineWidth(2);
  hProtonTrackPlane1CaloPResidual->SetLineStyle(5);
  y_max = max( y_max, hProtonTrackPlane1CaloPResidual->GetMaximum() );

  TH1 *hProtonTrackPlane2CaloPResidual = sProtonTrackPlane2CaloPResidual.ToTH1(TargetPOT);
  hProtonTrackPlane2CaloPResidual->SetName("hProtonTrackPlane2CaloPResidual");
  hProtonTrackPlane2CaloPResidual->SetLineColor(kViolet);
  hProtonTrackPlane2CaloPResidual->SetLineWidth(2);
  hProtonTrackPlane2CaloPResidual->SetLineStyle(7);
  y_max = max( y_max, hProtonTrackPlane2CaloPResidual->GetMaximum() );

  hProtonTrackRangePResidual->Draw("histsame");
  hProtonTrackMCSPResidual->Draw("histsame");
  hProtonTrackCombinedPResidual->Draw("histsame");
  hProtonTrackCaloPResidual->Draw("histsame");
  hProtonTrackPlane0CaloPResidual->Draw("histsame");
  hProtonTrackPlane1CaloPResidual->Draw("histsame");
  hProtonTrackPlane2CaloPResidual->Draw("histsame");

  TLegend *lg_ProtonPResidual = new TLegend(0.35, 0.70, 0.90, 0.90);
  lg_ProtonPResidual->SetBorderSize(0);
  lg_ProtonPResidual->SetFillStyle(0);
  lg_ProtonPResidual->AddEntry(hProtonTrackRangePResidual, "Range", "l");
  lg_ProtonPResidual->AddEntry(hProtonTrackMCSPResidual, "MCS", "l");
  lg_ProtonPResidual->AddEntry(hProtonTrackCombinedPResidual, "Range (Contained) + MCS (Exiting)", "l");
  lg_ProtonPResidual->AddEntry(hProtonTrackCaloPResidual, "Bestplane calo", "l");
  lg_ProtonPResidual->AddEntry(hProtonTrackPlane0CaloPResidual, "Ind. 1 calo", "l");
  lg_ProtonPResidual->AddEntry(hProtonTrackPlane1CaloPResidual, "Ind. 2 calo", "l");
  lg_ProtonPResidual->AddEntry(hProtonTrackPlane2CaloPResidual, "Coll. calo", "l");

  lg_ProtonPResidual->Draw();

  latex_POT.DrawLatex(0.20, 0.96, str_TargetPOT);

  hist_dummy_ProtonPResidual->GetYaxis()->SetRangeUser(0., y_max*1.1);
  c_ProtonPResidual->SaveAs(outputDir+"ProtonPResidual.pdf");

  outputfile->cd();
  c_ProtonPResidual->Write();

  c_ProtonPResidual->Close();

  //====   ResidualFraction

  y_max = -999.;
  TCanvas *c_ProtonPResidualFraction = new TCanvas("c_ProtonPResidualFraction", "", 800, 800);
  canvas_margin(c_ProtonPResidualFraction);
  c_ProtonPResidualFraction->cd();

  TH1D *hist_dummy_ProtonPResidualFraction = new TH1D("hist_dummy_ProtonPResidualFraction", "", 600, -3., 3.);
  hist_axis(hist_dummy_ProtonPResidualFraction);
  hist_dummy_ProtonPResidualFraction->GetXaxis()->SetTitle("Proton momentum (Reco-Truth)/Truth");
  hist_dummy_ProtonPResidualFraction->Draw("hist");

  TH1 *hProtonTrackRangePResidualFraction = sProtonTrackRangePResidualFraction.ToTH1(TargetPOT);
  hProtonTrackRangePResidualFraction->SetName("hProtonTrackRangePResidualFraction");
  hProtonTrackRangePResidualFraction->SetLineColor(kGreen);
  hProtonTrackRangePResidualFraction->SetLineWidth(2);
  y_max = max( y_max, hProtonTrackRangePResidualFraction->GetMaximum() );

  TH1 *hProtonTrackMCSPResidualFraction = sProtonTrackMCSPResidualFraction.ToTH1(TargetPOT);
  hProtonTrackMCSPResidualFraction->SetName("hProtonTrackMCSPResidualFraction");
  hProtonTrackMCSPResidualFraction->SetLineColor(kBlue);
  hProtonTrackMCSPResidualFraction->SetLineWidth(2);
  y_max = max( y_max, hProtonTrackMCSPResidualFraction->GetMaximum() );

  TH1 *hProtonTrackCombinedPResidualFraction = sProtonTrackCombinedPResidualFraction.ToTH1(TargetPOT);
  hProtonTrackCombinedPResidualFraction->SetName("hProtonTrackCombinedPResidualFraction");
  hProtonTrackCombinedPResidualFraction->SetLineColor(kBlack);
  hProtonTrackCombinedPResidualFraction->SetLineWidth(2);
  y_max = max( y_max, hProtonTrackCombinedPResidualFraction->GetMaximum() );

  TH1 *hProtonTrackCaloPResidualFraction = sProtonTrackCaloPResidualFraction.ToTH1(TargetPOT);
  hProtonTrackCaloPResidualFraction->SetName("hProtonTrackCaloPResidualFraction");
  hProtonTrackCaloPResidualFraction->SetLineColor(kViolet);
  hProtonTrackCaloPResidualFraction->SetLineWidth(2);
  y_max = max( y_max, hProtonTrackCaloPResidualFraction->GetMaximum() );

  TH1 *hProtonTrackPlane0CaloPResidualFraction = sProtonTrackPlane0CaloPResidualFraction.ToTH1(TargetPOT);
  hProtonTrackPlane0CaloPResidualFraction->SetName("hProtonTrackPlane0CaloPResidualFraction");
  hProtonTrackPlane0CaloPResidualFraction->SetLineColor(kViolet);
  hProtonTrackPlane0CaloPResidualFraction->SetLineWidth(2);
  hProtonTrackPlane0CaloPResidualFraction->SetLineStyle(3);
  y_max = max( y_max, hProtonTrackPlane0CaloPResidualFraction->GetMaximum() );

  TH1 *hProtonTrackPlane1CaloPResidualFraction = sProtonTrackPlane1CaloPResidualFraction.ToTH1(TargetPOT);
  hProtonTrackPlane1CaloPResidualFraction->SetName("hProtonTrackPlane1CaloPResidualFraction");
  hProtonTrackPlane1CaloPResidualFraction->SetLineColor(kViolet);
  hProtonTrackPlane1CaloPResidualFraction->SetLineWidth(2);
  hProtonTrackPlane1CaloPResidualFraction->SetLineStyle(5);
  y_max = max( y_max, hProtonTrackPlane1CaloPResidualFraction->GetMaximum() );

  TH1 *hProtonTrackPlane2CaloPResidualFraction = sProtonTrackPlane2CaloPResidualFraction.ToTH1(TargetPOT);
  hProtonTrackPlane2CaloPResidualFraction->SetName("hProtonTrackPlane2CaloPResidualFraction");
  hProtonTrackPlane2CaloPResidualFraction->SetLineColor(kViolet);
  hProtonTrackPlane2CaloPResidualFraction->SetLineWidth(2);
  hProtonTrackPlane2CaloPResidualFraction->SetLineStyle(7);
  y_max = max( y_max, hProtonTrackPlane2CaloPResidualFraction->GetMaximum() );

  hProtonTrackRangePResidualFraction->Draw("histsame");
  hProtonTrackMCSPResidualFraction->Draw("histsame");
  hProtonTrackCombinedPResidualFraction->Draw("histsame");
  hProtonTrackCaloPResidualFraction->Draw("histsame");
  hProtonTrackPlane0CaloPResidualFraction->Draw("histsame");
  hProtonTrackPlane1CaloPResidualFraction->Draw("histsame");
  hProtonTrackPlane2CaloPResidualFraction->Draw("histsame");

  TLegend *lg_ProtonPResidualFraction = new TLegend(0.35, 0.70, 0.90, 0.90);
  lg_ProtonPResidualFraction->SetBorderSize(0);
  lg_ProtonPResidualFraction->SetFillStyle(0);
  lg_ProtonPResidualFraction->AddEntry(hProtonTrackRangePResidualFraction, "Range", "l");
  lg_ProtonPResidualFraction->AddEntry(hProtonTrackMCSPResidualFraction, "MCS", "l");
  lg_ProtonPResidualFraction->AddEntry(hProtonTrackCombinedPResidualFraction, "Range (Contained) + MCS (Exiting)", "l");
  lg_ProtonPResidualFraction->AddEntry(hProtonTrackCaloPResidualFraction, "Bestplane calo", "l");
  lg_ProtonPResidualFraction->AddEntry(hProtonTrackPlane0CaloPResidualFraction, "Ind. 1 calo", "l");
  lg_ProtonPResidualFraction->AddEntry(hProtonTrackPlane1CaloPResidualFraction, "Ind. 2 calo", "l");
  lg_ProtonPResidualFraction->AddEntry(hProtonTrackPlane2CaloPResidualFraction, "Coll. calo", "l");

  lg_ProtonPResidualFraction->Draw();

  latex_POT.DrawLatex(0.20, 0.96, str_TargetPOT);

  hist_dummy_ProtonPResidualFraction->GetYaxis()->SetRangeUser(0., y_max*1.1);
  c_ProtonPResidualFraction->SaveAs(outputDir+"ProtonPResidualFraction.pdf");

  outputfile->cd();
  c_ProtonPResidualFraction->Write();

  c_ProtonPResidualFraction->Close();


  //==== Neutrnio

  //====   Neutrino energy

  y_max = -999.;
  TCanvas *c_NeutrinoE = new TCanvas("c_NeutrinoE", "", 800, 800);
  canvas_margin(c_NeutrinoE);
  c_NeutrinoE->cd();

  TH1D *hist_dummy_NeutrinoE = new TH1D("hist_dummy_NeutrinoE", "", int( (xEnergyMax-xEnergyMin)/dxEnergy ), xEnergyMin, xEnergyMax);
  hist_axis(hist_dummy_NeutrinoE);
  hist_dummy_NeutrinoE->GetXaxis()->SetTitle("Neutrino energy (GeV)");
  hist_dummy_NeutrinoE->GetXaxis()->SetRangeUser(0.,5.0);
  hist_dummy_NeutrinoE->Draw("hist");

  TH1 *hNeutrinoCaloEnergy = sNeutrinoCaloEnergy.ToTH1(TargetPOT);
  hNeutrinoCaloEnergy->SetName("hNeutrinoCaloEnergy");
  hNeutrinoCaloEnergy->SetLineColor(kGreen);
  hNeutrinoCaloEnergy->SetLineWidth(2);
  y_max = max( y_max, hNeutrinoCaloEnergy->GetMaximum() );

  TH1 *hNeutrinoCombinedEnergy = sNeutrinoCombinedEnergy.ToTH1(TargetPOT);
  hNeutrinoCombinedEnergy->SetName("hNeutrinoCombinedEnergy");
  hNeutrinoCombinedEnergy->SetLineColor(kBlue);
  hNeutrinoCombinedEnergy->SetLineWidth(2);
  y_max = max( y_max, hNeutrinoCombinedEnergy->GetMaximum() );

  TH1 *hNeutrinoQE = sNeutrinoQE.ToTH1(TargetPOT);
  hNeutrinoQE->SetName("hNeutrinoQE");
  hNeutrinoQE->SetLineColor(kBlack);
  hNeutrinoQE->SetLineWidth(2);
  y_max = max( y_max, hNeutrinoQE->GetMaximum() );

  TH1 *hNeutrinoFakeRecoEnergy = sNeutrinoFakeRecoEnergy.ToTH1(TargetPOT);
  hNeutrinoFakeRecoEnergy->SetName("hNeutrinoFakeRecoEnergy");
  hNeutrinoFakeRecoEnergy->SetLineColor(kRed);
  hNeutrinoFakeRecoEnergy->SetLineWidth(2);
  y_max = max( y_max, hNeutrinoFakeRecoEnergy->GetMaximum() );

  TH1 *hNeutrinoTruthE = sNeutrinoTruthE.ToTH1(TargetPOT);
  hNeutrinoTruthE->SetName("hNeutrinoTruthE");
  hNeutrinoTruthE->SetLineColor(kRed);
  hNeutrinoTruthE->SetLineWidth(2);
  hNeutrinoTruthE->SetLineStyle(3);
  y_max = max( y_max, hNeutrinoTruthE->GetMaximum() );

  hNeutrinoCaloEnergy->Draw("histsame");
  hNeutrinoCombinedEnergy->Draw("histsame");
  hNeutrinoQE->Draw("histsame");
  hNeutrinoFakeRecoEnergy->Draw("histsame");
  hNeutrinoTruthE->Draw("histsame");

  TLegend *lg_NeutrinoE = new TLegend(0.35, 0.70, 0.90, 0.90);
  lg_NeutrinoE->SetBorderSize(0);
  lg_NeutrinoE->SetFillStyle(0);
  lg_NeutrinoE->AddEntry(hNeutrinoCaloEnergy, "Muon (calo) + Proton (calo)", "l");
  lg_NeutrinoE->AddEntry(hNeutrinoCombinedEnergy, "Muon (combined) + Proton (calo)", "l");
  lg_NeutrinoE->AddEntry(hNeutrinoQE, "E_{QE}", "l");
  lg_NeutrinoE->AddEntry(hNeutrinoFakeRecoEnergy, "Muon (Truth) + Proton (Truth)", "l");
  lg_NeutrinoE->AddEntry(hNeutrinoTruthE, "Truth", "l");

  lg_NeutrinoE->Draw();

  latex_POT.DrawLatex(0.20, 0.96, str_TargetPOT);

  hist_dummy_NeutrinoE->GetYaxis()->SetRangeUser(0., y_max*1.1);
  c_NeutrinoE->SaveAs(outputDir+"NeutrinoE.pdf");

  outputfile->cd();
  c_NeutrinoE->Write();

  c_NeutrinoE->Close();

  //====   Truth vs Reco 2D : only save TH2, not canvas
  TH2 *hNeutrino_TruthE_vs_CombinedEnergy = sNeutrino_TruthE_vs_CombinedEnergy.ToTH2(TargetPOT);
  hNeutrino_TruthE_vs_CombinedEnergy->SetName("hNeutrino_TruthE_vs_CombinedEnergy");
  outputfile->cd();
  hNeutrino_TruthE_vs_CombinedEnergy->Write();
  //====   Truth vs Reco 2D : only save TH2, not canvas
  TH2 *hNeutrino_TruthE_vs_QE = sNeutrino_TruthE_vs_QE.ToTH2(TargetPOT);
  hNeutrino_TruthE_vs_QE->SetName("hNeutrino_TruthE_vs_QE");
  outputfile->cd();
  hNeutrino_TruthE_vs_QE->Write();

  //====   Residual

  y_max = -999.;
  TCanvas *c_NeutrinoEResidual = new TCanvas("c_NeutrinoEResidual", "", 800, 800);
  canvas_margin(c_NeutrinoEResidual);
  c_NeutrinoEResidual->cd();

  TH1D *hist_dummy_NeutrinoEResidual = new TH1D("hist_dummy_NeutrinoEResidual", "", 600, -3., 3.);
  hist_axis(hist_dummy_NeutrinoEResidual);
  hist_dummy_NeutrinoEResidual->GetXaxis()->SetTitle("Neutrino momentum (Reco-Truth) (GeV)");
  hist_dummy_NeutrinoEResidual->Draw("hist");

  TH1 *hNeutrinoCaloEnergyResidual = sNeutrinoCaloEnergyResidual.ToTH1(TargetPOT);
  hNeutrinoCaloEnergyResidual->SetName("hNeutrinoCaloEnergyResidual");
  hNeutrinoCaloEnergyResidual->SetLineColor(kGreen);
  hNeutrinoCaloEnergyResidual->SetLineWidth(2);
  y_max = max( y_max, hNeutrinoCaloEnergyResidual->GetMaximum() );

  TH1 *hNeutrinoCombinedEnergyResidual = sNeutrinoCombinedEnergyResidual.ToTH1(TargetPOT);
  hNeutrinoCombinedEnergyResidual->SetName("hNeutrinoCombinedEnergyResidual");
  hNeutrinoCombinedEnergyResidual->SetLineColor(kBlue);
  hNeutrinoCombinedEnergyResidual->SetLineWidth(2);
  y_max = max( y_max, hNeutrinoCombinedEnergyResidual->GetMaximum() );

  TH1 *hNeutrinoQEResidual = sNeutrinoQEResidual.ToTH1(TargetPOT);
  hNeutrinoQEResidual->SetName("hNeutrinoQEResidual");
  hNeutrinoQEResidual->SetLineColor(kBlack);
  hNeutrinoQEResidual->SetLineWidth(2);
  y_max = max( y_max, hNeutrinoQEResidual->GetMaximum() );

  hNeutrinoCaloEnergyResidual->Draw("histsame");
  hNeutrinoCombinedEnergyResidual->Draw("histsame");
  hNeutrinoQEResidual->Draw("histsame");

  TLegend *lg_NeutrinoEResidual = new TLegend(0.35, 0.70, 0.90, 0.90);
  lg_NeutrinoEResidual->SetBorderSize(0);
  lg_NeutrinoEResidual->SetFillStyle(0);
  lg_NeutrinoEResidual->AddEntry(hNeutrinoCaloEnergyResidual, "Muon (calo) + Proton (calo)", "l");
  lg_NeutrinoEResidual->AddEntry(hNeutrinoCombinedEnergyResidual, "Muon (combined) + Proton (calo)", "l");
  lg_NeutrinoEResidual->AddEntry(hNeutrinoQEResidual, "E_{QE}", "l");

  lg_NeutrinoEResidual->Draw();

  latex_POT.DrawLatex(0.20, 0.96, str_TargetPOT);

  hist_dummy_NeutrinoEResidual->GetYaxis()->SetRangeUser(0., y_max*1.1);
  c_NeutrinoEResidual->SaveAs(outputDir+"NeutrinoEResidual.pdf");

  outputfile->cd();
  c_NeutrinoEResidual->Write();

  c_NeutrinoEResidual->Close();

  //====   ResidualFraction

  y_max = -999.;
  TCanvas *c_NeutrinoEResidualFraction = new TCanvas("c_NeutrinoEResidualFraction", "", 800, 800);
  canvas_margin(c_NeutrinoEResidualFraction);
  c_NeutrinoEResidualFraction->cd();

  TH1D *hist_dummy_NeutrinoEResidualFraction = new TH1D("hist_dummy_NeutrinoEResidualFraction", "", 600, -3., 3.);
  hist_axis(hist_dummy_NeutrinoEResidualFraction);
  hist_dummy_NeutrinoEResidualFraction->GetXaxis()->SetTitle("Neutrino momentum (Reco-Truth)/Truth (GeV)");
  hist_dummy_NeutrinoEResidualFraction->Draw("hist");

  TH1 *hNeutrinoCaloEnergyResidualFraction = sNeutrinoCaloEnergyResidualFraction.ToTH1(TargetPOT);
  hNeutrinoCaloEnergyResidualFraction->SetName("hNeutrinoCaloEnergyResidualFraction");
  hNeutrinoCaloEnergyResidualFraction->SetLineColor(kGreen);
  hNeutrinoCaloEnergyResidualFraction->SetLineWidth(2);
  y_max = max( y_max, hNeutrinoCaloEnergyResidualFraction->GetMaximum() );

  TH1 *hNeutrinoCombinedEnergyResidualFraction = sNeutrinoCombinedEnergyResidualFraction.ToTH1(TargetPOT);
  hNeutrinoCombinedEnergyResidualFraction->SetName("hNeutrinoCombinedEnergyResidualFraction");
  hNeutrinoCombinedEnergyResidualFraction->SetLineColor(kBlue);
  hNeutrinoCombinedEnergyResidualFraction->SetLineWidth(2);
  y_max = max( y_max, hNeutrinoCombinedEnergyResidualFraction->GetMaximum() );

  TH1 *hNeutrinoQEResidualFraction = sNeutrinoQEResidualFraction.ToTH1(TargetPOT);
  hNeutrinoQEResidualFraction->SetName("hNeutrinoQEResidualFraction");
  hNeutrinoQEResidualFraction->SetLineColor(kBlack);
  hNeutrinoQEResidualFraction->SetLineWidth(2);
  y_max = max( y_max, hNeutrinoQEResidualFraction->GetMaximum() );

  hNeutrinoCaloEnergyResidualFraction->Draw("histsame");
  hNeutrinoCombinedEnergyResidualFraction->Draw("histsame");
  hNeutrinoQEResidualFraction->Draw("histsame");

  TLegend *lg_NeutrinoEResidualFraction = new TLegend(0.35, 0.70, 0.90, 0.90);
  lg_NeutrinoEResidualFraction->SetBorderSize(0);
  lg_NeutrinoEResidualFraction->SetFillStyle(0);
  lg_NeutrinoEResidualFraction->AddEntry(hNeutrinoCaloEnergyResidualFraction, "Muon (calo) + Proton (calo)", "l");
  lg_NeutrinoEResidualFraction->AddEntry(hNeutrinoCombinedEnergyResidualFraction, "Muon (combined) + Proton (calo)", "l");
  lg_NeutrinoEResidualFraction->AddEntry(hNeutrinoQEResidualFraction, "E_{QE}", "l");

  lg_NeutrinoEResidualFraction->Draw();

  latex_POT.DrawLatex(0.20, 0.96, str_TargetPOT);

  hist_dummy_NeutrinoEResidualFraction->GetYaxis()->SetRangeUser(0., y_max*1.1);
  c_NeutrinoEResidualFraction->SaveAs(outputDir+"NeutrinoEResidualFraction.pdf");

  outputfile->cd();
  c_NeutrinoEResidualFraction->Write();

  c_NeutrinoEResidualFraction->Close();


  //==== close file

  outputfile->Close();

}

