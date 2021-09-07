#include "myEfficiency.h"

myEfficiency::myEfficiency(){

  cout << "[myEfficiency::myEfficiency] called" << endl;
  TargetPOT = 6.6e20;
  str_TargetPOT = "6.6e20 POT";
  outputName = "output.root";
  vec_Spectrums.clear();
  cout << "[myEfficiency::myEfficiency] Finished" << endl;

}

void myEfficiency::initialize(){

  cout << "[myEfficiency::initialize] outputDir = " << outputDir << endl;
  gSystem->mkdir(outputDir, kTRUE);
  TString outputName = "output.root";
  outputfile = new TFile(outputDir+outputName,"RECREATE");

}

void myEfficiency::bookTruth(SpectrumLoader& loader, Cut cut){

  cout << "[myEfficiency::bookTruth] called" << endl;

  //==== some binnings
  double xEnergyMin = 0.; // GeV
  double xEnergyMax = 10.;
  double dxEnergy = 0.1;
  const Binning binsEnergy = Binning::Simple( int( (xEnergyMax-xEnergyMin)/dxEnergy ), xEnergyMin, xEnergyMax );
  double xChi2MuonMin = 0.;
  double xChi2MuonMax = 200.;
  double dxChi2Muon = 1.;
  const Binning binsChi2Muon = Binning::Simple( int( (xChi2MuonMax-xChi2MuonMin)/dxChi2Muon ), xChi2MuonMin, xChi2MuonMax );
  double xChi2ProtonMin = 0.;
  double xChi2ProtonMax = 600.;
  double dxChi2Proton = 1.;
  const Binning binsChi2Proton = Binning::Simple( int( (xChi2ProtonMax-xChi2ProtonMin)/dxChi2Proton ), xChi2ProtonMin, xChi2ProtonMax );
  const Binning binsNormalizedChi2 = Binning::Simple( 100, 0., 10. );

  cout << "[myEfficiency::bookTruth] Delcaring HistAxis" << endl;

  //==== HistAxis
  const HistAxis axMuonTrackMatchedTruthP("MuonTrackMatchedTruthP", binsEnergy, varMuonTrackMatchedTruthP);
  const HistAxis axProtonTrackMatchedTruthP("ProtonTrackMatchedTruthP", binsEnergy, varProtonTrackMatchedTruthP);
  const HistAxis axNeutrinoTruthE("NeutrinoTruthE", binsEnergy, varNeutrinoTruthE);
  const HistAxis axMuonTruthP("MuonTruthP", binsEnergy, varMuonTruthP);
  const HistAxis axProtonTruthP("ProtonTruthP", binsEnergy, varProtonTruthP);
  //=== matching..
  const HistAxis axTruthMuonMatchedTrackChi2Proton("TruthMuonMatchedTrackChi2Proton", binsChi2Proton, varTruthMuonMatchedTrackChi2Proton);
  const HistAxis axTruthMuonMatchedTrackChi2Muon("TruthMuonMatchedTrackChi2Muon", binsChi2Muon, varTruthMuonMatchedTrackChi2Muon);
  const HistAxis axTruthProtonMatchedTrackChi2Proton("TruthProtonMatchedTrackChi2Proton", binsChi2Proton, varTruthProtonMatchedTrackChi2Proton);
  const HistAxis axTruthProtonMatchedTrackChi2Muon("TruthProtonMatchedTrackChi2Muon", binsChi2Muon, varTruthProtonMatchedTrackChi2Muon);
  const HistAxis axTruthMuonMatchedTrackNormalizedChi2Proton("TruthMuonMatchedTrackNormalizedChi2Proton", binsNormalizedChi2, varTruthMuonMatchedTrackNormalizedChi2Proton);
  const HistAxis axTruthMuonMatchedTrackNormalizedChi2Muon("TruthMuonMatchedTrackNormalizedChi2Muon", binsNormalizedChi2, varTruthMuonMatchedTrackNormalizedChi2Muon);
  const HistAxis axTruthProtonMatchedTrackNormalizedChi2Proton("TruthProtonMatchedTrackNormalizedChi2Proton", binsNormalizedChi2, varTruthProtonMatchedTrackNormalizedChi2Proton);
  const HistAxis axTruthProtonMatchedTrackNormalizedChi2Muon("TruthProtonMatchedTrackNormalizedChi2Muon", binsNormalizedChi2, varTruthProtonMatchedTrackNormalizedChi2Muon);

  cout << "[myEfficiency::bookTruth] Delcaring Spectrum" << endl;

  //==== Spectrum
  Spectrum *sMuonTrackMatchedTruthP = new Spectrum(loader, axMuonTrackMatchedTruthP, kNoSpillCut, cut);
  Spectrum *sProtonTrackMatchedTruthP = new Spectrum(loader, axProtonTrackMatchedTruthP, kNoSpillCut, cut);
  Spectrum *sNeutrinoTruthE = new Spectrum(loader, axNeutrinoTruthE, kNoSpillCut, cut);
  Spectrum *sMuonTruthP = new Spectrum(loader, axMuonTruthP, kNoSpillCut, cut);
  Spectrum *sProtonTruthP = new Spectrum(loader, axProtonTruthP, kNoSpillCut, cut);
  Spectrum *sTruthMuonMatchedTrackChi2Proton = new Spectrum(loader, axTruthMuonMatchedTrackChi2Proton, kNoSpillCut, cut);
  Spectrum *sTruthMuonMatchedTrackChi2Muon = new Spectrum(loader, axTruthMuonMatchedTrackChi2Muon, kNoSpillCut, cut);
  Spectrum *sTruthMuonMatchedTrackChi2_Muon_vs_Proton = new Spectrum(loader, axTruthMuonMatchedTrackChi2Muon, axTruthMuonMatchedTrackChi2Proton, kNoSpillCut, cut);
  Spectrum *sTruthProtonMatchedTrackChi2Proton = new Spectrum(loader, axTruthProtonMatchedTrackChi2Proton, kNoSpillCut, cut);
  Spectrum *sTruthProtonMatchedTrackChi2Muon = new Spectrum(loader, axTruthProtonMatchedTrackChi2Muon, kNoSpillCut, cut);
  Spectrum *sTruthProtonMatchedTrackChi2_Muon_vs_Proton = new Spectrum(loader, axTruthProtonMatchedTrackChi2Muon, axTruthProtonMatchedTrackChi2Proton, kNoSpillCut, cut);
  Spectrum *sTruthMuonMatchedTrackNormalizedChi2Proton = new Spectrum(loader, axTruthMuonMatchedTrackNormalizedChi2Proton, kNoSpillCut, cut);
  Spectrum *sTruthMuonMatchedTrackNormalizedChi2Muon = new Spectrum(loader, axTruthMuonMatchedTrackNormalizedChi2Muon, kNoSpillCut, cut);
  Spectrum *sTruthMuonMatchedTrackNormalizedChi2_Muon_vs_Proton = new Spectrum(loader, axTruthMuonMatchedTrackNormalizedChi2Muon, axTruthMuonMatchedTrackNormalizedChi2Proton, kNoSpillCut, cut);
  Spectrum *sTruthProtonMatchedTrackNormalizedChi2Proton = new Spectrum(loader, axTruthProtonMatchedTrackNormalizedChi2Proton, kNoSpillCut, cut);
  Spectrum *sTruthProtonMatchedTrackNormalizedChi2Muon = new Spectrum(loader, axTruthProtonMatchedTrackNormalizedChi2Muon, kNoSpillCut, cut);
  Spectrum *sTruthProtonMatchedTrackNormalizedChi2_Muon_vs_Proton = new Spectrum(loader, axTruthProtonMatchedTrackNormalizedChi2Muon, axTruthProtonMatchedTrackNormalizedChi2Proton, kNoSpillCut, cut);

  cout << "[myEfficiency::bookTruth] Saving spectrums" << endl;

  vec_Spectrums.push_back(sMuonTrackMatchedTruthP);
  vec_Spectrums.push_back(sProtonTrackMatchedTruthP);
  vec_Spectrums.push_back(sNeutrinoTruthE);
  vec_Spectrums.push_back(sMuonTruthP);
  vec_Spectrums.push_back(sProtonTruthP);
  vec_Spectrums.push_back(sTruthMuonMatchedTrackChi2Proton);
  vec_Spectrums.push_back(sTruthMuonMatchedTrackChi2Muon);
  vec_Spectrums.push_back(sTruthMuonMatchedTrackChi2_Muon_vs_Proton);
  vec_Spectrums.push_back(sTruthProtonMatchedTrackChi2Proton);
  vec_Spectrums.push_back(sTruthProtonMatchedTrackChi2Muon);
  vec_Spectrums.push_back(sTruthProtonMatchedTrackChi2_Muon_vs_Proton);
  vec_Spectrums.push_back(sTruthMuonMatchedTrackNormalizedChi2Proton);
  vec_Spectrums.push_back(sTruthMuonMatchedTrackNormalizedChi2Muon);
  vec_Spectrums.push_back(sTruthMuonMatchedTrackNormalizedChi2_Muon_vs_Proton);
  vec_Spectrums.push_back(sTruthProtonMatchedTrackNormalizedChi2Proton);
  vec_Spectrums.push_back(sTruthProtonMatchedTrackNormalizedChi2Muon);
  vec_Spectrums.push_back(sTruthProtonMatchedTrackNormalizedChi2_Muon_vs_Proton);

  cout << "[myEfficiency::bookTruth] Finished" << endl;

}

void myEfficiency::bookRecoMuon(SpectrumLoader& loader, Cut cut){

  cout << "[myEfficiency::bookRecoMuon] called" << endl;

  //==== some binnings
  double xEnergyMin = 0.; // GeV
  double xEnergyMax = 10.;
  double dxEnergy = 0.1;
  const Binning binsEnergy = Binning::Simple( int( (xEnergyMax-xEnergyMin)/dxEnergy ), xEnergyMin, xEnergyMax );

  cout << "[myEfficiency::bookRecoMuon] Delcaring HistAxis" << endl;

  //==== HistAxis
  const HistAxis axMuonTrackRangeP("MuonTrackRangeP", binsEnergy, varMuonTrackRangeP);
  const HistAxis axMuonTrackMCSP("MuonTrackMCSP", binsEnergy, varMuonTrackMCSP);
  const HistAxis axMuonTrackCombinedP("MuonTrackCombinedP", binsEnergy, varMuonTrackCombinedP);
  const HistAxis axMuonTrackCaloP("MuonTrackCaloP", binsEnergy, varMuonTrackCaloP);

  cout << "[myEfficiency::bookRecoMuon] Delcaring Spectrum" << endl;

  //==== Spectrum
  Spectrum *sMuonTrackRangeP = new Spectrum(loader, axMuonTrackRangeP, kNoSpillCut, cut);
  Spectrum *sMuonTrackMCSP = new Spectrum(loader, axMuonTrackMCSP, kNoSpillCut, cut);
  Spectrum *sMuonTrackCombinedP = new Spectrum(loader, axMuonTrackCombinedP, kNoSpillCut, cut);
  Spectrum *sMuonTrackCaloP = new Spectrum(loader, axMuonTrackCaloP, kNoSpillCut, cut);

  cout << "[myEfficiency::bookRecoMuon] Saving spectrums" << endl;

  vec_Spectrums.push_back(sMuonTrackRangeP);
  vec_Spectrums.push_back(sMuonTrackMCSP);
  vec_Spectrums.push_back(sMuonTrackCombinedP);
  vec_Spectrums.push_back(sMuonTrackCaloP);

  cout << "[myEfficiency::bookRecoMuon] Finished" << endl;

}

void myEfficiency::bookRecoProton(SpectrumLoader& loader, Cut cut){

  cout << "[myEfficiency::bookRecoProton] called" << endl;

  //==== some binnings
  double xEnergyMin = 0.; // GeV
  double xEnergyMax = 10.;
  double dxEnergy = 0.1;
  const Binning binsEnergy = Binning::Simple( int( (xEnergyMax-xEnergyMin)/dxEnergy ), xEnergyMin, xEnergyMax );

  cout << "[myEfficiency::bookRecoProton] Delcaring HistAxis" << endl;

  //==== HistAxis
  const HistAxis axProtonTrackRangeP("ProtonTrackRangeP", binsEnergy, varProtonTrackRangeP);
  const HistAxis axProtonTrackMCSP("ProtonTrackMCSP", binsEnergy, varProtonTrackMCSP);
  const HistAxis axProtonTrackCombinedP("ProtonTrackCombinedP", binsEnergy, varProtonTrackCombinedP);
  const HistAxis axProtonTrackCaloP("ProtonTrackCaloP", binsEnergy, varProtonTrackCaloP);

  cout << "[myEfficiency::bookRecoProton] Delcaring Spectrum" << endl;

  //==== Spectrum
  Spectrum *sProtonTrackRangeP = new Spectrum(loader, axProtonTrackRangeP, kNoSpillCut, cut);
  Spectrum *sProtonTrackMCSP = new Spectrum(loader, axProtonTrackMCSP, kNoSpillCut, cut);
  Spectrum *sProtonTrackCombinedP = new Spectrum(loader, axProtonTrackCombinedP, kNoSpillCut, cut);
  Spectrum *sProtonTrackCaloP = new Spectrum(loader, axProtonTrackCaloP, kNoSpillCut, cut);

  cout << "[myEfficiency::bookRecoProton] Saving spectrums" << endl;

  vec_Spectrums.push_back(sProtonTrackRangeP);
  vec_Spectrums.push_back(sProtonTrackMCSP);
  vec_Spectrums.push_back(sProtonTrackCombinedP);
  vec_Spectrums.push_back(sProtonTrackCaloP);

  cout << "[myEfficiency::bookRecoProton] Finished" << endl;

}

void myEfficiency::bookRecoNeutrino(SpectrumLoader& loader, Cut cut){

  cout << "[myEfficiency::bookRecoNeutrino] called" << endl;

  //==== some binnings
  double xEnergyMin = 0.; // GeV
  double xEnergyMax = 10.;
  double dxEnergy = 0.1;
  const Binning binsEnergy = Binning::Simple( int( (xEnergyMax-xEnergyMin)/dxEnergy ), xEnergyMin, xEnergyMax );

  cout << "[myEfficiency::bookRecoNeutrino] Delcaring HistAxis" << endl;

  //==== HistAxis
  const HistAxis axNeutrinoCaloEnergy("NeutrinoCaloEnergy", binsEnergy, varNeutrinoCaloEnergy);
  const HistAxis axNeutrinoCombinedEnergy("NeutrinoCombinedEnergy", binsEnergy, varNeutrinoCombinedEnergy);
  const HistAxis axNeutrinoQE("NeutrinoQE", binsEnergy, varNeutrinoQE);
  const HistAxis axNeutrinoFakeRecoEnergy("NeutrinoFakeRecoEnergy", binsEnergy, varNeutrinoFakeRecoEnergy);

  cout << "[myEfficiency::bookRecoNeutrino] Delcaring Spectrum" << endl;

  //==== Spectrum
  Spectrum *sNeutrinoCaloEnergy = new Spectrum(loader, axNeutrinoCaloEnergy, kNoSpillCut, cut);
  Spectrum *sNeutrinoCombinedEnergy = new Spectrum(loader, axNeutrinoCombinedEnergy, kNoSpillCut, cut);
  Spectrum *sNeutrinoQE = new Spectrum(loader, axNeutrinoQE, kNoSpillCut, cut);
  Spectrum *sNeutrinoFakeRecoEnergy = new Spectrum(loader, axNeutrinoFakeRecoEnergy, kNoSpillCut, cut);

  cout << "[myEfficiency::bookRecoNeutrino] Saving spectrums" << endl;

  vec_Spectrums.push_back(sNeutrinoCaloEnergy);
  vec_Spectrums.push_back(sNeutrinoCombinedEnergy);
  vec_Spectrums.push_back(sNeutrinoQE);
  vec_Spectrums.push_back(sNeutrinoFakeRecoEnergy);

  cout << "[myEfficiency::bookRecoNeutrino] Finished" << endl;

}

void myEfficiency::saveHistograms(){

  outputfile->cd();
  for(unsigned int i=0; i<vec_Spectrums.size(); i++){

    if(vec_Spectrums.at(i)->GetBinnings().size()==1){
      TH1 *h = vec_Spectrums.at(i)->ToTH1(TargetPOT);
      TString hName = vec_Spectrums.at(i)->GetLabels()[0];
      cout << "[myEfficiency::saveHistograms] Writing TH1, \"" << hName << "\"" << endl;
      h->SetName(hName); 
      h->Write();
    }
    else if(vec_Spectrums.at(i)->GetBinnings().size()==2){
      TH2 *h = vec_Spectrums.at(i)->ToTH2(TargetPOT);
      TString hName = vec_Spectrums.at(i)->GetLabels()[0];
      cout << "[myEfficiency::saveHistograms] Writing TH2, \"" << hName << "\"" << endl;
      h->SetName(hName);
      h->Write();
    }


  }

}

myEfficiency::~myEfficiency(){

  cout << "[myEfficiency::~myEfficiency] called" << endl;

  outputfile->Close();
  vec_Spectrums.clear();

  cout << "[myEfficiency::~myEfficiency] Finished" << endl;

}

