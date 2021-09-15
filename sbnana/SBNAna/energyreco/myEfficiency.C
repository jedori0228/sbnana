#include "sbnana/SBNAna/energyreco/myEfficiency.h"
#include "sbnana/SBNAna/energyreco/myTempNuMIVariables.h"

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

void myEfficiency::bookTruth(SpectrumLoader& loader, Cut cut, SpillCut spillCut){

  cout << "[myEfficiency::bookTruth] called" << endl;
  //==== some binnings
  double xEnergyMin = 0.; // GeV
  double xEnergyMax = 10.;
  double dxEnergy = 0.1;
  const Binning binsEnergy = Binning::Simple( int( (xEnergyMax-xEnergyMin)/dxEnergy ), xEnergyMin, xEnergyMax );
  double xCosineThetaMin = -1;
  double xCosineThetaMax = 1.;
  double dxCosineTheta = 0.05;
  const Binning binsCosineTheta = Binning::Simple( int( (xCosineThetaMax-xCosineThetaMin)/dxCosineTheta ), xCosineThetaMin, xCosineThetaMax );

  cout << "[myEfficiency::bookTruth] Delcaring HistAxis" << endl;

  //==== HistAxis
  const HistAxis axNeutrinoTruthE("NeutrinoTruthE", binsEnergy, varNeutrinoTruthE);
  const HistAxis axMuonTruthP("MuonTruthP", binsEnergy, varMuonTruthP);
  const HistAxis axMuonTruthCosineTheta("MuonTruthCosineTheta", binsCosineTheta, varMuonTruthCosineTheta);
  const HistAxis axProtonTruthP("ProtonTruthP", binsEnergy, varProtonTruthP);
  const HistAxis axProtonTruthCosineTheta("ProtonTruthCosineTheta", binsEnergy, varProtonTruthCosineTheta);
  //==== matching
  const HistAxis axTruthMuonMatchedTrackRangeP("TruthMuonMatchedTrackRangeP", binsEnergy, varTruthMuonMatchedTrackRangeP);
  const HistAxis axTruthMuonMatchedTrackMCSP("TruthMuonMatchedTrackMCSP", binsEnergy, varTruthMuonMatchedTrackMCSP);
  const HistAxis axTruthMuonMatchedTrackCombinedP("TruthMuonMatchedTrackCombinedP", binsEnergy, varTruthMuonMatchedTrackCombinedP);
  const HistAxis axTruthProtonMatchedTrackRangeP("TruthProtonMatchedTrackRangeP", binsEnergy, varTruthProtonMatchedTrackRangeP);
  const HistAxis axTruthProtonMatchedTrackMCSP("TruthProtonMatchedTrackMCSP", binsEnergy, varTruthProtonMatchedTrackMCSP);
  const HistAxis axTruthProtonMatchedTrackCombinedP("TruthProtonMatchedTrackCombinedP", binsEnergy, varTruthProtonMatchedTrackCombinedP);

  cout << "[myEfficiency::bookTruth] Delcaring Spectrum" << endl;

  //==== Spectrum
  Spectrum *sNeutrinoTruthE = new Spectrum(loader, axNeutrinoTruthE, spillCut, cut);
  Spectrum *sMuonTruthP = new Spectrum(loader, axMuonTruthP, spillCut, cut);
  Spectrum *sMuonTruthCosineTheta = new Spectrum(loader, axMuonTruthCosineTheta, spillCut, cut);
  Spectrum *sProtonTruthP = new Spectrum(loader, axProtonTruthP, spillCut, cut);
  Spectrum *sProtonTruthCosineTheta = new Spectrum(loader, axProtonTruthCosineTheta, spillCut, cut);
  Spectrum *sTruthMuonMatchedTrackRangeP = new Spectrum(loader, axTruthMuonMatchedTrackRangeP, spillCut, cut);
  Spectrum *sTruthMuonMatchedTrackMCSP = new Spectrum(loader, axTruthMuonMatchedTrackMCSP, spillCut, cut);
  Spectrum *sTruthMuonMatchedTrackCombinedP = new Spectrum(loader, axTruthMuonMatchedTrackCombinedP, spillCut, cut);
  Spectrum *sTruthProtonMatchedTrackRangeP = new Spectrum(loader, axTruthProtonMatchedTrackRangeP, spillCut, cut);
  //Spectrum *sTruthProtonMatchedTrackMCSP = new Spectrum(loader, axTruthProtonMatchedTrackMCSP, spillCut, cut);
  Spectrum *sTruthProtonMatchedTrackCombinedP = new Spectrum(loader, axTruthProtonMatchedTrackCombinedP, spillCut, cut);

  cout << "[myEfficiency::bookTruth] Saving spectrums" << endl;

  vec_Spectrums.push_back(sNeutrinoTruthE);
  vec_Spectrums.push_back(sMuonTruthP);
  vec_Spectrums.push_back(sMuonTruthCosineTheta);
  vec_Spectrums.push_back(sProtonTruthP);
  vec_Spectrums.push_back(sProtonTruthCosineTheta);
  vec_Spectrums.push_back(sTruthMuonMatchedTrackRangeP);
  vec_Spectrums.push_back(sTruthMuonMatchedTrackMCSP);
  vec_Spectrums.push_back(sTruthMuonMatchedTrackCombinedP);
  vec_Spectrums.push_back(sTruthProtonMatchedTrackRangeP);
  //vec_Spectrums.push_back(sTruthProtonMatchedTrackMCSP);
  vec_Spectrums.push_back(sTruthProtonMatchedTrackCombinedP);
  cout << "[myEfficiency::bookTruth] Finished" << endl;

}

void myEfficiency::bookRecoMuon(SpectrumLoader& loader, Cut cut, SpillCut spillCut){

  cout << "[myEfficiency::bookRecoMuon] called" << endl;

  //==== some binnings
  double xEnergyMin = 0.; // GeV
  double xEnergyMax = 10.;
  double dxEnergy = 0.1;
  const Binning binsEnergy = Binning::Simple( int( (xEnergyMax-xEnergyMin)/dxEnergy ), xEnergyMin, xEnergyMax );
  double xCosineThetaMin = -1;
  double xCosineThetaMax = 1.;
  double dxCosineTheta = 0.05;
  const Binning binsCosineTheta = Binning::Simple( int( (xCosineThetaMax-xCosineThetaMin)/dxCosineTheta ), xCosineThetaMin, xCosineThetaMax );
  const Binning binsEResidual = Binning::Simple(600, -3., 3.);
  const Binning binsEResidualFraction = Binning::Simple(600, -3., 3.);
  const Binning binsCosineThetaResidualFraction = Binning::Simple(600, -0.3, 0.3);
  const Binning binsPDG = Binning::Simple(100000, -5000., 5000.);

  cout << "[myEfficiency::bookRecoMuon] Delcaring HistAxis" << endl;

  //==== HistAxis
  const HistAxis axTempNuMIMuonRecoP("TempNuMIMuonRecoP", binsEnergy, varTempNuMIMuonRecoP);
  const HistAxis axTempNuMIMuonTrueP("TempNuMIMuonTrueP", binsEnergy, varTempNuMIMuonTrueP);
  const HistAxis axTempNuMIMuonTruePDG("TempNuMIMuonTruePDG", binsPDG, varTempNuMIMuonTruePDG);
  const HistAxis axTempNuMIMuonPResidual("TempNuMIMuonPResidual", binsEResidual, varTempNuMIMuonPResidual);
  const HistAxis axTempNuMIMuonPResidualFraction("TempNuMIMuonPResidualFraction", binsEResidualFraction, varTempNuMIMuonPResidualFraction);
  const HistAxis axTempNuMIMuonRecoCosineTheta("TempNuMIMuonRecoCosineTheta", binsCosineTheta, varTempNuMIMuonRecoCosineTheta);
  const HistAxis axTempNuMIMuonTrueCosineTheta("TempNuMIMuonTrueCosineTheta", binsCosineTheta, varTempNuMIMuonTrueCosineTheta);
  const HistAxis axTempNuMIMuonCosineThetaResidual("TempNuMIMuonCosineThetaResidual", binsEResidual, varTempNuMIMuonCosineThetaResidual);
  const HistAxis axTempNuMIMuonCosineThetaResidualFraction("TempNuMIMuonCosineThetaResidualFraction", binsCosineThetaResidualFraction, varTempNuMIMuonCosineThetaResidualFraction);

  cout << "[myEfficiency::bookRecoMuon] Delcaring Spectrum" << endl;

  //==== Spectrum
  Spectrum *sTempNuMIMuonRecoP = new Spectrum(loader, axTempNuMIMuonRecoP, spillCut, cut);
  Spectrum *sTempNuMIMuonTrueP = new Spectrum(loader, axTempNuMIMuonTrueP, spillCut, cut);
  Spectrum *sTempNuMIMuonTruePDG = new Spectrum(loader, axTempNuMIMuonTruePDG, spillCut, cut);
  Spectrum *sTempNuMIMuon_TrueP_vs_RecoP = new Spectrum(loader, axTempNuMIMuonTrueP, axTempNuMIMuonRecoP, spillCut, cut);
  Spectrum *sTempNuMIMuonPResidual = new Spectrum(loader, axTempNuMIMuonPResidual, spillCut, cut);
  Spectrum *sTempNuMIMuonPResidualFraction = new Spectrum(loader, axTempNuMIMuonPResidualFraction, spillCut, cut);
  Spectrum *sTempNuMIMuonRecoCosineTheta = new Spectrum(loader, axTempNuMIMuonRecoCosineTheta, spillCut, cut);
  Spectrum *sTempNuMIMuonTrueCosineTheta = new Spectrum(loader, axTempNuMIMuonTrueCosineTheta, spillCut, cut);
  Spectrum *sTempNuMIMuon_TrueCosineTheta_vs_RecoCosineTheta = new Spectrum(loader, axTempNuMIMuonTrueCosineTheta, axTempNuMIMuonRecoCosineTheta, spillCut, cut);
  Spectrum *sTempNuMIMuonCosineThetaResidual = new Spectrum(loader, axTempNuMIMuonCosineThetaResidual, spillCut, cut);
  Spectrum *sTempNuMIMuonCosineThetaResidualFraction = new Spectrum(loader, axTempNuMIMuonCosineThetaResidualFraction, spillCut, cut);
  Spectrum *sTempNuMIMuon_True_P_vs_CosineTheta = new Spectrum(loader, axTempNuMIMuonTrueP, axTempNuMIMuonTrueCosineTheta, spillCut, cut);
  Spectrum *sTempNuMIMuon_Reco_P_vs_CosineTheta = new Spectrum(loader, axTempNuMIMuonRecoP, axTempNuMIMuonRecoCosineTheta, spillCut, cut);
  //==== x : cosine theta, y : reco-true / true
  Spectrum *sTempNuMIMuon_TrueCosineTheta_vs_TrueCosineThetaResidualFraction = new Spectrum(loader, axTempNuMIMuonTrueCosineTheta, axTempNuMIMuonCosineThetaResidualFraction, spillCut, cut);

  cout << "[myEfficiency::bookRecoMuon] Saving spectrums" << endl;

  vec_Spectrums.push_back(sTempNuMIMuonRecoP);
  vec_Spectrums.push_back(sTempNuMIMuonTrueP);
  vec_Spectrums.push_back(sTempNuMIMuonTruePDG);
  vec_Spectrums.push_back(sTempNuMIMuon_TrueP_vs_RecoP);
  vec_Spectrums.push_back(sTempNuMIMuonPResidual);
  vec_Spectrums.push_back(sTempNuMIMuonPResidualFraction);
  vec_Spectrums.push_back(sTempNuMIMuonRecoCosineTheta);
  vec_Spectrums.push_back(sTempNuMIMuonTrueCosineTheta);
  vec_Spectrums.push_back(sTempNuMIMuon_TrueCosineTheta_vs_RecoCosineTheta);
  vec_Spectrums.push_back(sTempNuMIMuonCosineThetaResidual);
  vec_Spectrums.push_back(sTempNuMIMuonCosineThetaResidualFraction);
  vec_Spectrums.push_back(sTempNuMIMuon_True_P_vs_CosineTheta);
  vec_Spectrums.push_back(sTempNuMIMuon_Reco_P_vs_CosineTheta);
  vec_Spectrums.push_back(sTempNuMIMuon_TrueCosineTheta_vs_TrueCosineThetaResidualFraction);

  cout << "[myEfficiency::bookRecoMuon] Finished" << endl;

}

void myEfficiency::bookRecoProton(SpectrumLoader& loader, Cut cut, SpillCut spillCut){

  cout << "[myEfficiency::bookRecoProton] called" << endl;

  //==== some binnings
  double xEnergyMin = 0.; // GeV
  double xEnergyMax = 10.;
  double dxEnergy = 0.1;
  const Binning binsEnergy = Binning::Simple( int( (xEnergyMax-xEnergyMin)/dxEnergy ), xEnergyMin, xEnergyMax );
  double xCosineThetaMin = -1;
  double xCosineThetaMax = 1.;
  double dxCosineTheta = 0.05;
  const Binning binsCosineTheta = Binning::Simple( int( (xCosineThetaMax-xCosineThetaMin)/dxCosineTheta ), xCosineThetaMin, xCosineThetaMax );
  const Binning binsEResidual = Binning::Simple(600, -3., 3.);
  const Binning binsEResidualFraction = Binning::Simple(600, -3., 3.);
  const Binning binsPDG = Binning::Simple(100000, -5000., 5000.);
  const Binning binsN = Binning::Simple(5, 0., 5.);

  cout << "[myEfficiency::bookRecoProton] Delcaring HistAxis" << endl;

  //==== HistAxis
  const HistAxis axTempNuMIProtonRecoP("TempNuMIProtonRecoP", binsEnergy, varTempNuMIProtonRecoP);
  const HistAxis axTempNuMIProtonCaloP("TempNuMIProtonCaloP", binsEnergy, varTempNuMIProtonCaloP);
  const HistAxis axTempNuMIProtonTrueP("TempNuMIProtonTrueP", binsEnergy, varTempNuMIProtonTrueP);
  const HistAxis axTempNuMIProtonTruePDG("TempNuMIProtonTruePDG", binsPDG, varTempNuMIProtonTruePDG);
  const HistAxis axTempNuMIProtonPResidual("TempNuMIProtonPResidual", binsEResidual, varTempNuMIProtonPResidual);
  const HistAxis axTempNuMIProtonPResidualFraction("TempNuMIProtonPResidualFraction", binsEResidualFraction, varTempNuMIProtonPResidualFraction);
  const HistAxis axTempNuMIProtonRecoCosineTheta("TempNuMIProtonRecoCosineTheta", binsCosineTheta, varTempNuMIProtonRecoCosineTheta);
  const HistAxis axTempNuMIProtonTrueCosineTheta("TempNuMIProtonTrueCosineTheta", binsCosineTheta, varTempNuMIProtonTrueCosineTheta);
  const HistAxis axTempNuMIProtonCosineThetaResidual("TempNuMIProtonCosineThetaResidual", binsEResidual, varTempNuMIProtonCosineThetaResidual);
  const HistAxis axTempNuMIProtonCosineThetaResidualFraction("TempNuMIProtonCosineThetaResidualFraction", binsEResidualFraction, varTempNuMIProtonCosineThetaResidualFraction);
  //==== stub
  const HistAxis axTempNuMINStub("TempNuMINStub", binsN, varTempNuMINStub);

  cout << "[myEfficiency::bookRecoProton] Delcaring Spectrum" << endl;

  //==== Spectrum
  Spectrum *sTempNuMIProtonRecoP = new Spectrum(loader, axTempNuMIProtonRecoP, spillCut, cut);
  Spectrum *sTempNuMIProtonCaloP = new Spectrum(loader, axTempNuMIProtonCaloP, spillCut, cut);
  Spectrum *sTempNuMIProtonTrueP = new Spectrum(loader, axTempNuMIProtonTrueP, spillCut, cut);
  Spectrum *sTempNuMIProtonTruePDG = new Spectrum(loader, axTempNuMIProtonTruePDG, spillCut, cut);
  Spectrum *sTempNuMIProton_TrueP_vs_RecoP = new Spectrum(loader, axTempNuMIProtonTrueP, axTempNuMIProtonRecoP, spillCut, cut);
  Spectrum *sTempNuMIProtonPResidual = new Spectrum(loader, axTempNuMIProtonPResidual, spillCut, cut);
  Spectrum *sTempNuMIProtonPResidualFraction = new Spectrum(loader, axTempNuMIProtonPResidualFraction, spillCut, cut);
  Spectrum *sTempNuMIProtonRecoCosineTheta = new Spectrum(loader, axTempNuMIProtonRecoCosineTheta, spillCut, cut);
  Spectrum *sTempNuMIProtonTrueCosineTheta = new Spectrum(loader, axTempNuMIProtonTrueCosineTheta, spillCut, cut);
  Spectrum *sTempNuMIProton_TrueCosineTheta_vs_RecoCosineTheta = new Spectrum(loader, axTempNuMIProtonTrueCosineTheta, axTempNuMIProtonRecoCosineTheta, spillCut, cut);
  Spectrum *sTempNuMIProtonCosineThetaResidual = new Spectrum(loader, axTempNuMIProtonCosineThetaResidual, spillCut, cut);
  Spectrum *sTempNuMIProtonCosineThetaResidualFraction = new Spectrum(loader, axTempNuMIProtonCosineThetaResidualFraction, spillCut, cut);
  Spectrum *sTempNuMIProton_True_P_vs_CosineTheta = new Spectrum(loader, axTempNuMIProtonTrueP, axTempNuMIProtonTrueCosineTheta, spillCut, cut);
  Spectrum *sTempNuMIProton_Reco_P_vs_CosineTheta = new Spectrum(loader, axTempNuMIProtonRecoP, axTempNuMIProtonRecoCosineTheta, spillCut, cut);
  Spectrum *sTempNuMINStub = new Spectrum(loader, axTempNuMINStub, spillCut, cut);

  cout << "[myEfficiency::bookRecoProton] Saving spectrums" << endl;

  vec_Spectrums.push_back(sTempNuMIProtonRecoP);
  vec_Spectrums.push_back(sTempNuMIProtonCaloP);
  vec_Spectrums.push_back(sTempNuMIProtonTrueP);
  vec_Spectrums.push_back(sTempNuMIProtonTruePDG);
  vec_Spectrums.push_back(sTempNuMIProton_TrueP_vs_RecoP);
  vec_Spectrums.push_back(sTempNuMIProtonPResidual);
  vec_Spectrums.push_back(sTempNuMIProtonPResidualFraction);
  vec_Spectrums.push_back(sTempNuMIProtonRecoCosineTheta);
  vec_Spectrums.push_back(sTempNuMIProtonTrueCosineTheta);
  vec_Spectrums.push_back(sTempNuMIProton_TrueCosineTheta_vs_RecoCosineTheta);
  vec_Spectrums.push_back(sTempNuMIProtonCosineThetaResidual);
  vec_Spectrums.push_back(sTempNuMIProtonCosineThetaResidualFraction);
  vec_Spectrums.push_back(sTempNuMIProton_True_P_vs_CosineTheta);
  vec_Spectrums.push_back(sTempNuMIProton_Reco_P_vs_CosineTheta);
  vec_Spectrums.push_back(sTempNuMINStub);

  cout << "[myEfficiency::bookRecoProton] Finished" << endl;

}

void myEfficiency::bookRecoNeutrino(SpectrumLoader& loader, Cut cut, SpillCut spillCut){

  cout << "[myEfficiency::bookRecoNeutrino] called" << endl;

  //==== some binnings
  double xEnergyMin = 0.; // GeV
  double xEnergyMax = 10.;
  double dxEnergy = 0.1;
  const Binning binsEnergy = Binning::Simple( int( (xEnergyMax-xEnergyMin)/dxEnergy ), xEnergyMin, xEnergyMax );
  const Binning binsEResidual = Binning::Simple(600, -3., 3.);
  const Binning binsEResidualFraction = Binning::Simple(600, -3., 3.);

  cout << "[myEfficiency::bookRecoNeutrino] Delcaring HistAxis" << endl;

  //==== HistAxis
  const HistAxis axTempNuMINeutrinoCombinedEnergy("TempNuMINeutrinoCombinedEnergy", binsEnergy, varTempNuMINeutrinoCombinedEnergy);
  const HistAxis axTempNuMINeutrinoQE("TempNuMINeutrinoQE", binsEnergy, varTempNuMINeutrinoQE);
  const HistAxis axTempNuMINeutrinoCombinedEnergyResidual("TempNuMINeutrinoCombinedEnergyResidual", binsEResidual, varTempNuMINeutrinoCombinedEnergyResidual);
  const HistAxis axTempNuMINeutrinoQEResidual("TempNuMINeutrinoQEResidual", binsEResidual, varTempNuMINeutrinoQEResidual);
  const HistAxis axTempNuMINeutrinoCombinedEnergyResidualFraction("TempNuMINeutrinoCombinedEnergyResidualFraction", binsEResidualFraction, varTempNuMINeutrinoCombinedEnergyResidualFraction);
  const HistAxis axTempNuMINeutrinoQEResidualFraction("TempNuMINeutrinoQEResidualFraction", binsEResidualFraction, varTempNuMINeutrinoQEResidualFraction);

  cout << "[myEfficiency::bookRecoNeutrino] Delcaring Spectrum" << endl;

  //==== Spectrum
  Spectrum *sTempNuMINeutrinoCombinedEnergy = new Spectrum(loader, axTempNuMINeutrinoCombinedEnergy, spillCut, cut);
  Spectrum *sTempNuMINeutrinoQE = new Spectrum(loader, axTempNuMINeutrinoQE, spillCut, cut);
  Spectrum *sTempNuMINeutrinoCombinedEnergyResidual = new Spectrum(loader, axTempNuMINeutrinoCombinedEnergyResidual, spillCut, cut);
  Spectrum *sTempNuMINeutrinoQEResidual = new Spectrum(loader, axTempNuMINeutrinoQEResidual, spillCut, cut);
  Spectrum *sTempNuMINeutrinoCombinedEnergyResidualFraction = new Spectrum(loader, axTempNuMINeutrinoCombinedEnergyResidualFraction, spillCut, cut);
  Spectrum *sTempNuMINeutrinoQEResidualFraction = new Spectrum(loader, axTempNuMINeutrinoQEResidualFraction, spillCut, cut);

  //==== Systematics
  for(unsigned int i=0; i<systWs.size(); i++){
    Spectrum *sTempNuMINeutrinoCombinedEnergySyst = new Spectrum(loader, axTempNuMINeutrinoCombinedEnergy, spillCut, cut, kNoShift, systWs.at(i));
    std::string baseLabel = sTempNuMINeutrinoCombinedEnergySyst->GetLabels()[0];
    std::string newLabel = baseLabel+"_"+ISysts.at(i)->ShortName();
    sTempNuMINeutrinoCombinedEnergySyst->SetLabel(0, newLabel);
    cout << "[myEfficiency::bookRecoNeutrino] newLabel = " << newLabel << endl;
    vec_Spectrums.push_back(sTempNuMINeutrinoCombinedEnergySyst);
  }


  cout << "[myEfficiency::bookRecoNeutrino] Saving spectrums" << endl;

  vec_Spectrums.push_back(sTempNuMINeutrinoCombinedEnergy);
  vec_Spectrums.push_back(sTempNuMINeutrinoQE);
  vec_Spectrums.push_back(sTempNuMINeutrinoCombinedEnergyResidual);
  vec_Spectrums.push_back(sTempNuMINeutrinoQEResidual);
  vec_Spectrums.push_back(sTempNuMINeutrinoCombinedEnergyResidualFraction);
  vec_Spectrums.push_back(sTempNuMINeutrinoQEResidualFraction);

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

void myEfficiency::setSystematicWeights(){

  ISysts = GetSBNGenieWeightSysts();
  for(unsigned int i=0; i<ISysts.size(); ++i){
    systWs.push_back(GetUniverseWeight(ISysts, i));
  }

}

myEfficiency::~myEfficiency(){

  cout << "[myEfficiency::~myEfficiency] called" << endl;

  outputfile->Close();
  vec_Spectrums.clear();

  cout << "[myEfficiency::~myEfficiency] Finished" << endl;

}

