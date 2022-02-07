#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_HistoProducer.h"

using namespace ICARUSNumuXsec;

HistoProducer::HistoProducer(){

  cout << "[HistoProducer::HistoProducer] called" << endl;
  TargetPOT = 6.6e20;
  str_TargetPOT = "6.6e20 POT";
  outputName = "output.root";
  currentCutName = "DefaultCutName";
  vec_cutNames.clear();
  map_cutName_to_vec_Spectrums.clear();
  map_cutName_to_vec_SystEnsembleSpectrumPairs.clear();

  doTruthMatch = false;
  doPerformanceStudy = false;
  fillBeamInfo = true;
  fillNominal = true;

  cout << "[HistoProducer::HistoProducer] Finished" << endl;

}

void HistoProducer::initialize(){

  cout << "[HistoProducer::initialize] outputDir = " << outputDir << endl;
  gSystem->mkdir(outputDir, kTRUE);
  outputfile = new TFile(outputDir+"/"+outputName,"RECREATE");

}

bool HistoProducer::setCut(TString cutName){

  currentCutName = cutName;
  if( std::find(vec_cutNames.begin(), vec_cutNames.end(), cutName) != vec_cutNames.end() ){
    cout << "[HistoProducer::setCut] cutName = " << cutName << " already exist.." << endl;
    return false;
  }
  else{
    cout << "[HistoProducer::setCut] adding cutName = " << cutName << endl;
    vec_cutNames.push_back(cutName);
    return true;
  }

}

void HistoProducer::bookSpectrums(SpectrumLoader& loader, SpillCut spillCut, Cut cut){

  cout << "[HistoProducer::bookSpectrums] called" << endl;

  //==== Binning
  double xFMScoreMin = -2.;
  double xFMScoreMax = 150.;
  double dxFMScore = 1.;
  const Binning binsFMScore = Binning::Simple( int( (xFMScoreMax-xFMScoreMin)/dxFMScore ), xFMScoreMin, xFMScoreMax );
  double xFMTimeMin = -62;
  double xFMTimeMax = 60.;
  double dxFMTime = 1.;
  const Binning binsFMTime = Binning::Simple( int( (xFMTimeMax-xFMTimeMin)/dxFMTime ), xFMTimeMin, xFMTimeMax );
  double xNuScoreMin = -0.1;
  double xNuScoreMax = 1.;
  double dxNuScore = 0.05;
  const Binning binsNuScore = Binning::Simple( int( (xNuScoreMax-xNuScoreMin)/dxNuScore ), xNuScoreMin, xNuScoreMax );
  double xEnergyMin = 0.; // GeV
  double xEnergyMax = 10.;
  double dxEnergy = 0.1;
  const Binning binsEnergy = Binning::Simple( int( (xEnergyMax-xEnergyMin)/dxEnergy ), xEnergyMin, xEnergyMax );
  double xMomentumMin = 0.; // GeV
  double xMomentumMax = 5.;
  double dxMomentum = 0.1;
  const Binning binsMomentum = Binning::Simple( int( (xMomentumMax-xMomentumMin)/dxMomentum ), xMomentumMin, xMomentumMax );
  double xCosineThetaMin = -1;
  double xCosineThetaMax = 1.;
  double dxCosineTheta = 0.1;
  const Binning binsCosineTheta = Binning::Simple( int( (xCosineThetaMax-xCosineThetaMin)/dxCosineTheta ), xCosineThetaMin, xCosineThetaMax );
  const Binning binsOne = Binning::Simple( 1, 0., 1.);
  const Binning binsLength = Binning::Simple(2000, 0., 2000.);
  const Binning binsNhits = Binning::Simple(200, 0., 2000.);
  const Binning binsCharge = Binning::Simple(1000, 0., 10000.);
  const Binning binsXPosition = Binning::Simple(1000, -500., 500.);
  const Binning binsYPosition = Binning::Simple(400, -200., 200.);
  const Binning binsZPosition = Binning::Simple(2000, -1000., 1000.);
  const Binning binsChi2 = Binning::Simple(400, 0., 400);

  //==== HistAxis
  //====   Slice variables
  const HistAxis axFMScore("FMScore", binsFMScore, varFMScore);
  const HistAxis axFMTime("FMTime", binsFMTime, varFMTime);
  const HistAxis axNuScore("NuScore", binsNuScore, varNuScore);
  const HistAxis axCountSlice("CountSlice", binsOne, varCountSlice);
  const HistAxis axSliceTrackNhitsPlane0("SliceTrackNhitsPlane0", binsNhits, varSliceTrackNhitsPlane0);
  const HistAxis axSliceTrackNhitsPlane1("SliceTrackNhitsPlane1", binsNhits, varSliceTrackNhitsPlane1);
  const HistAxis axSliceTrackNhitsPlane2("SliceTrackNhitsPlane2", binsNhits, varSliceTrackNhitsPlane2);
  const HistAxis axSliceShowerNhitsPlane0("SliceShowerNhitsPlane0", binsNhits, varSliceShowerNhitsPlane0);
  const HistAxis axSliceShowerNhitsPlane1("SliceShowerNhitsPlane1", binsNhits, varSliceShowerNhitsPlane1);
  const HistAxis axSliceShowerNhitsPlane2("SliceShowerNhitsPlane2", binsNhits, varSliceShowerNhitsPlane2);
  const HistAxis axSliceTrackChargePlane0("SliceTrackChargePlane0", binsCharge, varSliceTrackChargePlane0);
  const HistAxis axSliceTrackChargePlane1("SliceTrackChargePlane1", binsCharge, varSliceTrackChargePlane1);
  const HistAxis axSliceTrackChargePlane2("SliceTrackChargePlane2", binsCharge, varSliceTrackChargePlane2);
  //====   Truth variables
  //====     Neutrinos
  const HistAxis axNeutrinoTruthE("NeutrinoTruthE", binsEnergy, varNeutrinoTruthE);
  const HistAxis axTruthQ2("TruthQ2", binsEnergy, varTruthQ2);
  const HistAxis axTruthq0_lab("Truthq0_lab", binsEnergy, varTruthq0_lab);
  const HistAxis axTruthmodq_lab("Truthmodq_lab", binsEnergy, varTruthmodq_lab);
  //====     Muon
  const HistAxis axMuonTruthP("MuonTruthP", binsMomentum, varMuonTruthP);
  const HistAxis axMuonTruthCosineTheta("MuonTruthCosineTheta", binsCosineTheta, varMuonTruthCosineTheta);
  const HistAxis axMuonTruthNuMICosineTheta("MuonTruthNuMICosineTheta", binsCosineTheta, varMuonTruthNuMICosineTheta);
  const HistAxis axTruthMuonMatchedTrackChi2Muon("TruthMuonMatchedTrackChi2Muon", binsChi2, varTruthMuonMatchedTrackChi2Muon);
  const HistAxis axTruthMuonMatchedTrackChi2Proton("TruthMuonMatchedTrackChi2Proton", binsChi2, varTruthMuonMatchedTrackChi2Proton);
  //====     Proton
  const HistAxis axProtonTruthP("ProtonTruthP", binsMomentum, varProtonTruthP);
  const HistAxis axProtonTruthCosineTheta("ProtonTruthCosineTheta", binsCosineTheta, varProtonTruthCosineTheta);
  const HistAxis axProtonTruthNuMICosineTheta("ProtonTruthNuMICosineTheta", binsCosineTheta, varProtonTruthNuMICosineTheta);
  const HistAxis axTruthProtonMatchedTrackChi2Proton("TruthProtonMatchedTrackChi2Proton", binsChi2, varTruthProtonMatchedTrackChi2Proton);
  //====     Muon+Proton
  const HistAxis axTruthMuonProtonCosineTheta("TruthMuonProtonCosineTheta", binsCosineTheta, varTruthMuonProtonCosineTheta);
  //====   Reco variables
  //====     Muon
  const HistAxis axMuonRecoP("MuonRecoP", binsMomentum, varMuonRecoP);
  const HistAxis axMuonRecoCosineTheta("MuonRecoCosineTheta", binsCosineTheta, varMuonRecoCosineTheta);
  const HistAxis axMuonRecoNuMICosineTheta("MuonRecoNuMICosineTheta", binsCosineTheta, varMuonRecoNuMICosineTheta);
  const HistAxis axMuonLength("MuonLength", binsLength, varMuonLength);
  const HistAxis axMuonChi2Muon("MuonChi2Muon", binsChi2, varMuonChi2Muon);
  const HistAxis axMuonChi2Proton("MuonChi2Proton", binsChi2, varMuonChi2Proton);
  //====     Proton
  const HistAxis axProtonRecoP("ProtonRecoP", binsMomentum, varProtonRecoP);
  const HistAxis axProtonRecoCosineTheta("ProtonRecoCosineTheta", binsCosineTheta, varProtonRecoCosineTheta);
  const HistAxis axProtonRecoNuMICosineTheta("ProtonRecoNuMICosineTheta", binsCosineTheta, varProtonRecoNuMICosineTheta);
  const HistAxis axProtonLength("ProtonLength", binsLength, varProtonLength);
  const HistAxis axProtonChi2Proton("ProtonChi2Proton", binsChi2, varMuonChi2Proton);
  //====     Neutrino
  const HistAxis axNeutrinoCombinedEnergy("NeutrinoCombinedEnergy", binsEnergy, varNeutrinoCombinedEnergy);
  const HistAxis axNeutrinoQE("NeutrinoQE", binsEnergy, varNeutrinoQE);

  //==== Spectrum
  //====   Slice variables
  Spectrum *sFMScore = new Spectrum(loader, axFMScore, spillCut, cut);
  Spectrum *sFMTime = new Spectrum(loader, axFMTime, spillCut, cut);
  Spectrum *sNuScore = new Spectrum(loader, axNuScore, spillCut, cut);
  Spectrum *sCountSlice = new Spectrum(loader, axCountSlice, spillCut, cut);
  Spectrum *sSliceTrackNhitsPlane0 = new Spectrum(loader, axSliceTrackNhitsPlane0, spillCut, cut);
  Spectrum *sSliceTrackNhitsPlane1 = new Spectrum(loader, axSliceTrackNhitsPlane1, spillCut, cut);
  Spectrum *sSliceTrackNhitsPlane2 = new Spectrum(loader, axSliceTrackNhitsPlane2, spillCut, cut);
  Spectrum *sSliceShowerNhitsPlane0 = new Spectrum(loader, axSliceShowerNhitsPlane0, spillCut, cut);
  Spectrum *sSliceShowerNhitsPlane1 = new Spectrum(loader, axSliceShowerNhitsPlane1, spillCut, cut);
  Spectrum *sSliceShowerNhitsPlane2 = new Spectrum(loader, axSliceShowerNhitsPlane2, spillCut, cut);
  Spectrum *sSliceTrackChargePlane0 = new Spectrum(loader, axSliceTrackChargePlane0, spillCut, cut);
  Spectrum *sSliceTrackChargePlane1 = new Spectrum(loader, axSliceTrackChargePlane1, spillCut, cut);
  Spectrum *sSliceTrackChargePlane2 = new Spectrum(loader, axSliceTrackChargePlane2, spillCut, cut);
  Spectrum *sAllTrackStartPositionX = new Spectrum("AllTrackStartPositionX", binsXPosition, loader, varAllTrackStartPositionX, spillCut, cut);
  Spectrum *sAllTrackStartPositionY = new Spectrum("AllTrackStartPositionY", binsYPosition, loader, varAllTrackStartPositionY, spillCut, cut);
  Spectrum *sAllTrackStartPositionZ = new Spectrum("AllTrackStartPositionZ", binsZPosition, loader, varAllTrackStartPositionZ, spillCut, cut);
  //====   Truth variables
  //====     Neutrinos
  Spectrum *sNeutrinoTruthE = new Spectrum(loader, axNeutrinoTruthE, spillCut, cut);
  Spectrum *sTruthQ2 = new Spectrum(loader, axTruthQ2, spillCut, cut);
  Spectrum *sTruthq0_lab = new Spectrum(loader, axTruthq0_lab, spillCut, cut);
  Spectrum *sTruthmodq_lab = new Spectrum(loader, axTruthmodq_lab, spillCut, cut);
  //====     Muon
  Spectrum *sMuonTruthP = new Spectrum(loader, axMuonTruthP, spillCut, cut);
  Spectrum *sMuonTruthCosineTheta = new Spectrum(loader, axMuonTruthCosineTheta, spillCut, cut);
  Spectrum *sMuonTruthNuMICosineTheta = new Spectrum(loader, axMuonTruthNuMICosineTheta, spillCut, cut);
  Spectrum *sTruthMuonMatchedTrackChi2Muon = new Spectrum(loader, axTruthMuonMatchedTrackChi2Muon, spillCut, cut);
  Spectrum *sTruthMuonMatchedTrackChi2Proton = new Spectrum(loader, axTruthMuonMatchedTrackChi2Proton, spillCut, cut);
  //====     Proton
  Spectrum *sProtonTruthP = new Spectrum(loader, axProtonTruthP, spillCut, cut);
  Spectrum *sProtonTruthCosineTheta = new Spectrum(loader, axProtonTruthCosineTheta, spillCut, cut);
  Spectrum *sProtonTruthNuMICosineTheta = new Spectrum(loader, axProtonTruthNuMICosineTheta, spillCut, cut);
  Spectrum *sTruthProtonMatchedTrackChi2Proton = new Spectrum(loader, axTruthProtonMatchedTrackChi2Proton, spillCut, cut);
  //====     Muon+Proton
  Spectrum *sTruthMuonProtonCosineTheta = new Spectrum(loader, axTruthMuonProtonCosineTheta, spillCut, cut);
  //====   Reco variables
  //====     Muon
  Spectrum *sMuonRecoP = new Spectrum(loader, axMuonRecoP, spillCut, cut);
  Spectrum *sMuonRecoCosineTheta = new Spectrum(loader, axMuonRecoCosineTheta, spillCut, cut);
  Spectrum *sMuonRecoNuMICosineTheta = new Spectrum(loader, axMuonRecoNuMICosineTheta, spillCut, cut);
  Spectrum *sMuonLength = new Spectrum(loader, axMuonLength, spillCut, cut);
  Spectrum *sMuonChi2Muon = new Spectrum(loader, axMuonChi2Muon, spillCut, cut);
  Spectrum *sMuonChi2Proton = new Spectrum(loader, axMuonChi2Proton, spillCut, cut);
  //====     Proton
  Spectrum *sProtonRecoP = new Spectrum(loader, axProtonRecoP, spillCut, cut);
  Spectrum *sProtonRecoCosineTheta = new Spectrum(loader, axProtonRecoCosineTheta, spillCut, cut);
  Spectrum *sProtonRecoNuMICosineTheta = new Spectrum(loader, axProtonRecoNuMICosineTheta, spillCut, cut);
  Spectrum *sProtonLength = new Spectrum(loader, axProtonLength, spillCut, cut);
  Spectrum *sProtonChi2Proton = new Spectrum(loader, axProtonChi2Proton, spillCut, cut);
  //====     Neutrino
  Spectrum *sNeutrinoCombinedEnergy = new Spectrum(loader, axNeutrinoCombinedEnergy, spillCut, cut);
  Spectrum *sNeutrinoQE = new Spectrum(loader, axNeutrinoQE, spillCut, cut);
  //====   Reco vs Truth
  Spectrum *sMuonRecoP_vs_MuonTruthP = new Spectrum(loader, axMuonRecoP, axMuonTruthP, spillCut, cut);
  Spectrum *sMuonRecoCosineTheta_vs_MuonTruthCosineTheta = new Spectrum(loader, axMuonRecoCosineTheta, axMuonTruthCosineTheta, spillCut, cut);
  Spectrum *sMuonRecoNuMICosineTheta_vs_MuonTruthNuMICosineTheta = new Spectrum(loader, axMuonRecoNuMICosineTheta, axMuonTruthNuMICosineTheta, spillCut, cut);
  Spectrum *sNeutrinoCombinedEnergy_vs_NeutrinoTruthE = new Spectrum(loader, axNeutrinoCombinedEnergy, axNeutrinoTruthE, spillCut, cut);

  if(fillNominal){
    map_cutName_to_vec_Spectrums[currentCutName].push_back(sFMScore);
    map_cutName_to_vec_Spectrums[currentCutName].push_back(sFMTime);
    map_cutName_to_vec_Spectrums[currentCutName].push_back(sNuScore);
    map_cutName_to_vec_Spectrums[currentCutName].push_back(sCountSlice);
    map_cutName_to_vec_Spectrums[currentCutName].push_back(sSliceTrackNhitsPlane0);
    map_cutName_to_vec_Spectrums[currentCutName].push_back(sSliceTrackNhitsPlane1);
    map_cutName_to_vec_Spectrums[currentCutName].push_back(sSliceTrackNhitsPlane2);
    map_cutName_to_vec_Spectrums[currentCutName].push_back(sSliceShowerNhitsPlane0);
    map_cutName_to_vec_Spectrums[currentCutName].push_back(sSliceShowerNhitsPlane1);
    map_cutName_to_vec_Spectrums[currentCutName].push_back(sSliceShowerNhitsPlane2);
    map_cutName_to_vec_Spectrums[currentCutName].push_back(sSliceTrackChargePlane0);
    map_cutName_to_vec_Spectrums[currentCutName].push_back(sSliceTrackChargePlane1);
    map_cutName_to_vec_Spectrums[currentCutName].push_back(sSliceTrackChargePlane2);
    map_cutName_to_vec_Spectrums[currentCutName].push_back(sAllTrackStartPositionX);
    map_cutName_to_vec_Spectrums[currentCutName].push_back(sAllTrackStartPositionY);
    map_cutName_to_vec_Spectrums[currentCutName].push_back(sAllTrackStartPositionZ);
    map_cutName_to_vec_Spectrums[currentCutName].push_back(sNeutrinoTruthE);
    map_cutName_to_vec_Spectrums[currentCutName].push_back(sTruthQ2);
    map_cutName_to_vec_Spectrums[currentCutName].push_back(sTruthq0_lab);
    map_cutName_to_vec_Spectrums[currentCutName].push_back(sTruthmodq_lab);
    map_cutName_to_vec_Spectrums[currentCutName].push_back(sMuonTruthP);
    map_cutName_to_vec_Spectrums[currentCutName].push_back(sMuonTruthCosineTheta);
    map_cutName_to_vec_Spectrums[currentCutName].push_back(sMuonTruthNuMICosineTheta);
    map_cutName_to_vec_Spectrums[currentCutName].push_back(sTruthMuonMatchedTrackChi2Muon);
    map_cutName_to_vec_Spectrums[currentCutName].push_back(sTruthMuonMatchedTrackChi2Proton);
    map_cutName_to_vec_Spectrums[currentCutName].push_back(sProtonTruthP);
    map_cutName_to_vec_Spectrums[currentCutName].push_back(sProtonTruthCosineTheta);
    map_cutName_to_vec_Spectrums[currentCutName].push_back(sProtonTruthNuMICosineTheta);
    map_cutName_to_vec_Spectrums[currentCutName].push_back(sTruthProtonMatchedTrackChi2Proton);
    map_cutName_to_vec_Spectrums[currentCutName].push_back(sTruthMuonProtonCosineTheta);
    map_cutName_to_vec_Spectrums[currentCutName].push_back(sMuonRecoP);
    map_cutName_to_vec_Spectrums[currentCutName].push_back(sMuonRecoCosineTheta);
    map_cutName_to_vec_Spectrums[currentCutName].push_back(sMuonRecoNuMICosineTheta);
    map_cutName_to_vec_Spectrums[currentCutName].push_back(sMuonLength);
    map_cutName_to_vec_Spectrums[currentCutName].push_back(sMuonChi2Muon);
    map_cutName_to_vec_Spectrums[currentCutName].push_back(sMuonChi2Proton);
    map_cutName_to_vec_Spectrums[currentCutName].push_back(sProtonRecoP);
    map_cutName_to_vec_Spectrums[currentCutName].push_back(sProtonRecoCosineTheta);
    map_cutName_to_vec_Spectrums[currentCutName].push_back(sProtonRecoNuMICosineTheta);
    map_cutName_to_vec_Spectrums[currentCutName].push_back(sProtonLength);
    map_cutName_to_vec_Spectrums[currentCutName].push_back(sProtonChi2Proton);
    map_cutName_to_vec_Spectrums[currentCutName].push_back(sNeutrinoCombinedEnergy);
    map_cutName_to_vec_Spectrums[currentCutName].push_back(sNeutrinoQE);
    map_cutName_to_vec_Spectrums[currentCutName].push_back(sMuonRecoP_vs_MuonTruthP);
    map_cutName_to_vec_Spectrums[currentCutName].push_back(sMuonRecoCosineTheta_vs_MuonTruthCosineTheta);
    map_cutName_to_vec_Spectrums[currentCutName].push_back(sMuonRecoNuMICosineTheta_vs_MuonTruthNuMICosineTheta);
    map_cutName_to_vec_Spectrums[currentCutName].push_back(sNeutrinoCombinedEnergy_vs_NeutrinoTruthE);
  }

  //==== Systematic
  for(unsigned int i_syst=0; i_syst<IGENIESysts.size(); i_syst++){
    const ISyst* s = IGENIESysts.at(i_syst);

    addUpDownSystematic(loader, axNeutrinoTruthE, spillCut, cut, currentCutName, s);
    addUpDownSystematic(loader, axMuonTruthP, spillCut, cut, currentCutName, s);
    addUpDownSystematic(loader, axMuonTruthCosineTheta, spillCut, cut, currentCutName, s);
    addUpDownSystematic(loader, axMuonTruthNuMICosineTheta, spillCut, cut, currentCutName, s);
    addUpDownSystematic(loader, axMuonRecoP, axMuonTruthP, spillCut, cut, currentCutName, s);
    addUpDownSystematic(loader, axMuonRecoCosineTheta, axMuonTruthCosineTheta, spillCut, cut, currentCutName, s);
    addUpDownSystematic(loader, axMuonRecoNuMICosineTheta, axMuonTruthNuMICosineTheta, spillCut, cut, currentCutName, s);
    addUpDownSystematic(loader, axNeutrinoCombinedEnergy, axNeutrinoTruthE, spillCut, cut, currentCutName, s);

  }

  for(unsigned int i_syst=0; i_syst<IFluxSysts.size(); i_syst++){
    const ISyst* s = IFluxSysts.at(i_syst);

    addUpDownSystematic(loader, axNeutrinoTruthE, spillCut, cut, currentCutName, s);
    addUpDownSystematic(loader, axMuonTruthP, spillCut, cut, currentCutName, s);
    addUpDownSystematic(loader, axMuonTruthCosineTheta, spillCut, cut, currentCutName, s);
    addUpDownSystematic(loader, axMuonTruthNuMICosineTheta, spillCut, cut, currentCutName, s);
    addUpDownSystematic(loader, axMuonRecoP, axMuonTruthP, spillCut, cut, currentCutName, s);
    addUpDownSystematic(loader, axMuonRecoCosineTheta, axMuonTruthCosineTheta, spillCut, cut, currentCutName, s);
    addUpDownSystematic(loader, axMuonRecoNuMICosineTheta, axMuonTruthNuMICosineTheta, spillCut, cut, currentCutName, s);
    addUpDownSystematic(loader, axNeutrinoCombinedEnergy, axNeutrinoTruthE, spillCut, cut, currentCutName, s);

  }

  for(unsigned int i_syst=0; i_syst<IDetectorSysts.size(); i_syst++){
    const ISyst* s = IDetectorSysts.at(i_syst);

    addUpDownSystematic(loader, axNeutrinoTruthE, spillCut, cut, currentCutName, s);
    addUpDownSystematic(loader, axMuonTruthP, spillCut, cut, currentCutName, s);
    addUpDownSystematic(loader, axMuonTruthCosineTheta, spillCut, cut, currentCutName, s);
    addUpDownSystematic(loader, axMuonTruthNuMICosineTheta, spillCut, cut, currentCutName, s);
    addUpDownSystematic(loader, axMuonRecoP, axMuonTruthP, spillCut, cut, currentCutName, s);
    addUpDownSystematic(loader, axMuonRecoCosineTheta, axMuonTruthCosineTheta, spillCut, cut, currentCutName, s);
    addUpDownSystematic(loader, axMuonRecoNuMICosineTheta, axMuonTruthNuMICosineTheta, spillCut, cut, currentCutName, s);
    addUpDownSystematic(loader, axNeutrinoCombinedEnergy, axNeutrinoTruthE, spillCut, cut, currentCutName, s);


  }

}

void HistoProducer::saveHistograms(){

  outputfile->cd();

  cout << "[HistoProducer::saveHistograms] Number of cuts = " << vec_cutNames.size() << endl;
  const unsigned int nCutName = vec_cutNames.size();
  if(nCutName>0){

    TH1D *hist_CutNames = new TH1D("hist_CutNames", "", nCutName, 0., 1.*nCutName);

    for(unsigned int ic=0; ic<nCutName; ic++){

      TString cutName = vec_cutNames.at(ic);
      //TString dirName = "cutName"+TString::Itoa(ic,10);
      TString dirName = cutName; // I'd rather make the cutName simpler and use it as dirName too
      TDirectory *dir = outputfile->GetDirectory(dirName);
      if(!dir){
        outputfile->mkdir(dirName);
        dir = outputfile->GetDirectory(dirName);
      }
      outputfile->cd(dirName);

      vector<Spectrum *> vec_Spectrums = map_cutName_to_vec_Spectrums[cutName];
      vector< pair<TString, EnsembleSpectrum *> > vec_SystEnsembleSpectrumPairs = map_cutName_to_vec_SystEnsembleSpectrumPairs[cutName];
      vector< pair<TString, Spectrum *> > vec_SystSpectrumPairs = map_cutName_to_vec_SystSpectrumPairs[cutName];

      //==== BinContent = number of Spectrum
      hist_CutNames->SetBinContent(ic+1, vec_Spectrums.size());
      hist_CutNames->GetXaxis()->SetBinLabel(ic+1, cutName);

      cout << "[HistoProducer::saveHistograms] cutName = " << cutName << endl;
      cout << "[HistoProducer::saveHistograms]   Directory name = " << dirName << endl;

      //==== Spectrum
      cout << "[HistoProducer::saveHistograms]   Number of Spectrum = " << vec_Spectrums.size() << endl;
      for(unsigned int i=0; i<vec_Spectrums.size(); i++){

        cout << "[HistoProducer::saveHistograms]   " << i << "-th Spectrum.." << endl;

        if(i==0 && fillBeamInfo){
          cout << "[HistoProducer::saveHistograms]     POT = " << vec_Spectrums.at(i)->POT() << endl;
          cout << "[HistoProducer::saveHistograms]     Livetime = " << vec_Spectrums.at(i)->Livetime() << endl;
          TH1D *hist_POT = new TH1D("POT_"+dirName, "POT", 1, 0., 1.);
          hist_POT->SetBinContent(1, vec_Spectrums.at(i)->POT());
          TH1D *hist_Livetime = new TH1D("Livetime_"+dirName, "Livetime", 1, 0., 1.);
          hist_Livetime->SetBinContent(1, vec_Spectrums.at(i)->Livetime());
          hist_POT->Write();
          hist_Livetime->Write();
        }

        TString hName = vec_Spectrums.at(i)->GetLabels()[0];

        if(vec_Spectrums.at(i)->GetBinnings().size()==1){
          TH1 *h = vec_Spectrums.at(i)->ToTH1(TargetPOT);
          cout << "[HistoProducer::saveHistograms]     Writing TH1, \"" << hName << "\"" << endl;
          h->SetName(hName+"_"+dirName);
          h->Write();
        }
        else if(vec_Spectrums.at(i)->GetBinnings().size()==2){
          TH2 *h = vec_Spectrums.at(i)->ToTH2(TargetPOT);
          cout << "[HistoProducer::saveHistograms]     Writing TH2, \"" << hName << "\"" << endl;
          h->SetName(hName+"_"+dirName);
          h->Write();
        }

      }

      //==== EnsembleSpectrum
      cout << "[HistoProducer::saveHistograms]   Number of SystematicEnsembleSpectrum = " << vec_SystEnsembleSpectrumPairs.size() << endl;
      for(unsigned int i=0; i<vec_SystEnsembleSpectrumPairs.size(); i++){

        cout << "[HistoProducer::saveHistograms]   " << i << "-th EnsembleSpectrum.." << endl;

        TString systematicName = vec_SystEnsembleSpectrumPairs.at(i).first;
        Spectrum sNominal = vec_SystEnsembleSpectrumPairs.at(i).second->Nominal();
        TString baseLabel = sNominal.GetLabels()[0];
        dir->mkdir(baseLabel+"_"+systematicName+"_"+dirName);
        dir->cd(baseLabel+"_"+systematicName+"_"+dirName);

        if(sNominal.GetBinnings().size()==1){
          TString newLabel = baseLabel+"_NominalFromES";
          TH1 *h = sNominal.ToTH1(TargetPOT);
          h->SetName(newLabel+"_"+dirName);
          h->Write();
          cout << "[HistoProducer::saveHistograms]     Nominal histogram = " << baseLabel << endl;
          for(unsigned int iu=0; iu<vec_SystEnsembleSpectrumPairs.at(i).second->NUniverses(); ++iu){
            TH1 *hU = vec_SystEnsembleSpectrumPairs.at(i).second->Universe(iu).ToTH1(TargetPOT);
            //TString systName = IGENIESysts.at(iu)->ShortName();
            hU->SetName(baseLabel+"_"+systematicName+"_Univ"+TString::Itoa(iu,10)+"_"+dirName);
            hU->Write();
          }
          TGraphAsymmErrors *grErrorBand = vec_SystEnsembleSpectrumPairs.at(i).second->ErrorBand(TargetPOT);
          grErrorBand->SetName(baseLabel+"_"+systematicName+"_ErrorBandFromES_"+dirName);
          grErrorBand->Write();
        }
        else if(sNominal.GetBinnings().size()==2){
          for(unsigned int iu=0; iu<vec_SystEnsembleSpectrumPairs.at(i).second->NUniverses(); ++iu){
            TH2 *hU = vec_SystEnsembleSpectrumPairs.at(i).second->Universe(iu).ToTH2(TargetPOT);
            //TString systName = IGENIESysts.at(iu)->ShortName();
            hU->SetName(baseLabel+"_"+systematicName+"_Univ"+TString::Itoa(iu,10)+"_"+dirName);
            hU->Write();
          }
        }

        dir->cd();

      }

      //==== Systematic
      cout << "[HistoProducer::saveHistograms]   Number of SystematicSpectrum = " << vec_SystSpectrumPairs.size() << endl;

      for(unsigned int i=0; i<vec_SystSpectrumPairs.size(); i++){

        cout << "[HistoProducer::saveHistograms]   " << i << "-th Systematic Spectrum.." << endl;

        TString systematicName = vec_SystSpectrumPairs.at(i).first;

        if(vec_SystSpectrumPairs.at(i).second->GetBinnings().size()==1){
          TH1 *h = vec_SystSpectrumPairs.at(i).second->ToTH1(TargetPOT);
          TString baseLabel = vec_SystSpectrumPairs.at(i).second->GetLabels()[0];
          TString newLabel = baseLabel+"_"+systematicName;
          h->SetName(newLabel+"_"+dirName);
          h->Write();
          cout << "[HistoProducer::saveHistograms]     --> Done" << endl;
        }
        else if(vec_SystSpectrumPairs.at(i).second->GetBinnings().size()==2){
          TH2 *h = vec_SystSpectrumPairs.at(i).second->ToTH2(TargetPOT);
          TString baseLabel = vec_SystSpectrumPairs.at(i).second->GetLabels()[0];
          TString newLabel = baseLabel+"_"+systematicName;
          h->SetName(newLabel+"_"+dirName);
          h->Write();
          cout << "[HistoProducer::saveHistograms]     --> Done" << endl;
        }

      }

      outputfile->cd();

    } // END loop cutname

    outputfile->cd();
    hist_CutNames->Write();

  } // END if nCutName>0

  if(fillNominal){
    outputfile->cd();
    TH1D *hist_TargetPOT = new TH1D("hist_TargetPOT", "", 1, 0., 1.);
    hist_TargetPOT->SetBinContent(1, TargetPOT);
    hist_TargetPOT->Write();
  }

}

void HistoProducer::setSystematicWeights(){

  cout << "[HistoProducer::setSystematicWeights] Setting GENIE systematics" << endl;
  IGENIESysts = GetSBNGenieWeightSysts();
  //==== For EnsembleSpectrum
  for(unsigned int i=0; i<IGENIESysts.size(); ++i){
    cout << "[HistoProducer::setSystematicWeights] " << i << " : " << IGENIESysts.at(i)->ShortName() << endl;

    vector<const ISyst*> this_GENIESyst;
    this_GENIESyst.clear();
    this_GENIESyst.push_back(IGENIESysts.at(i));

    vector<Var> weis;
    for(unsigned iu=0; iu<100; iu++){
      weis.push_back( GetUniverseWeight(this_GENIESyst, iu) );
    }
    vec_UniverseWeightsForEachGENIESource.push_back( weis );

  }
  cout << "[HistoProducer::setSystematicWeights] Setting flux systematics" << endl;
  IFluxSysts = GetAllNuMIFluxSysts(5);
  for(unsigned int i=0; i<IFluxSysts.size(); i++){
    cout << "[HistoProducer::setSystematicWeights] Syst = " << IFluxSysts.at(i)->ShortName() << endl;
  }

  cout << "[HistoProducer::setSystematicWeights] Setting detector systematics" << endl;
  IDetectorSysts.push_back( new MuonMomentumScaleSyst(0.02) );

}

void HistoProducer::addEnsembleSpectrum(SpectrumLoader& loader, const HistAxis& ax, SpillCut spillCut, Cut cut, TString currentCutName, vector<Var> varWeights, TString systName){

  map_cutName_to_vec_SystEnsembleSpectrumPairs[currentCutName].push_back( 
    std::make_pair( systName, new EnsembleSpectrum(loader, ax, spillCut, cut, varWeights) )
  );

}

void HistoProducer::addEnsembleSpectrum(SpectrumLoader& loader, const HistAxis& axX, const HistAxis& axY, SpillCut spillCut, Cut cut, TString currentCutName, vector<Var> varWeights, TString systName){

  map_cutName_to_vec_SystEnsembleSpectrumPairs[currentCutName].push_back(
    std::make_pair( systName, new EnsembleSpectrum(loader, axX, axY, spillCut, cut, varWeights) )
  );

}

void HistoProducer::addUpDownSystematic(SpectrumLoader& loader, const HistAxis& ax, SpillCut spillCut, Cut cut, TString currentCutName, const ISyst* s){

  map_cutName_to_vec_SystSpectrumPairs[currentCutName].push_back(
    std::make_pair( s->ShortName()+"_Up", new Spectrum(loader, ax, spillCut, cut, SystShifts(s, +1)) )
  );
  map_cutName_to_vec_SystSpectrumPairs[currentCutName].push_back(
    std::make_pair( s->ShortName()+"_Down", new Spectrum(loader, ax, spillCut, cut, SystShifts(s, -1)) )
  );

}

void HistoProducer::addUpDownSystematic(SpectrumLoader& loader, const HistAxis& axX, const HistAxis& axY, SpillCut spillCut, Cut cut, TString currentCutName, const ISyst* s){

  map_cutName_to_vec_SystSpectrumPairs[currentCutName].push_back(
    std::make_pair( s->ShortName()+"_Up", new Spectrum(loader, axX, axY, spillCut, cut, SystShifts(s, +1)) )
  );
  map_cutName_to_vec_SystSpectrumPairs[currentCutName].push_back(
    std::make_pair( s->ShortName()+"_Down", new Spectrum(loader, axX, axY, spillCut, cut, SystShifts(s, -1)) )
  );

}

HistoProducer::~HistoProducer(){

  cout << "[HistoProducer::~HistoProducer] called" << endl;

  outputfile->Close();

  cout << "[HistoProducer::~HistoProducer] output file : " << outputDir+"/"+outputName << endl;

  cout << "[HistoProducer::~HistoProducer] Finished" << endl;

}

