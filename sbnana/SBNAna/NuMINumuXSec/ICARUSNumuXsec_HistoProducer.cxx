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
  double xFMScoreMin = 0.;
  double xFMScoreMax = 150.;
  double dxFMScore = 1.;
  const Binning binsFMScore = Binning::Simple( int( (xFMScoreMax-xFMScoreMin)/dxFMScore ), xFMScoreMin, xFMScoreMax );
  double xFMTimeMin = -60;
  double xFMTimeMax = 60.;
  double dxFMTime = 1.;
  const Binning binsFMTime = Binning::Simple( int( (xFMTimeMax-xFMTimeMin)/dxFMTime ), xFMTimeMin, xFMTimeMax );
  double xNuScoreMin = 0.;
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

  //==== HistAxis
  //====   Slice variables
  const HistAxis axFMScore("FMScore", binsFMScore, varFMScore);
  const HistAxis axFMTime("FMTime", binsFMTime, varFMTime);
  const HistAxis axNuScore("NuScore", binsNuScore, varNuScore);
  const HistAxis axCountSlice("CountSlice", binsOne, varCountSlice);
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
  //====     Proton
  const HistAxis axProtonTruthP("ProtonTruthP", binsMomentum, varProtonTruthP);
  const HistAxis axProtonTruthCosineTheta("ProtonTruthCosineTheta", binsCosineTheta, varProtonTruthCosineTheta);
  //====     Muon+Proton
  const HistAxis axTruthMuonProtonCosineTheta("TruthMuonProtonCosineTheta", binsCosineTheta, varTruthMuonProtonCosineTheta);
  //====   Reco variables
  //====     Muon
  const HistAxis axMuonRecoP("MuonRecoP", binsMomentum, varMuonRecoP);
  const HistAxis axMuonRecoCosineTheta("MuonRecoCosineTheta", binsCosineTheta, varMuonRecoCosineTheta);
  const HistAxis axMuonRecoNuMICosineTheta("MuonRecoNuMICosineTheta", binsCosineTheta, varMuonRecoNuMICosineTheta);
  //====     Proton
  const HistAxis axProtonRecoP("ProtonRecoP", binsMomentum, varProtonRecoP);
  const HistAxis axProtonRecoCosineTheta("ProtonRecoCosineTheta", binsCosineTheta, varProtonRecoCosineTheta);
  //====     Neutrino
  const HistAxis axNeutrinoCombinedEnergy("NeutrinoCombinedEnergy", binsEnergy, varNeutrinoCombinedEnergy);
  const HistAxis axNeutrinoQE("NeutrinoQE", binsEnergy, varNeutrinoQE);

  //==== Spectrum
  //====   Slice variables
  Spectrum *sFMScore = new Spectrum(loader, axFMScore, spillCut, cut);
  Spectrum *sFMTime = new Spectrum(loader, axFMTime, spillCut, cut);
  Spectrum *sNuScore = new Spectrum(loader, axNuScore, spillCut, cut);
  Spectrum *sCountSlice = new Spectrum(loader, axCountSlice, spillCut, cut);
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
  //====     Proton
  Spectrum *sProtonTruthP = new Spectrum(loader, axProtonTruthP, spillCut, cut);
  Spectrum *sProtonTruthCosineTheta = new Spectrum(loader, axProtonTruthCosineTheta, spillCut, cut);
  //====     Muon+Proton
  Spectrum *sTruthMuonProtonCosineTheta = new Spectrum(loader, axTruthMuonProtonCosineTheta, spillCut, cut);
  //====   Reco variables
  //====     Muon
  Spectrum *sMuonRecoP = new Spectrum(loader, axMuonRecoP, spillCut, cut);
  Spectrum *sMuonRecoCosineTheta = new Spectrum(loader, axMuonRecoCosineTheta, spillCut, cut);
  Spectrum *sMuonRecoNuMICosineTheta = new Spectrum(loader, axMuonRecoNuMICosineTheta, spillCut, cut);
  //====     Proton
  Spectrum *sProtonRecoP = new Spectrum(loader, axProtonRecoP, spillCut, cut);
  Spectrum *sProtonRecoCosineTheta = new Spectrum(loader, axProtonRecoCosineTheta, spillCut, cut);
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
    map_cutName_to_vec_Spectrums[currentCutName].push_back(sNeutrinoTruthE);
    map_cutName_to_vec_Spectrums[currentCutName].push_back(sTruthQ2);
    map_cutName_to_vec_Spectrums[currentCutName].push_back(sTruthq0_lab);
    map_cutName_to_vec_Spectrums[currentCutName].push_back(sTruthmodq_lab);
    map_cutName_to_vec_Spectrums[currentCutName].push_back(sMuonTruthP);
    map_cutName_to_vec_Spectrums[currentCutName].push_back(sMuonTruthCosineTheta);
    map_cutName_to_vec_Spectrums[currentCutName].push_back(sMuonTruthNuMICosineTheta);
    map_cutName_to_vec_Spectrums[currentCutName].push_back(sProtonTruthP);
    map_cutName_to_vec_Spectrums[currentCutName].push_back(sProtonTruthCosineTheta);
    map_cutName_to_vec_Spectrums[currentCutName].push_back(sTruthMuonProtonCosineTheta);
    map_cutName_to_vec_Spectrums[currentCutName].push_back(sMuonRecoP);
    map_cutName_to_vec_Spectrums[currentCutName].push_back(sMuonRecoCosineTheta);
    map_cutName_to_vec_Spectrums[currentCutName].push_back(sMuonRecoNuMICosineTheta);
    map_cutName_to_vec_Spectrums[currentCutName].push_back(sProtonRecoP);
    map_cutName_to_vec_Spectrums[currentCutName].push_back(sProtonRecoCosineTheta);
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

void HistoProducer::bookSlice(SpectrumLoader& loader, SpillCut spillCut, Cut cut){

  cout << "[HistoProducer::bookSlice] called" << endl;
  //==== some binnings
  double xFMScoreMin = 0.;
  double xFMScoreMax = 150.;
  double dxFMScore = 1.;
  const Binning binsFMScore = Binning::Simple( int( (xFMScoreMax-xFMScoreMin)/dxFMScore ), xFMScoreMin, xFMScoreMax );
  double xFMTimeMin = -50;
  double xFMTimeMax = 50.;
  double dxFMTime = 1.;
  const Binning binsFMTime = Binning::Simple( int( (xFMTimeMax-xFMTimeMin)/dxFMTime ), xFMTimeMin, xFMTimeMax );
  const Binning binsOne = Binning::Simple( 1, 0., 1.);
  const Binning binsIntCode = Binning::Simple( 17, -2, 15 );

  cout << "[HistoProducer::bookSlice] Delcaring HistAxis" << endl;

  //==== HistAxis
  const HistAxis axCountSlice("CountSlice", binsOne, varCountSlice);
  const HistAxis axFMScore("FMScore", binsFMScore, varFMScore);
  const HistAxis axFMTime("FMTime", binsFMTime, varFMTime);
  const HistAxis axGENIEIntCode("GENIEIntCode", binsIntCode, varGENIEIntCode);

  cout << "[HistoProducer::bookSlice] Delcaring Spectrum" << endl;

  //==== Spectrum
  Spectrum *sCountSlice = new Spectrum(loader, axCountSlice, spillCut, cut);
  Spectrum *sFMScore = new Spectrum(loader, axFMScore, spillCut, cut);
  Spectrum *sFMTime = new Spectrum(loader, axFMTime, spillCut, cut);
  //Spectrum *sFMScore_vs_FMTime = new Spectrum(loader, axFMScore, axFMTime, spillCut, cut);
  Spectrum *sGENIEIntCode = new Spectrum(loader, axGENIEIntCode, spillCut, cut);

  cout << "[HistoProducer::bookTruth] Saving spectrums" << endl;

  map_cutName_to_vec_Spectrums[currentCutName].push_back(sCountSlice);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sFMScore);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sFMTime);
  //map_cutName_to_vec_Spectrums[currentCutName].push_back(sFMScore_vs_FMTime);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sGENIEIntCode);

}

void HistoProducer::bookTruth(SpectrumLoader& loader, SpillCut spillCut, Cut cut){

  cout << "[HistoProducer::bookTruth] called" << endl;
  //==== some binnings
  double xEnergyMin = 0.; // GeV
  double xEnergyMax = 10.;
  double dxEnergy = 0.1;
  const Binning binsEnergy = Binning::Simple( int( (xEnergyMax-xEnergyMin)/dxEnergy ), xEnergyMin, xEnergyMax );
  double xCosineThetaMin = -1;
  double xCosineThetaMax = 1.;
  double dxCosineTheta = 0.05;
  const Binning binsCosineTheta = Binning::Simple( int( (xCosineThetaMax-xCosineThetaMin)/dxCosineTheta ), xCosineThetaMin, xCosineThetaMax );
  const Binning binsN = Binning::Simple(5, 0., 5.);
  const Binning binsLength = Binning::Simple(2000, 0., 2000.);
  const Binning binsChi2B = Binning::Simple(1500, 0., 150);
  const Binning binsReducedChi2B = Binning::Simple(1000, 0., 10);

  cout << "[HistoProducer::bookTruth] Delcaring HistAxis" << endl;

  //==== HistAxis
  const HistAxis axNeutrinoTruthE("NeutrinoTruthE", binsEnergy, varNeutrinoTruthE);
  const HistAxis axTruthQ2("TruthQ2", binsEnergy, varTruthQ2);
  const HistAxis axTruthq0_lab("Truthq0_lab", binsEnergy, varTruthq0_lab);
  const HistAxis axTruthmodq_lab("Truthmodq_lab", binsEnergy, varTruthmodq_lab);
  //==== Muon
  const HistAxis axMuonTruthP("MuonTruthP", binsEnergy, varMuonTruthP);
  const HistAxis axMuonTruthCosineTheta("MuonTruthCosineTheta", binsCosineTheta, varMuonTruthCosineTheta);
  //==== Proton
  const HistAxis axProtonTruthP("ProtonTruthP", binsEnergy, varProtonTruthP);
  const HistAxis axProtonTruthCosineTheta("ProtonTruthCosineTheta", binsCosineTheta, varProtonTruthCosineTheta);
  //==== Muon+Proton
  const HistAxis axTruthMuonProtonCosineTheta("TruthMuonProtonCosineTheta", binsCosineTheta, varTruthMuonProtonCosineTheta);
  //==== Truth match
  const HistAxis axTruthMuonMatchedTrackRangeP("TruthMuonMatchedTrackRangeP", binsEnergy, varTruthMuonMatchedTrackRangeP);
  const HistAxis axTruthMuonMatchedTrackMCSP("TruthMuonMatchedTrackMCSP", binsEnergy, varTruthMuonMatchedTrackMCSP);
  const HistAxis axTruthMuonMatchedTrackCombinedP("TruthMuonMatchedTrackCombinedP", binsEnergy, varTruthMuonMatchedTrackCombinedP);
  const HistAxis axTruthProtonMatchedTrackRangeP("TruthProtonMatchedTrackRangeP", binsEnergy, varTruthProtonMatchedTrackRangeP);
  const HistAxis axTruthProtonMatchedTrackMCSP("TruthProtonMatchedTrackMCSP", binsEnergy, varTruthProtonMatchedTrackMCSP);
  const HistAxis axTruthProtonMatchedTrackCombinedP("TruthProtonMatchedTrackCombinedP", binsEnergy, varTruthProtonMatchedTrackCombinedP);
  const HistAxis axTruthProtonMatchedTrackLength("TruthProtonMatchedTrackLength", binsLength, varTruthProtonMatchedTrackLength);
  const HistAxis axTruthProtonMatchedTrackChi2Proton("TruthProtonMatchedTrackChi2Proton", binsChi2B, varTruthProtonMatchedTrackChi2Proton);
  const HistAxis axTruthProtonMatchedTrackReducedChi2Proton("TruthProtonMatchedTrackReducedChi2Proton", binsReducedChi2B, varTruthProtonMatchedTrackReducedChi2Proton);
  //==== stub
  //const HistAxis axNStub("NStub", binsN, varNStub);

  cout << "[HistoProducer::bookTruth] Delcaring Spectrum" << endl;

  //==== Spectrum
  Spectrum *sNeutrinoTruthE = new Spectrum(loader, axNeutrinoTruthE, spillCut, cut);
  Spectrum *sTruthQ2 = new Spectrum(loader, axTruthQ2, spillCut, cut);
  Spectrum *sTruthq0_lab = new Spectrum(loader, axTruthq0_lab, spillCut, cut);
  Spectrum *sTruthmodq_lab = new Spectrum(loader, axTruthmodq_lab, spillCut, cut);
  Spectrum *sTruth_modq_lab_vs_q0_lab = new Spectrum(loader, axTruthmodq_lab, axTruthq0_lab, spillCut, cut);
  //==== Muon
  Spectrum *sMuonTruthP = new Spectrum(loader, axMuonTruthP, spillCut, cut);
  Spectrum *sMuonTruthCosineTheta = new Spectrum(loader, axMuonTruthCosineTheta, spillCut, cut);
  Spectrum *sMuonTruth_P_vs_CosineTheta = new Spectrum(loader, axMuonTruthP, axMuonTruthCosineTheta, spillCut, cut);
  //==== Proton
  Spectrum *sProtonTruthP = new Spectrum(loader, axProtonTruthP, spillCut, cut);
  Spectrum *sProtonTruthCosineTheta = new Spectrum(loader, axProtonTruthCosineTheta, spillCut, cut);
  //==== Muon+Proton
  Spectrum *sTruthMuonProtonCosineTheta = new Spectrum(loader, axTruthMuonProtonCosineTheta, spillCut, cut);
  //Spectrum *sNStub = new Spectrum(loader, axNStub, spillCut, cut);

  cout << "[HistoProducer::bookTruth] Saving spectrums" << endl;

  map_cutName_to_vec_Spectrums[currentCutName].push_back(sNeutrinoTruthE);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sTruthQ2);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sTruthq0_lab);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sTruthmodq_lab);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sTruth_modq_lab_vs_q0_lab);
  //==== Muon
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sMuonTruthP);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sMuonTruthCosineTheta);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sMuonTruth_P_vs_CosineTheta);
  //==== Proton
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sProtonTruthP);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sProtonTruthCosineTheta);
  //==== Muon+Proton
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sTruthMuonProtonCosineTheta);
  //map_cutName_to_vec_Spectrums[currentCutName].push_back(sNStub);

  for(unsigned int i_syst=0; i_syst<IGENIESysts.size(); i_syst++){
    const ISyst* s = IGENIESysts.at(i_syst);

    //addEnsembleSpectrum(loader, axNeutrinoTruthE, spillCut, cut, currentCutName, vec_UniverseWeightsForEachGENIESource.at(i_syst), s->ShortName());

    addUpDownSystematic(loader, axNeutrinoTruthE, spillCut, cut, currentCutName, s);
    addUpDownSystematic(loader, axMuonTruthP, spillCut, cut, currentCutName, s);
    addUpDownSystematic(loader, axMuonTruthCosineTheta, spillCut, cut, currentCutName, s);
    addUpDownSystematic(loader, axProtonTruthP, spillCut, cut, currentCutName, s);

  }

  for(unsigned int i_syst=0; i_syst<IFluxSysts.size(); i_syst++){
    const ISyst* s = IFluxSysts.at(i_syst);

    addUpDownSystematic(loader, axNeutrinoTruthE, spillCut, cut, currentCutName, s);
    addUpDownSystematic(loader, axMuonTruthP, spillCut, cut, currentCutName, s);
    addUpDownSystematic(loader, axMuonTruthCosineTheta, spillCut, cut, currentCutName, s);
    addUpDownSystematic(loader, axProtonTruthP, spillCut, cut, currentCutName, s);

  }

  if(doTruthMatch){

    //==== Truth match
    Spectrum *sTruthMuonMatchedTrackRangeP = new Spectrum(loader, axTruthMuonMatchedTrackRangeP, spillCut, cut);
    Spectrum *sTruthMuonMatchedTrackMCSP = new Spectrum(loader, axTruthMuonMatchedTrackMCSP, spillCut, cut);
    Spectrum *sTruthMuonMatchedTrackCombinedP = new Spectrum(loader, axTruthMuonMatchedTrackCombinedP, spillCut, cut);
    Spectrum *sTruthProtonMatchedTrackRangeP = new Spectrum(loader, axTruthProtonMatchedTrackRangeP, spillCut, cut);
    Spectrum *sTruthProtonMatchedTrackMCSP = new Spectrum(loader, axTruthProtonMatchedTrackMCSP, spillCut, cut);
    Spectrum *sTruthProtonMatchedTrackCombinedP = new Spectrum(loader, axTruthProtonMatchedTrackCombinedP, spillCut, cut);
    Spectrum *sTruthProtonMatchedTrackLength = new Spectrum(loader, axTruthProtonMatchedTrackLength, spillCut, cut);
    Spectrum *sTruthProtonMatchedTrackChi2Proton = new Spectrum(loader, axTruthProtonMatchedTrackChi2Proton, spillCut, cut);
    Spectrum *sTruthProtonMatchedTrackReducedChi2Proton = new Spectrum(loader, axTruthProtonMatchedTrackReducedChi2Proton, spillCut, cut);
    //Spectrum *sTruthProtonMatchedTrack_Length_vs_Chi2Proton = new Spectrum(loader, axTruthProtonMatchedTrackLength, axTruthProtonMatchedTrackChi2Proton, spillCut, cut);
    //Spectrum *sTruthProtonMatchedTrack_Length_vs_ReducedChi2Proton = new Spectrum(loader, axTruthProtonMatchedTrackLength, axTruthProtonMatchedTrackReducedChi2Proton, spillCut, cut);

    //==== Truth match
    map_cutName_to_vec_Spectrums[currentCutName].push_back(sTruthMuonMatchedTrackRangeP);
    map_cutName_to_vec_Spectrums[currentCutName].push_back(sTruthMuonMatchedTrackMCSP);
    map_cutName_to_vec_Spectrums[currentCutName].push_back(sTruthMuonMatchedTrackCombinedP);
    map_cutName_to_vec_Spectrums[currentCutName].push_back(sTruthProtonMatchedTrackRangeP);
    map_cutName_to_vec_Spectrums[currentCutName].push_back(sTruthProtonMatchedTrackMCSP);
    map_cutName_to_vec_Spectrums[currentCutName].push_back(sTruthProtonMatchedTrackCombinedP);
    map_cutName_to_vec_Spectrums[currentCutName].push_back(sTruthProtonMatchedTrackLength);
    map_cutName_to_vec_Spectrums[currentCutName].push_back(sTruthProtonMatchedTrackChi2Proton);
    map_cutName_to_vec_Spectrums[currentCutName].push_back(sTruthProtonMatchedTrackReducedChi2Proton);
    //map_cutName_to_vec_Spectrums[currentCutName].push_back(sTruthProtonMatchedTrack_Length_vs_Chi2Proton);
    //map_cutName_to_vec_Spectrums[currentCutName].push_back(sTruthProtonMatchedTrack_Length_vs_ReducedChi2Proton);

  }

  cout << "[HistoProducer::bookTruth] Finished" << endl;

}

void HistoProducer::bookRecoMuon(SpectrumLoader& loader, SpillCut spillCut, Cut cut){

  cout << "[HistoProducer::bookRecoMuon] called" << endl;

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
  const Binning binsLength = Binning::Simple(2000, 0., 2000.);
  const Binning binsChi2 = Binning::Simple(200, 0., 200);
  const Binning binsReducedChi2 = Binning::Simple(100, 0., 10);

  cout << "[HistoProducer::bookRecoMuon] Delcaring HistAxis" << endl;

  //==== HistAxis
  const HistAxis axMuonRecoP("MuonRecoP", binsEnergy, varMuonRecoP);
  const HistAxis axMuonCaloPlane0P("MuonCaloPlane0P", binsEnergy, varMuonCaloPlane0P);
  const HistAxis axMuonCaloPlane1P("MuonCaloPlane1P", binsEnergy, varMuonCaloPlane1P);
  const HistAxis axMuonCaloPlane2P("MuonCaloPlane2P", binsEnergy, varMuonCaloPlane2P);
  const HistAxis axMuonBestmatchP("MuonBestmatchP", binsEnergy, varMuonBestmatchP);
  const HistAxis axMuonBestmatchPDG("MuonBestmatchPDG", binsPDG, varMuonBestmatchPDG);
  const HistAxis axMuonPResidual("MuonPResidual", binsEResidual, varMuonPResidual);
  const HistAxis axMuonPResidualFraction("MuonPResidualFraction", binsEResidualFraction, varMuonPResidualFraction);
  const HistAxis axMuonCaloPlane0PResidualFraction("MuonCaloPlane0PResidualFraction", binsEResidualFraction, varMuonCaloPlane0PResidualFraction);
  const HistAxis axMuonCaloPlane1PResidualFraction("MuonCaloPlane1PResidualFraction", binsEResidualFraction, varMuonCaloPlane1PResidualFraction);
  const HistAxis axMuonCaloPlane2PResidualFraction("MuonCaloPlane2PResidualFraction", binsEResidualFraction, varMuonCaloPlane2PResidualFraction);
  const HistAxis axMuonRecoCosineTheta("MuonRecoCosineTheta", binsCosineTheta, varMuonRecoCosineTheta);
  const HistAxis axMuonBestmatchCosineTheta("MuonBestmatchCosineTheta", binsCosineTheta, varMuonBestmatchCosineTheta);
  const HistAxis axMuonCosineThetaResidual("MuonCosineThetaResidual", binsEResidual, varMuonCosineThetaResidual);
  const HistAxis axMuonLength("MuonLength", binsLength, varMuonLength);
  const HistAxis axMuonChi2Muon("MuonChi2Muon", binsChi2, varMuonChi2Muon);
  const HistAxis axMuonReducedChi2Muon("MuonReducedChi2Muon", binsReducedChi2, varMuonReducedChi2Muon);

  cout << "[HistoProducer::bookRecoMuon] Delcaring Spectrum" << endl;

  //==== Spectrum
  Spectrum *sMuonRecoP = new Spectrum(loader, axMuonRecoP, spillCut, cut);
  Spectrum *sMuonBestmatchP = new Spectrum(loader, axMuonBestmatchP, spillCut, cut);
  Spectrum *sMuon_TrueP_vs_RecoP = new Spectrum(loader, axMuonBestmatchP, axMuonRecoP, spillCut, cut);
  Spectrum *sMuon_RecoP_vs_TrueP = new Spectrum(loader, axMuonRecoP, axMuonBestmatchP, spillCut, cut);
  Spectrum *sMuonRecoCosineTheta = new Spectrum(loader, axMuonRecoCosineTheta, spillCut, cut);
  Spectrum *sMuonBestmatchCosineTheta = new Spectrum(loader, axMuonBestmatchCosineTheta, spillCut, cut);
  Spectrum *sMuon_TrueCosineTheta_vs_RecoCosineTheta = new Spectrum(loader, axMuonBestmatchCosineTheta, axMuonRecoCosineTheta, spillCut, cut);
  Spectrum *sMuon_RecoCosineTheta_vs_TrueCosineTheta = new Spectrum(loader, axMuonRecoCosineTheta, axMuonBestmatchCosineTheta, spillCut, cut);
  Spectrum *sMuon_True_P_vs_CosineTheta = new Spectrum(loader, axMuonBestmatchP, axMuonBestmatchCosineTheta, spillCut, cut);
  Spectrum *sMuon_Reco_P_vs_CosineTheta = new Spectrum(loader, axMuonRecoP, axMuonRecoCosineTheta, spillCut, cut);
  Spectrum *sMuonLength = new Spectrum(loader, axMuonLength, spillCut, cut);

/*
  //==== Systematics
  for(unsigned int i=0; i<systGENIEWs.size(); i++){
    Spectrum *sMuonRecoPSyst = new Spectrum(loader, axMuonRecoP, spillCut, cut, kNoShift, systGENIEWs.at(i));
    std::string baseLabel = sMuonRecoPSyst->GetLabels()[0];
    std::string newLabel = baseLabel+"_"+IGENIESysts.at(i)->ShortName();
    sMuonRecoPSyst->SetLabel(0, newLabel);
    map_cutName_to_vec_Spectrums[currentCutName].push_back(sMuonRecoPSyst);
  }
*/

  cout << "[HistoProducer::bookRecoMuon] Saving spectrums" << endl;

  map_cutName_to_vec_Spectrums[currentCutName].push_back(sMuonRecoP);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sMuonBestmatchP);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sMuon_TrueP_vs_RecoP);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sMuon_RecoP_vs_TrueP);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sMuonRecoCosineTheta);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sMuonBestmatchCosineTheta);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sMuon_TrueCosineTheta_vs_RecoCosineTheta);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sMuon_RecoCosineTheta_vs_TrueCosineTheta);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sMuon_True_P_vs_CosineTheta);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sMuon_Reco_P_vs_CosineTheta);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sMuonLength);

  for(unsigned int i_syst=0; i_syst<IGENIESysts.size(); i_syst++){
    const ISyst* s = IGENIESysts.at(i_syst);
/*
    addEnsembleSpectrum(loader, axMuonRecoP, spillCut, cut, currentCutName, vec_UniverseWeightsForEachGENIESource.at(i_syst), s->ShortName());
    addEnsembleSpectrum(loader, axMuonBestmatchP, spillCut, cut, currentCutName, vec_UniverseWeightsForEachGENIESource.at(i_syst), s->ShortName());
    addEnsembleSpectrum(loader, axMuonRecoCosineTheta, spillCut, cut, currentCutName, vec_UniverseWeightsForEachGENIESource.at(i_syst), s->ShortName());
    addEnsembleSpectrum(loader, axMuonBestmatchCosineTheta, spillCut, cut, currentCutName, vec_UniverseWeightsForEachGENIESource.at(i_syst), s->ShortName());

    addEnsembleSpectrum(loader, axMuonRecoP, axMuonBestmatchP, spillCut, cut, currentCutName, vec_UniverseWeightsForEachGENIESource.at(i_syst), s->ShortName());
*/

    addUpDownSystematic(loader, axMuonRecoP, axMuonBestmatchP, spillCut, cut, currentCutName, s);
    addUpDownSystematic(loader, axMuonRecoP, axMuonBestmatchP, spillCut, cut, currentCutName, s);

  }


  for(unsigned int i_syst=0; i_syst<IFluxSysts.size(); i_syst++){
    const ISyst* s = IFluxSysts.at(i_syst);

    addUpDownSystematic(loader, axMuonRecoP, axMuonBestmatchP, spillCut, cut, currentCutName, s);

  }



  if(doPerformanceStudy){

    Spectrum *sMuonCaloPlane0P = new Spectrum(loader, axMuonCaloPlane0P, spillCut, cut);
    Spectrum *sMuonCaloPlane1P = new Spectrum(loader, axMuonCaloPlane1P, spillCut, cut);
    Spectrum *sMuonCaloPlane2P = new Spectrum(loader, axMuonCaloPlane2P, spillCut, cut);
    Spectrum *sMuonBestmatchPDG = new Spectrum(loader, axMuonBestmatchPDG, spillCut, cut);
    Spectrum *sMuonPResidual = new Spectrum(loader, axMuonPResidual, spillCut, cut);
    Spectrum *sMuonPResidualFraction = new Spectrum(loader, axMuonPResidualFraction, spillCut, cut);
    Spectrum *sMuonCaloPlane0PResidualFraction = new Spectrum(loader, axMuonCaloPlane0PResidualFraction, spillCut, cut);
    Spectrum *sMuonCaloPlane1PResidualFraction = new Spectrum(loader, axMuonCaloPlane1PResidualFraction, spillCut, cut);
    Spectrum *sMuonCaloPlane2PResidualFraction = new Spectrum(loader, axMuonCaloPlane2PResidualFraction, spillCut, cut);
    Spectrum *sMuonCosineThetaResidual = new Spectrum(loader, axMuonCosineThetaResidual, spillCut, cut);
    Spectrum *sMuon_TrueCosineTheta_vs_TrueCosineThetaResidual = new Spectrum(loader, axMuonBestmatchCosineTheta, axMuonCosineThetaResidual, spillCut, cut);
    Spectrum *sMuonChi2Muon = new Spectrum(loader, axMuonChi2Muon, spillCut, cut);
    Spectrum *sMuonReducedChi2Muon = new Spectrum(loader, axMuonReducedChi2Muon, spillCut, cut);
    //Spectrum *sMuon_Length_vs_Chi2Muon = new Spectrum(loader, axMuonLength, axMuonChi2Muon, spillCut, cut);
    //Spectrum *sMuon_Length_vs_ReducedChi2Muon = new Spectrum(loader, axMuonLength, axMuonReducedChi2Muon, spillCut, cut);

    map_cutName_to_vec_Spectrums[currentCutName].push_back(sMuonCaloPlane0P);
    map_cutName_to_vec_Spectrums[currentCutName].push_back(sMuonCaloPlane1P);
    map_cutName_to_vec_Spectrums[currentCutName].push_back(sMuonCaloPlane2P);
    map_cutName_to_vec_Spectrums[currentCutName].push_back(sMuonBestmatchPDG);
    map_cutName_to_vec_Spectrums[currentCutName].push_back(sMuonPResidual);
    map_cutName_to_vec_Spectrums[currentCutName].push_back(sMuonPResidualFraction);
    map_cutName_to_vec_Spectrums[currentCutName].push_back(sMuonCaloPlane0PResidualFraction);
    map_cutName_to_vec_Spectrums[currentCutName].push_back(sMuonCaloPlane1PResidualFraction);
    map_cutName_to_vec_Spectrums[currentCutName].push_back(sMuonCaloPlane2PResidualFraction);
    map_cutName_to_vec_Spectrums[currentCutName].push_back(sMuonCosineThetaResidual);
    map_cutName_to_vec_Spectrums[currentCutName].push_back(sMuon_TrueCosineTheta_vs_TrueCosineThetaResidual);
    map_cutName_to_vec_Spectrums[currentCutName].push_back(sMuonChi2Muon);
    map_cutName_to_vec_Spectrums[currentCutName].push_back(sMuonReducedChi2Muon);
    //map_cutName_to_vec_Spectrums[currentCutName].push_back(sMuon_Length_vs_Chi2Muon);
    //map_cutName_to_vec_Spectrums[currentCutName].push_back(sMuon_Length_vs_ReducedChi2Muon);

  }

  cout << "[HistoProducer::bookRecoMuon] Finished" << endl;

}

void HistoProducer::bookMuonPerformance(SpectrumLoader& loader, SpillCut spillCut, Cut cut){

  cout << "[HistoProducer::bookMuonPerformance] called" << endl;

  //==== some binnings
  double xEnergyMin = 0.; // GeV
  double xEnergyMax = 5.;
  double dxEnergy = 0.1;
  const Binning binsEnergy = Binning::Simple( int( (xEnergyMax-xEnergyMin)/dxEnergy ), xEnergyMin, xEnergyMax );
  double xCosineThetaMin = -1;
  double xCosineThetaMax = 1.;
  double dxCosineTheta = 0.05;
  const Binning binsCosineTheta = Binning::Simple( int( (xCosineThetaMax-xCosineThetaMin)/dxCosineTheta ), xCosineThetaMin, xCosineThetaMax );
  const Binning binsEResidual = Binning::Simple(600, -3., 3.);
  const Binning binsEResidualFraction = Binning::Simple(600, -3., 3.);
  const Binning binsPDG = Binning::Simple(100000, -5000., 5000.);
  const Binning binsLength = Binning::Simple(2000, 0., 2000.);
  const Binning binsChi2 = Binning::Simple(200, 0., 200);
  const Binning binsReducedChi2 = Binning::Simple(100, 0., 10);

  cout << "[HistoProducer::bookMuonPerformance] Delcaring HistAxis" << endl;

  //==== HistAxis
  //====   Reco
  const HistAxis axMuonRecoP("MuonRecoP", binsEnergy, varMuonRecoP);
  const HistAxis axMuonCaloPlane0P("MuonCaloPlane0P", binsEnergy, varMuonCaloPlane0P);
  const HistAxis axMuonCaloPlane1P("MuonCaloPlane1P", binsEnergy, varMuonCaloPlane1P);
  const HistAxis axMuonCaloPlane2P("MuonCaloPlane2P", binsEnergy, varMuonCaloPlane2P);
  const HistAxis axMuonRecoCosineTheta("MuonRecoCosineTheta", binsCosineTheta, varMuonRecoCosineTheta);
  const HistAxis axMuonLength("MuonLength", binsLength, varMuonLength);
  const HistAxis axMuonChi2Muon("MuonChi2Muon", binsChi2, varMuonChi2Muon);
  const HistAxis axMuonReducedChi2Muon("MuonReducedChi2Muon", binsReducedChi2, varMuonReducedChi2Muon);
  //====   Truth muon
  const HistAxis axMuonTruthP("MuonTruthP", binsEnergy, varMuonTruthP);
  const HistAxis axMuonTruthCosineTheta("MuonTruthCosineTheta", binsCosineTheta, varMuonTruthCosineTheta);
  //====   Bestmatch truth (may not be a muon)
  const HistAxis axMuonBestmatchPDG("MuonBestmatchPDG", binsPDG, varMuonBestmatchPDG);
  const HistAxis axMuonBestmatchP("MuonBestmatchP", binsEnergy, varMuonBestmatchP);
  const HistAxis axMuonBestmatchCosineTheta("MuonBestmatchCosineTheta", binsCosineTheta, varMuonBestmatchCosineTheta);
  //====   Residual from bestmatch truth
  const HistAxis axMuonPResidual("MuonPResidual", binsEResidual, varMuonPResidual);
  const HistAxis axMuonPResidualFraction("MuonPResidualFraction", binsEResidualFraction, varMuonPResidualFraction);
  const HistAxis axMuonCaloPlane0PResidualFraction("MuonCaloPlane0PResidualFraction", binsEResidualFraction, varMuonCaloPlane0PResidualFraction);
  const HistAxis axMuonCaloPlane1PResidualFraction("MuonCaloPlane1PResidualFraction", binsEResidualFraction, varMuonCaloPlane1PResidualFraction);
  const HistAxis axMuonCaloPlane2PResidualFraction("MuonCaloPlane2PResidualFraction", binsEResidualFraction, varMuonCaloPlane2PResidualFraction);
  const HistAxis axMuonCosineThetaResidual("MuonCosineThetaResidual", binsEResidual, varMuonCosineThetaResidual);

  cout << "[HistoProducer::bookMuonPerformance] Delcaring Spectrum" << endl;

  //==== Spectrum
  //====   Reco
  Spectrum *sMuonRecoP = new Spectrum(loader, axMuonRecoP, spillCut, cut);
  Spectrum *sMuonCaloPlane0P = new Spectrum(loader, axMuonCaloPlane0P, spillCut, cut);
  Spectrum *sMuonCaloPlane1P = new Spectrum(loader, axMuonCaloPlane1P, spillCut, cut);
  Spectrum *sMuonCaloPlane2P = new Spectrum(loader, axMuonCaloPlane2P, spillCut, cut);
  Spectrum *sMuonRecoCosineTheta = new Spectrum(loader, axMuonRecoCosineTheta, spillCut, cut);
  Spectrum *sMuonLength = new Spectrum(loader, axMuonLength, spillCut, cut);
  Spectrum *sMuonChi2Muon = new Spectrum(loader, axMuonChi2Muon, spillCut, cut);
  Spectrum *sMuonReducedChi2Muon = new Spectrum(loader, axMuonReducedChi2Muon, spillCut, cut);
  //====   Truth muon
  Spectrum *sMuonTruthP = new Spectrum(loader, axMuonTruthP, spillCut, cut);
  Spectrum *sMuonTruthCosineTheta = new Spectrum(loader, axMuonTruthCosineTheta, spillCut, cut);
  //====   Bestmatch truth (may not be a muon)
  Spectrum *sMuonBestmatchPDG = new Spectrum(loader, axMuonBestmatchPDG, spillCut, cut);
  Spectrum *sMuonBestmatchP = new Spectrum(loader, axMuonBestmatchP, spillCut, cut);
  Spectrum *sMuonBestmatchCosineTheta = new Spectrum(loader, axMuonBestmatchCosineTheta, spillCut, cut);
  //====   Residual from bestmatch truth
  Spectrum *sMuonPResidual = new Spectrum(loader, axMuonPResidual, spillCut, cut);
  Spectrum *sMuonPResidualFraction = new Spectrum(loader, axMuonPResidualFraction, spillCut, cut);
  Spectrum *sMuonCaloPlane0PResidualFraction = new Spectrum(loader, axMuonCaloPlane0PResidualFraction, spillCut, cut);
  Spectrum *sMuonCaloPlane1PResidualFraction = new Spectrum(loader, axMuonCaloPlane1PResidualFraction, spillCut, cut);
  Spectrum *sMuonCaloPlane2PResidualFraction = new Spectrum(loader, axMuonCaloPlane2PResidualFraction, spillCut, cut);
  Spectrum *sMuonCosineThetaResidual = new Spectrum(loader, axMuonCosineThetaResidual, spillCut, cut);

  cout << "[HistoProducer::bookMuonPerformance] Saving spectrums" << endl;

  //====   Reco
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sMuonRecoP);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sMuonCaloPlane0P);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sMuonCaloPlane1P);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sMuonCaloPlane2P);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sMuonRecoCosineTheta);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sMuonLength);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sMuonChi2Muon);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sMuonReducedChi2Muon);
  //====   Truth muon
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sMuonTruthP);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sMuonTruthCosineTheta);
  //====   Bestmatch truth (may not be a muon)
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sMuonBestmatchPDG);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sMuonBestmatchP);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sMuonBestmatchCosineTheta);
  //====   Residual from bestmatch truth
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sMuonPResidual);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sMuonPResidualFraction);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sMuonCaloPlane0PResidualFraction);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sMuonCaloPlane1PResidualFraction);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sMuonCaloPlane2PResidualFraction);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sMuonCosineThetaResidual);

  cout << "[HistoProducer::bookRecoMuon] Finished" << endl;

}

void HistoProducer::bookRecoProton(SpectrumLoader& loader, SpillCut spillCut, Cut cut){

  cout << "[HistoProducer::bookRecoProton] called" << endl;

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
  const Binning binsLength = Binning::Simple(2000, 0., 200.);
  const Binning binsChi2 = Binning::Simple(1500, 0., 150);
  const Binning binsReducedChi2 = Binning::Simple(1000, 0., 10);

  cout << "[HistoProducer::bookRecoProton] Delcaring HistAxis" << endl;

  //==== HistAxis
  const HistAxis axProtonRecoP("ProtonRecoP", binsEnergy, varProtonRecoP);
  const HistAxis axProtonCaloP("ProtonCaloP", binsEnergy, varProtonCaloP);
  const HistAxis axProtonBestmatchP("ProtonBestmatchP", binsEnergy, varProtonBestmatchP);
  const HistAxis axProtonBestmatchPDG("ProtonBestmatchPDG", binsPDG, varProtonBestmatchPDG);
  const HistAxis axProtonPResidual("ProtonPResidual", binsEResidual, varProtonPResidual);
  const HistAxis axProtonPResidualFraction("ProtonPResidualFraction", binsEResidualFraction, varProtonPResidualFraction);
  const HistAxis axProtonRecoCosineTheta("ProtonRecoCosineTheta", binsCosineTheta, varProtonRecoCosineTheta);
  const HistAxis axProtonBestmatchCosineTheta("ProtonBestmatchCosineTheta", binsCosineTheta, varProtonBestmatchCosineTheta);
  const HistAxis axProtonCosineThetaResidual("ProtonCosineThetaResidual", binsEResidual, varProtonCosineThetaResidual);
  const HistAxis axProtonLength("ProtonLength", binsLength, varProtonLength);
  const HistAxis axProtonChi2Proton("ProtonChi2Proton", binsChi2, varProtonChi2Proton);
  const HistAxis axProtonReducedChi2Proton("ProtonReducedChi2Proton", binsReducedChi2, varProtonReducedChi2Proton);

  cout << "[HistoProducer::bookRecoProton] Delcaring Spectrum" << endl;

  //==== Spectrum
  Spectrum *sProtonRecoP = new Spectrum(loader, axProtonRecoP, spillCut, cut);
  Spectrum *sProtonBestmatchP = new Spectrum(loader, axProtonBestmatchP, spillCut, cut);
  Spectrum *sProton_TrueP_vs_RecoP = new Spectrum(loader, axProtonBestmatchP, axProtonRecoP, spillCut, cut);
  Spectrum *sProton_RecoP_vs_TrueP = new Spectrum(loader, axProtonRecoP, axProtonBestmatchP, spillCut, cut);
  Spectrum *sProtonRecoCosineTheta = new Spectrum(loader, axProtonRecoCosineTheta, spillCut, cut);
  Spectrum *sProtonBestmatchCosineTheta = new Spectrum(loader, axProtonBestmatchCosineTheta, spillCut, cut);
  Spectrum *sProton_TrueCosineTheta_vs_RecoCosineTheta = new Spectrum(loader, axProtonBestmatchCosineTheta, axProtonRecoCosineTheta, spillCut, cut);
  Spectrum *sProton_RecoCosineTheta_vs_TrueCosineTheta = new Spectrum(loader, axProtonRecoCosineTheta, axProtonBestmatchCosineTheta, spillCut, cut);
  Spectrum *sProton_True_P_vs_CosineTheta = new Spectrum(loader, axProtonBestmatchP, axProtonBestmatchCosineTheta, spillCut, cut);
  Spectrum *sProton_Reco_P_vs_CosineTheta = new Spectrum(loader, axProtonRecoP, axProtonRecoCosineTheta, spillCut, cut);
  Spectrum *sProtonLength = new Spectrum(loader, axProtonLength, spillCut, cut);

  cout << "[HistoProducer::bookRecoProton] Saving spectrums" << endl;

  map_cutName_to_vec_Spectrums[currentCutName].push_back(sProtonRecoP);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sProtonBestmatchP);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sProton_TrueP_vs_RecoP);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sProton_RecoP_vs_TrueP);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sProtonRecoCosineTheta);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sProtonBestmatchCosineTheta);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sProton_TrueCosineTheta_vs_RecoCosineTheta);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sProton_RecoCosineTheta_vs_TrueCosineTheta);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sProton_True_P_vs_CosineTheta);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sProton_Reco_P_vs_CosineTheta);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sProtonLength);

/*
  for(unsigned int i_syst=0; i_syst<IGENIESysts.size(); i_syst++){
    const ISyst* s = IGENIESysts.at(i_syst);

    addEnsembleSpectrum(loader, axProtonRecoP, spillCut, cut, currentCutName, vec_UniverseWeightsForEachGENIESource.at(i_syst), s->ShortName());
    addEnsembleSpectrum(loader, axProtonBestmatchP, spillCut, cut, currentCutName, vec_UniverseWeightsForEachGENIESource.at(i_syst), s->ShortName());

  }
*/

  if(doPerformanceStudy){

    Spectrum *sProtonCaloP = new Spectrum(loader, axProtonCaloP, spillCut, cut);
    Spectrum *sProtonBestmatchPDG = new Spectrum(loader, axProtonBestmatchPDG, spillCut, cut);
    Spectrum *sProtonPResidual = new Spectrum(loader, axProtonPResidual, spillCut, cut);
    Spectrum *sProtonPResidualFraction = new Spectrum(loader, axProtonPResidualFraction, spillCut, cut);
    Spectrum *sProtonCosineThetaResidual = new Spectrum(loader, axProtonCosineThetaResidual, spillCut, cut);
    Spectrum *sProtonChi2Proton = new Spectrum(loader, axProtonChi2Proton, spillCut, cut);
    Spectrum *sProtonReducedChi2Proton = new Spectrum(loader, axProtonReducedChi2Proton, spillCut, cut);
    //Spectrum *sProton_Length_vs_Chi2Proton = new Spectrum(loader, axProtonLength, axProtonChi2Proton, spillCut, cut);
    //Spectrum *sProton_Length_vs_ReducedChi2Proton = new Spectrum(loader, axProtonLength, axProtonReducedChi2Proton, spillCut, cut);

    map_cutName_to_vec_Spectrums[currentCutName].push_back(sProtonCaloP);
    map_cutName_to_vec_Spectrums[currentCutName].push_back(sProtonBestmatchPDG);
    map_cutName_to_vec_Spectrums[currentCutName].push_back(sProtonPResidual);
    map_cutName_to_vec_Spectrums[currentCutName].push_back(sProtonPResidualFraction);
    map_cutName_to_vec_Spectrums[currentCutName].push_back(sProtonCosineThetaResidual);
    map_cutName_to_vec_Spectrums[currentCutName].push_back(sProtonChi2Proton);
    map_cutName_to_vec_Spectrums[currentCutName].push_back(sProtonReducedChi2Proton);
    //map_cutName_to_vec_Spectrums[currentCutName].push_back(sProton_Length_vs_Chi2Proton);
    //map_cutName_to_vec_Spectrums[currentCutName].push_back(sProton_Length_vs_ReducedChi2Proton);

  }

  cout << "[HistoProducer::bookRecoProton] Finished" << endl;

}

void HistoProducer::bookRecoNeutrino(SpectrumLoader& loader, SpillCut spillCut, Cut cut){

  cout << "[HistoProducer::bookRecoNeutrino] called" << endl;

  //==== some binnings
  double xEnergyMin = 0.; // GeV
  double xEnergyMax = 10.;
  double dxEnergy = 0.1;
  const Binning binsEnergy = Binning::Simple( int( (xEnergyMax-xEnergyMin)/dxEnergy ), xEnergyMin, xEnergyMax );
  const Binning binsEResidual = Binning::Simple(600, -3., 3.);
  const Binning binsEResidualFraction = Binning::Simple(600, -3., 3.);

  cout << "[HistoProducer::bookRecoNeutrino] Delcaring HistAxis" << endl;

  //==== HistAxis
  const HistAxis axNeutrinoCombinedEnergy("NeutrinoCombinedEnergy", binsEnergy, varNeutrinoCombinedEnergy);
  const HistAxis axNeutrinoTestEnergy("NeutrinoTestEnergy", binsEnergy, varNeutrinoTestEnergy);
  const HistAxis axNeutrinoQE("NeutrinoQE", binsEnergy, varNeutrinoQE);
  const HistAxis axNeutrinoCombinedEnergyResidual("NeutrinoCombinedEnergyResidual", binsEResidual, varNeutrinoCombinedEnergyResidual);
  const HistAxis axNeutrinoQEResidual("NeutrinoQEResidual", binsEResidual, varNeutrinoQEResidual);
  const HistAxis axNeutrinoCombinedEnergyResidualFraction("NeutrinoCombinedEnergyResidualFraction", binsEResidualFraction, varNeutrinoCombinedEnergyResidualFraction);
  const HistAxis axNeutrinoQEResidualFraction("NeutrinoQEResidualFraction", binsEResidualFraction, varNeutrinoQEResidualFraction);
  const HistAxis axNeutrinoTruthE("NeutrinoTruthE", binsEnergy, varNeutrinoTruthE);

  cout << "[HistoProducer::bookRecoNeutrino] Delcaring Spectrum" << endl;

  //==== Spectrum
  Spectrum *sNeutrinoCombinedEnergy = new Spectrum(loader, axNeutrinoCombinedEnergy, spillCut, cut);
  Spectrum *sNeutrinoTestEnergy = new Spectrum(loader, axNeutrinoTestEnergy, spillCut, cut);
  Spectrum *sNeutrinoQE = new Spectrum(loader, axNeutrinoQE, spillCut, cut);
  Spectrum *sNeutrino_RecoCombinedE_vs_TruthE = new Spectrum(loader, axNeutrinoCombinedEnergy, axNeutrinoTruthE, spillCut, cut);

  cout << "[HistoProducer::bookRecoNeutrino] Saving spectrums" << endl;

  map_cutName_to_vec_Spectrums[currentCutName].push_back(sNeutrinoCombinedEnergy);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sNeutrinoTestEnergy);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sNeutrinoQE);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sNeutrino_RecoCombinedE_vs_TruthE);

  for(unsigned int i_syst=0; i_syst<IGENIESysts.size(); i_syst++){
    const ISyst* s = IGENIESysts.at(i_syst);

    addEnsembleSpectrum(loader, axNeutrinoCombinedEnergy, spillCut, cut, currentCutName, vec_UniverseWeightsForEachGENIESource.at(i_syst), s->ShortName());
    addEnsembleSpectrum(loader, axNeutrinoQE, spillCut, cut, currentCutName, vec_UniverseWeightsForEachGENIESource.at(i_syst), s->ShortName());

  }

  for(unsigned int i_syst=0; i_syst<IFluxSysts.size(); i_syst++){
    const ISyst* s = IFluxSysts.at(i_syst);

    addUpDownSystematic(loader, axNeutrinoCombinedEnergy, axNeutrinoTruthE, spillCut, cut, currentCutName, s);

  }

  if(doPerformanceStudy){

    Spectrum *sNeutrinoCombinedEnergyResidual = new Spectrum(loader, axNeutrinoCombinedEnergyResidual, spillCut, cut);
    Spectrum *sNeutrinoQEResidual = new Spectrum(loader, axNeutrinoQEResidual, spillCut, cut);
    Spectrum *sNeutrinoCombinedEnergyResidualFraction = new Spectrum(loader, axNeutrinoCombinedEnergyResidualFraction, spillCut, cut);
    Spectrum *sNeutrinoQEResidualFraction = new Spectrum(loader, axNeutrinoQEResidualFraction, spillCut, cut);

    map_cutName_to_vec_Spectrums[currentCutName].push_back(sNeutrinoCombinedEnergyResidual);
    map_cutName_to_vec_Spectrums[currentCutName].push_back(sNeutrinoQEResidual);
    map_cutName_to_vec_Spectrums[currentCutName].push_back(sNeutrinoCombinedEnergyResidualFraction);
    map_cutName_to_vec_Spectrums[currentCutName].push_back(sNeutrinoQEResidualFraction);

  }

  cout << "[HistoProducer::bookRecoNeutrino] Finished" << endl;

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

  outputfile->cd();
  TH1D *hist_TargetPOT = new TH1D("hist_TargetPOT", "", 1, 0., 1.);
  hist_TargetPOT->SetBinContent(1, TargetPOT);
  hist_TargetPOT->Write();

}

void HistoProducer::setSystematicWeights(){

  cout << "[HistoProducer::setSystematicWeights] Setting GENIE systematics" << endl;
  IGENIESysts = GetSBNGenieWeightSysts();
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

