#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_HistoProducer.h"

HistoProducer::HistoProducer(){

  cout << "[HistoProducer::HistoProducer] called" << endl;
  SystProviderPrefix = "GENIEReWeight_ICARUS_v2";
  TrueTreeFilled = false;
  NNuMIFluxPCA = 12;
  TargetPOT = 6.0e20;
  str_TargetPOT = "6.0e20 POT";
  outputName = "output.root";
  currentCutName = "DefaultCutName";
  vec_cutNames.clear();
  map_cutName_to_vec_Spectrums.clear();
  map_cutName_to_vec_SystEnsembleSpectrumPairs.clear();

  ApplyNuMIPPFXCVWeight = false;

  FillMetaData = true;
  FillSystematics = false;

  FillGENIESyst = true;
  FillFlux = true;

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

void HistoProducer::FillMetadata(SpectrumLoader& loader){

  map_cutName_to_vec_Spectrums["Metadata"].push_back( new Spectrum("TriggerInfoTriggerType", Binning::Simple(11, -1., 10.), loader, TriggerInfoTriggerType, kNoSpillCut) );
  map_cutName_to_vec_Spectrums["Metadata"].push_back( new Spectrum("TriggerInfoSourceType", Binning::Simple(11, -1., 10.), loader, TriggerInfoSourceType, kNoSpillCut) );

}

void HistoProducer::FillFlashMatching(SpectrumLoader& loader, SpillCut spillCut, Cut cut){

  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("FMScore", Binning::Simple(102, -2., 100.), loader, varFMScore, spillCut, cut) );
  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("FMTime", Binning::Simple(122, -62., 60.), loader, varFMTime, spillCut, cut) );

}

void HistoProducer::FillLongestTrack(SpectrumLoader& loader, SpillCut spillCut, Cut cut){

  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("LongestTrackLength", Binning::Simple(100, 0.,500.), loader, LongestTrackLength, spillCut, cut) );
  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("LongestTrackDirectionX", Binning::Simple(40, -1., 1.), loader, LongestTrackDirectionX, spillCut, cut) );
  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("LongestTrackDirectionY", Binning::Simple(40, -1., 1.), loader, LongestTrackDirectionY, spillCut, cut) );
  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("LongestTrackDirectionZ", Binning::Simple(40, -1., 1.), loader, LongestTrackDirectionZ, spillCut, cut) );
  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("LongestTrackDirectionXZ", Binning::Simple(40, -1., 1.), loader, LongestTrackDirectionXZ, spillCut, cut) );
  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("LongestTrackForceDownDirectionX", Binning::Simple(40, -1., 1.), loader, LongestTrackForceDownDirectionX, spillCut, cut) );
  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("LongestTrackForceDownDirectionY", Binning::Simple(40, -1., 1.), loader, LongestTrackForceDownDirectionY, spillCut, cut) );
  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("LongestTrackForceDownDirectionZ", Binning::Simple(40, -1., 1.), loader, LongestTrackForceDownDirectionZ, spillCut, cut) );

}

void HistoProducer::FillSpill(SpectrumLoader& loader, SpillCut spillCut, Cut cut){

  FillCVSpectrum(loader, "CountSpill", spillvarCountSpill, Binning::Simple(1, 0.,1.), spillCut);
  FillCVSpectrum(loader, "TriggerWithinGate", TriggerWithinGate, Binning::Simple(300, -10.,20.), spillCut);
  //FillCVSpectrum(loader, "TriggerWithinGateWithCut", TriggerWithinGate, Binning::Simple(300, -10.,20.), spillCut, cut);

}
void HistoProducer::FillSlice(SpectrumLoader& loader, SpillCut spillCut, Cut cut){

  FillCVSpectrum(loader, "CountSlice", varCountSlice, Binning::Simple(1, 0.,1.), spillCut, cut);
  FillCVSpectrum(loader, "IsClearCosmic", varIsClearCosmic, Binning::Simple(2, 0., 2.), spillCut, cut);
  FillCVSpectrum(loader, "CRLongestTrackDirY", varCRLongestTrackDirY, Binning::Simple(20, -1., 1.), spillCut, cut);
  FillCVSpectrum(loader, "GENIEIntCode", kNuMITrueMode, Binning::Simple(17, -2., 15.), spillCut, cut);

  // stub test
  FillCVSpectrum(loader, "NStubs", NStubs, Binning::Simple(10, 0., 10.), spillCut, cut);

  FillCVSpectrum(loader, "StubCollectionCharges", ICARUSNumuXsec::StubCollectionCharges, Binning::Simple(200., 0., 20.), spillCut, cut);
  FillCVSpectrum(loader, "NNuPFP", NNuPFP, Binning::Simple(30, 0., 30.), spillCut, cut);
  FillCVSpectrum(loader, "NCosmicPFP", NCosmicPFP, Binning::Simple(10, 0., 10.), spillCut, cut);

  FillCVSpectrum(loader, "NNuPFP_vs_NCosmicPFP", NNuPFP, Binning::Simple(30, 0., 30.), NCosmicPFP, Binning::Simple(10, 0., 10.), spillCut, cut);

}

// - NuMu event selection
void HistoProducer::NuMIXSec(SpectrumLoader& loader, SpillCut spillCut, Cut cut){
  FillSlice(loader, spillCut, cut);
  FillSpill(loader, spillCut, cut);
  TruthStudy(loader, spillCut, cut);
  TwoTrackAnalysis(loader, spillCut, cut);
}

void HistoProducer::TruthStudy(SpectrumLoader& loader, SpillCut spillCut, Cut cut){

  using namespace ICARUSNumuXsec::TruthMatch;

  // Scattering

  FillCVandSystSpectrum(loader, "TruthE", kNuMITrueNuE, Binning::Simple(50, 0., 10.), spillCut, cut);

  FillCVSpectrum(loader, "TruthQ2", kNuMITrueQ2, Binning::Simple(60, 0., 6.), spillCut, cut);
  FillCVSpectrum(loader, "Truthq0_lab", kNuMITrueq0, Binning::Simple(60, 0., 6.), spillCut, cut);
  FillCVSpectrum(loader, "Truthmodq_lab", kNuMITrueq3, Binning::Simple(60, 0., 6.), spillCut, cut);
  FillCVSpectrum(loader, "TruthW", kNuMITruew, Binning::Simple(60, 0., 6.), spillCut, cut);

  // Muon

  FillCVSpectrum(loader, "TruthMuonKE", kNuMIMuonTrueKE, Binning::Simple(300, 0., 3.), spillCut, cut);
  FillCVSpectrum(loader, "TruthMuonP", kNuMIMuonTrueP, Binning::Simple(50, 0., 5.), spillCut, cut);
  FillCVSpectrum(loader, "TruthMuonPt", kNuMITrueMuonPt, Binning::Simple(20, 0., 2.), spillCut, cut);
  FillCVSpectrum(loader, "TruthMuonLength", TruthMuonLength, Binning::Simple(500, 0., 500), spillCut, cut);
  FillCVSpectrum(loader, "TruthMuonNuCosineTheta", kNuMIMuonNuCosineTheta, Binning::Simple(40, -1., 1.), spillCut, cut);
  FillCVSpectrum(loader, "TruthMuonMatchedTrackCustomChi2InelasticPionCollection", TruthMuonMatchedTrackCustomChi2InelasticPionCollection, Binning::Simple(80, 0., 80), spillCut, cut);

  FillCVSpectrum(loader, "TruthMuonMichelStartProcess", TruthMuonMichelStartProcess, Binning::Simple(65, 0., 65.), spillCut, cut);
  FillCVSpectrum(loader, "TruthMuonMichelKE", TruthMuonMichelKE, Binning::Simple(100, 0., 0.1), spillCut, cut);
  FillCVSpectrum(loader, "TruthMuonMichelMatchedShowerKE", TruthMuonMichelMatchedShowerKE, Binning::Simple(100, 0., 0.1), spillCut, cut);
  FillCVSpectrum(loader, "TruthMuonMichelMatchedShowerDistanceFromMuonEnd", TruthMuonMichelMatchedShowerDistanceFromMuonEnd, Binning::Simple(300, 0., 30.), spillCut, cut);

  // Proton

  FillCVSpectrum(loader, "TruthProtonKE", kNuMIProtonTrueKE, Binning::Simple(100, 0., 1.0), spillCut, cut);
  FillCVSpectrum(loader, "TruthProtonP", kNuMIProtonTrueP, Binning::Simple(100, 0., 1.0), spillCut, cut);
  FillCVSpectrum(loader, "TruthProtonLength", TruthProtonLength, Binning::Simple(100, 0., 100.), spillCut, cut);
  FillCVSpectrum(loader, "TruthProtonNuCosineTheta", kNuMIProtonNuCosineTheta, Binning::Simple(40, -1., 1.), spillCut, cut);
  FillCVSpectrum(loader, "TruthProtonMatchedTrackCustomChi2InelasticPionCollection", TruthProtonMatchedTrackCustomChi2InelasticPionCollection, Binning::Simple(80, 0., 80), spillCut, cut);

  // Muon+Proton

  FillCVSpectrum(loader, "TruthMuonProtonCosineTheta", kNuMITrueCosThMuP, Binning::Simple(40, -1., 1.), spillCut, cut);
  FillCVSpectrum(loader, "TruthdeltaPT", kNuMITruedeltaPT, Binning::Simple(40, 0., 2.0), spillCut, cut);
  FillCVSpectrum(loader, "TruthdeltaPTx", kNuMITruedeltaPTx, Binning::Simple(60, -1.5, 1.5), spillCut, cut);
  FillCVSpectrum(loader, "TruthdeltaPTy", kNuMITruedeltaPTy, Binning::Simple(60, -1.5, 1.5), spillCut, cut);

  // Charged Pion

  FillCVSpectrum(loader, "TruthChargedPionKE", TruthChargedPionKE, Binning::Simple(100, 0., 1.0), spillCut, cut);
  FillCVSpectrum(loader, "TruthChargedPionLength", TruthChargedPionLength, Binning::Simple(500, 0., 500.), spillCut, cut);
  FillCVSpectrum(loader, "TruthChargedPionHasReco", TruthChargedPionHasReco, Binning::Simple(6, -1., 5.), spillCut, cut);
  FillCVSpectrum(loader, "TruthChargedPionMatchedTrackScore", TruthChargedPionMatchedTrackScore, Binning::Simple(100, 0., 1.), spillCut, cut);
  FillCVSpectrum(loader, "TruthChargedPionMatchedTrackEndProcess", TruthChargedPionMatchedTrackEndProcess, Binning::Simple(65, 0., 65.), spillCut, cut);
  FillCVSpectrum(loader, "TruthChargedPionMatchedTrackCustomChi2MuonCollection", TruthChargedPionMatchedTrackCustomChi2MuonCollection, Binning::Simple(100, 0., 100.), spillCut, cut);
  FillCVSpectrum(loader, "TruthChargedPionMatchedTrackCustomChi2ProtonCollection", TruthChargedPionMatchedTrackCustomChi2ProtonCollection, Binning::Simple(200, 0., 200.), spillCut, cut);
  FillCVSpectrum(loader, "TruthChargedPionMatchedTrackCustomChi2InelasticPionCollection", TruthChargedPionMatchedTrackCustomChi2InelasticPionCollection, Binning::Simple(80, 0., 80), spillCut, cut);
  // ChargedPion michel
  FillCVSpectrum(loader, "TruthChargedPionMichelExistence", TruthChargedPionMichelExistence, Binning::Simple(3, -1., 2.), spillCut, cut);
  FillCVSpectrum(loader, "TruthChargedPionMichelStartProcess", TruthChargedPionMichelStartProcess, Binning::Simple(65, 0., 65.), spillCut, cut);
  FillCVSpectrum(loader, "TruthChargedPionMichelKE", TruthChargedPionMichelKE, Binning::Simple(100, 0., 0.1), spillCut, cut);
  FillCVSpectrum(loader, "TruthChargedPionMichelMatchedShowerKE", TruthChargedPionMichelMatchedShowerKE, Binning::Simple(100, 0., 0.1), spillCut, cut);
  FillCVSpectrum(loader, "TruthChargedPionMichelMatchedShowerDistanceFromMuonEnd", TruthChargedPionMichelMatchedShowerDistanceFromMuonEnd, Binning::Simple(300, 0., 30.), spillCut, cut);
  FillCVSpectrum(loader, "TruthChargedPionMichelHasReco", TruthChargedPionMichelHasReco, Binning::Simple(6, -1., 5.), spillCut, cut);
  FillCVSpectrum(loader, "TruthChargedPionMichelMatchedSlice", ICARUSNumuXsec::TruthMatch::TruthChargedPionMichelMatchedSlice, Binning::Simple(8, -3., 5.), spillCut);

  // Neutral Pion
  FillCVSpectrum(loader, "TruthNeutralPionKE", TruthNeutralPionKE, Binning::Simple(100, 0., 1.0), spillCut, cut);
  FillCVSpectrum(loader, "TruthNeutralPionNMatchedShower", TruthNeutralPionNMatchedShower, Binning::Simple(10, 0., 10.), spillCut, cut);

  // Number of particles
  FillCVSpectrum(loader, "TruthNProton", kNuMITrueNProton, Binning::Simple(10, 0., 10.), spillCut, cut);
  FillCVSpectrum(loader, "TruthNNeutron", kNuMITrueNNeutron, Binning::Simple(10, 0., 10.), spillCut, cut);
  FillCVSpectrum(loader, "TruthNPip", kNuMITrueNpip, Binning::Simple(10, 0., 10.), spillCut, cut);
  FillCVSpectrum(loader, "TruthNPim", kNuMITrueNpim, Binning::Simple(10, 0., 10.), spillCut, cut);
  FillCVSpectrum(loader, "TruthNPi0", kNuMITrueNpi0, Binning::Simple(10, 0., 10.), spillCut, cut);

}

void HistoProducer::TwoTrackAnalysis(SpectrumLoader& loader, SpillCut spillCut, Cut cut){

  using namespace ICARUSNumuXsec::TwoTrack;
  using namespace ICARUSNumuXsec::TruthMatch;

  FillCVSpectrum(loader, "NPrimaryTracks", NPrimaryTracks, Binning::Simple(10, 0., 10.), spillCut, cut);
  FillCVSpectrum(loader, "MuonTrackLength", MuonTrackLength, Binning::Simple(500, 0., 500), spillCut, cut);
  FillCVandSystSpectrum(loader, "MuonTrackP", kNuMIMuonCandidateRecoP, Binning::Simple(50, 0., 5.), spillCut, cut);
  FillCVandSystSpectrum(loader, "MuonTrackPt", kNuMIRecoMuonPt, Binning::Simple(20, 0., 2.), spillCut, cut);
  FillCVSpectrum(loader, "MuonTrackDirX", MuonTrackDirX, Binning::Simple(40, -1., 1), spillCut, cut);
  FillCVSpectrum(loader, "MuonTrackDirY", MuonTrackDirY, Binning::Simple(40, -1., 1), spillCut, cut);
  FillCVSpectrum(loader, "MuonTrackDirZ", MuonTrackDirZ, Binning::Simple(40, -1., 1), spillCut, cut);
  FillCVandSystSpectrum(loader, "MuonTrackNuMICosineTheta", kNuMIRecoCosThBeam, Binning::Simple(40, -1., 1.), spillCut, cut);
  FillCVandSystSpectrum(loader, "MuonTrackNuMIToVtxCosineTheta", kNuMIRecoCosThVtx, Binning::Simple(40, -1., 1.), spillCut, cut);
  // truth matching
  FillCVSpectrum(loader, "MuonTrackTruthPDG", MuonTrackTruthPDG, Binning::Simple(6000, -3000, 3000.), spillCut, cut);
  FillCVSpectrum(loader, "MuonTrackTruthMatchedPrimaryType", MuonTrackTruthMatchedPrimaryType, Binning::Simple(6, -1., 5.), spillCut, cut);
  // true primary
  FillCVSpectrum(loader, "TruthMuonP_vs_MuonTrackP", kNuMIMuonTrueP, Binning::Simple(50, 0., 5.), kNuMIMuonCandidateRecoP, Binning::Simple(50, 0., 5.), spillCut, cut);
  FillCVSpectrum(loader, "TruthMuonPt_vs_MuonTrackPt", kNuMITrueMuonPt, Binning::Simple(20, 0., 2.), kNuMIRecoMuonPt, Binning::Simple(20, 0., 2.), spillCut, cut);
  FillCVSpectrum(loader, "TruthMuonNuCosineTheta_vs_MuonTrackNuMIToVtxCosineTheta", kNuMITrueCosThVtx, Binning::Simple(40, -1., 1.), kNuMIRecoCosThVtx, Binning::Simple(40, -1., 1.), spillCut, cut);

  FillCVSpectrum(loader, "ProtonTrackP", kNuMIProtonCandidateRecoP, Binning::Simple(20, 0., 2.), spillCut, cut);
  FillCVSpectrum(loader, "ProtonTrackDirX", ProtonTrackDirX, Binning::Simple(40, -1., 1), spillCut, cut);
  FillCVSpectrum(loader, "ProtonTrackDirY", ProtonTrackDirY, Binning::Simple(40, -1., 1), spillCut, cut);
  FillCVSpectrum(loader, "ProtonTrackDirZ", ProtonTrackDirZ, Binning::Simple(40, -1., 1), spillCut, cut);
  FillCVSpectrum(loader, "ProtonTrackLength", ProtonTrackLength, Binning::Simple(100, 0., 100.), spillCut, cut);
  FillCVSpectrum(loader, "ProtonTrackNuMICosineTheta", kNuMIProtonRecoCosThBeam, Binning::Simple(40, -1., 1.), spillCut, cut);
  FillCVSpectrum(loader, "ProtonTrackNuMIToVtxCosineTheta", kNuMIProtonTrueCosThBeam, Binning::Simple(40, -1., 1.), spillCut, cut);
  FillCVSpectrum(loader, "ProtonTrackNHitsCollection", ProtonTrackNHitsCollection, Binning::Simple(200, 0., 200.), spillCut, cut);
  FillCVSpectrum(loader, "ProtonTrackChi2MuonCollection", ProtonTrackChi2MuonCollection, Binning::Simple(150, 0., 150.), spillCut, cut);
  FillCVSpectrum(loader, "ProtonTrackChi2ProtonCollection", ProtonTrackChi2ProtonCollection, Binning::Simple(400, 0., 400.), spillCut, cut);
  // truth match
  FillCVSpectrum(loader, "ProtonTrackTruthPDG", ProtonTrackTruthPDG, Binning::Simple(6000, -3000, 3000.), spillCut, cut);
  FillCVSpectrum(loader, "ProtonTrackTruthMatchedPrimaryType", ProtonTrackTruthMatchedPrimaryType, Binning::Simple(6, -1., 5.), spillCut, cut);
  // true primary
  FillCVSpectrum(loader, "TruthProtonP_vs_ProtonTrackP", kNuMIProtonTrueP, Binning::Simple(20, 0., 2.), kNuMIProtonCandidateRecoP, Binning::Simple(20, 0., 2.), spillCut, cut);
  FillCVSpectrum(loader, "TruthProtonNuCosineTheta_vs_ProtonTrackNuMIToVtxCosineTheta", kNuMIProtonTrueCosThVtx, Binning::Simple(40, -1., 1.), kNuMIProtonRecoCosThVtx, Binning::Simple(40, -1., 1.), spillCut, cut);


  FillCVandSystSpectrum(loader, "MuonProtonCosineTheta", kNuMIRecoCosThMuP, Binning::Simple(40, -1., 1.), spillCut, cut);
  FillCVSpectrum(loader, "TruthMuonProtonCosineTheta_vs_MuonProtonCosineTheta", kNuMITrueCosThMuP, Binning::Simple(40, -1., 1.), kNuMIRecoCosThMuP, Binning::Simple(40, -1., 1.), spillCut, cut);

  FillCVandSystSpectrum(loader, "deltaPT", kNuMIRecodeltaPT, Binning::Simple(40, 0., 2.0), spillCut, cut);
  FillCVandSystSpectrum(loader, "deltaPTx", kNuMIRecodeltaPTx, Binning::Simple(60, -3.0, 3.0), spillCut, cut);
  FillCVandSystSpectrum(loader, "deltaPTy", kNuMIRecodeltaPTy, Binning::Simple(60, -3.0, 3.0), spillCut, cut);

  FillCVSpectrum(loader, "TruthdeltaPT_vs_deltaPT", kNuMITruedeltaPT, Binning::Simple(40, 0., 2.0), kNuMIRecodeltaPT, Binning::Simple(40, 0., 2.0), spillCut, cut);
  FillCVSpectrum(loader, "TruthdeltaPTx_vs_deltaPTx", kNuMITruedeltaPTx, Binning::Simple(30, -1.5, 1.5), kNuMIRecodeltaPTx, Binning::Simple(30, -1.5, 1.5), spillCut, cut);
  FillCVSpectrum(loader, "TruthdeltaPTy_vs_deltaPTy", kNuMITruedeltaPTy, Binning::Simple(30, -1.5, 1.5), kNuMIRecodeltaPTy, Binning::Simple(30, -1.5, 1.5), spillCut, cut);

  // pions for side band
  // - stopping charged pion
  FillCVSpectrum(loader, "StoppedChargedPionTrackLength", StoppedChargedPionTrackLength, Binning::Simple(300, 0., 300.), spillCut, cut);
  FillCVandSystSpectrum(loader, "StoppedChargedPionTrackP", StoppedChargedPionTrackP, Binning::Simple(200, 0., 2.), spillCut, cut);
  FillCVSpectrum(loader, "StoppedChargedPionTrackNuMICosineTheta", StoppedChargedPionTrackNuMICosineTheta, Binning::Simple(40, -1., 1.), spillCut, cut);
  FillCVSpectrum(loader, "StoppedChargedPionTrackNuMIToVtxCosineTheta", StoppedChargedPionTrackNuMIToVtxCosineTheta, Binning::Simple(40, -1., 1.), spillCut, cut);
  FillCVSpectrum(loader, "StoppedChargedPionTrackChi2MuonCollection", StoppedChargedPionTrackChi2MuonCollection, Binning::Simple(100, 0., 100.), spillCut, cut);
  FillCVSpectrum(loader, "StoppedChargedPionTrackChi2ProtonCollection", StoppedChargedPionTrackChi2ProtonCollection, Binning::Simple(400, 0., 400.), spillCut, cut);
  // - inelastic pion
  FillCVSpectrum(loader, "InelasticChargedPionTrackLength", InelasticChargedPionTrackLength, Binning::Simple(200, 0., 200.), spillCut, cut);
  FillCVSpectrum(loader, "InelasticChargedPionTrackNuMICosineTheta", InelasticChargedPionTrackNuMICosineTheta, Binning::Simple(40, -1., 1.), spillCut, cut);
  FillCVSpectrum(loader, "InelasticChargedPionTrackNuMIToVtxCosineTheta", InelasticChargedPionTrackNuMIToVtxCosineTheta, Binning::Simple(40, -1., 1.), spillCut, cut);
  FillCVSpectrum(loader, "InelasticChargedPionTrackChi2MIPCollection", InelasticChargedPionTrackChi2MIPCollection, Binning::Simple(100, 0., 10.), spillCut, cut);
  // - neutral pion
  FillCVSpectrum(loader, "NNeutralPionPhotonShower", NNeutralPionPhotonShower, Binning::Simple(10, 0., 10.), spillCut, cut);
  FillCVandSystSpectrum(loader, "NeutralPionPhotonShowerSumEnergy", NeutralPionPhotonShowerSumEnergy, Binning::Simple(20, 0., 2.), spillCut, cut);
  FillCVSpectrum(loader, "NeutralPionPhotonShowerSumInvariantMass", NeutralPionPhotonShowerSumInvariantMass, Binning::Simple(60, 0., 0.3), spillCut, cut);
  FillCVSpectrum(loader, "NeutralPionPhotonShowerConvGaps", NeutralPionPhotonShowerConvGaps, Binning::Simple(120, 0., 30.), spillCut, cut);
  FillCVSpectrum(loader, "NeutralPionPhotonShowerLengths", NeutralPionPhotonShowerLengths, Binning::Simple(50, 0., 50.), spillCut, cut);
  FillCVSpectrum(loader, "NeutralPionPhotonShowerEnergies", NeutralPionPhotonShowerEnergies, Binning::Simple(20, 0., 1.), spillCut, cut);
  FillCVSpectrum(loader, "NeutralPionPhotonShowerdEdxs", NeutralPionPhotonShowerdEdxs, Binning::Simple(50, 0., 10.), spillCut, cut);

}

void HistoProducer::TrackPIDStudy(SpectrumLoader& loader, SpillCut spillCut, Cut cut){

  FillSpill(loader, spillCut, cut);
  FillSlice(loader, spillCut, cut);

  using namespace ICARUSNumuXsec::TwoTrack;

  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "RelaxedMuonTrackTrackScore",
      Binning::Simple(20, 0., 1.), loader,
      ICARUSNumuXsec::TwoTrack::Aux::RelaxedMuonTrackTrackScore,
      spillCut, cut
    )
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "RelaxedMuonTrackLength",
      Binning::Simple(500, 0., 500), loader,
      ICARUSNumuXsec::TwoTrack::Aux::RelaxedMuonTrackLength,
      spillCut, cut
    )
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "RelaxedMuonTrackP",
      Binning::Simple(50, 0., 5), loader,
      ICARUSNumuXsec::TwoTrack::Aux::RelaxedMuonTrackP,
      spillCut, cut
    )
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "RelaxedMuonTrackNuMICosineTheta",
      Binning::Simple(40, -1., 1), loader,
      ICARUSNumuXsec::TwoTrack::Aux::RelaxedMuonTrackNuMICosineTheta,
      spillCut, cut
    )
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "RelaxedMuonTrackPosAbsX",
      Binning::Simple(400, 0., 400.), loader,
      ICARUSNumuXsec::TwoTrack::Aux::RelaxedMuonTrackPosAbsX,
      spillCut, cut
    )
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "RelaxedMuonTrackDirX",
      Binning::Simple(40, -1., 1), loader,
      ICARUSNumuXsec::TwoTrack::Aux::RelaxedMuonTrackDirX,
      spillCut, cut
    )
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "RelaxedMuonTrackDirY",
      Binning::Simple(40, -1., 1), loader,
      ICARUSNumuXsec::TwoTrack::Aux::RelaxedMuonTrackDirY,
      spillCut, cut
    )
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "RelaxedMuonTrackDirZ",
      Binning::Simple(40, -1., 1), loader,
      ICARUSNumuXsec::TwoTrack::Aux::RelaxedMuonTrackDirZ,
      spillCut, cut
    )
  );

  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "RelaxedMuonTrackChi2Muon",
      Binning::Simple(150, 0., 150.), loader,
      ICARUSNumuXsec::TwoTrack::Aux::RelaxedMuonTrackChi2Muon,
      spillCut, cut
    )
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "RelaxedMuonTrackChi2Proton",
      Binning::Simple(400, 0., 400.), loader,
      ICARUSNumuXsec::TwoTrack::Aux::RelaxedMuonTrackChi2Proton,
      spillCut, cut
    )
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "RelaxedMuonTrackChi2MuonCollection",
      Binning::Simple(150, 0., 150.), loader,
      ICARUSNumuXsec::TwoTrack::Aux::RelaxedMuonTrackChi2MuonCollection,
      spillCut, cut
    )
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "RelaxedMuonTrackChi2ProtonCollection",
      Binning::Simple(400, 0., 400.), loader,
      ICARUSNumuXsec::TwoTrack::Aux::RelaxedMuonTrackChi2ProtonCollection,
      spillCut, cut
    )
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "RelaxedMuonTrackCustomChi2MuonCollection",
      Binning::Simple(150, 0., 150.), loader,
      ICARUSNumuXsec::TwoTrack::Aux::RelaxedMuonTrackCustomChi2MuonCollection,
      spillCut, cut
    )
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "RelaxedMuonTrackCustomChi2ProtonCollection",
      Binning::Simple(400, 0., 400.), loader,
      ICARUSNumuXsec::TwoTrack::Aux::RelaxedMuonTrackCustomChi2ProtonCollection,
      spillCut, cut
    )
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "RelaxedMuonTrackTruthStartProcess",
      Binning::Simple(65, 0., 65.), loader,
      ICARUSNumuXsec::TwoTrack::Aux::RelaxedMuonTrackTruthStartProcess,
      spillCut, cut
    )
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "RelaxedMuonTrackTruthEndProcess",
      Binning::Simple(65, 0., 65.), loader,
      ICARUSNumuXsec::TwoTrack::Aux::RelaxedMuonTrackTruthEndProcess,
      spillCut, cut
    )
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "RelaxedMuonTrackCollectionRR_vs_RelaxedMuonTrackCollectiondEdX", loader,
      Binning::Simple(50, 0., 25.), ICARUSNumuXsec::TwoTrack::Aux::RelaxedMuonTrackCollectionRR,
      Binning::Simple(50, 0., 5.), ICARUSNumuXsec::TwoTrack::Aux::RelaxedMuonTrackCollectiondEdX,
      spillCut, cut
    )
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "RelaxedMuonTrackTruthLength_vs_RelaxedMuonTrackLength", loader,
      Binning::Simple(50, 0., 500.), ICARUSNumuXsec::TwoTrack::Aux::RelaxedMuonTrackTruthLength,
      Binning::Simple(50, 0., 500.), ICARUSNumuXsec::TwoTrack::Aux::RelaxedMuonTrackLength,
      spillCut, cut
    )
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "RelaxedMuonTrackTruthP_vs_RelaxedMuonTrackP", loader,
      Binning::Simple(50, 0., 5.), ICARUSNumuXsec::TwoTrack::Aux::RelaxedMuonTrackTruthP,
      Binning::Simple(50, 0., 5.), ICARUSNumuXsec::TwoTrack::Aux::RelaxedMuonTrackP,
      spillCut, cut
    )
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "RelaxedMuonTrackTruthP_vs_RelaxedMuonTrackTruthPResFrac", loader,
      Binning::Simple(50, 0., 5.), ICARUSNumuXsec::TwoTrack::Aux::RelaxedMuonTrackTruthP,
      Binning::Simple(100, -1, 1.), ICARUSNumuXsec::TwoTrack::Aux::RelaxedMuonTrackTruthPResFrac,
      spillCut, cut
    )
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "RelaxedMuonTrackTruthOneOverP_vs_RelaxedMuonTrackTruthOneOverPResFrac", loader,
      Binning::Simple(50, 0., 5.), ICARUSNumuXsec::TwoTrack::Aux::RelaxedMuonTrackTruthOneOverP,
      Binning::Simple(100, -1, 1.), ICARUSNumuXsec::TwoTrack::Aux::RelaxedMuonTrackTruthOneOverPResFrac,
      spillCut, cut
    )
  );

  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "RelaxedMuonTrackTruthNuMICosineTheta_vs_RelaxedMuonTrackNuMICosineTheta", loader,
      Binning::Simple(40, -1., 1), ICARUSNumuXsec::TwoTrack::Aux::RelaxedMuonTrackTruthNuMICosineTheta,
      Binning::Simple(40, -1., 1), ICARUSNumuXsec::TwoTrack::Aux::RelaxedMuonTrackNuMICosineTheta,
      spillCut, cut
    )
  );


  // proton

  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "RelaxedProtonTrackLength",
      Binning::Simple(100, 0., 100), loader,
      ICARUSNumuXsec::TwoTrack::Aux::RelaxedProtonTrackLength,
      spillCut, cut
    )
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "RelaxedProtonTrackP",
      Binning::Simple(50, 0., 5), loader,
      ICARUSNumuXsec::TwoTrack::Aux::RelaxedProtonTrackP,
      spillCut, cut
    )
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "RelaxedProtonTrackNuMICosineTheta",
      Binning::Simple(40, -1., 1), loader,
      ICARUSNumuXsec::TwoTrack::Aux::RelaxedProtonTrackNuMICosineTheta,
      spillCut, cut
    )
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "RelaxedProtonTrackDirX",
      Binning::Simple(40, -1., 1), loader,
      ICARUSNumuXsec::TwoTrack::Aux::RelaxedProtonTrackDirX,
      spillCut, cut
    )
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "RelaxedProtonTrackDirY",
      Binning::Simple(40, -1., 1), loader,
      ICARUSNumuXsec::TwoTrack::Aux::RelaxedProtonTrackDirY,
      spillCut, cut
    )
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "RelaxedProtonTrackDirZ",
      Binning::Simple(40, -1., 1), loader,
      ICARUSNumuXsec::TwoTrack::Aux::RelaxedProtonTrackDirZ,
      spillCut, cut
    )
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "RelaxedProtonTrackTrackScore",
      Binning::Simple(20, 0., 1.), loader,
      ICARUSNumuXsec::TwoTrack::Aux::RelaxedProtonTrackTrackScore,
      spillCut, cut
    )
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "RelaxedProtonTrackChi2Muon",
      Binning::Simple(150, 0., 150.), loader,
      ICARUSNumuXsec::TwoTrack::Aux::RelaxedProtonTrackChi2Muon,
      spillCut, cut
    )
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "RelaxedProtonTrackChi2Proton",
      Binning::Simple(400, 0., 400.), loader,
      ICARUSNumuXsec::TwoTrack::Aux::RelaxedProtonTrackChi2Proton,
      spillCut, cut
    )
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "RelaxedProtonTrackChi2MuonCollection",
      Binning::Simple(150, 0., 150.), loader,
      ICARUSNumuXsec::TwoTrack::Aux::RelaxedProtonTrackChi2MuonCollection,
      spillCut, cut
    )
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "RelaxedProtonTrackChi2ProtonCollection",
      Binning::Simple(400, 0., 400.), loader,
      ICARUSNumuXsec::TwoTrack::Aux::RelaxedProtonTrackChi2ProtonCollection,
      spillCut, cut
    )
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "RelaxedProtonTrackCustomChi2MuonCollection",
      Binning::Simple(150, 0., 150.), loader,
      ICARUSNumuXsec::TwoTrack::Aux::RelaxedProtonTrackCustomChi2MuonCollection,
      spillCut, cut
    )
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "RelaxedProtonTrackCustomChi2ProtonCollection",
      Binning::Simple(400, 0., 400.), loader,
      ICARUSNumuXsec::TwoTrack::Aux::RelaxedProtonTrackCustomChi2ProtonCollection,
      spillCut, cut
    )
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "RelaxedProtonTrackNHitCollection",
      Binning::Simple(400, 0., 400.), loader,
      ICARUSNumuXsec::TwoTrack::Aux::RelaxedProtonTrackNHitCollection,
      spillCut, cut
    )
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "RelaxedProtonTrackCollectionRR_vs_RelaxedProtonTrackCollectiondEdX", loader,
      Binning::Simple(50, 0., 25.), ICARUSNumuXsec::TwoTrack::Aux::RelaxedProtonTrackCollectionRR,
      Binning::Simple(60, 0., 30.), ICARUSNumuXsec::TwoTrack::Aux::RelaxedProtonTrackCollectiondEdX,
      spillCut, cut
    )
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "RelaxedProtonTrackTruthStartProcess",
      Binning::Simple(65, 0., 65.), loader,
      ICARUSNumuXsec::TwoTrack::Aux::RelaxedProtonTrackTruthStartProcess,
      spillCut, cut
    )
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "RelaxedProtonTrackTruthEndProcess",
      Binning::Simple(65, 0., 65.), loader,
      ICARUSNumuXsec::TwoTrack::Aux::RelaxedProtonTrackTruthEndProcess,
      spillCut, cut
    )
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "RelaxedProtonTrackLength_vs_RelaxedProtonTrackChi2MuonCollection", loader,
      Binning::Simple(25, 0., 50.), ICARUSNumuXsec::TwoTrack::Aux::RelaxedProtonTrackLength,
      Binning::Simple(70, 0., 70.), ICARUSNumuXsec::TwoTrack::Aux::RelaxedProtonTrackChi2MuonCollection,
      spillCut, cut
    )
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "RelaxedProtonTrackTruthLength_vs_RelaxedProtonTrackLength", loader,
      Binning::Simple(100, 0., 100.), ICARUSNumuXsec::TwoTrack::Aux::RelaxedProtonTrackTruthLength,
      Binning::Simple(100, 0., 100.), ICARUSNumuXsec::TwoTrack::Aux::RelaxedProtonTrackLength,
      spillCut, cut
    )
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "RelaxedProtonTrackTruthP_vs_RelaxedProtonTrackP", loader,
      Binning::Simple(50, 0., 1.), ICARUSNumuXsec::TwoTrack::Aux::RelaxedProtonTrackTruthP,
      Binning::Simple(50, 0., 1.), ICARUSNumuXsec::TwoTrack::Aux::RelaxedProtonTrackP,
      spillCut, cut
    )
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "RelaxedProtonTrackTruthP_vs_RelaxedProtonTrackTruthPResFrac", loader,
      Binning::Simple(50, 0., 2.), ICARUSNumuXsec::TwoTrack::Aux::RelaxedProtonTrackTruthP,
      Binning::Simple(100, -1, 1.), ICARUSNumuXsec::TwoTrack::Aux::RelaxedProtonTrackTruthPResFrac,
      spillCut, cut
    )
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "RelaxedProtonTrackTruthNuMICosineTheta_vs_RelaxedProtonTrackNuMICosineTheta", loader,
      Binning::Simple(40, -1., 1), ICARUSNumuXsec::TwoTrack::Aux::RelaxedProtonTrackTruthNuMICosineTheta,
      Binning::Simple(40, -1., 1), ICARUSNumuXsec::TwoTrack::Aux::RelaxedProtonTrackNuMICosineTheta,
      spillCut, cut
    )
  );

  // ChargedPion

  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "RelaxedChargedPionTrackTrackScore",
      Binning::Simple(20, 0., 1.), loader,
      ICARUSNumuXsec::TwoTrack::Aux::RelaxedChargedPionTrackTrackScore,
      spillCut, cut
    )
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "RelaxedChargedPionTrackLength",
      Binning::Simple(500, 0., 500), loader,
      ICARUSNumuXsec::TwoTrack::Aux::RelaxedChargedPionTrackLength,
      spillCut, cut
    )
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "RelaxedChargedPionTrackP",
      Binning::Simple(20, 0., 2), loader,
      ICARUSNumuXsec::TwoTrack::Aux::RelaxedChargedPionTrackP,
      spillCut, cut
    )
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "RelaxedChargedPionTrackNuMICosineTheta",
      Binning::Simple(40, -1., 1), loader,
      ICARUSNumuXsec::TwoTrack::Aux::RelaxedChargedPionTrackNuMICosineTheta,
      spillCut, cut
    )
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "RelaxedChargedPionTrackCustomChi2MuonCollection",
      Binning::Simple(150, 0., 150.), loader,
      ICARUSNumuXsec::TwoTrack::Aux::RelaxedChargedPionTrackCustomChi2MuonCollection,
      spillCut, cut
    )
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "RelaxedChargedPionTrackCustomChi2ProtonCollection",
      Binning::Simple(400, 0., 400.), loader,
      ICARUSNumuXsec::TwoTrack::Aux::RelaxedChargedPionTrackCustomChi2ProtonCollection,
      spillCut, cut
    )
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "RelaxedChargedPionTrackCustomChi2InelasticPionCollection",
      Binning::Simple(80, 0., 80.), loader,
      ICARUSNumuXsec::TwoTrack::Aux::RelaxedChargedPionTrackCustomChi2InelasticPionCollection,
      spillCut, cut
    )
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "RelaxedChargedPionTrackFromVertex",
      Binning::Simple(100, 0., 10.), loader,
      ICARUSNumuXsec::TwoTrack::Aux::RelaxedChargedPionTrackFromVertex,
      spillCut, cut
    )
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "RelaxedChargedPionTrackNDaughter",
      Binning::Simple(20, 0., 20.), loader,
      ICARUSNumuXsec::TwoTrack::Aux::RelaxedChargedPionTrackNDaughter,
      spillCut, cut
    )
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "RelaxedChargedPionTrackTruthStartProcess",
      Binning::Simple(65, 0., 65.), loader,
      ICARUSNumuXsec::TwoTrack::Aux::RelaxedChargedPionTrackTruthStartProcess,
      spillCut, cut
    )
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "RelaxedChargedPionTrackTruthPDG",
      Binning::Simple(6000, -3000, 3000.), loader,
      ICARUSNumuXsec::TwoTrack::Aux::RelaxedChargedPionTrackTruthPDG,
      spillCut, cut
    )
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "RelaxedChargedPionTrackTruthEndProcess",
      Binning::Simple(65, 0., 65.), loader,
      ICARUSNumuXsec::TwoTrack::Aux::RelaxedChargedPionTrackTruthEndProcess,
      spillCut, cut
    )
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "RelaxedChargedPionTrackCollectionRR_vs_RelaxedChargedPionTrackCollectiondEdX", loader,
      Binning::Simple(50, 0., 25.), ICARUSNumuXsec::TwoTrack::Aux::RelaxedChargedPionTrackCollectionRR,
      Binning::Simple(50, 0., 5.), ICARUSNumuXsec::TwoTrack::Aux::RelaxedChargedPionTrackCollectiondEdX,
      spillCut, cut
    )
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "RelaxedChargedPionTrackCollectionFrontRR_vs_RelaxedChargedPionTrackCollectiondEdX", loader,
      Binning::Simple(50, 0., 25.), ICARUSNumuXsec::TwoTrack::Aux::RelaxedChargedPionTrackCollectionFrontRR,
      Binning::Simple(50, 0., 5.), ICARUSNumuXsec::TwoTrack::Aux::RelaxedChargedPionTrackCollectiondEdX,
      spillCut, cut
    )
  );

  // tagged muon
  FillCVSpectrum(loader, "MuonTrackLength", MuonTrackLength, Binning::Simple(500, 0., 500), spillCut, cut);

}

// - 230418_StubStudy
void HistoProducer::StubStudy(SpectrumLoader& loader, SpillCut spillCut, Cut cut){

  FillCVandSystSpectrum(loader, "NStubs", NStubs, Binning::Simple(10, 0., 10.), spillCut, cut);

  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum("StubCollectionCharges", Binning::Simple(200., 0., 20.), loader, ICARUSNumuXsec::StubCollectionCharges, spillCut, cut)
  );

  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "MuonTrackLength",
      Binning::Simple(500, 0., 500), loader,
      ICARUSNumuXsec::TwoTrack::MuonTrackLength,
      spillCut, cut
    )
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "MuonTrackP",
      Binning::Simple(50, 0., 5), loader,
      kNuMIMuonCandidateRecoP,
      spillCut, cut
    )
  );


}

// - 230517_TriggerEffStudy
void HistoProducer::TriggerEffStudy(SpectrumLoader& loader, SpillCut spillCut, Cut cut){

  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum("TriggerWithinGate", Binning::Simple(300, -10.,20.), loader, TriggerWithinGate, spillCut)
  );

  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum("NuMuSliceLongestTrackLenForTriggerEff", Binning::Simple(50, 0., 500.), loader, NuMuSliceLongestTrackLenForTriggerEff, spillCut)
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum("InTimeCosmicSliceLongestTrackLenForTriggerEff", Binning::Simple(50, 0., 500.), loader, InTimeCosmicSliceLongestTrackLenForTriggerEff, spillCut)
  );

}

// - 230524_MichelStudy
void HistoProducer::MichelStudy(SpectrumLoader& loader, SpillCut spillCut, Cut cut){

  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("TruthMuonMichelIndex", Binning::Simple(1, 0., 1.0), loader, ICARUSNumuXsec::TruthMatch::TruthMuonMichelIndex, spillCut, cut) );

/*
  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("TruthChargedPionWithMichel_TruthChargedPionKEs", Binning::Simple(100, 0., 1.), loader, ICARUSNumuXsec::MichelStudy::TruthChargedPionWithMichel_TruthChargedPionKEs, spillCut) );

  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("TruthChargedPionWithMichel_TruthMichelElectronKEs", Binning::Simple(100, 0., 0.1), loader, ICARUSNumuXsec::MichelStudy::TruthChargedPionWithMichel_TruthMichelElectronKEs, spillCut) );

  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("TruthChargedPionWithMichel_TruthMichelElectronMatchedRecoShowerEnergies", Binning::Simple(100, 0., 0.1), loader, ICARUSNumuXsec::MichelStudy::TruthChargedPionWithMichel_TruthMichelElectronMatchedRecoShowerEnergies, spillCut) );

  map_cutName_to_vec_Spectrums[currentCutName].push_back( 
    new Spectrum(
      "TruthChargedPionWithMichel_TruthMichelElectronKEs_vs_TruthMichelElectronMatchedRecoShowerEnergies", loader,
      Binning::Simple(100, 0., 0.1), ICARUSNumuXsec::MichelStudy::TruthChargedPionWithMichel_TruthMichelElectronKEs,
      Binning::Simple(100, 0., 0.1), ICARUSNumuXsec::MichelStudy::TruthChargedPionWithMichel_TruthMichelElectronMatchedRecoShowerEnergies,
      spillCut
    )
  );
  */

  //map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("MichelTest", Binning::Simple(1, 0., 1.), loader, ICARUSNumuXsec::MichelStudy::MichelTest, spillCut) );

/*

  using namespace ICARUSNumuXsec::TruthMatch;

  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("TruthChargedPionLength", Binning::Simple(500, 0., 500.), loader, TruthChargedPionLength, spillCut, cut) );
  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("TruthChargedPionKE", Binning::Simple(100, 0., 1.0), loader, TruthChargedPionKE, spillCut, cut) );

  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("TruthChargedPionMatchedTrackContainedness", Binning::Simple(2, 0., 2.), loader, TruthChargedPionMatchedTrackContainedness, spillCut, cut) );

  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("TruthChargedPionMatchedTrackEndProcess", Binning::Simple(65, 0., 65.), loader, TruthChargedPionMatchedTrackEndProcess, spillCut, cut) );

  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("TruthChargedPionNDaughters", Binning::Simple(10, 0., 10), loader, TruthChargedPionNDaughters, spillCut, cut) );

  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("TruthChargedPionDaughterPDGs", Binning::Simple(6000, -3000, 3000.), loader, TruthChargedPionDaughterPDGs, spillCut, cut) );

  // 2D
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "TruthChargedPionKE_vs_TruthChargedPionMatchedTrackEndProcess", loader,
      Binning::Simple(100, 0., 1.0), TruthChargedPionKE,
      Binning::Simple(65, 0., 65.), TruthChargedPionMatchedTrackEndProcess,
      spillCut, cut
    )
  );

  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "TruthChargedPionLength_vs_TruthChargedPionMatchedTrackEndProcess", loader,
      Binning::Simple(500, 0., 500.), TruthChargedPionLength,
      Binning::Simple(65, 0., 65.), TruthChargedPionMatchedTrackEndProcess,
      spillCut, cut
    )
  );
*/
}

// - 230814_MakeTree
void HistoProducer::MakeTrueTree(SpectrumLoader& loader){
  if(TrueTreeFilled) return;

  using namespace ICARUSNumuXsec::TwoTrack;
  using namespace ICARUSNumuXsec::TruthMatch;

  map_cutName_to_vec_Trees[currentCutName].push_back(

    new ana::Tree(
      "trueEvents", GetNuMITrueTreeLabels(),
      loader,
      GetNuMITrueTreeVars(),
      kNoSpillCut,
      false
    )

  );

  TrueTreeFilled = true;
}
void HistoProducer::MakeTree(SpectrumLoader& loader, SpillCut spillCut, Cut cut){

  MakeTrueTree(loader);

  using namespace ICARUSNumuXsec::TwoTrack;
  using namespace ICARUSNumuXsec::TruthMatch;

  map_cutName_to_vec_Trees[currentCutName].push_back(

    new ana::Tree(
      "selectedEvents", GetNuMIRecoTreeLabels(),
      loader,
      GetNuMIRecoTreeVars(),
      spillCut, cut,
      kNoShift, true, true
    )

  );

  std::vector<std::string> this_PsetNames;
  std::vector<const ISyst*> this_SigmaSysts;

  for(const std::string& name: ICARUSNumuXsec::GetGENIEMultisigmaKnobNames()){
    if(name=="FormZone"){
     continue;
    }
    std::string psetname = SystProviderPrefix+"_multisigma_"+name;
    this_PsetNames.push_back( name );
    this_SigmaSysts.push_back( new SBNWeightSyst(psetname) );
  }
  std::vector<const ISyst*> this_IFluxSysts = GetAllNuMIFluxSysts(NNuMIFluxPCA);
  for(unsigned int i=0; i<this_IFluxSysts.size(); i++){
    this_PsetNames.push_back( this_IFluxSysts.at(i)->ShortName() );
    this_SigmaSysts.push_back( this_IFluxSysts.at(i) );
  }


  map_cutName_to_vec_NSigmasTrees[currentCutName].push_back(

    new ana::NSigmasTree(
      (currentCutName+"_Sigma").Data(),
      this_PsetNames,
      loader,
      this_SigmaSysts,
      spillCut, cut,
      kNoShift, 3, false, false
    )

  );

}

void HistoProducer::Test(SpectrumLoader& loader, SpillCut spillCut, Cut cut){

  //map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("spillvarTest", Binning::Simple(1, 0., 1.), loader, spillvarTest, spillCut) );

  //map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("spillvarTest", Binning::Simple(1, 0., 1.), loader, ICARUSNumuXsec::TwoTrack::Aux::TestSpillVar, spillCut) );

  using namespace ICARUSNumuXsec::TwoTrack;
  using namespace ICARUSNumuXsec::TruthMatch;

  //map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("spillvarTest", Binning::Simple(1, 0., 1.), loader, TestVar, spillCut) );

  //map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("Test", Binning::Simple(1, 0., 1.), loader, ICARUSNumuXsec::TruthMatch::TruthChargedPionMichelIndex, spillCut, cut) );

  //map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("Test", Binning::Simple(8, -3., 5.), loader, ICARUSNumuXsec::TruthMatch::TruthChargedPionMichelMatchedSlice, spillCut) );

  //FillCVandSystSpectrum(loader, "TruthE", kNuMITrueNuE, Binning::Simple(50, 0., 10.), spillCut, cut);

/*
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "TruthProtonP_vs_TruthProtonLength", loader,
      Binning::Simple(100, 0., 1.0), ICARUSNumuXsec::TruthMatch::TruthProtonP,
      Binning::Simple(100, 0., 100.), ICARUSNumuXsec::TruthMatch::TruthProtonLength,
      spillCut, cut
    )
  );
*/

}

void HistoProducer::saveHistograms(){

  outputfile->cd();

  cout << "[HistoProducer::saveHistograms] Number of cuts = " << vec_cutNames.size() << endl;
  const unsigned int nCutName = vec_cutNames.size();
  double inputSamplePOT(-1.);
  double inputSampleLiveTime(-1.);
  if(nCutName>0){

    for(unsigned int ic=0; ic<nCutName; ic++){

      const TString cutName = vec_cutNames.at(ic);
      const TString dirName = cutName;
      TDirectory *dir = outputfile->GetDirectory(dirName);
      if(!dir){
        outputfile->mkdir(dirName);
        dir = outputfile->GetDirectory(dirName);
      }
      outputfile->cd(dirName);

      vector<Spectrum *> vec_Spectrums = map_cutName_to_vec_Spectrums[cutName];
      vector< pair<TString, EnsembleSpectrum *> > vec_SystEnsembleSpectrumPairs = map_cutName_to_vec_SystEnsembleSpectrumPairs[cutName];
      vector< pair<TString, Spectrum *> > vec_SystSpectrumPairs = map_cutName_to_vec_SystSpectrumPairs[cutName];
      vector<ana::Tree *> vec_Trees = map_cutName_to_vec_Trees[cutName];
      vector<ana::NSigmasTree *> vec_NSigmasTrees = map_cutName_to_vec_NSigmasTrees[cutName];

      cout << "[HistoProducer::saveHistograms] cutName = " << cutName << endl;
      cout << "[HistoProducer::saveHistograms]   Directory name = " << dirName << endl;

      // Spectrum
      cout << "[HistoProducer::saveHistograms]   Number of Spectrum = " << vec_Spectrums.size() << endl;
      for(unsigned int i=0; i<vec_Spectrums.size(); i++){

        cout << "[HistoProducer::saveHistograms]   " << i << "-th Spectrum.." << endl;

        if(ic==0 && i==0 && FillMetaData){
          cout << "[HistoProducer::saveHistograms]     POT = " << vec_Spectrums.at(i)->POT() << endl;
          cout << "[HistoProducer::saveHistograms]     Livetime = " << vec_Spectrums.at(i)->Livetime() << endl;
          inputSamplePOT = vec_Spectrums.at(i)->POT();
          inputSampleLiveTime = vec_Spectrums.at(i)->Livetime();
        }

        TString hName = vec_Spectrums.at(i)->GetLabels()[0];

        if(vec_Spectrums.at(i)->GetBinnings().size()==1){
          TH1 *h = vec_Spectrums.at(i)->ToTH1( vec_Spectrums.at(i)->POT() );
          cout << "[HistoProducer::saveHistograms]     Writing TH1, \"" << hName << "\"" << endl;
          h->SetName(hName+"_"+dirName);
          h->Write();
        }
        else if(vec_Spectrums.at(i)->GetBinnings().size()==2){
          TH2 *h = vec_Spectrums.at(i)->ToTH2( vec_Spectrums.at(i)->POT() );
          cout << "[HistoProducer::saveHistograms]     Writing TH2, \"" << hName << "\"" << endl;
          h->SetName(hName+"_"+dirName);
          h->Write();
        }
        else if(vec_Spectrums.at(i)->GetBinnings().size()==3){
          TH3 *h = vec_Spectrums.at(i)->ToTH3( vec_Spectrums.at(i)->POT() );
          cout << "[HistoProducer::saveHistograms]     Writing TH3, \"" << hName << "\"" << endl;
          h->SetName(hName+"_"+dirName);
          h->Write();
        }

      }

      // EnsembleSpectrum
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
          TH1 *h = sNominal.ToTH1( sNominal.POT() );
          h->SetName(newLabel+"_"+dirName);
          h->Write();
          cout << "[HistoProducer::saveHistograms]     Nominal histogram = " << baseLabel << endl;
          for(unsigned int iu=0; iu<vec_SystEnsembleSpectrumPairs.at(i).second->NUniverses(); ++iu){
            TH1 *hU = vec_SystEnsembleSpectrumPairs.at(i).second->Universe(iu).ToTH1( vec_SystEnsembleSpectrumPairs.at(i).second->POT() );
            hU->SetName(baseLabel+"_"+systematicName+"_Univ"+TString::Itoa(iu,10)+"_"+dirName);
            hU->Write();
          }
          TGraphAsymmErrors *grErrorBand = vec_SystEnsembleSpectrumPairs.at(i).second->ErrorBand(1., vec_SystEnsembleSpectrumPairs.at(i).second->POT() );
          grErrorBand->SetName(baseLabel+"_"+systematicName+"_ErrorBandFromES_"+dirName);
          grErrorBand->Write();
        }
        else if(sNominal.GetBinnings().size()==2){
          for(unsigned int iu=0; iu<vec_SystEnsembleSpectrumPairs.at(i).second->NUniverses(); ++iu){
            TH2 *hU = vec_SystEnsembleSpectrumPairs.at(i).second->Universe(iu).ToTH2( vec_SystEnsembleSpectrumPairs.at(i).second->POT() );
            hU->SetName(baseLabel+"_"+systematicName+"_Univ"+TString::Itoa(iu,10)+"_"+dirName);
            hU->Write();
          }
        }

        dir->cd();

      }

      // ISyst-based Systematic
      cout << "[HistoProducer::saveHistograms]   Number of SystematicSpectrum = " << vec_SystSpectrumPairs.size() << endl;

      for(unsigned int i=0; i<vec_SystSpectrumPairs.size(); i++){

        TString systematicName = vec_SystSpectrumPairs.at(i).first;

        cout << "[HistoProducer::saveHistograms]   " << i << "-th Systematic Spectrum; " << systematicName << endl;

        if(vec_SystSpectrumPairs.at(i).second->GetBinnings().size()==1){
          TH1 *h = vec_SystSpectrumPairs.at(i).second->ToTH1( vec_SystSpectrumPairs.at(i).second->POT()  );
          TString baseLabel = vec_SystSpectrumPairs.at(i).second->GetLabels()[0];
          TString newLabel = baseLabel+"_"+systematicName;
          h->SetName(newLabel+"_"+dirName);
          h->Write();
          cout << "[HistoProducer::saveHistograms]     --> Done" << endl;
        }
        else if(vec_SystSpectrumPairs.at(i).second->GetBinnings().size()==2){
          TH2 *h = vec_SystSpectrumPairs.at(i).second->ToTH2( vec_SystSpectrumPairs.at(i).second->POT() );
          TString baseLabel = vec_SystSpectrumPairs.at(i).second->GetLabels()[0];
          TString newLabel = baseLabel+"_"+systematicName;
          h->SetName(newLabel+"_"+dirName);
          h->Write();
          cout << "[HistoProducer::saveHistograms]     --> Done" << endl;
        }

      }

      // Tree
      cout << "[HistoProducer::saveHistograms]   Number of Trees = " << vec_Trees.size() << endl;
      for(unsigned int i=0; i<vec_Trees.size(); i++){

        cout << "[HistoProducer::saveHistograms]   " << i << "-th Tree.." << endl;

        vec_Trees.at(i)->SaveTo(dir);

      }
      cout << "[HistoProducer::saveHistograms]   Number of NSigmasTrees = " << vec_NSigmasTrees.size() << endl;
      for(unsigned int i=0; i<vec_NSigmasTrees.size(); i++){

        cout << "[HistoProducer::saveHistograms]   " << i << "-th NSigmasTree.." << endl;

        vec_NSigmasTrees.at(i)->SaveTo(dir);

      }

      outputfile->cd();

    } // END loop cutname

    outputfile->cd();

  } // END if nCutName>0

  // When we have too many histograms, we have to run
  // Nominal and Universes separeately, and hadd later
  // For this case, some histograms should not be hadd-ed
  if(FillMetaData){

    outputfile->cd();
    outputfile->mkdir("JobInfo");
    outputfile->cd("JobInfo");
    TH1D *hist_TargetPOT = new TH1D("hist_TargetPOT", "", 1, 0., 1.);
    hist_TargetPOT->SetBinContent(1, TargetPOT);
    hist_TargetPOT->Write();

    outputfile->cd();
    outputfile->mkdir("BeamInfo");
    outputfile->cd("BeamInfo");
    TH1D *hist_POT = new TH1D("POT_BeamInfo", "POT", 1, 0., 1.);
    hist_POT->SetBinContent(1, inputSamplePOT);
    TH1D *hist_Livetime = new TH1D("Livetime_BeamInfo", "Livetime", 1, 0., 1.);
    hist_Livetime->SetBinContent(1, inputSampleLiveTime);
    hist_POT->Write();
    hist_Livetime->Write();
    outputfile->cd();
  }

}

void HistoProducer::setSystematicWeights(){

  FillSystematics = true;

  if(FillGENIESyst){

    cout << "[HistoProducer::setSystematicWeights] Setting GENIE systematics" << endl;

    // Multisigma
    cout << "[HistoProducer::setSystematicWeights] - Adding multisigma" << std::endl;
    const std::vector<std::string> genieMultisigmaKnobNames = ICARUSNumuXsec::GetGENIEMultisigmaKnobNames();
    for(const std::string& name: genieMultisigmaKnobNames){
      if(name=="FormZone"){
       std::cout << "[HistoProducer::setSystematicWeights] Skipping FormZone" << std::endl;
       continue;
      }
      std::string psetname = SystProviderPrefix+"_multisigma_"+name;
      std::cout << "[HistoProducer::setSystematicWeights] Multisigma, " << name << " (psetname = " << psetname << ")" << std::endl;
      IGENIESysts.push_back( new SBNWeightSyst(psetname) );
    }

    // Shape morphs
    cout << "[HistoProducer::setSystematicWeights] - Adding morphing dials" << std::endl;
    const std::vector<std::string> genieMorphKnobNames = ICARUSNumuXsec::GetGENIEMorphKnobNames();
    for(const std::string& name: genieMorphKnobNames){
      if(name=="FormZone"){
       std::cout << "[HistoProducer::setSystematicWeights] Skipping FormZone" << std::endl;
       continue;
      }
      std::string psetname = SystProviderPrefix+"_multisigma_"+name;
      std::cout << "[HistoProducer::setSystematicWeights] Multisigma, " << name << " (psetname = " << psetname << ")" << std::endl;
      IGENIEMorphSysts.push_back( new SBNWeightSyst(psetname) );
    }

    // Multisim for dependent dials
    cout << "[HistoProducer::setSystematicWeights] - Adding dependent dials by multisim" << std::endl;
    const std::vector<std::string> genieDependentKnobNames = ICARUSNumuXsec::GetGENIEDependentKnobNames();
    for(const std::string& name: genieDependentKnobNames){
      map_DepDialName_to_UniverseWeights[name] = {};
      map_DepDialName_to_UniverseSpillWeights[name] = {};
      std::string psetname = SystProviderPrefix+"_multisim_"+name;
      std::cout << "[HistoProducer::setSystematicWeights] Dependent dial, " << name << " (psetname = " << psetname << ")" << std::endl;
      for(int u=0; u<100; u++){
        map_DepDialName_to_UniverseWeights[name].push_back( GetUniverseWeight(psetname, u) );
        map_DepDialName_to_UniverseSpillWeights[name].push_back( GetUniverseFirstNeutrinoWeight(psetname, u) );
      }
    }

  }

  if(FillFlux){
    cout << "[HistoProducer::setSystematicWeights] Setting flux systematics" << endl;
    //IFluxSysts = GetNuMIPCAFluxSysts(NNuMIFluxPCA);
    IFluxSysts = GetAllNuMIFluxSysts(NNuMIFluxPCA);
    for(unsigned int i=0; i<IFluxSysts.size(); i++){
      cout << "[HistoProducer::setSystematicWeights] Syst = " << IFluxSysts.at(i)->ShortName() << endl;
    }
  }

  //cout << "[HistoProducer::setSystematicWeights] Setting detector systematics" << endl;
  //IDetectorSysts.push_back( new MuonMomentumScaleSyst(0.02) );

}

// slice
template<class T> 
void HistoProducer::FillCVSpectrum(SpectrumLoader& loader, const std::string& label, const T& var, const Binning& binning, SpillCut spillCut, Cut cut, bool ForceFill){

  if(!ForceFill && FillSystematics) return;

  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum(label, binning, loader, var, spillCut, cut, kNoShift, GetGlobalWeight() ) );

}
template<class T>
void HistoProducer::FillCVSpectrum(SpectrumLoader& loader, const std::string& label, const T& var1, const Binning& binning1, const T& var2, const Binning& binning2, SpillCut spillCut, Cut cut, bool ForceFill){

  if(!ForceFill && FillSystematics) return;

  map_cutName_to_vec_Spectrums[currentCutName].push_back( 
    new Spectrum(
      label, loader, 
      binning1, var1,
      binning2, var2,
      spillCut, cut,
      kNoShift,
      GetGlobalWeight()
    )
  );
}
// spill
template<class T>
void HistoProducer::FillCVSpectrum(SpectrumLoader& loader, const std::string& label, const T& spillvar, const Binning& binning, SpillCut spillCut, bool ForceFill){

  if(!ForceFill && FillSystematics) return;

  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum(label, binning, loader, spillvar, spillCut) );

}
template<class T>
void HistoProducer::FillCVSpectrum(SpectrumLoader& loader, const std::string& label, const T& spillvar1, const Binning& binning1, const T& spillvar2, const Binning& binning2, SpillCut spillCut, bool ForceFill){

  if(!ForceFill && FillSystematics) return;

  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      label, loader,
      binning1, spillvar1,
      binning2, spillvar2,
      spillCut
    )
  );
}

void HistoProducer::FillSystSpectrum(SpectrumLoader& loader, const std::string& label, const Var& var, const Binning& binning, SpillCut spillCut, Cut cut){

  // ensemble spectrum
  if(map_DepDialName_to_UniverseWeights.size()>0){
    const HistAxis ax(label, binning, var);
    for(const auto& it: map_DepDialName_to_UniverseWeights){
      map_cutName_to_vec_SystEnsembleSpectrumPairs[currentCutName].push_back(
        std::make_pair( it.first, new EnsembleSpectrum(loader, ax, spillCut, cut, it.second, GetGlobalWeight() ) )
      );
    }
  }

  // multisigma
  for(const auto& s: IGENIESysts){
    map_cutName_to_vec_SystSpectrumPairs[currentCutName].push_back(
      std::make_pair( s->ShortName()+"_Up", new Spectrum(label, loader, binning, var, spillCut, cut, SystShifts(s, +1), GetGlobalWeight() ) )
    );
    map_cutName_to_vec_SystSpectrumPairs[currentCutName].push_back(
      std::make_pair( s->ShortName()+"_Down", new Spectrum(label, loader, binning, var, spillCut, cut, SystShifts(s, -1), GetGlobalWeight() ) )
    );
  }

  // morph
  for(const auto& s: IGENIEMorphSysts){
    // +1: Morphed to another model
    map_cutName_to_vec_SystSpectrumPairs[currentCutName].push_back(
      std::make_pair( s->ShortName()+"_Up", new Spectrum(label, loader, binning, var, spillCut, cut, SystShifts(s, +1), GetGlobalWeight() ) )
    );
    // 0: CV, so we don't have _Down. Just saving nominal
    map_cutName_to_vec_SystSpectrumPairs[currentCutName].push_back(
      std::make_pair( s->ShortName()+"_Down", new Spectrum(label, loader, binning, var, spillCut, cut, kNoShift, GetGlobalWeight() ) )
    );
  }

  //TEST
  // zexp
  Var wZExpMultiplied_Up = GetUniverseWeight("GENIEReWeight_ICARUS_v2_multisigma_ZExpA1CCQE", 1) 
                           * GetUniverseWeight("GENIEReWeight_ICARUS_v2_multisigma_ZExpA2CCQE", 1)
                           * GetUniverseWeight("GENIEReWeight_ICARUS_v2_multisigma_ZExpA3CCQE", 1)
                           * GetUniverseWeight("GENIEReWeight_ICARUS_v2_multisigma_ZExpA4CCQE", 1);
  map_cutName_to_vec_SystSpectrumPairs[currentCutName].push_back(
    std::make_pair( "GENIEReWeight_ICARUS_v2_multisigma_ZExpMultiplied_Up", new Spectrum(label, loader, binning, var, spillCut, cut, kNoShift, GetGlobalWeight()*wZExpMultiplied_Up ) )
  );  

  Var wZExpMultiplied_Down = GetUniverseWeight("GENIEReWeight_ICARUS_v2_multisigma_ZExpA1CCQE", 0)
                             * GetUniverseWeight("GENIEReWeight_ICARUS_v2_multisigma_ZExpA2CCQE", 0)
                             * GetUniverseWeight("GENIEReWeight_ICARUS_v2_multisigma_ZExpA3CCQE", 0)
                             * GetUniverseWeight("GENIEReWeight_ICARUS_v2_multisigma_ZExpA4CCQE", 0);
  map_cutName_to_vec_SystSpectrumPairs[currentCutName].push_back(
    std::make_pair( "GENIEReWeight_ICARUS_v2_multisigma_ZExpMultiplied_Down", new Spectrum(label, loader, binning, var, spillCut, cut, kNoShift, GetGlobalWeight()*wZExpMultiplied_Down ) )
  );

  // fsi_pi
  Var wFSIPiMultiplied_Up = GetUniverseWeight("GENIEReWeight_ICARUS_v2_multisigma_MFP_pi", 1)
                            * GetUniverseWeight("GENIEReWeight_ICARUS_v2_multisigma_FrCEx_pi", 1)
                            * GetUniverseWeight("GENIEReWeight_ICARUS_v2_multisigma_FrInel_pi", 1)
                            * GetUniverseWeight("GENIEReWeight_ICARUS_v2_multisigma_FrAbs_pi", 1)
                            * GetUniverseWeight("GENIEReWeight_ICARUS_v2_multisigma_FrPiProd_pi", 1);
  map_cutName_to_vec_SystSpectrumPairs[currentCutName].push_back(
    std::make_pair( "GENIEReWeight_ICARUS_v2_multisigma_FSIPiMultiplied_Up", new Spectrum(label, loader, binning, var, spillCut, cut, kNoShift, GetGlobalWeight()*wFSIPiMultiplied_Up ) )
  );

  Var wFSIPiMultiplied_Down = GetUniverseWeight("GENIEReWeight_ICARUS_v2_multisigma_MFP_pi", 0)
                              * GetUniverseWeight("GENIEReWeight_ICARUS_v2_multisigma_FrCEx_pi", 0)
                              * GetUniverseWeight("GENIEReWeight_ICARUS_v2_multisigma_FrInel_pi", 0)
                              * GetUniverseWeight("GENIEReWeight_ICARUS_v2_multisigma_FrAbs_pi", 0)
                              * GetUniverseWeight("GENIEReWeight_ICARUS_v2_multisigma_FrPiProd_pi", 0);
  map_cutName_to_vec_SystSpectrumPairs[currentCutName].push_back(
    std::make_pair( "GENIEReWeight_ICARUS_v2_multisigma_FSIPiMultiplied_Down", new Spectrum(label, loader, binning, var, spillCut, cut, kNoShift, GetGlobalWeight()*wFSIPiMultiplied_Down ) )
  );

  // fsi_N
  Var wFSINMultiplied_Up = GetUniverseWeight("GENIEReWeight_ICARUS_v2_multisigma_MFP_N", 1)
                            * GetUniverseWeight("GENIEReWeight_ICARUS_v2_multisigma_FrCEx_N", 1)
                            * GetUniverseWeight("GENIEReWeight_ICARUS_v2_multisigma_FrInel_N", 1)
                            * GetUniverseWeight("GENIEReWeight_ICARUS_v2_multisigma_FrAbs_N", 1)
                            * GetUniverseWeight("GENIEReWeight_ICARUS_v2_multisigma_FrPiProd_N", 1);
  map_cutName_to_vec_SystSpectrumPairs[currentCutName].push_back(
    std::make_pair( "GENIEReWeight_ICARUS_v2_multisigma_FSINMultiplied_Up", new Spectrum(label, loader, binning, var, spillCut, cut, kNoShift, GetGlobalWeight()*wFSINMultiplied_Up ) )
  );

  Var wFSINMultiplied_Down = GetUniverseWeight("GENIEReWeight_ICARUS_v2_multisigma_MFP_N", 0)
                              * GetUniverseWeight("GENIEReWeight_ICARUS_v2_multisigma_FrCEx_N", 0)
                              * GetUniverseWeight("GENIEReWeight_ICARUS_v2_multisigma_FrInel_N", 0)
                              * GetUniverseWeight("GENIEReWeight_ICARUS_v2_multisigma_FrAbs_N", 0)
                              * GetUniverseWeight("GENIEReWeight_ICARUS_v2_multisigma_FrPiProd_N", 0);
  map_cutName_to_vec_SystSpectrumPairs[currentCutName].push_back(
    std::make_pair( "GENIEReWeight_ICARUS_v2_multisigma_FSINMultiplied_Down", new Spectrum(label, loader, binning, var, spillCut, cut, kNoShift, GetGlobalWeight()*wFSINMultiplied_Down ) )
  );

  for(const auto& s: IFluxSysts){
    map_cutName_to_vec_SystSpectrumPairs[currentCutName].push_back(
      std::make_pair( s->ShortName()+"_Up", new Spectrum(label, loader, binning, var, spillCut, cut, SystShifts(s, +1), GetGlobalWeight() ) )
    );
    map_cutName_to_vec_SystSpectrumPairs[currentCutName].push_back(
      std::make_pair( s->ShortName()+"_Down", new Spectrum(label, loader, binning, var, spillCut, cut, SystShifts(s, -1), GetGlobalWeight() ) )
    );
  }

}

void HistoProducer::FillSystSpectrum(SpectrumLoader& loader, const std::string& label, const SpillVar& var, const Binning& binning, SpillCut spillCut){

  // ensemble spectrum
  if(map_DepDialName_to_UniverseSpillWeights.size()>0){
    for(const auto& it: map_DepDialName_to_UniverseSpillWeights){
      map_cutName_to_vec_SystEnsembleSpectrumPairs[currentCutName].push_back(
        std::make_pair( it.first, new EnsembleSpectrum(label, binning, loader, var, spillCut, it.second, kSpillUnweighted ) )
      );
    }
  }

  static const std::vector<std::string> genieMultisigmaKnobNames = ICARUSNumuXsec::GetGENIEMultisigmaKnobNames();
  for(const std::string& name: genieMultisigmaKnobNames){
    if(name=="FormZone"){
     continue;
    }
    std::string psetname = SystProviderPrefix+"_multisigma_"+name;

    SpillVar wUp = GetUniverseFirstNeutrinoWeight(psetname, 1);
    SpillVar wDown = GetUniverseFirstNeutrinoWeight(psetname, 0);

    map_cutName_to_vec_SystSpectrumPairs[currentCutName].push_back(
      std::make_pair( psetname+"_Up", new Spectrum(label, binning, loader, var, spillCut, wUp ) )
    );
    map_cutName_to_vec_SystSpectrumPairs[currentCutName].push_back(
      std::make_pair( psetname+"_Down", new Spectrum(label, binning, loader, var, spillCut, wUp ) )
    );

  }

}

void HistoProducer::FillCVandSystSpectrum(SpectrumLoader& loader, const std::string& label, const Var& var, const Binning& binning, SpillCut spillCut, Cut cut){

  FillCVSpectrum(loader, label, var, binning, spillCut, cut, true);
  if(FillSystematics) FillSystSpectrum(loader, label, var, binning, spillCut, cut);

}

void HistoProducer::FillCVandSystSpectrum(SpectrumLoader& loader, const std::string& label, const SpillVar& var, const Binning& binning, SpillCut spillCut){

  FillCVSpectrum(loader, label, var, binning, spillCut, true);
  if(FillSystematics) FillSystSpectrum(loader, label, var, binning, spillCut);

}

const Var HistoProducer::GetGlobalWeight(){

  Var GlobalWeight = kUnweighted;
  if(ApplyNuMIPPFXCVWeight) GlobalWeight = GlobalWeight * ana::kGetNuMIFluxWeight;

  return GlobalWeight;

}

HistoProducer::~HistoProducer(){

  cout << "[HistoProducer::~HistoProducer] called" << endl;

  outputfile->Close();

  cout << "[HistoProducer::~HistoProducer] output file : " << outputDir+"/"+outputName << endl;

  cout << "[HistoProducer::~HistoProducer] Finished" << endl;

}

