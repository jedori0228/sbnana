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

  nSigmasSaveMode = kVector;
  MakeGUNDAMTree = false;

  ApplyNuMIPPFXCVWeight = false;

  IsData = false;

  FillMetaData = true;
  FillSystematics = false;

  FillGENIESyst = true;
  FillNuSyst = false;
  FillFlux = true;
  FillGEANT4 = true;
  ApplyTrackSplit = false;
  FillFakeDataWeights = false;

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

  FillCVSpectrum(loader, "TruthMuonKE", kNuMITrueMuonKE, Binning::Simple(300, 0., 3.), spillCut, cut);
  FillCVSpectrum(loader, "TruthMuonP", kNuMIMuonTrueP, Binning::Simple(50, 0., 5.), spillCut, cut);
  FillCVSpectrum(loader, "TruthMuonPt", kNuMITrueMuonPt, Binning::Simple(20, 0., 2.), spillCut, cut);
  FillCVSpectrum(loader, "TruthMuonLength", TruthMuonLength, Binning::Simple(500, 0., 500), spillCut, cut);
  FillCVSpectrum(loader, "TruthMuonNuCosineTheta", kNuMITrueMuonNuCosineTheta, Binning::Simple(40, -1., 1.), spillCut, cut);
  FillCVSpectrum(loader, "TruthMuonMatchedTrackCustomChi2InelasticPionCollection", TruthMuonMatchedTrackCustomChi2InelasticPionCollection, Binning::Simple(80, 0., 80), spillCut, cut);

  FillCVSpectrum(loader, "TruthMuonMichelStartProcess", TruthMuonMichelStartProcess, Binning::Simple(65, 0., 65.), spillCut, cut);
  FillCVSpectrum(loader, "TruthMuonMichelKE", TruthMuonMichelKE, Binning::Simple(100, 0., 0.1), spillCut, cut);
  FillCVSpectrum(loader, "TruthMuonMichelMatchedShowerKE", TruthMuonMichelMatchedShowerKE, Binning::Simple(100, 0., 0.1), spillCut, cut);
  FillCVSpectrum(loader, "TruthMuonMichelMatchedShowerDistanceFromMuonEnd", TruthMuonMichelMatchedShowerDistanceFromMuonEnd, Binning::Simple(300, 0., 30.), spillCut, cut);

  // Proton

  FillCVSpectrum(loader, "TruthProtonKE", kNuMITrueProtonKE, Binning::Simple(100, 0., 1.0), spillCut, cut);
  FillCVSpectrum(loader, "TruthProtonP", kNuMIProtonTrueP, Binning::Simple(100, 0., 1.0), spillCut, cut);
  FillCVSpectrum(loader, "TruthProtonLength", TruthProtonLength, Binning::Simple(100, 0., 100.), spillCut, cut);
  FillCVSpectrum(loader, "TruthProtonNuCosineTheta", kNuMITrueProtonNuCosineTheta, Binning::Simple(40, -1., 1.), spillCut, cut);
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

  FillCVSpectrum(loader, "CountSlice", varCountSlice, Binning::Simple(1, 0.,1.), spillCut, cut);
  FillCVSpectrum(loader, "MuonChi2MuonOriginal", kNuMIRecoMuonTrackChi2Muon, Binning::Simple(300, 0., 30.), spillCut, cut);
  FillCVSpectrum(loader, "MuonFloatChi2MuonPlusMichelDebug", kNuMIRecoMuonFloatChi2MuonPlusMichelDebug, Binning::Simple(300, 0., 30.), spillCut, cut);

  FillCVSpectrum(loader, "MuonChi2MIP5cm", kNuMIRecoMuonChi2MIP5cm, Binning::Simple(200, 0., 20.), spillCut, cut);
  FillCVSpectrum(loader, "MuonChi2MIP10cm", kNuMIRecoMuonChi2MIP10cm, Binning::Simple(200, 0., 20.), spillCut, cut);
  FillCVSpectrum(loader, "MuonChi2MuonPlusMichel5cmShift", kNuMIRecoMuonChi2MuonPlusMichel5cmShift, Binning::Simple(1000, 0., 100.), spillCut, cut);
  FillCVSpectrum(loader, "MuonChi2MuonPlusMichel5cmNoShift", kNuMIRecoMuonChi2MuonPlusMichel5cmNoShift, Binning::Simple(1000, 0., 100.), spillCut, cut);
  FillCVSpectrum(loader, "MuonChi2MuonPlusMichel10cmShift", kNuMIRecoMuonChi2MuonPlusMichel10cmShift, Binning::Simple(1000, 0., 100.), spillCut, cut);
  FillCVSpectrum(loader, "MuonChi2MuonPlusMichel10cmNoShift", kNuMIRecoMuonChi2MuonPlusMichel10cmNoShift, Binning::Simple(1000, 0., 100.), spillCut, cut);

  FillCVSpectrum(loader, "MuonFloatChi2MuonPlusMichel1cm", kNuMIRecoMuonFloatChi2MuonPlusMichel1cm, Binning::Simple(1000, 0., 100.), spillCut, cut);
  FillCVSpectrum(loader, "MuonFloatChi2MuonPlusMichel2cm", kNuMIRecoMuonFloatChi2MuonPlusMichel2cm, Binning::Simple(1000, 0., 100.), spillCut, cut);
  FillCVSpectrum(loader, "MuonFloatChi2MuonPlusMichel3cm", kNuMIRecoMuonFloatChi2MuonPlusMichel3cm, Binning::Simple(1000, 0., 100.), spillCut, cut);
  FillCVSpectrum(loader, "MuonFloatChi2MuonPlusMichel4cm", kNuMIRecoMuonFloatChi2MuonPlusMichel4cm, Binning::Simple(1000, 0., 100.), spillCut, cut);
  FillCVSpectrum(loader, "MuonFloatChi2MuonPlusMichel5cm", kNuMIRecoMuonFloatChi2MuonPlusMichel5cm, Binning::Simple(1000, 0., 100.), spillCut, cut);
  FillCVSpectrum(loader, "MuonFloatChi2MuonPlusMichel10cm", kNuMIRecoMuonFloatChi2MuonPlusMichel10cm, Binning::Simple(1000, 0., 100.), spillCut, cut);
  FillCVSpectrum(loader, "MuonFloatChi2MuonPlusMichel15cm", kNuMIRecoMuonFloatChi2MuonPlusMichel15cm, Binning::Simple(1000, 0., 100.), spillCut, cut);

  FillCVSpectrum(loader, "MuonFloatChi2MuonPlusMichel15cmDelta", kNuMIRecoMuonFloatChi2MuonPlusMichel15cmDelta, Binning::Simple(3000, 0., 30.), spillCut, cut);
  FillCVSpectrum(loader, "MuonFloatChi2MuonPlusMichelDebugDelta", kNuMIRecoMuonFloatChi2MuonPlusMichelDebugDelta, Binning::Simple(3000, 0., 30.), spillCut, cut);

  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "MuonRR_vs_MuondEdX", loader,
      Binning::Simple(30, 0., 30.), kNuMIMuonCandidateRR,
      Binning::Simple(100., 0., 10.), kNuMIMuonCandidatedEdX,
      spillCut, cut
    )
  );

  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "Chi2Muon_Original_vs_Chi2Muon_Fitted10cm", loader,
      Binning::Simple(300., 0., 30.), kNuMIRecoMuonTrackChi2Muon,
      Binning::Simple(300, 0., 30.), kNuMIRecoMuonFloatChi2MuonPlusMichel10cm,
      spillCut, cut
    )
  );


}

void HistoProducer::MakeTree(SpectrumLoader& loader, SpillCut spillCut, Cut cut){

  // Some spectrums

  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "OpFlashFirstTime", Binning::Simple(500, -20., 30.),
      loader,
      OpFlashFirstTime,
      kNoSpillCut
    )
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "OpFlashTime", Binning::Simple(500, -20., 30.),
      loader,
      OpFlashTime,
      kNoSpillCut
    )
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "TriggerTime", Binning::Simple(500, -20., 30.),
      loader,
      kNuMISpillTriggerTime,
      kNoSpillCut
    )
  );

  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "OpFlashFirstTime_ValidTrig", Binning::Simple(500, -20., 30.),
      loader,
      OpFlashFirstTime,
      kNuMIValidTrigger
    )
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "OpFlashTime_ValidTrig", Binning::Simple(500, -20., 30.),
      loader,
      OpFlashTime,
      kNuMIValidTrigger
    )
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "TriggerTime_ValidTrig", Binning::Simple(500, -20., 30.),
      loader,
      kNuMISpillTriggerTime,
      kNuMIValidTrigger
    )
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "TriggerTime_ValidTrig_WithG3ChaseWeight", Binning::Simple(500, -20., 30.),
      loader,
      kNuMISpillTriggerTime,
      kNuMIValidTrigger,
      kNuMIG3ChaseSpillWeightByClosesetNu
    )
  );

  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "OpFlashTimeAfterSignalSelection", Binning::Simple(500, -20., 30.),
      loader,
      OpFlashTimeAfterSignalSelection,
      kNoSpillCut
    )
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "OpFlashTimeAfterSignalSelection_ValidTrig", Binning::Simple(500, -20., 30.),
      loader,
      OpFlashTimeAfterSignalSelection,
      kNuMIValidTrigger
    )
  );

  // base variables

  std::vector<std::string> this_reco_labels = GetNuMIRecoTreeLabels();
  std::vector<Var> this_reco_vars = GetNuMIRecoTreeVars();

  // NSigma

  std::vector<std::string> this_NSigmasPsetNames;
  std::vector<const ISyst*> this_NSigmasISysts;
  std::vector<std::pair<int,int>> this_NSigmasPairs;
  std::vector<std::vector<double>> this_NSigmas;

  // NUniverses

  std::vector<std::string> this_NUniversesPsetNames;
  std::vector<std::vector<Var>> this_NUniversesVarVectors;
  std::vector<std::vector<TruthVar>> this_NUniversesTruthVarVectors;
  std::vector<unsigned int> this_NUniversesNUnivs;

  if(FillFakeDataWeights){

    std::vector<std::string> this_lowq2rec_study_labels = GetFakeDataWeightsLabels();
    std::vector<Var> this_lowq2rec_study_vars = GetFakeDataWeightsVars();

    this_reco_labels.insert(
      this_reco_labels.end(),
      this_lowq2rec_study_labels.begin(),
      this_lowq2rec_study_labels.end()
    );

    this_reco_vars.insert(
      this_reco_vars.end(),
      this_lowq2rec_study_vars.begin(),
      this_lowq2rec_study_vars.end()
    );

  }

  if(FillSystematics){

    // NSigma

    for(unsigned int i=0; i<genieMultisigmaKnobNames.size(); i++){
      this_NSigmasPsetNames.push_back( genieMultisigmaKnobNames.at(i) );
      this_NSigmasISysts.push_back( IGENIESysts.at(i) );
      this_NSigmasPairs.push_back( std::make_pair(-3, 3) );
      this_NSigmas.push_back( {-3, -2, -1, 0, 1, 2, 3} );
    }
    for(unsigned int i=0; i<genieMorphKnobNames.size(); i++){
      this_NSigmasPsetNames.push_back( genieMorphKnobNames.at(i) );
      this_NSigmasISysts.push_back( IGENIEMorphSysts.at(i) );
      this_NSigmasPairs.push_back( std::make_pair(0, 1) );
      this_NSigmas.push_back( {-1, -0.5, 0, 0.5, 1} );
    }
    for(unsigned int i=0; i<IFluxSysts.size(); i++){
      this_NSigmasPsetNames.push_back( IFluxSysts.at(i)->ShortName() );
      this_NSigmasISysts.push_back( IFluxSysts.at(i) );
      this_NSigmasPairs.push_back( std::make_pair(-3, 3) );
      this_NSigmas.push_back( {-3, -2, -1, 0, 1, 2, 3} );
    }
    for(unsigned int i=0; i<IDetectorSysts.size(); i++){
      this_NSigmasPsetNames.push_back( IDetectorSysts.at(i)->ShortName() );
      this_NSigmasISysts.push_back( IDetectorSysts.at(i) );
      this_NSigmasPairs.push_back( std::make_pair(-3, 3) );
      this_NSigmas.push_back( {-3, -2, -1, 0, 1, 2, 3} );
    }

    // NUniverses

    for(const std::string& name: genieDependentKnobNames){
      this_NUniversesPsetNames.push_back( name );
      this_NUniversesVarVectors.push_back( map_DepDialName_to_UniverseWeights[name] );
      this_NUniversesTruthVarVectors.push_back( map_DepDialName_to_TruthUniverseWeights[name] );
      this_NUniversesNUnivs.push_back( 100 );
    }
    for(const std::string& name: geant4DependentKnobNames){
      this_NUniversesPsetNames.push_back( name );
      this_NUniversesVarVectors.push_back( map_DepDialName_to_UniverseWeights[name] );
      this_NUniversesTruthVarVectors.push_back( map_DepDialName_to_TruthUniverseWeights[name] );
      this_NUniversesNUnivs.push_back( 1000 );
    }

  }

  if(IsData){

    // Data CV only

    this_reco_labels.push_back( "IsData/i" );
    this_reco_vars.push_back(
      Var([](const caf::SRSliceProxy* slc) -> int {
        return 1;
      })
    );

    // CV

    map_cutName_to_vec_Trees[currentCutName].push_back(
      new ana::Tree(
        "selectedEvents", this_reco_labels,
        loader,
        this_reco_vars,
        spillCut, cut,
        kNoShift,
        true, true
      )
    );

    if(FillSystematics){

      // NSigmasTree

      map_cutName_to_vec_NSigmasTrees[currentCutName].push_back(
        new ana::NSigmasTree(
          "selectedEvents_NSigmas",
          this_NSigmasPsetNames,
          loader,
          this_NSigmasISysts,
          //this_NSigmasPairs,
          this_NSigmas,
          spillCut, cut,
          kNoShift,
          true, true
        )
      );

      // NUniversesTree

      map_cutName_to_vec_NUniversesTrees[currentCutName].push_back(

        new ana::NUniversesTree(
          "selectedEvents_NUniverses",
          this_NUniversesPsetNames,
          loader,
          this_NUniversesVarVectors,
          this_NUniversesNUnivs,
          spillCut, cut,
          kNoShift,
          true, true
        )

      );


    }

  }
  else{

    // MC

    this_reco_labels.push_back( "IsData/i" );
    this_reco_vars.push_back( 
      Var([](const caf::SRSliceProxy* slc) -> int {
        return 0;
      })
    );

    // CV

    map_cutName_to_vec_Trees[currentCutName].push_back(
      new ana::Tree(
        "selectedEvents", this_reco_labels,
        loader,
        this_reco_vars,
        spillCut, cut,
        ApplyTrackSplit ? SystShifts(&kTrackSplittingSyst, +1.) : kNoShift,
        true, true
      )
    );

    if(FillSystematics){

      // Shifted tree
      // 1) Calo dE/dX up
      map_cutName_to_vec_Trees[currentCutName].push_back(
        new ana::Tree(
          "selectedEvents_CalodEdXShiftUp", this_reco_labels,
          loader,
          this_reco_vars,
          spillCut, cut,
          ApplyTrackSplit ? SystShifts( {{&kTrackSplittingSyst, +1.}, {&kCalodEdXShiftSyst, +1.}} ) : SystShifts(&kCalodEdXShiftSyst, +1.),
          true, true
        )
      );
      // 2) Calo dE/dX down
      map_cutName_to_vec_Trees[currentCutName].push_back(
        new ana::Tree(
          "selectedEvents_CalodEdXShiftDown", this_reco_labels,
          loader,
          this_reco_vars,
          spillCut, cut,
          ApplyTrackSplit ? SystShifts( {{&kTrackSplittingSyst, +1.}, {&kCalodEdXShiftSyst, -1.}} ) : SystShifts(&kCalodEdXShiftSyst, -1.),
          true, true
        )
      );
      // 3) Calo gain up
      map_cutName_to_vec_Trees[currentCutName].push_back(
        new ana::Tree(
          "selectedEvents_CaloGainShiftUp", this_reco_labels,
          loader,
          this_reco_vars,
          spillCut, cut,
          ApplyTrackSplit ? SystShifts( {{&kTrackSplittingSyst, +1.}, {&kCaloGainShiftSyst, +1.}} ) : SystShifts(&kCaloGainShiftSyst, +1.),
          true, true
        )
      );
      // 4) Calo gain down
      map_cutName_to_vec_Trees[currentCutName].push_back(
        new ana::Tree(
          "selectedEvents_CaloGainShiftDown", this_reco_labels,
          loader,
          this_reco_vars,
          spillCut, cut,
          ApplyTrackSplit ? SystShifts( {{&kTrackSplittingSyst, +1.}, {&kCaloGainShiftSyst, -1.}} ) : SystShifts(&kCaloGainShiftSyst, -1.),
          true, true
        )
      );

      // NSigmasTree

      map_cutName_to_vec_NSigmasTrees[currentCutName].push_back(
        new ana::NSigmasTree(
          "selectedEvents_NSigmas",
          this_NSigmasPsetNames,
          loader,
          this_NSigmasISysts,
          //this_NSigmasPairs,
          this_NSigmas,
          spillCut, cut,
          ApplyTrackSplit ? SystShifts(&kTrackSplittingSyst, +1.) : kNoShift,
          true, true
        )
      );

      // NUniversesTree

      map_cutName_to_vec_NUniversesTrees[currentCutName].push_back(

        new ana::NUniversesTree(
          "selectedEvents_NUniverses",
          this_NUniversesPsetNames,
          loader,
          this_NUniversesVarVectors,
          this_NUniversesNUnivs,
          spillCut, cut,
          ApplyTrackSplit ? SystShifts(&kTrackSplittingSyst, +1.) : kNoShift,
          true, true
        )

      );

    }

    // TrueTree
    if(TrueTreeFilled) return;

    std::vector< std::pair<std::string, Cut> > RecoCutsForEffs = {
      std::make_pair("", kNuMISelection_1muNp0pi),
      std::make_pair("_Contained", kNuMISelection_1muNp0pi && kNuMIMuonCandidateContained),
      std::make_pair("_Exiting", kNuMISelection_1muNp0pi && !kNuMIMuonCandidateContained),
    };

    for(unsigned i_Cut=0; i_Cut<RecoCutsForEffs.size(); i_Cut++){

      map_cutName_to_vec_Trees[currentCutName].push_back(
        new ana::Tree(
          "trueEvents"+RecoCutsForEffs[i_Cut].first, 
          GetNuMITrueTreeLabels(), loader, GetNuMITrueTreeVars(), kNuMIValidTrigger, kTruthCut_IsSignal,
          RecoCutsForEffs[i_Cut].second,
          ApplyTrackSplit ? SystShifts(&kTrackSplittingSyst, +1.) : kNoShift,
          true
        )
      );

      // Fill systematics

      if(FillSystematics){

        map_cutName_to_vec_Trees[currentCutName].push_back(
          new ana::Tree(
            "trueEvents"+RecoCutsForEffs[i_Cut].first+"_CalodEdXShiftUp",
            GetNuMITrueTreeLabels(), loader, GetNuMITrueTreeVars(), kNuMIValidTrigger, kTruthCut_IsSignal,
            RecoCutsForEffs[i_Cut].second,
            ApplyTrackSplit ? SystShifts( {{&kTrackSplittingSyst, +1.}, {&kCalodEdXShiftSyst, +1.}} ) : SystShifts(&kCalodEdXShiftSyst, +1.),
            true
          )
        );

        map_cutName_to_vec_Trees[currentCutName].push_back(
          new ana::Tree(
            "trueEvents"+RecoCutsForEffs[i_Cut].first+"_CalodEdXShiftDown",
            GetNuMITrueTreeLabels(), loader, GetNuMITrueTreeVars(), kNuMIValidTrigger, kTruthCut_IsSignal,
            RecoCutsForEffs[i_Cut].second,
            ApplyTrackSplit ? SystShifts( {{&kTrackSplittingSyst, +1.}, {&kCalodEdXShiftSyst, -1.}} ) : SystShifts(&kCalodEdXShiftSyst, -1.),
            true
          )
        );

        map_cutName_to_vec_Trees[currentCutName].push_back(
          new ana::Tree(
            "trueEvents"+RecoCutsForEffs[i_Cut].first+"_CaloGainShiftUp",
            GetNuMITrueTreeLabels(), loader, GetNuMITrueTreeVars(), kNuMIValidTrigger, kTruthCut_IsSignal,
            RecoCutsForEffs[i_Cut].second,
            ApplyTrackSplit ? SystShifts( {{&kTrackSplittingSyst, +1.}, {&kCaloGainShiftSyst, +1.}} ) : SystShifts(&kCaloGainShiftSyst, +1.),
            true
          )
        );

        map_cutName_to_vec_Trees[currentCutName].push_back(
          new ana::Tree(
            "trueEvents"+RecoCutsForEffs[i_Cut].first+"_CaloGainShiftDown",
            GetNuMITrueTreeLabels(), loader, GetNuMITrueTreeVars(), kNuMIValidTrigger, kTruthCut_IsSignal,
            RecoCutsForEffs[i_Cut].second,
            ApplyTrackSplit ? SystShifts( {{&kTrackSplittingSyst, +1.}, {&kCaloGainShiftSyst, -1.}} ) : SystShifts(&kCaloGainShiftSyst, -1.),
            true
          )
        );

        for(unsigned int i=0; i<IDetectorSysts.size(); i++){
          map_cutName_to_vec_Trees[currentCutName].push_back(
            new ana::Tree(
              "trueEvents"+RecoCutsForEffs[i_Cut].first+"_"+IDetectorSysts.at(i)->ShortName()+"Up",
              {"Dummy"}, loader, {DummyTruthVar}, kNuMIValidTrigger, kTruthCut_IsSignal,
              RecoCutsForEffs[i_Cut].second,
              ApplyTrackSplit ? SystShifts( {{&kTrackSplittingSyst, +1.}, {IDetectorSysts.at(i), +1.}} ) : SystShifts(IDetectorSysts.at(i), +1.),
              true
            )
          );
          map_cutName_to_vec_Trees[currentCutName].push_back(
            new ana::Tree(
              "trueEvents"+RecoCutsForEffs[i_Cut].first+"_"+IDetectorSysts.at(i)->ShortName()+"Down",
              {"Dummy"}, loader, {DummyTruthVar}, kNuMIValidTrigger, kTruthCut_IsSignal,
              RecoCutsForEffs[i_Cut].second,
              ApplyTrackSplit ? SystShifts( {{&kTrackSplittingSyst, +1.}, {IDetectorSysts.at(i), -1.}} ) : SystShifts(IDetectorSysts.at(i), -1.),
              true
            )
          );

        } // END IDetectorSysts loop

      } // END if FillSystematics

    } // END RecoCutsForEffs loop

    if(FillSystematics){

      map_cutName_to_vec_NSigmasTrees[currentCutName].push_back(
        new ana::NSigmasTree(
          "trueEvents_NSigmas",
          this_NSigmasPsetNames,
          loader,
          this_NSigmasISysts,
          //this_NSigmasPairs,
          this_NSigmas,
          kTruthCut_IsSignal,
          ApplyTrackSplit ? SystShifts(&kTrackSplittingSyst, +1.) : kNoShift,
          true
        )
      );

      map_cutName_to_vec_NUniversesTrees[currentCutName].push_back(
        new ana::NUniversesTree(
          "trueEvents_NUniverses",
          this_NUniversesPsetNames,
          loader,
          this_NUniversesTruthVarVectors,
          this_NUniversesNUnivs,
          kTruthCut_IsSignal,
          ApplyTrackSplit ? SystShifts(&kTrackSplittingSyst, +1.) : kNoShift,
          true
        )
      );

    }

    TrueTreeFilled = true;

  } // IF MC



}

void HistoProducer::NuMIXSecBkgdStudy(SpectrumLoader& loader, SpillCut spillCut, Cut cut){

  FillCVSpectrum(loader, "ChargedPionMichelMatchedPfpScores", kNuMIChargedPionMichelMatchedPfpScores, Binning::Simple(40, 0., 1.), spillCut, cut);
  FillCVSpectrum(loader, "ChargedPionMichelMatchedPfpTrackLengths", kNuMIChargedPionMichelMatchedPfpTrackLengths, Binning::Simple(50, 0., 50.), spillCut, cut);
  FillCVSpectrum(loader, "ChargedPionMichelMatchedPfpTrackDistances", kNuMIChargedPionMichelMatchedPfpTrackDistances, Binning::Simple(200, 0., 200.), spillCut, cut);
  FillCVSpectrum(loader, "ChargedPionMichelMatchedPfpShowerEnergies", kNuMIChargedPionMichelMatchedPfpShowerEnergies, Binning::Simple(109, 0., 0.1), spillCut, cut);
  FillCVSpectrum(loader, "ChargedPionMichelMatchedPfpShowerLengths", kNuMIChargedPionMichelMatchedPfpShowerLengths, Binning::Simple(50, 0., 50.), spillCut, cut);
  FillCVSpectrum(loader, "ChargedPionMichelMatchedPfpShowerGaps", kNuMIChargedPionMichelMatchedPfpShowerGaps, Binning::Simple(50, 0., 50.), spillCut, cut);
  FillCVSpectrum(loader, "ChargedPionMichelMatchedPfpShowerEnergySum", kNuMIChargedPionMichelMatchedPfpShowerEnergySum, Binning::Simple(200, 0., 0.2), spillCut, cut);

  FillCVSpectrum(loader, "TrueChargedPionMichelEnergy", kNuMITrueChargedPionMichelEnergy, Binning::Simple(109, 0., 0.1), spillCut, cut);


  FillCVSpectrum(loader, "ChargedPionMichelMatchedPfpTrackIsPrimaries", kNuMIChargedPionMichelMatchedPfpTrackIsPrimaries, Binning::Simple(2, 0., 2.), spillCut, cut);
  FillCVSpectrum(loader, "ChargedPionMichelMatchedPfpTrackHitPurities", kNuMIChargedPionMichelMatchedPfpTrackHitPurities, Binning::Simple(40., 0., 1.), spillCut, cut);
  FillCVSpectrum(loader, "ChargedPionMichelMatchedPfpTrackHitCompletenesses", kNuMIChargedPionMichelMatchedPfpTrackHitCompletenesses, Binning::Simple(40., 0., 1.), spillCut, cut);

  FillCVSpectrum(loader, "ChargedPionMichelMatchedPfpShowerOpeningAngles", kNuMIChargedPionMichelMatchedPfpShowerOpeningAngles, Binning::Simple(100, -1., M_PI), spillCut, cut);

  //FillCVSpectrum(loader, "MichelCandidateShowerEnergySum", kNuMIMichelCandidateShowerEnergySum, Binning::Simple(200, 0., 0.2), spillCut, cut);


/*
  FillCVSpectrum(loader, "NChargedPionMatchedTracks", kNuMINChargedPionMatchedTracks, Binning::Simple(10, 0., 10.), spillCut, cut);
  FillCVSpectrum(loader, "ChargedPionMatchedTrackScores", kNuMIChargedPionMatchedTrackScores, Binning::Simple(10, 0., 1.), spillCut, cut);
  FillCVSpectrum(loader, "ChargedPionEndProcess", kNuMITrueChargedPionEndProcess, Binning::Simple(65, 0., 65.), spillCut, cut);
  FillCVSpectrum(loader, "ChargedPionMatchedTrackLengths", kNuMIChargedPionMatchedTrackLengths, Binning::Simple(20, 0., 200.), spillCut, cut);
  FillCVSpectrum(loader, "ChargedPionMatchedTrackChi2Muons", kNuMIChargedPionMatchedTrackChi2Muons, Binning::Simple(100, 0., 100.), spillCut, cut);
  FillCVSpectrum(loader, "ChargedPionMatchedTrackChi2Protons", kNuMIChargedPionMatchedTrackChi2Protons, Binning::Simple(200, 0., 200.), spillCut, cut);

  FillCVSpectrum(loader, "NChargedPionMatchedShowers", kNuMINChargedPionMatchedShowers, Binning::Simple(10, 0., 10.), spillCut, cut);
  FillCVSpectrum(loader, "ChargedPionMatchedShowerScores", kNuMIChargedPionMatchedShowerScores, Binning::Simple(10, 0., 1.), spillCut, cut);
  FillCVSpectrum(loader, "ChargedPionMatchedShowerGaps", kNuMIChargedPionMatchedShowerGaps, Binning::Simple(30, 0., 30.), spillCut, cut);
  FillCVSpectrum(loader, "ChargedPionMatchedShowerIsPrimaries", kNuMIChargedPionMatchedShowerIsPrimaries, Binning::Simple(2, 0., 2.), spillCut, cut);
  FillCVSpectrum(loader, "ChargedPionMatchedShowerEnergies", kNuMIChargedPionMatchedShowerEnergies, Binning::Simple(100, 0., 1.), spillCut, cut);

  FillCVSpectrum(loader, "NChargedPionShowerCandidates", kNuMINChargedPionShowerCandidates, Binning::Simple(10, 0., 10.), spillCut, cut);
*/

/*
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "TEST",
      Binning::Simple(1, 0., 1.),
      loader,
      kTruth_ChargedPionMichelIndexTEST,
      kNoTruthCut,
      kNoSpillCut
    )
  );
*/

}

// - 230908_PIDStudy
void HistoProducer::MakePIDStudyTree(SpectrumLoader& loader, SpillCut spillCut, Cut cut){

  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "OpFlashFirstTime", Binning::Simple(500, -20., 30.),
      loader,
      OpFlashFirstTime,
      kNoSpillCut
    )
  );

  std::vector<std::string> labels = {
    // IsCosmicSlice, for cosmic slice reweight
    "IsCosmicSlice/i",
    // Weight
    "FluxWeight",
    "FluxWeightWithG3Chase",
    "FluxWeightWithG4Updated",
    "SPPCVCorrection",
    // Muon
    "MuonSelection/i",
    "MuonMatchType/i",
    "MuonLength",
    "MuonChi2Muon",
    "MuonChi2Proton",
    "MuonTrackScore",
    "MuonMatchedTruthPDG/i",
    "MuonMatchedTruthIntID/i",
    "MuonMatchedTruthContained/i",
    // TagMuon
    "TagMuonLength",
    // Proton
    "ProtonSelection/i",
    "ProtonMatchType/i",
    "ProtonP",
    "ProtonLength",
    "ProtonChi2Muon",
    "ProtonChi2Proton",
    "ProtonTrackScore",
    "ProtonMatchedTruthPDG/i",
    "ProtonMatchedTruthIntID/i",
    "ProtonMatchedTruthContained/i",
    // ChargedPion
    "ChargedPionSelection/i",
    "ChargedPionMatchType/i",
    "ChargedPionLength",
    "ChargedPionChi2Muon",
    "ChargedPionChi2Proton",
    "ChargedPionTrackScore",
    "ChargedPionMatchedTruthPDG/i",
    "ChargedPionMatchedTruthIntID/i",
    "ChargedPionMatchedTruthContained/i",
  };
  std::vector<Var> varlists = {
    // IsCosmicSlice, for cosmic slice reweight
    IsCosmicSlice,
    // Weight
    kGetNuMIFluxWeight,
    kGetNuMIFluxWeightG3Chase,
    kGetNuMIFluxWeightUpdated,
    kNuMISPPCVCorrection,
    // Muon
    kNuMIIsRelaxedMuonSelection,
    kNuMIRelaxedMuonTrackMatchType,
    kNuMIRelaxedMuonTrackLength,
    kNuMIRelaxedMuonTrackChi2Muon,
    kNuMIRelaxedMuonTrackChi2Proton,
    kNuMIRelaxedMuonTrackScore,
    kNuMIRelaxedMuonTrackMatchedTruthPDG,
    kNuMIRelaxedMuonTrackMatchedTruthIntID,
    kNuMIRelaxedMuonTrackMatchedTruthContained,
    // TagMuon
    kNuMITagMuonLength,
    // Proton
    kNuMIIsRelaxedProtonSelection,
    kNuMIRelaxedProtonTrackMatchType,
    kNuMIRelaxedProtonTrackP,
    kNuMIRelaxedProtonTrackLength,
    kNuMIRelaxedProtonTrackChi2Muon,
    kNuMIRelaxedProtonTrackChi2Proton,
    kNuMIRelaxedProtonTrackScore,
    kNuMIRelaxedProtonTrackMatchedTruthPDG,
    kNuMIRelaxedProtonTrackMatchedTruthIntID,
    kNuMIRelaxedProtonTrackMatchedTruthContained,
    // ChargedPion
    kNuMIIsRelaxedChargedPionSelection,
    kNuMIRelaxedChargedPionTrackMatchType,
    kNuMIRelaxedChargedPionTrackLength,
    kNuMIRelaxedChargedPionTrackChi2Muon,
    kNuMIRelaxedChargedPionTrackChi2Proton,
    kNuMIRelaxedChargedPionTrackScore,
    kNuMIRelaxedChargedPionTrackMatchedTruthPDG,
    kNuMIRelaxedChargedPionTrackMatchedTruthIntID,
    kNuMIRelaxedChargedPionTrackMatchedTruthContained,
  };

  if(IsData){

    // Data

    // CV

    map_cutName_to_vec_Trees[currentCutName].push_back(
      new ana::Tree(
        ("PIDStudyTree_"+currentCutName).Data(), labels,
        loader,
        varlists,
        spillCut, cut,
        kNoShift,
        true, true
      )
    );

  }
  else{

    // MC

    // CV

    map_cutName_to_vec_Trees[currentCutName].push_back(
      new ana::Tree(
        ("PIDStudyTree_"+currentCutName).Data(), labels,
        loader,
        varlists,
        spillCut, cut,
        ApplyTrackSplit ? SystShifts(&kTrackSplittingSyst, +1.) : kNoShift,
        true, true
      )
    );

    // Shifts

    map_cutName_to_vec_Trees[currentCutName].push_back(
      new ana::Tree(
        ("PIDStudyTree_CalodEdXShiftUp_"+currentCutName).Data(), labels,
        loader,
        varlists,
        spillCut, cut,
        ApplyTrackSplit ? SystShifts( {{&kTrackSplittingSyst, +1.}, {&kCalodEdXShiftSyst, +1.}} ) : SystShifts(&kCalodEdXShiftSyst, +1.),
        true, true
      )
    );
    map_cutName_to_vec_Trees[currentCutName].push_back(
      new ana::Tree(
        ("PIDStudyTree_CalodEdXShiftDown_"+currentCutName).Data(), labels,
        loader,
        varlists,
        spillCut, cut,
        ApplyTrackSplit ? SystShifts( {{&kTrackSplittingSyst, +1.}, {&kCalodEdXShiftSyst, -1.}} ) : SystShifts(&kCalodEdXShiftSyst, -1.),
        true, true
      )
    );
    map_cutName_to_vec_Trees[currentCutName].push_back(
      new ana::Tree(
        ("PIDStudyTree_CaloGainShiftUp_"+currentCutName).Data(), labels,
        loader,
        varlists,
        spillCut, cut,
        ApplyTrackSplit ? SystShifts( {{&kTrackSplittingSyst, +1.}, {&kCaloGainShiftSyst, +1.}} ) : SystShifts(&kCaloGainShiftSyst, +1.),
        true, true
      )
    );
    map_cutName_to_vec_Trees[currentCutName].push_back(
      new ana::Tree(
        ("PIDStudyTree_CaloGainShiftDown_"+currentCutName).Data(), labels,
        loader,
        varlists,
        spillCut, cut,
        ApplyTrackSplit ? SystShifts( {{&kTrackSplittingSyst, +1.}, {&kCaloGainShiftSyst, -1.}} ) : SystShifts(&kCaloGainShiftSyst, -1.),
        true, true
      )
    );

    if(!FillSystematics) return;

    // NSigmas
    std::vector<std::string> this_NSigmasPsetNames;
    std::vector<const ISyst*> this_NSigmasISysts;
    std::vector<std::pair<int,int>> this_NSigmasPairs;
    std::vector<std::vector<double>> this_NSigmas;
    for(unsigned int i=0; i<genieMultisigmaKnobNames.size(); i++){
      this_NSigmasPsetNames.push_back( genieMultisigmaKnobNames.at(i) );
      this_NSigmasISysts.push_back( IGENIESysts.at(i) );
      this_NSigmasPairs.push_back( std::make_pair(-3, 3) );
      this_NSigmas.push_back( {-3, -2, -1, 0, 1, 2, 3} );
    }
    for(unsigned int i=0; i<genieMorphKnobNames.size(); i++){
      this_NSigmasPsetNames.push_back( genieMorphKnobNames.at(i) );
      this_NSigmasISysts.push_back( IGENIEMorphSysts.at(i) );
      this_NSigmasPairs.push_back( std::make_pair(0, 1) );
      this_NSigmas.push_back( {-1, -0.5, 0, 0.5, 1} );
    }
    for(unsigned int i=0; i<IFluxSysts.size(); i++){
      this_NSigmasPsetNames.push_back( IFluxSysts.at(i)->ShortName() );
      this_NSigmasISysts.push_back( IFluxSysts.at(i) );
      this_NSigmasPairs.push_back( std::make_pair(-3, 3) );
      this_NSigmas.push_back( {-3, -2, -1, 0, 1, 2, 3} );
    }
    for(unsigned int i=0; i<IDetectorSysts.size(); i++){
      this_NSigmasPsetNames.push_back( IDetectorSysts.at(i)->ShortName() );
      this_NSigmasISysts.push_back( IDetectorSysts.at(i) );
      this_NSigmasPairs.push_back( std::make_pair(-3, 3) );
      this_NSigmas.push_back( {-3, -2, -1, 0, 1, 2, 3} );
    }

    map_cutName_to_vec_NSigmasTrees[currentCutName].push_back(
      new ana::NSigmasTree(
        ("PIDStudyTree_NSigmas_"+currentCutName).Data(),
        this_NSigmasPsetNames,
        loader,
        this_NSigmasISysts,
        this_NSigmas,
        spillCut, cut,
        ApplyTrackSplit ? SystShifts(&kTrackSplittingSyst, +1.) : kNoShift,
        true, true
      )
    );

    // NUniverses
    std::vector<std::string> this_NUniversesPsetNames;
    std::vector<std::vector<Var>> this_NUniversesVarVectors;
    std::vector<std::vector<TruthVar>> this_NUniversesTruthVarVectors;
    std::vector<unsigned int> this_NUniversesNUnivs;
    for(const std::string& name: genieDependentKnobNames){
      this_NUniversesPsetNames.push_back( name );
      this_NUniversesVarVectors.push_back( map_DepDialName_to_UniverseWeights[name] );
      this_NUniversesTruthVarVectors.push_back( map_DepDialName_to_TruthUniverseWeights[name] );
      this_NUniversesNUnivs.push_back( 100 );
    }
    for(const std::string& name: geant4DependentKnobNames){
      this_NUniversesPsetNames.push_back( name );
      this_NUniversesVarVectors.push_back( map_DepDialName_to_UniverseWeights[name] );
      this_NUniversesTruthVarVectors.push_back( map_DepDialName_to_TruthUniverseWeights[name] );
      this_NUniversesNUnivs.push_back( 1000 );
    }
    map_cutName_to_vec_NUniversesTrees[currentCutName].push_back(
      new ana::NUniversesTree(
        ("PIDStudyTree_NUniverses_"+currentCutName).Data(),
        this_NUniversesPsetNames,
        loader,
        this_NUniversesVarVectors,
        this_NUniversesNUnivs,
        spillCut, cut,
        ApplyTrackSplit ? SystShifts(&kTrackSplittingSyst, +1.) : kNoShift,
        true, true
      )
    );

  }


}

// - 230918_DetSyst
void HistoProducer::MakeDetSystStudyTree(SpectrumLoader& loader, SpillCut spillCut, Cut cut){

  std::vector<std::string> labels = {
    "dedx",
    "dqdx",
    "integral",
    "sumadc",
    "t",
  };

  std::vector<MultiVar> varlists = {
    Hit_dedx,
    Hit_dqdx,
    Hit_integral,
    Hit_sumadc,
    Hit_t,
  };

  map_cutName_to_vec_Trees[currentCutName].push_back(
    new ana::Tree(
      ("DetStudyStudyTree_"+currentCutName).Data(), labels,
      loader,
      varlists,
      spillCut, cut,
      kNoShift, true, true
    )
  );

}

// - 231015_CutFlow
void HistoProducer::MakeCutFlow(SpectrumLoader& loader, SpillCut spillCut, Cut cut){

  if(currentCutName!="AllSamples_NoCut"){

    map_cutName_to_vec_Spectrums[currentCutName].push_back(
      new Spectrum(
        "MuonP", Binning::Simple(300, 0., 3.), loader,
        kTruth_MuonP,
        kTruthCut_IsSignal, kNuMIValidTrigger, cut, kNoShift, kGetTruthNuMIFluxWeight
      )
    );
    map_cutName_to_vec_Spectrums[currentCutName].push_back(
      new Spectrum(
        "MuonP_WithoutPhaseSpaceCut", Binning::Simple(300, 0., 3.), loader,
        kTruth_MuonP,
        kTruthCut_IsSignalWithoutPhaseSpaceCut, kNuMIValidTrigger, cut, kNoShift, kGetTruthNuMIFluxWeight
      )
    );

    map_cutName_to_vec_Spectrums[currentCutName].push_back(
      new Spectrum(
        "ProtonP", Binning::Simple(200, 0., 2.), loader,
        kTruth_ProtonP,
        kTruthCut_IsSignal, kNuMIValidTrigger, cut, kNoShift, kGetTruthNuMIFluxWeight
      )
    );
    map_cutName_to_vec_Spectrums[currentCutName].push_back(
      new Spectrum(
        "ProtonP_WithoutPhaseSpaceCut", Binning::Simple(200, 0., 2.), loader,
        kTruth_ProtonP,
        kTruthCut_IsSignalWithoutPhaseSpaceCut, kNuMIValidTrigger, cut, kNoShift, kGetTruthNuMIFluxWeight
      )
    );

    map_cutName_to_vec_Spectrums[currentCutName].push_back(
      new Spectrum(
        "deltaPT", Binning::Simple(120, 0., 1.2), loader,
        kTruth_deltaPT,
        kTruthCut_IsSignal, kNuMIValidTrigger, cut, kNoShift, kGetTruthNuMIFluxWeight
      )
    );

  }

  else{

    map_cutName_to_vec_Spectrums[currentCutName].push_back(
      new Spectrum(
        "MuonP", Binning::Simple(300, 0., 3.), loader,
        kTruth_MuonP,
        kTruthCut_IsSignal, kNuMIValidTrigger, kNoShift, kGetTruthNuMIFluxWeight
      )
    );
    map_cutName_to_vec_Spectrums[currentCutName].push_back(
      new Spectrum(
        "MuonP_WithoutPhaseSpaceCut", Binning::Simple(300, 0., 3.), loader,
        kTruth_MuonP,
        kTruthCut_IsSignalWithoutPhaseSpaceCut, kNuMIValidTrigger, kNoShift, kGetTruthNuMIFluxWeight
      )
    );

    map_cutName_to_vec_Spectrums[currentCutName].push_back(
      new Spectrum(
        "ProtonP", Binning::Simple(200, 0., 2.), loader,
        kTruth_ProtonP,
        kTruthCut_IsSignal, kNuMIValidTrigger, kNoShift, kGetTruthNuMIFluxWeight
      )
    );
    map_cutName_to_vec_Spectrums[currentCutName].push_back(
      new Spectrum(
        "ProtonP_WithoutPhaseSpaceCut", Binning::Simple(200, 0., 2.), loader,
        kTruth_ProtonP,
        kTruthCut_IsSignalWithoutPhaseSpaceCut, kNuMIValidTrigger, kNoShift, kGetTruthNuMIFluxWeight
      )
    );

    map_cutName_to_vec_Spectrums[currentCutName].push_back(
      new Spectrum(
        "deltaPT", Binning::Simple(120, 0., 1.2), loader,
        kTruth_deltaPT,
        kTruthCut_IsSignal, kNuMIValidTrigger, kNoShift, kGetTruthNuMIFluxWeight
      )
    );

  }


}
void HistoProducer::MakeCutFlowTree(SpectrumLoader& loader, SpillCut spillCut, Cut cut){

  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "OpFlashFirstTime", Binning::Simple(500, -20., 30.),
      loader,
      OpFlashFirstTime,
      kNoSpillCut
    )
  );

  std::vector<std::string> true_labels = {
    // Weight
    "FluxWeight",
    "FluxWeightWithG3Chase",
    "FluxWeightWithG4Updated",
    // Nu E
    "TrueE",
    // Muon
    "TrueMuonP",
    "TrueMuonCos",
    "TrueMuonLength",
    // Proton
    "TrueProtonP",
    "TrueProtonCos",
    "TrueProtonLength",
    // Muon+Proton
    "TrueMuonProtonCos",
    // TKI
    "TruedeltaPT",
    "TruedeltaPTx",
    "TruedeltaPTy",
    "TruedeltaalphaT",
    "TruedeltaphiT",
  };

  std::vector<TruthVar> true_vars = {
    // Weight
    kGetTruthNuMIFluxWeight,
    kGetTruthNuMIFluxWeightG3Chase,
    kGetTruthNuMIFluxWeightUpdated,
    // Nu E
    kTruth_NeutrinoE,
    // Muon
    kTruth_MuonP,
    kTruth_MuonNuCosineTheta,
    kTruth_MuonLength,
    // Proton
    kTruth_ProtonP,
    kTruth_ProtonNuCosineTheta,
    kTruth_ProtonLength,
    // Muon+Proton
    kTruth_CosThMuonProton,
    // TKI
    kTruth_deltaPT,
    kTruth_deltaPTx,
    kTruth_deltaPTy,
    kTruth_deltaalphaT,
    kTruth_deltaphiT,
  };

  map_cutName_to_vec_Trees[currentCutName].push_back(
    new ana::Tree(
      ("trueEvents_WithPhaseSpaceCut_"+currentCutName).Data(),
      true_labels, loader, true_vars, kNuMIValidTrigger, kTruthCut_IsSignal,
      cut,
      ApplyTrackSplit ? SystShifts(&kTrackSplittingSyst, +1.) : kNoShift,
      true
    )
  );

  map_cutName_to_vec_Trees[currentCutName].push_back(
    new ana::Tree(
      ("trueEvents_WithoutPhaseSpaceCut_"+currentCutName).Data(),
      true_labels, loader, true_vars, kNuMIValidTrigger, kTruthCut_IsSignalWithoutPhaseSpaceCut,
      cut,
      ApplyTrackSplit ? SystShifts(&kTrackSplittingSyst, +1.) : kNoShift,
      true
    )
  );

  std::cout << "currentCutName = " << currentCutName << std::endl;
  if( currentCutName=="AllSamples_NoCut" ){

    std::vector<Var> reco_vars;
    reco_vars.push_back(kGetNuMIFluxWeightG3Chase);
    reco_vars.push_back(kNuMISPPQ2RW);
    reco_vars.push_back(kNuMISPPTpiMINERvAFittedReweight);
    reco_vars.push_back(kNuMISplitTrackCVCorrection);
    reco_vars.push_back(kNuMISliceSignalType);
    reco_vars.push_back(kNuMISliceSignalTypeWithoutOOPS);

    std::vector<std::string> reco_labels;
    reco_labels.push_back( "FluxWeightWithG3Chase" );
    reco_labels.push_back( "SPPQ2RW" );
    reco_labels.push_back( "SPPTpiMINERvAFittedReweight" );
    reco_labels.push_back( "TrackSplitRW" );
    reco_labels.push_back( "IsSignal/I" );
    reco_labels.push_back( "IsSignalWithoutOOPS/I" );


    reco_vars.push_back( Var([](const caf::SRSliceProxy* slc) -> int {return kNuMIVertexInFV(slc) ? 1 : 0;}) );
    reco_vars.push_back( Var([](const caf::SRSliceProxy* slc) -> int {return kNuMINotClearCosmic(slc) ? 1 : 0;}) );
    reco_vars.push_back( Var([](const caf::SRSliceProxy* slc) -> int {return kNuMIHasMuonCandidate(slc) ? 1 : 0;}) );
    reco_vars.push_back( Var([](const caf::SRSliceProxy* slc) -> int {return kNuMIHasProtonCandidate(slc) ? 1 : 0;}) );
    reco_vars.push_back( Var([](const caf::SRSliceProxy* slc) -> int {return kNuMIProtonCandidateRecoPTreshold(slc) ? 1 : 0;}) );
    reco_vars.push_back( Var([](const caf::SRSliceProxy* slc) -> int {return kNuMIAllPrimaryHadronsContained(slc) ? 1 : 0;}) );
    reco_vars.push_back( Var([](const caf::SRSliceProxy* slc) -> int {return kNuMINoSecondPrimaryMuonlikeTracks(slc) ? 1 : 0;}) );
    reco_vars.push_back( Var([](const caf::SRSliceProxy* slc) -> int {return kNuMICutPhotons(slc) ? 1 : 0;}) );
    reco_vars.push_back( Var([](const caf::SRSliceProxy* slc) -> int {return kNuMIMuonCandidateContained(slc) ? 1 : 0;}) );

    reco_labels.push_back( "Pass_VtxInFV/I" );
    reco_labels.push_back( "Pass_NotClearCosmic/I" );
    reco_labels.push_back( "Pass_HasMuon/I" );
    reco_labels.push_back( "Pass_HasProton/I" );
    reco_labels.push_back( "Pass_ProtonPCut/I" );
    reco_labels.push_back( "Pass_PrimaryHadronContained/I" );
    reco_labels.push_back( "Pass_NoChargedPionTrack/I" );
    reco_labels.push_back( "Pass_NoNeutralPionShower/I" );
    reco_labels.push_back( "Pass_MuonContained/I" );

    map_cutName_to_vec_Trees[currentCutName].push_back(
      new ana::Tree(
        ("selectedEvents_"+currentCutName).Data(), reco_labels,
        loader,
        reco_vars,
        spillCut, cut,
        ApplyTrackSplit ? SystShifts(&kTrackSplittingSyst, +1.) : kNoShift,
        true, true
      )
    );

  }

}

// - 231217_CountNuMINu
void HistoProducer::MakeNuMINuCountTree(SpectrumLoader& loader, SpillCut spillCut, Cut cut){

  std::vector<std::string> labels = {
    "IsAV/I",
    "IsFV/I",
    "NuPDG/I",
    "NuMode/I",
    "NuCC/I",
    "NuE",
    "MuonP",
    "MuonLength",
    // Weight
    "FluxWeight",
    "FluxWeightWithG3Chase",
    "FluxWeightWithG4Updated",
    // SPP RW for res events
    "IsSPP/i",
    "SPPCVCorrection",
  };

  std::vector<TruthVar> varlists = {
    kTruth_IsVertexInAV,
    kTruth_IsVertexInFV,
    kTruth_NeutrinoPDG,
    kTruth_NeutrinoMode,
    kTruth_IsCC,
    kTruth_NeutrinoE,
    kTruth_MuonP,
    kTruth_MuonLength,
    // Weight
    kGetTruthNuMIFluxWeight,
    kGetTruthNuMIFluxWeightG3Chase,
    kGetTruthNuMIFluxWeightUpdated,
    // SPP RW for res events
    kTruth_IsSPP,
    kTruth_NuMISPPCVCorrection,
  };

  map_cutName_to_vec_Trees[currentCutName].push_back(
    new ana::Tree(
      "NuMINuCount",
      labels,
      loader,
      varlists,
      kNoSpillCut,
      kNoTruthCut,
      kNoCut,
      kNoShift,
      true
    )
  );

}
// - 240307_BNBFlux
void HistoProducer::MakeBNBFluxTree(SpectrumLoader& loader, SpillCut spillCut, Cut cut){

  std::vector<std::string> truthvar_labels = {
    "TrueE",
    "TrueIsFHC/i",
    "TruePDG/i",
    "Weight",
  };
  std::vector<TruthVar> truthvar_vars = {
    kTruth_NeutrinoE,
    kTruth_IsFHC,
    kTruth_NeutrinoPDG,
    kTruth_BNBDefaultWeight,
  };

  map_cutName_to_vec_Trees[currentCutName].push_back(
    new ana::Tree(
      "trueEvents",
      truthvar_labels, loader, truthvar_vars,
      kNoSpillCut, kNoTruthCut,
      kNoCut,
      kNoShift,
      true
    )
  );

  std::vector<std::string> systNames = {
"expskin",
"horncurrent",
"kminus",
"kplus",
"kzero",
"nucleoninexsec",
"nucleonqexsec",
"nucleontotxsec",
"piminus",
"pioninexsec",
"piontotxsec",
  };
  std::vector<std::vector<TruthVar>> systTruthVarVectors;
  std::vector<unsigned int> systNUnivs;

  for(auto& systName: systNames){
    std::vector<TruthVar> vec_VarVector;
    for(int u=0; u<1000; u++){
      std::string psetname = systName+"_Flux";
      vec_VarVector.push_back( GetTruthUniverseWeight(psetname, u) );
    }
    systTruthVarVectors.push_back( vec_VarVector );
    systNUnivs.push_back( 1000 );
  }


  map_cutName_to_vec_NUniversesTrees[currentCutName].push_back(
    new ana::NUniversesTree(
      "trueEvents_NUniverses",
      systNames,
      loader,
      systTruthVarVectors,
      systNUnivs,
      kNoTruthCut,
      kNoShift, true
    )
  );


}

// - 240315_FSICovTree
void HistoProducer::MakeFSICovTree(SpectrumLoader& loader, SpillCut spillCut, Cut cut){

  std::vector<std::string> truthvar_labels = {
    // Weight
    "FluxWeight",
    "FluxWeightWithG3Chase",
    "FluxWeightWithG4Updated",
    "SPPCVCorrection",
    // Muon
    "TrueMuonCos",
    // Proton
    "TrueProtonCos",
    // Muon+Proton
    "TrueMuonProtonCos",
    // TKI
    "TruedeltaPT",
    "TruedeltaPTx",
    "TruedeltaPTy",
    "TruedeltaalphaT",
    "TruedeltaphiT",
  };
  std::vector<TruthVar> truthvar_vars = {
    // Weight
    kGetTruthNuMIFluxWeight,
    kGetTruthNuMIFluxWeightG3Chase,
    kGetTruthNuMIFluxWeightUpdated,
    kTruth_NuMISPPCVCorrection,
    // Muon
    kTruth_MuonNuCosineTheta,
    // Proton
    kTruth_ProtonNuCosineTheta,
    // Muon+Proton
    kTruth_CosThMuonProton,
    // TKI
    kTruth_deltaPT,
    kTruth_deltaPTx,
    kTruth_deltaPTy,
    kTruth_deltaalphaT,
    kTruth_deltaphiT,
  };

  map_cutName_to_vec_Trees[currentCutName].push_back(
    new ana::Tree(
      "trueEvents",
      truthvar_labels, loader, truthvar_vars,
      kNoSpillCut, kNuMITrueNuMuCCInFV,
      kNoCut,
      kNoShift,
      true
    )
  );

  // Multisigma

  std::vector<std::string> MultisigmaSystNames = {
    "MFP_pi",
    "FrCEx_pi",
    "FrInel_pi",
    "FrAbs_pi",
    "FrPiProd_pi",
    "MFP_N",
    "FrCEx_N",
    "FrInel_N",
    "FrAbs_N",
    "FrPiProd_N",
  };

  std::vector<const ISyst*> this_NSigmasISysts;
  std::vector<std::vector<double>> this_NSigmas;
  for(auto& name: MultisigmaSystNames){
    std::string psetname = SystProviderPrefix+"_multisigma_"+name;
    this_NSigmasISysts.push_back( new SBNWeightSyst(psetname) );
    this_NSigmas.push_back( {-3, -2, -1, 0, 1, 2, 3} );
  }

  map_cutName_to_vec_NSigmasTrees[currentCutName].push_back(
    new ana::NSigmasTree(
      "trueEvents_NSigmas",
      MultisigmaSystNames,
      loader,
      this_NSigmasISysts,
      this_NSigmas,
      kNuMITrueNuMuCCInFV,
      kNoShift, true
    )
  );

  // Multisim

  std::vector<std::string> MultisimSystNames = {
    "FSI_pi_VariationResponse",
    "FSI_N_VariationResponse",
  };
  std::vector<std::vector<TruthVar>> MultisimSystTruthVarVectors;
  std::vector<unsigned int> MultisimSystNUnivs;

  for(auto& systName: MultisimSystNames){
    std::vector<TruthVar> vec_VarVector;
    for(int u=0; u<100; u++){
      std::string psetname = SystProviderPrefix+"_multisim_"+systName;
      vec_VarVector.push_back( GetTruthUniverseWeight(psetname, u) );
    }
    MultisimSystTruthVarVectors.push_back( vec_VarVector );
    MultisimSystNUnivs.push_back( 100 );
  }


  map_cutName_to_vec_NUniversesTrees[currentCutName].push_back(
    new ana::NUniversesTree(
      "trueEvents_NUniverses",
      MultisimSystNames,
      loader,
      MultisimSystTruthVarVectors,
      MultisimSystNUnivs,
      kNuMITrueNuMuCCInFV,
      kNoShift, true
    )
  );

}

// - 240712_DupeCheck
void HistoProducer::MakeDupeCheckTree(SpectrumLoader& loader, SpillCut spillCut, Cut cut){

  std::vector<std::string> this_reco_labels;
  std::vector<Var> this_reco_vars;

  this_reco_labels.push_back( "Dummy/i" );
  this_reco_vars.push_back(
    Var([](const caf::SRSliceProxy* slc) -> int {
      return 1;
    })
  );

  map_cutName_to_vec_Trees[currentCutName].push_back(
    new ana::Tree(
      "selectedEvents", this_reco_labels,
      loader,
      this_reco_vars,
      spillCut, cut,
      kNoShift, true, true
    )
  );


}

// - 240805_EventListTree
void HistoProducer::MakeEventListTree(SpectrumLoader& loader, SpillCut spillCut, Cut cut){

  std::vector<std::string> this_labels;
  std::vector<SpillVar> this_vars;

  this_labels.push_back( "Dummy/i" );
  this_vars.push_back(
    SpillVar([](const caf::SRSpillProxy* slc) -> int {
      return 1;
    })
  );

  map_cutName_to_vec_Trees[currentCutName].push_back(
    new ana::Tree(
      "EvenvList", this_labels,
      loader,
      this_vars,
      spillCut,
      true
    )
  );


}

// - 241021_BeamQualTree
void HistoProducer::MakeBeamQualTree(SpectrumLoader& loader, SpillCut spillCut, Cut cut){


  // All 
  std::vector<std::string> SpillMultiVar_labels = {
    "TRTGTDAll",
    "TR101DAll",
    "HornCurrentAll",
    "POTInSpillAll",
    "BeamHPTGT",
    "BeamVPTGT",
    "BeamPosHAll",
    "BeamPosVAll",
    "BeamWidthHAll",
    "BeamWidthVAll",
    "EventAll/I",
  };
  std::vector<SpillMultiVar> SpillMultiVar_vars = {
    kTRTGTDAll,
    kTR101DAll,
    kHornCurrentAll,
    kPOTInSpillAll,
    kBeamHPTGT,
    kBeamVPTGT,
    kBeamPosHAll,
    kBeamPosVAll,
    kBeamWidthHAll,
    kBeamWidthVAll,
    kEventAll,
  };


/*
  std::vector<std::string> SpillMultiVar_labels = {
"Test",
  };
  std::vector<SpillMultiVar> SpillMultiVar_vars = {
spillvarTest
  };
*/

  map_cutName_to_vec_Trees[currentCutName].push_back(
    new ana::Tree(
      "BeamQualTreeAll",
      SpillMultiVar_labels,
      loader,
      SpillMultiVar_vars,
      spillCut,
      true
    )
  );


}

// - 241021_InTimeCosmicOverlapTree
void HistoProducer::MakeInTimeCosmicOverlapTree(SpectrumLoader& loader, SpillCut spillCut, Cut cut){

  std::vector<std::string> SpillVar_labels = {
"IntimeNeutrinoTime",
"HasIntimeCosmicParticle/I",
"HasIntimeCosmicParticleWithCut/I",
"IntimeCosmicParticleLongestLength",
  };
  std::vector<SpillVar> SpillVar_vars = {
kNuMIIntimeNeutrinoTime,
kNuMIHasIntimeCosmicParticle,
kNuMIHasIntimeCosmicParticleWithCut,
kNuMIIntimeCosmicParticleLongestLength,
  };


  map_cutName_to_vec_Trees[currentCutName].push_back(
    new ana::Tree(
      "InTimeCosmicOverlapTree",
      SpillVar_labels,
      loader,
      SpillVar_vars,
      spillCut,
      true
    )
  );

}

// - 241031_LifetimeVariationTree
void HistoProducer::MakeLifetimeVariationTree(SpectrumLoader& loader, SpillCut spillCut, Cut cut){

  std::vector<std::string> Var_labels = {
"MuonMatchedTrack_TrackScore",
"MuonMatchedTrack_Chi2Muon",
"MuonMatchedTrack_Chi2Proton",
"ProtonMatchedTrack_TrackScore",
"ProtonMatchedTrack_Chi2Muon",
"ProtonMatchedTrack_Chi2Proton",
  };
  std::vector<Var> Var_vars = {
kNuMI_MuonMatchedTrack_TrackScore,
kNuMI_MuonMatchedTrack_Chi2Muon,
kNuMI_MuonMatchedTrack_Chi2Proton,
kNuMI_ProtonMatchedTrack_TrackScore,
kNuMI_ProtonMatchedTrack_Chi2Muon,
kNuMI_ProtonMatchedTrack_Chi2Proton,
  };


  map_cutName_to_vec_Trees[currentCutName].push_back(
    new ana::Tree(
      "LifetimeVariationTree",
      Var_labels,
      loader,
      Var_vars,
      spillCut, cut,
      ApplyTrackSplit ? SystShifts(&kTrackSplittingSyst, +1.) : kNoShift,
      true, true
    )
  );

}

void HistoProducer::Test(SpectrumLoader& loader, SpillCut spillCut, Cut cut){

  // base variables

  std::vector<std::string> this_reco_labels = GetNuMIRecoTreeLabels();
  std::vector<Var> this_reco_vars = GetNuMIRecoTreeVars();

    // MC

    this_reco_labels.push_back( "IsData/i" );
    this_reco_vars.push_back(
      Var([](const caf::SRSliceProxy* slc) -> int {
        return 0;
      })
    );

    // CV

    map_cutName_to_vec_Trees[currentCutName].push_back(
      new ana::Tree(
        "selectedEvents", this_reco_labels,
        loader,
        this_reco_vars,
        spillCut, cut,
        //kNoShift,
        SystShifts(&kTrackSplittingSyst, +1.),
        true, true
      )
    );


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
      vector<ana::NUniversesTree *> vec_NUniversesTrees = map_cutName_to_vec_NUniversesTrees[cutName];

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

      // If no Tree, continue;
      if(vec_Trees.size()==0){
        continue;
      }

      if(MakeGUNDAMTree){

        for(unsigned int i=0; i<vec_Trees.size(); i++){
          if(vec_Trees.at(i)->Name()=="selectedEvents"){
            cout << "[HistoProducer::saveHistograms]   Skipping Tree:" << vec_Trees.at(i)->Name() << "; will be merged to NSigmasTree" << std::endl;
            continue;
          }
          if(vec_Trees.at(i)->Name()=="trueEvents"){
            cout << "[HistoProducer::saveHistograms]   Skipping Tree:" << vec_Trees.at(i)->Name() << "; will be merged to NSigmasTree" << std::endl;
            continue;
          }
          cout << "[HistoProducer::saveHistograms]   Writing Tree:" << vec_Trees.at(i)->Name() << std::endl;
          vec_Trees.at(i)->SaveTo(dir);
        }
        for(unsigned int i=0; i<vec_NSigmasTrees.size(); i++){
          if(vec_NSigmasTrees.at(i)->Name()=="selectedEvents_NSigmas"){
            for(unsigned int j=0; j<vec_Trees.size(); j++){
              if(vec_Trees.at(j)->Name()=="selectedEvents"){
                cout << "[HistoProducer::saveHistograms]   Merging " << vec_Trees.at(j)->Name() << " to " << vec_NSigmasTrees.at(i)->Name() << std::endl;
                vec_NSigmasTrees.at(i)->MergeTree( *vec_Trees.at(j) );
                break;
              }
            }
          }
          if(vec_NSigmasTrees.at(i)->Name()=="trueEvents_NSigmas"){
            for(unsigned int j=0; j<vec_Trees.size(); j++){
              if(vec_Trees.at(j)->Name()=="trueEvents"){
                cout << "[HistoProducer::saveHistograms]   Merging " << vec_Trees.at(j)->Name() << " to " << vec_NSigmasTrees.at(i)->Name() << std::endl;
                vec_NSigmasTrees.at(i)->MergeTree( *vec_Trees.at(j) );
                break;
              }
            }
          }
          cout << "[HistoProducer::saveHistograms]   Writing NSigmasTree: " << vec_NSigmasTrees.at(i)->Name() << " (save mode = " << nSigmasSaveMode << ")" << std::endl;
          if(nSigmasSaveMode==kVector) vec_NSigmasTrees.at(i)->SaveTo(dir);
          else if(nSigmasSaveMode==kSpline) vec_NSigmasTrees.at(i)->SaveToSplines(dir);
          else if(nSigmasSaveMode==kGraph) vec_NSigmasTrees.at(i)->SaveToGraphs(dir);
          else if(nSigmasSaveMode==kTClonesArrays) vec_NSigmasTrees.at(i)->SaveToTClonesArrays(dir);
        }
        for(unsigned int i=0; i<vec_NUniversesTrees.size(); i++){
          cout << "[HistoProducer::saveHistograms]   Writing UniversesTree: " << vec_NUniversesTrees.at(i)->Name() << std::endl;
          vec_NUniversesTrees.at(i)->SaveTo(dir);
        }

      }
      else{

        for(unsigned int i=0; i<vec_Trees.size(); i++){
           cout << "[HistoProducer::saveHistograms]   Writing Tree:" << vec_Trees.at(i)->Name() << std::endl;
          vec_Trees.at(i)->SaveTo(dir);
        }
        for(unsigned int i=0; i<vec_NSigmasTrees.size(); i++){
          cout << "[HistoProducer::saveHistograms]   Writing NSigmasTree: " << vec_NSigmasTrees.at(i)->Name() << " (save mode = " << nSigmasSaveMode << ")" << std::endl;
          if(nSigmasSaveMode==kVector) vec_NSigmasTrees.at(i)->SaveTo(dir);
          else if(nSigmasSaveMode==kSpline) vec_NSigmasTrees.at(i)->SaveToSplines(dir);
          else if(nSigmasSaveMode==kGraph) vec_NSigmasTrees.at(i)->SaveToGraphs(dir);
          else if(nSigmasSaveMode==kTClonesArrays) vec_NSigmasTrees.at(i)->SaveToTClonesArrays(dir);
        }
        for(unsigned int i=0; i<vec_NUniversesTrees.size(); i++){
          cout <<  "[HistoProducer::saveHistograms]   Writing UniversesTree: " << vec_NUniversesTrees.at(i)->Name() << std::endl;
          vec_NUniversesTrees.at(i)->SaveTo(dir);
        }

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
    genieMultisigmaKnobNames = ICARUSNumuXsec::GetGENIEMultisigmaKnobNames();
    for(const std::string& name: genieMultisigmaKnobNames){
      if(name=="FormZone"){
       std::cout << "[HistoProducer::setSystematicWeights] Skipping FormZone" << std::endl;
       continue;
      }
      std::string psetname = SystProviderPrefix+"_multisigma_"+name;
      std::cout << "[HistoProducer::setSystematicWeights] Multisigma, " << name << " (psetname = " << psetname << ")" << std::endl;
      IGENIESysts.push_back( new SBNWeightSyst(psetname) );
    }
    // Adding custom
    // 1) pi syst
    genieMultisigmaKnobNames.push_back( "CC1piTPi" );
    IGENIESysts.push_back( new NuMIXSecPiSyst("CC1piTPi", "CC1piTpi") );
    genieMultisigmaKnobNames.push_back( "LowQ2Suppression" );
    IGENIESysts.push_back( new NuMIXSecLowQ2Suppression("LowQ2Suppression", "LowQ2Suppression") );
    // 2) nusyst
    if(FillNuSyst){

      // PCAed z-exp parameters
      genieMultisigmaKnobNames.push_back( "ZExpPCAB1" );
      IGENIESysts.push_back( new SBNWeightSyst("ZExpPCAWeighter_SBNNuSyst_ZExpPCA_multisigma_b1") );
      genieMultisigmaKnobNames.push_back( "ZExpPCAB2" );
      IGENIESysts.push_back( new SBNWeightSyst("ZExpPCAWeighter_SBNNuSyst_ZExpPCA_multisigma_b2") );
      genieMultisigmaKnobNames.push_back( "ZExpPCAB3" );
      IGENIESysts.push_back( new SBNWeightSyst("ZExpPCAWeighter_SBNNuSyst_ZExpPCA_multisigma_b3") );
      genieMultisigmaKnobNames.push_back( "ZExpPCAB4" );
      IGENIESysts.push_back( new SBNWeightSyst("ZExpPCAWeighter_SBNNuSyst_ZExpPCA_multisigma_b4") );

      // Emiss

      std::string EmissPsetNamePrefix = "DIRT2_Emiss_SBNNuSyst_EMiss_multisigma_";

      std::vector<std::string> Emiss_tgts = {"Ar", "C"};
      std::vector<std::string> Emiss_nucs = {"p", "n"};

      for(const auto& Emiss_tgt: Emiss_tgts){
        for(const auto& Emis_nuc: Emiss_nucs){

          genieMultisigmaKnobNames.push_back( "Emiss_CorrTail_"+Emiss_tgt+"_"+Emis_nuc );
          IGENIESysts.push_back(
            new SBNWeightAbsVarSyst(
              EmissPsetNamePrefix+"Emiss_CorrTail_"+Emiss_tgt+"_"+Emis_nuc,
              {
                std::make_pair(0.25, -2.0),
                std::make_pair(0.5, -1.0),
                std::make_pair(1., 0.),
                std::make_pair(2., +1.0),
                std::make_pair(4., +2.0),
              }
            )
          );

          genieMultisigmaKnobNames.push_back( "Emiss_Linear_"+Emiss_tgt+"_"+Emis_nuc );
          IGENIESysts.push_back(
            new SBNWeightAbsVarSyst(
              EmissPsetNamePrefix+"Emiss_Linear_"+Emiss_tgt+"_"+Emis_nuc,
              { 
                std::make_pair(-30., -2.0),
                std::make_pair(-10., -1.0),
                std::make_pair(0., 0.),
                std::make_pair(10., +1.0),
                std::make_pair(30., +2.0),
              }
            )
          );

          genieMultisigmaKnobNames.push_back( "Emiss_ShiftPeak_"+Emiss_tgt+"_"+Emis_nuc );
          IGENIESysts.push_back(
            new SBNWeightAbsVarSyst(
              EmissPsetNamePrefix+"Emiss_ShiftPeak_"+Emiss_tgt+"_"+Emis_nuc,
              {
                std::make_pair(0., -2.0),
                std::make_pair(0.5, -1.0),
                std::make_pair(1., 0.),
                std::make_pair(1.5, +1.0),
                std::make_pair(2.0, +2.0),
              }
            )
          );


        } // END loop nucleon

      } // END loop target



    } // END if FillNuSyst


    // Morphs
    cout << "[HistoProducer::setSystematicWeights] - Adding morphing dials" << std::endl;
    genieMorphKnobNames = ICARUSNumuXsec::GetGENIEMorphKnobNames();
    for(const std::string& name: genieMorphKnobNames){
      if(name=="FormZone"){
       std::cout << "[HistoProducer::setSystematicWeights] Skipping FormZone" << std::endl;
       continue;
      }
      std::string psetname = SystProviderPrefix+"_multisigma_"+name;
      std::cout << "[HistoProducer::setSystematicWeights] Multisigma, " << name << " (psetname = " << psetname << ")" << std::endl;
      //IGENIEMorphSysts.push_back( new SBNWeightSyst(psetname) );
      IGENIEMorphSysts.push_back( new SBNWeightMirrorSyst(psetname) );
    }
    // Adding custom
    // 1) nusyst
    if(FillNuSyst){

      // Lars' 2p2h; independent ones

      std::vector<std::string> Lars2p2hDialNames = {
        "XSecShape_CCMEC",
        "XSecShape_CCMEC_Empirical",
        "XSecShape_CCMEC_Martini",
        "EnergyDependence_CCMEC",
      };

      for(const std::string& name: Lars2p2hDialNames){

        std::string dialname = "Lars2p2h_"+name;
        std::string psetname = "GENIEReWeight_SBNNuSyst_GENIE_multisigma_"+name;

        std::cout << "[HistoProducer::setSystematicWeights] Multisigma, " << name << " (psetname = " << psetname << ")" << std::endl;

        genieMorphKnobNames.push_back( dialname );
        IGENIEMorphSysts.push_back( new SBNWeightMirrorSyst(psetname) );

      }

      // Lars' 2p2h; dependent, dec ang mec

      genieMorphKnobNames.push_back( "Lars2p2p_DecayAngMEC_P1Variation_P2CV" );
      IGENIEMorphSysts.push_back( 
        new SBNWeightUnivSyst(
          "Lars2p2p_DecayAngMEC_P1Variation_P2CV",
          "GENIEReWeight_SBNNuSyst_GENIE_multisigma_DecayAngMECVariationResponse",
          {
            std::make_pair(20, -1.0),
            std::make_pair(10, -0.5),
            std::make_pair(0, 0.),
            std::make_pair(5, +0.5),
            std::make_pair(15, +1.0),
          }
        )
      );

      genieMorphKnobNames.push_back( "Lars2p2p_DecayAngMEC_P1Variation_P2p1sig" );
      IGENIEMorphSysts.push_back(
        new SBNWeightUnivSyst(
          "Lars2p2p_DecayAngMEC_P1Variation_P2p1sig",
          "GENIEReWeight_SBNNuSyst_GENIE_multisigma_DecayAngMECVariationResponse",
          { 
            std::make_pair(23, -1.0),
            std::make_pair(13, -0.5),
            std::make_pair(3, 0.),
            std::make_pair(8, +0.5),
            std::make_pair(18, +1.0),
          }
        )
      );

      genieMorphKnobNames.push_back( "Lars2p2p_DecayAngMEC_P1p1sig_P2Variation" );
      IGENIEMorphSysts.push_back( 
        new SBNWeightUnivSyst(
          "Lars2p2p_DecayAngMEC_P1p1sig_P2Variation",
          "GENIEReWeight_SBNNuSyst_GENIE_multisigma_DecayAngMECVariationResponse",
          {
            std::make_pair(19, -1.0),
            std::make_pair(17, -0.5),
            std::make_pair(15, 0.),
            std::make_pair(16, +0.5),
            std::make_pair(18, +1.0),
          }
        )
      );

/*
      // This is always 1..
      genieMorphKnobNames.push_back( "Lars2p2p_DecayAngMEC_P1CV_P2Variation" );
      IGENIEMorphSysts.push_back( 
        new SBNWeightUnivSyst(
          "Lars2p2p_DecayAngMEC_P1CV_P2Variation",
          "GENIEReWeight_SBNNuSyst_GENIE_multisigma_DecayAngMECVariationResponse",
          {
            std::make_pair(4, -1.0),
            std::make_pair(2, -0.5),
            std::make_pair(0, 0.),
            std::make_pair(1, +0.5),
            std::make_pair(3, +1.0),
          }
        )
      );
*/

      // CCQE RPA

      genieMorphKnobNames.push_back( "CCQERPAReweight" );
      IGENIEMorphSysts.push_back( new SBNWeightMirrorSyst("CCQERPAReweight_SBNNuSyst_CCQERPAReweight_multisigma_CCQERPAReweight") );

      // FSI; hA-to-hN

      genieMorphKnobNames.push_back( "FSIReweight_hN" );
      IGENIEMorphSysts.push_back( new SBNWeightMirrorSyst("FSIReweight_SBNNuSyst_FSI_hNReweight_multisigma_FSIReweight") );

      // FSI; hA-to-INCL

      genieMorphKnobNames.push_back( "FSIReweight_INCL" );
      IGENIEMorphSysts.push_back( new SBNWeightMirrorSyst("FSIReweight_SBNNuSyst_FSI_INCLReweight_multisigma_FSIReweight") );

      // FSI; hA-to-G4BC

      genieMorphKnobNames.push_back( "FSIReweight_G4BC" );
      IGENIEMorphSysts.push_back( new SBNWeightMirrorSyst("FSIReweight_SBNNuSyst_FSI_G4BCReweight_multisigma_FSIReweight") );

    }

    // Multisim for dependent dials
    cout << "[HistoProducer::setSystematicWeights] - Adding dependent dials by multisim" << std::endl;
    genieDependentKnobNames = ICARUSNumuXsec::GetGENIEDependentKnobNames();
    for(const std::string& name: genieDependentKnobNames){
      map_DepDialName_to_UniverseWeights[name] = {};
      map_DepDialName_to_TruthUniverseWeights[name] = {};
      std::string psetname = SystProviderPrefix+"_multisim_"+name;
      std::cout << "[HistoProducer::setSystematicWeights] Dependent dial, " << name << " (psetname = " << psetname << ")" << std::endl;
      for(int u=0; u<100; u++){
        map_DepDialName_to_UniverseWeights[name].push_back( GetUniverseWeight(psetname, u) );
        map_DepDialName_to_TruthUniverseWeights[name].push_back( GetTruthUniverseWeight(psetname, u) );
      }
    }
  }

  if(FillGEANT4){
    // GEANT4
    geant4DependentKnobNames = ICARUSNumuXsec::GetGEANT4MultisimKnobNames();
    for(const std::string& name: geant4DependentKnobNames){
      map_DepDialName_to_UniverseWeights[name] = {};
      map_DepDialName_to_TruthUniverseWeights[name] = {};
      std::string psetname = name;
      std::cout << "[HistoProducer::setSystematicWeights] GEANT4 multisim dial, " << name << " (psetname = " << psetname << ")" << std::endl;
      for(int u=0; u<1000; u++){
        map_DepDialName_to_UniverseWeights[name].push_back( GetUniverseWeight(psetname, u) );
        map_DepDialName_to_TruthUniverseWeights[name].push_back( GetTruthUniverseWeight(psetname, u) );
      }
    }
  }

  if(FillFlux){
    cout << "[HistoProducer::setSystematicWeights] Setting flux systematics" << endl;
    IFluxSysts = GetAllNuMIFluxSysts(NNuMIFluxPCA);
    IFluxSysts.push_back( GetNuMIBeamShiftSyst() );
    IFluxSysts.push_back( new NuMIBeamG3ChaseSyst("numi_beam_G3Chase", "numi_beam_G3Chase") );

    for(unsigned int i=0; i<IFluxSysts.size(); i++){
      cout << "[HistoProducer::setSystematicWeights] Syst = " << IFluxSysts.at(i)->ShortName() << endl;
    }
  }

  cout << "[HistoProducer::setSystematicWeights] Setting detector systematics" << endl;
  IDetectorSysts.push_back(
    new NuMIXSecDetectorSysts(
      NuMIXSecDetectorSysts::kFrontIndPlaneGain,
      "NuMIXSecFrontIndPlaneGainSyst",
      "Front ind. plane gain #pm10%"
    )
  );
  IDetectorSysts.push_back(
    new NuMIXSecDetectorSysts(
    NuMIXSecDetectorSysts::kFrontIndPlaneNoise,
    "NuMIXSecFrontIndPlaneNoiseSyst",
    "Front ind. plane noise +10%"
    )
  );
  IDetectorSysts.push_back(
    new NuMIXSecDetectorSysts(
    NuMIXSecDetectorSysts::kFrontIndPlaneSignalShape,
    "NuMIXSecFrontIndPlaneSignalShapeSyst",
    "Front ind. plane signal shape"
    )
  );

  IDetectorSysts.push_back(
    new NuMIXSecDetectorSysts(
    NuMIXSecDetectorSysts::kFrontIndPlaneSignalShapeFitted,
    "NuMIXSecFrontIndPlaneSignalShapeFittedSyst",
    "Front ind. plane signal shape, fitted"
    )
  );
  IDetectorSysts.push_back(
    new NuMIXSecDetectorSysts(
    NuMIXSecDetectorSysts::kMiddleIndPlaneTransparency,
    "NuMIXSecMiddleIndPlaneTransparencySyst",
    "Middle ind. plane transparency"
    )
  );

  IDetectorSysts.push_back(
    new NuMIXSecDetectorSysts(
    NuMIXSecDetectorSysts::kCaloGain,
    "NuMIXSecCaloGainSyst",
    "Calo gain up down"
    )
  );

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

