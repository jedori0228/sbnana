#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_HistoProducer.h"


using namespace ICARUSNumuXsec;

HistoProducer::HistoProducer(){

  cout << "[HistoProducer::HistoProducer] called" << endl;
  TargetPOT = 6.0e20;
  str_TargetPOT = "6.0e20 POT";
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
  double xNuScoreMin = -0.2;
  double xNuScoreMax = 1.;
  double dxNuScore = 0.01;
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
  const Binning binsFMChargeQ = Binning::Simple(502, -2, 500.);
  const Binning binsXPosition = Binning::Simple(1000, -500., 500.);
  const Binning binsYPosition = Binning::Simple(400, -200., 200.);
  const Binning binsZPosition = Binning::Simple(2000, -1000., 1000.);
  const Binning binsFMLightPE = Binning::Simple(1002, -2, 1000.);
  const Binning binsChi2 = Binning::Simple(400, 0., 400);
  const Binning binsPDG = Binning::Simple(10000, -5000., 5000.);

  //==== HistAxis
  //====   Slice variables
  const HistAxis axVertexRecoX("VertexRecoX", binsXPosition, varVertexRecoX);
  const HistAxis axVertexRecoY("VertexRecoY", binsYPosition, varVertexRecoY);
  const HistAxis axVertexRecoZ("VertexRecoZ", binsZPosition, varVertexRecoZ);
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
  //====   NuScore vars
  const HistAxis axSliceNuNFinalStatePfos("SliceNuNFinalStatePfos", Binning::Simple(20,0.,20.), varSliceNuNFinalStatePfos);
  const HistAxis axSliceNuNHitsTotal("SliceNuNHitsTotal", Binning::Simple(500,0.,5000.), varSliceNuNHitsTotal);
  const HistAxis axSliceNuVertexY("SliceNuVertexY", binsYPosition, varSliceNuVertexY);
  const HistAxis axSliceNuWeightedDirZ("SliceNuWeightedDirZ", Binning::Simple(100,-1,1), varSliceNuWeightedDirZ);
  const HistAxis axSliceNuNSpacePointsInSphere("SliceNuNSpacePointsInSphere", Binning::Simple(300,0.,300.), varSliceNuNSpacePointsInSphere);
  const HistAxis axSliceNuEigenRatioInSphere("SliceNuEigenRatioInSphere", Binning::Simple(100,0.,1.), varSliceNuEigenRatioInSphere);
  const HistAxis axSliceCRLongestTrackDirY("SliceCRLongestTrackDirY", Binning::Simple(100,-1,1), varSliceCRLongestTrackDirY);
  const HistAxis axSliceCRLongestTrackDeflection("SliceCRLongestTrackDeflection", Binning::Simple(40,0.,2.), varSliceCRLongestTrackDeflection);
  const HistAxis axSliceCRFracHitsInLongestTrack("SliceCRFracHitsInLongestTrack", Binning::Simple(50,0.,1.), varSliceCRFracHitsInLongestTrack);
  const HistAxis axSliceCRNHitsMax("SliceCRNHitsMax", Binning::Simple(1000,0.,1000.), varSliceCRNHitsMax);
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
  const HistAxis axTruthMuonMatchedTrackScore("TruthMuonMatchedTrackScore", Binning::Simple(100,0.,1.), varTruthMuonMatchedTrackScore);
  const HistAxis axTruthMuonMatchedTrackVertexDistance("TruthMuonMatchedTrackVertexDistance", Binning::Simple(300,0.,30), varTruthMuonMatchedTrackVertexDistance);
  const HistAxis axTruthMuonMatchedShowerScore("TruthMuonMatchedShowerScore", Binning::Simple(100,0.,1.), varTruthMuonMatchedShowerScore);
  //====     Proton
  const HistAxis axProtonTruthP("ProtonTruthP", binsMomentum, varProtonTruthP);
  const HistAxis axProtonTruthCosineTheta("ProtonTruthCosineTheta", binsCosineTheta, varProtonTruthCosineTheta);
  const HistAxis axProtonTruthNuMICosineTheta("ProtonTruthNuMICosineTheta", binsCosineTheta, varProtonTruthNuMICosineTheta);
  const HistAxis axTruthProtonMatchedTrackChi2Proton("TruthProtonMatchedTrackChi2Proton", binsChi2, varTruthProtonMatchedTrackChi2Proton);
  const HistAxis axTruthProtonMatchedTrackChi2ProtonCollection("TruthProtonMatchedTrackChi2ProtonCollection", binsChi2, varTruthProtonMatchedTrackChi2ProtonCollection);
  const HistAxis axTruthProtonMatchedTrackScore("TruthProtonMatchedTrackScore", Binning::Simple(100,0.,1.), varTruthProtonMatchedTrackScore);
  const HistAxis axTruthProtonMatchedTrackVertexDistance("TruthProtonMatchedTrackVertexDistance", Binning::Simple(300,0.,30.), varTruthProtonMatchedTrackVertexDistance);
  const HistAxis axTruthProtonMatchedShowerScore("TruthProtonMatchedShowerScore", Binning::Simple(100,0.,1.), varTruthProtonMatchedShowerScore);
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
  const HistAxis axMuonRecoStartX("MuonRecoStartX", binsXPosition, varMuonRecoStartX);
  const HistAxis axMuonRecoStartY("MuonRecoStartY", binsYPosition, varMuonRecoStartY);
  const HistAxis axMuonRecoStartZ("MuonRecoStartZ", binsZPosition, varMuonRecoStartZ);
  const HistAxis axMuonRecoDirX("MuonRecoDirX", Binning::Simple(100,-1,1), varMuonRecoDirX);
  const HistAxis axMuonRecoDirY("MuonRecoDirY", Binning::Simple(100,-1,1), varMuonRecoDirY);
  const HistAxis axMuonRecoDirZ("MuonRecoDirZ", Binning::Simple(100,-1,1), varMuonRecoDirZ);
  const HistAxis axMuonRecoForceDownDirX("MuonRecoForceDownDirX", Binning::Simple(100,-1,1), varMuonRecoForceDownDirX);
  const HistAxis axMuonRecoForceDownDirY("MuonRecoForceDownDirY", Binning::Simple(100,-1,1), varMuonRecoForceDownDirY);
  const HistAxis axMuonRecoForceDownDirZ("MuonRecoForceDownDirZ", Binning::Simple(100,-1,1), varMuonRecoForceDownDirZ);
  //====     Proton
  const HistAxis axNProtonCandTrack("NProtonCandTrack", Binning::Simple(10, 0., 10.), varNProtonCandTrack);
  const HistAxis axNProtonCandMatched("NProtonCandMatched", Binning::Simple(10, 0., 10.), varNProtonCandMatched);
  const HistAxis axProtonRecoP("ProtonRecoP", binsMomentum, varProtonRecoP);
  const HistAxis axProtonRecoCosineTheta("ProtonRecoCosineTheta", binsCosineTheta, varProtonRecoCosineTheta);
  const HistAxis axProtonRecoNuMICosineTheta("ProtonRecoNuMICosineTheta", binsCosineTheta, varProtonRecoNuMICosineTheta);
  const HistAxis axProtonLength("ProtonLength", binsLength, varProtonLength);
  const HistAxis axProtonChi2Proton("ProtonChi2Proton", binsChi2, varProtonChi2Proton);
  const HistAxis axProtonBestmatchPDG("ProtonBestmatchPDG", binsPDG, varProtonBestmatchPDG);
  //====     Neutrino
  const HistAxis axNeutrinoCombinedEnergy("NeutrinoCombinedEnergy", binsEnergy, varNeutrinoCombinedEnergy);
  const HistAxis axNeutrinoQE("NeutrinoQE", binsEnergy, varNeutrinoQE);

  if(fillNominal){
    //====   Slice variables
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axVertexRecoX, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axVertexRecoY, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axVertexRecoZ, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axFMScore, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axFMTime, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axNuScore, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axCountSlice, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axSliceTrackNhitsPlane0, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axSliceTrackNhitsPlane1, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axSliceTrackNhitsPlane2, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axSliceShowerNhitsPlane0, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axSliceShowerNhitsPlane1, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axSliceShowerNhitsPlane2, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axSliceTrackChargePlane0, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axSliceTrackChargePlane1, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axSliceTrackChargePlane2, spillCut, cut));
    //====   NuScore vars
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axSliceNuNFinalStatePfos, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axSliceNuNHitsTotal, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axSliceNuVertexY, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axSliceNuWeightedDirZ, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axSliceNuNSpacePointsInSphere, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axSliceNuEigenRatioInSphere, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axSliceCRLongestTrackDirY, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axSliceCRLongestTrackDeflection, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axSliceCRFracHitsInLongestTrack, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axSliceCRNHitsMax, spillCut, cut));
    //====   NuScore vs Something
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axNuScore, axSliceCRLongestTrackDirY, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum("AllTrackStartPositionX", binsXPosition, loader, varAllTrackStartPositionX, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum("AllTrackStartPositionY", binsYPosition, loader, varAllTrackStartPositionY, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum("AllTrackStartPositionZ", binsZPosition, loader, varAllTrackStartPositionZ, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum("FMChargeQ", binsFMChargeQ, loader, varFMChargeQ, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum("FMChargeCenterX", binsXPosition, loader, varFMChargeCenterX, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum("FMChargeCenterY", binsYPosition, loader, varFMChargeCenterY, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum("FMChargeCenterZ", binsZPosition, loader, varFMChargeCenterZ, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum("FMLightCenterX", binsXPosition, loader, varFMLightCenterX, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum("FMLightCenterY", binsYPosition, loader, varFMLightCenterY, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum("FMLightCenterZ", binsZPosition, loader, varFMLightCenterZ, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum("FMLightPE", binsFMLightPE, loader, varFMLightPE, spillCut, cut));
    //====   Truth variables
    //====     Neutrinos
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axNeutrinoTruthE, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axTruthQ2, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axTruthq0_lab, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axTruthmodq_lab, spillCut, cut));
    //====     Muon
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axMuonTruthP, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axMuonTruthCosineTheta, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axMuonTruthNuMICosineTheta, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axTruthMuonMatchedTrackChi2Muon, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axTruthMuonMatchedTrackChi2Proton, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axTruthMuonMatchedTrackScore, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axTruthMuonMatchedTrackVertexDistance, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axTruthMuonMatchedShowerScore, spillCut, cut));
    //====     Proton
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axProtonTruthP, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axProtonTruthCosineTheta, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axProtonTruthNuMICosineTheta, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axTruthProtonMatchedTrackChi2Proton, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axTruthProtonMatchedTrackChi2ProtonCollection, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axTruthProtonMatchedTrackScore, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axTruthProtonMatchedTrackVertexDistance, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axTruthProtonMatchedShowerScore, spillCut, cut));
    //====     Muon+Proton
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axTruthMuonProtonCosineTheta, spillCut, cut));
    //====   Reco variables
    //====     Muon
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axMuonRecoP, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axMuonRecoCosineTheta, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axMuonRecoNuMICosineTheta, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axMuonLength, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axMuonChi2Muon, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axMuonChi2Proton, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axMuonRecoStartX, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axMuonRecoStartY, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axMuonRecoStartZ, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axMuonRecoDirX, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axMuonRecoDirY, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axMuonRecoDirZ, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axMuonRecoForceDownDirX, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axMuonRecoForceDownDirY, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axMuonRecoForceDownDirZ, spillCut, cut));
    //====     Proton
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axNProtonCandTrack, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axNProtonCandMatched, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axProtonRecoP, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axProtonRecoCosineTheta, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axProtonRecoNuMICosineTheta, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axProtonLength, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axProtonChi2Proton, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axProtonBestmatchPDG, spillCut, cut));
    //====     Neutrino
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axNeutrinoCombinedEnergy, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axNeutrinoQE, spillCut, cut));
    //====   Reco vs Truth
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axMuonRecoP, axMuonTruthP, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axMuonRecoCosineTheta, axMuonTruthCosineTheta, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axMuonRecoNuMICosineTheta, axMuonTruthNuMICosineTheta, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axProtonRecoP, axProtonTruthP, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axNeutrinoCombinedEnergy, axNeutrinoTruthE, spillCut, cut));
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

void HistoProducer::bookRecoPerformance(SpectrumLoader& loader, SpillCut spillCut, Cut cut){

  //==== Binning
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
  const Binning binsResidualFraction = Binning::Simple(1000, -5., 5.);
  const Binning binsPDG = Binning::Simple(10000, -5000., 5000.);

  //==== Axis
  //====   Truth variables
  //====     Neutrinos
  const HistAxis axNeutrinoTruthE("NeutrinoTruthE", binsEnergy, varNeutrinoTruthE);
  const HistAxis axTruthQ2("TruthQ2", binsEnergy, varTruthQ2);
  const HistAxis axTruthq0_lab("Truthq0_lab", binsEnergy, varTruthq0_lab);
  const HistAxis axTruthmodq_lab("Truthmodq_lab", binsEnergy, varTruthmodq_lab);
  //====     Muon
  const HistAxis axMuonTruthP("MuonTruthP", binsMomentum, varMuonTruthP);
  const HistAxis axMuonTruthT("MuonTruthT", Binning::Simple(200., 0., 2.), varMuonTruthT);
  const HistAxis axMuonTruthCosineTheta("MuonTruthCosineTheta", binsCosineTheta, varMuonTruthCosineTheta);
  const HistAxis axMuonTruthNuMICosineTheta("MuonTruthNuMICosineTheta", binsCosineTheta, varMuonTruthNuMICosineTheta);
  const HistAxis axTruthMuonMatchedTrackChi2Muon("TruthMuonMatchedTrackChi2Muon", binsChi2, varTruthMuonMatchedTrackChi2Muon);
  const HistAxis axTruthMuonMatchedTrackChi2Proton("TruthMuonMatchedTrackChi2Proton", binsChi2, varTruthMuonMatchedTrackChi2Proton);
  const HistAxis axTruthMuonMatchedTrackLength("TruthMuonMatchedTrackLength", binsLength, varTruthMuonMatchedTrackLength);
  //====     Proton
  const HistAxis axProtonTruthP("ProtonTruthP", binsMomentum, varProtonTruthP);
  const HistAxis axProtonTruthT("ProtonTruthT", Binning::Simple(100., 0., 1.), varProtonTruthT);
  const HistAxis axProtonTruthCosineTheta("ProtonTruthCosineTheta", binsCosineTheta, varProtonTruthCosineTheta);
  const HistAxis axProtonTruthNuMICosineTheta("ProtonTruthNuMICosineTheta", binsCosineTheta, varProtonTruthNuMICosineTheta);
  const HistAxis axTruthProtonMatchedTrackChi2Proton("TruthProtonMatchedTrackChi2Proton", binsChi2, varTruthProtonMatchedTrackChi2Proton);
  const HistAxis axTruthProtonMatchedTrackLength("TruthProtonMatchedTrackLength", binsLength, varTruthProtonMatchedTrackLength);
  //const HistAxis axTruthProtonMatchedStubE("TruthProtonMatchedStubE", binsMomentum, varTruthProtonMatchedStubE);
  //====   Reco variables
  //====     Muon
  const HistAxis axMuonRecoP("MuonRecoP", binsMomentum, varMuonRecoP);
  const HistAxis axMuonBestmatchP("MuonBestmatchP", binsMomentum, varMuonBestmatchP);
  const HistAxis axMuonBestmatchPDG("MuonBestmatchPDG", binsPDG, varMuonBestmatchPDG);
  const HistAxis axMuonRecoCosineTheta("MuonRecoCosineTheta", binsCosineTheta, varMuonRecoCosineTheta);
  const HistAxis axMuonRecoNuMICosineTheta("MuonRecoNuMICosineTheta", binsCosineTheta, varMuonRecoNuMICosineTheta);
  const HistAxis axMuonLength("MuonLength", binsLength, varMuonLength);
  const HistAxis axMuonChi2Muon("MuonChi2Muon", binsChi2, varMuonChi2Muon);
  const HistAxis axMuonChi2Proton("MuonChi2Proton", binsChi2, varMuonChi2Proton);
  //====     Proton
  const HistAxis axProtonRecoP("ProtonRecoP", binsMomentum, varProtonRecoP);
  const HistAxis axProtonBestmatchP("ProtonBestmatchP", binsMomentum, varProtonBestmatchP);
  const HistAxis axProtonBestmatchPDG("ProtonBestmatchPDG", binsPDG, varProtonBestmatchPDG);
  const HistAxis axProtonRecoCosineTheta("ProtonRecoCosineTheta", binsCosineTheta, varProtonRecoCosineTheta);
  const HistAxis axProtonRecoNuMICosineTheta("ProtonRecoNuMICosineTheta", binsCosineTheta, varProtonRecoNuMICosineTheta);
  const HistAxis axProtonLength("ProtonLength", Binning::Simple(200, 0., 200.), varProtonLength);
  const HistAxis axProtonChi2Proton("ProtonChi2Proton", binsChi2, varProtonChi2Proton);
  //====     Neutrino
  const HistAxis axNeutrinoCombinedEnergy("NeutrinoCombinedEnergy", binsEnergy, varNeutrinoCombinedEnergy);
  const HistAxis axNeutrinoQE("NeutrinoQE", binsEnergy, varNeutrinoQE);
  //====   Residual
  //====     Muon
  const HistAxis axMuonPResidualFraction("MuonPResidualFraction", binsResidualFraction, varMuonPResidualFraction);
  //====     Proton
  const HistAxis axProtonPResidualFraction("ProtonPResidualFraction", binsResidualFraction, varProtonPResidualFraction);

  //==== Spectrum
  //====   Truth variables
  //====     Neutrinos
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axNeutrinoTruthE, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axTruthQ2, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axTruthq0_lab, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axTruthmodq_lab, spillCut, cut));
  //====     Muon
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axMuonTruthP, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axMuonTruthCosineTheta, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axMuonTruthNuMICosineTheta, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axTruthMuonMatchedTrackChi2Muon, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axTruthMuonMatchedTrackChi2Proton, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axTruthMuonMatchedTrackLength, spillCut, cut));
  //====     Proton
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axProtonTruthP, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axProtonTruthCosineTheta, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axProtonTruthNuMICosineTheta, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axTruthProtonMatchedTrackChi2Proton, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axTruthProtonMatchedTrackLength, spillCut, cut));
  //map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axTruthProtonMatchedStubE, spillCut, cut));
  //====   Reco variables
  //====     Muon
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axMuonRecoP, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axMuonBestmatchP, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axMuonBestmatchPDG, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axMuonRecoCosineTheta, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axMuonRecoNuMICosineTheta, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axMuonLength, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axMuonChi2Muon, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axMuonChi2Proton, spillCut, cut));
  //====     Proton
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axProtonRecoP, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axProtonBestmatchP, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axProtonBestmatchPDG, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axProtonRecoCosineTheta, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axProtonRecoNuMICosineTheta, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axProtonLength, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axProtonChi2Proton, spillCut, cut));
  //====     Neutrino
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axNeutrinoCombinedEnergy, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axNeutrinoQE, spillCut, cut));
  //====   Reco vs Truth
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axMuonRecoP, axMuonTruthP, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axMuonRecoCosineTheta, axMuonTruthCosineTheta, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axMuonRecoNuMICosineTheta, axMuonTruthNuMICosineTheta, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axNeutrinoCombinedEnergy, axNeutrinoTruthE, spillCut, cut));
  //====   Residual
  //====     Muon
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axMuonPResidualFraction, spillCut, cut));
  //====     Proton
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axProtonPResidualFraction, spillCut, cut));
  //====   Extra
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axMuonTruthP, axMuonBestmatchP, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axProtonTruthP, axProtonBestmatchP, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axMuonRecoP, axMuonPResidualFraction, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axMuonBestmatchP, axMuonPResidualFraction, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axMuonTruthT, axTruthMuonMatchedTrackLength, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axProtonTruthT, axTruthProtonMatchedTrackLength, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axMuonTruthT, axMuonLength, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axProtonTruthT, axProtonLength, spillCut, cut));

}

void HistoProducer::bookQuickCutEfficiency(SpectrumLoader& loader, SpillCut spillCut, Cut cut){

  const HistAxis axCountSlice("CountSlice", Binning::Simple( 1, 0., 1.), varCountSlice);
  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum(loader, axCountSlice, spillCut, cut) );

}

void HistoProducer::bookCRTStudy(SpectrumLoader& loader, SpillCut spillCut, Cut cut){

  map_cutName_to_vec_Spectrums[currentCutName].push_back( 
    new Spectrum("CRTHitT0", Binning::Simple(400, -10, 30), loader, spillvarCRTHitT0, spillCut)
  );


  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum("CRTHitXvsY", "",
      loader,
      Binning::Simple(120, -504, 504), spillvarCRTHitPosX,
      Binning::Simple(100, -395, 435), spillvarCRTHitPosY,
      spillCut)
  );

}

void HistoProducer::bookTEMP(SpectrumLoader& loader, SpillCut spillCut, Cut cut){

  const Binning binsChi2 = Binning::Simple(400, 0., 400);

  const HistAxis axMuonChi2Muon("MuonChi2Muon", binsChi2, varMuonChi2Muon);
  const HistAxis axMuonChi2Proton("MuonChi2Proton", binsChi2, varMuonChi2Proton);
  const HistAxis axMuonRecoTrackFromVertex("MuonRecoTrackFromVertex", Binning::Simple(300, 0., 30.), varMuonRecoTrackFromVertex);
  const HistAxis axSliceNuVertexY("SliceNuVertexY", Binning::Simple(400,-200,200), varSliceNuVertexY);

  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axMuonChi2Muon, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axMuonChi2Proton, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axMuonRecoTrackFromVertex, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axSliceNuVertexY, spillCut, cut));

}

void HistoProducer::saveHistograms(){

  outputfile->cd();

  cout << "[HistoProducer::saveHistograms] Number of cuts = " << vec_cutNames.size() << endl;
  const unsigned int nCutName = vec_cutNames.size();
  double inputSamplePOT(-1.);
  double inputSampleLiveTime(-1.);
  if(nCutName>0){

    TH1D *hist_CutNames = new TH1D("hist_CutNames", "", nCutName, 0., 1.*nCutName);

    for(unsigned int ic=0; ic<nCutName; ic++){

      TString cutName = vec_cutNames.at(ic);
      TString dirName = cutName;
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

        if(ic==0 && i==0 && fillBeamInfo){
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
          TH1 *h = sNominal.ToTH1( sNominal.POT() );
          h->SetName(newLabel+"_"+dirName);
          h->Write();
          cout << "[HistoProducer::saveHistograms]     Nominal histogram = " << baseLabel << endl;
          for(unsigned int iu=0; iu<vec_SystEnsembleSpectrumPairs.at(i).second->NUniverses(); ++iu){
            TH1 *hU = vec_SystEnsembleSpectrumPairs.at(i).second->Universe(iu).ToTH1( vec_SystEnsembleSpectrumPairs.at(i).second->POT() );
            //TString systName = IGENIESysts.at(iu)->ShortName();
            hU->SetName(baseLabel+"_"+systematicName+"_Univ"+TString::Itoa(iu,10)+"_"+dirName);
            hU->Write();
          }
          TGraphAsymmErrors *grErrorBand = vec_SystEnsembleSpectrumPairs.at(i).second->ErrorBand(1., vec_SystEnsembleSpectrumPairs.at(i).second->POT()  );
          grErrorBand->SetName(baseLabel+"_"+systematicName+"_ErrorBandFromES_"+dirName);
          grErrorBand->Write();
        }
        else if(sNominal.GetBinnings().size()==2){
          for(unsigned int iu=0; iu<vec_SystEnsembleSpectrumPairs.at(i).second->NUniverses(); ++iu){
            TH2 *hU = vec_SystEnsembleSpectrumPairs.at(i).second->Universe(iu).ToTH2( vec_SystEnsembleSpectrumPairs.at(i).second->POT() );
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

      outputfile->cd();

    } // END loop cutname

    outputfile->cd();
    hist_CutNames->Write();

  } // END if nCutName>0

  //==== This should not be added by hadd
  if(fillNominal){
    outputfile->cd();
    outputfile->mkdir("JobInfo");
    outputfile->cd("JobInfo");
    TH1D *hist_TargetPOT = new TH1D("hist_TargetPOT", "", 1, 0., 1.);
    hist_TargetPOT->SetBinContent(1, TargetPOT);
    hist_TargetPOT->Write();
    outputfile->cd();
  }
  if(fillBeamInfo){
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

  cout << "[HistoProducer::setSystematicWeights] Setting GENIE systematics" << endl;

  const std::vector<std::string> genieKnobNames = GetSBNGenieWeightNames();
  for(const std::string& name: genieKnobNames){
    std::string psetname(name);
    IGENIESysts.push_back( new SBNWeightSyst(name) );
  }
/*
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
*/
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

