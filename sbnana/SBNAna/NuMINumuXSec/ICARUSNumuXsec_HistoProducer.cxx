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
  const Binning binsCosineTheta2 = Binning::Simple( int( (xCosineThetaMax-xCosineThetaMin)/0.01 ), xCosineThetaMin, xCosineThetaMax );
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
  const HistAxis axFMTimeDataTemp("FMTime", binsFMTime, varFMTimeDataTemp);
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
  //====   Vertex
  const HistAxis axVtxResidualX("VtxResidualX", Binning::Simple(100, -50., 50.), varVtxResidualX);
  const HistAxis axVtxResidualY("VtxResidualY", Binning::Simple(100, -50., 50.), varVtxResidualY);
  const HistAxis axVtxResidualZ("VtxResidualZ", Binning::Simple(100, -50., 50.), varVtxResidualZ);
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
  //====     Proton
  const HistAxis axProtonTruthP("ProtonTruthP", binsMomentum, varProtonTruthP);
  const HistAxis axProtonTruthCosineTheta("ProtonTruthCosineTheta", binsCosineTheta, varProtonTruthCosineTheta);
  const HistAxis axProtonTruthNuMICosineTheta("ProtonTruthNuMICosineTheta", binsCosineTheta, varProtonTruthNuMICosineTheta);
  //====     Muon+Proton
  const HistAxis axTruthMuonProtonCosineTheta("TruthMuonProtonCosineTheta", binsCosineTheta2, varTruthMuonProtonCosineTheta);
  //====   Reco variables
  //====     Muon
  const HistAxis axNMuonCandTrack("NMuonCandTrack", Binning::Simple(5, 0., 5.), varNMuonCandTrack);
  const HistAxis axMuonRecoP("MuonRecoP", binsMomentum, varMuonRecoP);
  const HistAxis axMuonRecoCosineTheta("MuonRecoCosineTheta", binsCosineTheta, varMuonRecoCosineTheta);
  const HistAxis axMuonRecoNuMICosineTheta("MuonRecoNuMICosineTheta", binsCosineTheta, varMuonRecoNuMICosineTheta);
  const HistAxis axMuonLength("MuonLength", binsLength, varMuonLength);
  const HistAxis axMuonChi2Muon("MuonChi2Muon", binsChi2, varMuonChi2Muon);
  const HistAxis axMuonChi2Proton("MuonChi2Proton", binsChi2, varMuonChi2Proton);
  const HistAxis axMuonChi2Pion("MuonChi2Pion", binsChi2, varMuonChi2Pion);
  const HistAxis axMuonRecoStartX("MuonRecoStartX", binsXPosition, varMuonRecoStartX);
  const HistAxis axMuonRecoStartY("MuonRecoStartY", binsYPosition, varMuonRecoStartY);
  const HistAxis axMuonRecoStartZ("MuonRecoStartZ", binsZPosition, varMuonRecoStartZ);
  const HistAxis axMuonRecoDirectionX("MuonRecoDirectionX", Binning::Simple(100,-1,1), varMuonRecoDirectionX);
  const HistAxis axMuonRecoDirectionY("MuonRecoDirectionY", Binning::Simple(100,-1,1), varMuonRecoDirectionY);
  const HistAxis axMuonRecoDirectionZ("MuonRecoDirectionZ", Binning::Simple(100,-1,1), varMuonRecoDirectionZ);
  const HistAxis axMuonRecoForceDownDirectionX("MuonRecoForceDownDirectionX", Binning::Simple(100,-1,1), varMuonRecoForceDownDirectionX);
  const HistAxis axMuonRecoForceDownDirectionY("MuonRecoForceDownDirectionY", Binning::Simple(100,-1,1), varMuonRecoForceDownDirectionY);
  const HistAxis axMuonRecoForceDownDirectionZ("MuonRecoForceDownDirectionZ", Binning::Simple(100,-1,1), varMuonRecoForceDownDirectionZ);
  const HistAxis axMuonBestmatchPDG("MuonBestmatchPDG", binsPDG, varMuonBestmatchPDG);
  //====     Proton
  const HistAxis axNProtonCandTrack("NProtonCandTrack", Binning::Simple(10, 0., 10.), varNProtonCandTrack);
  const HistAxis axNProtonCandMatched("NProtonCandMatched", Binning::Simple(10, 0., 10.), varNProtonCandMatched);
  const HistAxis axProtonRecoP("ProtonRecoP", binsMomentum, varProtonRecoP);
  const HistAxis axProtonRecoCosineTheta("ProtonRecoCosineTheta", binsCosineTheta, varProtonRecoCosineTheta);
  const HistAxis axProtonRecoNuMICosineTheta("ProtonRecoNuMICosineTheta", binsCosineTheta, varProtonRecoNuMICosineTheta);
  const HistAxis axProtonLength("ProtonLength", binsLength, varProtonLength);
  const HistAxis axProtonChi2Proton("ProtonChi2Proton", binsChi2, varProtonChi2Proton);
  const HistAxis axProtonBestmatchPDG("ProtonBestmatchPDG", binsPDG, varProtonBestmatchPDG);
  //====     Muon+Proton
  const HistAxis axMuonProtonCosineTheta("MuonProtonCosineTheta", binsCosineTheta2, varMuonProtonCosineTheta);
  //====     Neutrino
  const HistAxis axNeutrinoCombinedEnergy("NeutrinoCombinedEnergy", binsEnergy, varNeutrinoCombinedEnergy);
  //====     Residuals
  const HistAxis axNeutrinoCombinedEnergyResidual("NeutrinoCombinedEnergyResidual", Binning::Simple(100, -5., 5.), varNeutrinoCombinedEnergyResidual);
  const HistAxis axNeutrinoCombinedEnergyResidualFraction("NeutrinoCombinedEnergyResidualFraction", Binning::Simple(1000, -2., 2.), varNeutrinoCombinedEnergyResidualFraction);

  const HistAxis axNeutrinoQEEnergy("NeutrinoQEEnergy", binsEnergy, varNeutrinoQEEnergy);

  if(fillNominal){
    //====   SpillVar
    map_cutName_to_vec_Spectrums[currentCutName].push_back(
      new Spectrum("SpillMCNeutrinoE", binsEnergy, loader, spillvarMCNeutrinoE, spillCut)
    );
    map_cutName_to_vec_Spectrums[currentCutName].push_back(
      new Spectrum("CountSpill", Binning::Simple(1, 0., 1.), loader, spillvarCountSpill, spillCut)
    );

    //====   Slice variables
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axVertexRecoX, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axVertexRecoY, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axVertexRecoZ, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axFMScore, spillCut, cut));
    if(IsData && sampleName!="NUMI:1" && sampleName!="NUMI:8515"){
      map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axFMTimeDataTemp, spillCut, cut));
    }
    else{
      map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axFMTime, spillCut, cut));
    }
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
/*
    //====   Vertex
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axVtxResidualX, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axVtxResidualY, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axVtxResidualZ, spillCut, cut));
*/
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
/*
    //====   All track
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum("AllTrackStartPositionX", binsXPosition, loader, varAllTrackStartPositionX, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum("AllTrackStartPositionY", binsYPosition, loader, varAllTrackStartPositionY, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum("AllTrackStartPositionZ", binsZPosition, loader, varAllTrackStartPositionZ, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum("AllTrackDirectionX", Binning::Simple(40, -1., 1.), loader, varAllTrackDirectionX, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum("AllTrackDirectionY", Binning::Simple(40, -1., 1.), loader, varAllTrackDirectionY, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum("AllTrackDirectionZ", Binning::Simple(40, -1., 1.), loader, varAllTrackDirectionZ, spillCut, cut));
*/
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum("AllTrackLength", Binning::Simple(90, 0., 900.), loader, varAllTrackLength, spillCut, cut));
/*
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum("AllTrackMatchedTruthDirectionX", Binning::Simple(40, -1., 1.), loader, varAllTrackMatchedTruthDirectionX, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum("AllTrackMatchedTruthDirectionY", Binning::Simple(40, -1., 1.), loader, varAllTrackMatchedTruthDirectionY, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum("AllTrackMatchedTruthDirectionZ", Binning::Simple(40, -1., 1.), loader, varAllTrackMatchedTruthDirectionZ, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum("AllTrackMatchedTruthLength", Binning::Simple(90, 0., 900.), loader, varAllTrackMatchedTruthLength, spillCut, cut));
    //==== TODO test
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum("TestSelectedTrackStartPositionX", binsXPosition, loader, varTestSelectedTrackStartPositionX, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum("TestSelectedTrackStartPositionY", binsYPosition, loader, varTestSelectedTrackStartPositionY, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum("TestSelectedTrackStartPositionZ", binsZPosition, loader, varTestSelectedTrackStartPositionZ, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(
      new Spectrum(
        "TestSelectedTrackStartPositionZvsY", loader,
        Binning::Simple(180, -900., 900.), varTestSelectedTrackStartPositionZ,
        Binning::Simple(36, -180., 180.), varTestSelectedTrackStartPositionY,
        spillCut, cut
      )
    );
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum("TestSelectedTrackEndPositionX", binsXPosition, loader, varTestSelectedTrackEndPositionX, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum("TestSelectedTrackEndPositionY", binsYPosition, loader, varTestSelectedTrackEndPositionY, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum("TestSelectedTrackEndPositionZ", binsZPosition, loader, varTestSelectedTrackEndPositionZ, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum("TestSelectedTrackDirectionX", Binning::Simple(40, -1., 1.), loader, varTestSelectedTrackDirectionX, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum("TestSelectedTrackDirectionY", Binning::Simple(40, -1., 1.), loader, varTestSelectedTrackDirectionY, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum("TestSelectedTrackDirectionZ", Binning::Simple(40, -1., 1.), loader, varTestSelectedTrackDirectionZ, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum("TestSelectedTrackLength", Binning::Simple(90, 0., 900.), loader, varTestSelectedTrackLength, spillCut, cut));
    //==== its matched version
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum("TestSelectedTrackMatchedTruthStartPositionX", binsXPosition, loader, varTestSelectedTrackMatchedTruthStartPositionX, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum("TestSelectedTrackMatchedTruthStartPositionY", binsYPosition, loader, varTestSelectedTrackMatchedTruthStartPositionY, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum("TestSelectedTrackMatchedTruthStartPositionZ", binsZPosition, loader, varTestSelectedTrackMatchedTruthStartPositionZ, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum("TestSelectedTrackMatchedTruthEndPositionX", binsXPosition, loader, varTestSelectedTrackMatchedTruthEndPositionX, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum("TestSelectedTrackMatchedTruthEndPositionY", binsYPosition, loader, varTestSelectedTrackMatchedTruthEndPositionY, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum("TestSelectedTrackMatchedTruthEndPositionZ", binsZPosition, loader, varTestSelectedTrackMatchedTruthEndPositionZ, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum("TestSelectedTrackMatchedTruthDirectionX", Binning::Simple(40, -1., 1.), loader, varTestSelectedTrackMatchedTruthDirectionX, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum("TestSelectedTrackMatchedTruthDirectionY", Binning::Simple(40, -1., 1.), loader, varTestSelectedTrackMatchedTruthDirectionY, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum("TestSelectedTrackMatchedTruthDirectionZ", Binning::Simple(40, -1., 1.), loader, varTestSelectedTrackMatchedTruthDirectionZ, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum("TestSelectedTrackMatchedTruthLength", Binning::Simple(90, 0., 900.), loader, varTestSelectedTrackMatchedTruthLength, spillCut, cut));
*/
    //==== Longest track
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum("LongestTrackDirectionX", Binning::Simple(40, -1., 1.), loader, varLongestTrackDirectionX, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum("LongestTrackDirectionY", Binning::Simple(40, -1., 1.), loader, varLongestTrackDirectionY, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum("LongestTrackDirectionZ", Binning::Simple(40, -1., 1.), loader, varLongestTrackDirectionZ, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum("LongestTrackDirectionXZ", Binning::Simple(40, 0., 1.), loader, varLongestTrackDirectionXZ, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum("LongestTrackLength", Binning::Simple(90, 0., 900.), loader, varLongestTrackLength, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum("LongestTrackChi2Muon", Binning::Simple(50, 0., 100.), loader, varLongestTrackChi2Muon, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum("LongestTrackChi2Proton", Binning::Simple(40, 0., 400.), loader, varLongestTrackChi2Proton, spillCut, cut));

    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum("LongestTrackForceDownDirectionX", Binning::Simple(40, -1., 1.), loader, varLongestTrackForceDownDirectionX, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum("LongestTrackForceDownDirectionY", Binning::Simple(40, -1., 1.), loader, varLongestTrackForceDownDirectionY, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum("LongestTrackForceDownDirectionZ", Binning::Simple(40, -1., 1.), loader, varLongestTrackForceDownDirectionZ, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum("LongestTrackForceDownStartPositionX", Binning::Simple(40, -1., 1.), loader, varLongestTrackForceDownStartPositionX, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum("LongestTrackForceDownStartPositionY", Binning::Simple(40, -1., 1.), loader, varLongestTrackForceDownStartPositionY, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum("LongestTrackForceDownStartPositionZ", Binning::Simple(40, -1., 1.), loader, varLongestTrackForceDownStartPositionZ, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum("LongestTrackForceDownEndPositionX", Binning::Simple(40, -1., 1.), loader, varLongestTrackForceDownEndPositionX, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum("LongestTrackForceDownEndPositionY", Binning::Simple(40, -1., 1.), loader, varLongestTrackForceDownEndPositionY, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum("LongestTrackForceDownEndPositionZ", Binning::Simple(40, -1., 1.), loader, varLongestTrackForceDownEndPositionZ, spillCut, cut));

  if(!IsData){
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum("LongestTrackMatchedTruthDirectionX", Binning::Simple(40, -1., 1.), loader, varLongestTrackMatchedTruthDirectionX, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum("LongestTrackMatchedTruthDirectionY", Binning::Simple(40, -1., 1.), loader, varLongestTrackMatchedTruthDirectionY, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum("LongestTrackMatchedTruthDirectionZ", Binning::Simple(40, -1., 1.), loader, varLongestTrackMatchedTruthDirectionZ, spillCut, cut));
  }
/*
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum("FMChargeQ", binsFMChargeQ, loader, varFMChargeQ, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum("FMChargeCenterX", binsXPosition, loader, varFMChargeCenterX, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum("FMChargeCenterY", binsYPosition, loader, varFMChargeCenterY, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum("FMChargeCenterZ", binsZPosition, loader, varFMChargeCenterZ, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum("FMLightCenterX", binsXPosition, loader, varFMLightCenterX, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum("FMLightCenterY", binsYPosition, loader, varFMLightCenterY, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum("FMLightCenterZ", binsZPosition, loader, varFMLightCenterZ, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum("FMLightPE", binsFMLightPE, loader, varFMLightPE, spillCut, cut));
*/

    //====   Truth variables
    if(!IsData){
      //====     Neutrinos
      map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axNeutrinoTruthE, spillCut, cut));
      map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axTruthQ2, spillCut, cut));
      map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axTruthq0_lab, spillCut, cut));
      map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axTruthmodq_lab, spillCut, cut));
      //====     Muon
      map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axMuonTruthP, spillCut, cut));
      map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axMuonTruthCosineTheta, spillCut, cut));
      map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axMuonTruthNuMICosineTheta, spillCut, cut));
      map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum("MuonTruthDirectionX", Binning::Simple(100, -1., 1.), loader, varMuonTruthDirectionX, spillCut, cut));
      map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum("MuonTruthDirectionY", Binning::Simple(100, -1., 1.), loader, varMuonTruthDirectionY, spillCut, cut));
      map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum("MuonTruthDirectionZ", Binning::Simple(100, -1., 1.), loader, varMuonTruthDirectionZ, spillCut, cut));
      //====     Proton
      map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axProtonTruthP, spillCut, cut));
      map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axProtonTruthCosineTheta, spillCut, cut));
      map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axProtonTruthNuMICosineTheta, spillCut, cut));
      //====     Muon+Proton
      map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axTruthMuonProtonCosineTheta, spillCut, cut));
    }
    //====   Reco variables
    //====     Muon
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axNMuonCandTrack, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axMuonRecoP, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axMuonRecoCosineTheta, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axMuonRecoNuMICosineTheta, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axMuonLength, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axMuonChi2Muon, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axMuonChi2Proton, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axMuonChi2Pion, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axMuonRecoStartX, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axMuonRecoStartY, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axMuonRecoStartZ, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axMuonRecoDirectionX, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axMuonRecoDirectionY, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axMuonRecoDirectionZ, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axMuonRecoForceDownDirectionX, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axMuonRecoForceDownDirectionY, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axMuonRecoForceDownDirectionZ, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axMuonBestmatchPDG, spillCut, cut));
    //====     Proton
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axNProtonCandTrack, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axNProtonCandMatched, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axProtonRecoP, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axProtonRecoCosineTheta, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axProtonRecoNuMICosineTheta, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axProtonLength, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axProtonChi2Proton, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axProtonBestmatchPDG, spillCut, cut));
    //====     Muon+Proton
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axMuonProtonCosineTheta, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axMuonLength, axMuonProtonCosineTheta, spillCut, cut));
    //====     Neutrino
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axNeutrinoCombinedEnergy, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axNeutrinoCombinedEnergyResidual, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axNeutrinoCombinedEnergyResidualFraction, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axNeutrinoQEEnergy, spillCut, cut));
/*
    //====   Reco vs Truth
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axMuonRecoP, axMuonTruthP, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axMuonRecoCosineTheta, axMuonTruthCosineTheta, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axMuonRecoNuMICosineTheta, axMuonTruthNuMICosineTheta, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axProtonRecoP, axProtonTruthP, spillCut, cut));
    map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axNeutrinoCombinedEnergy, axNeutrinoTruthE, spillCut, cut));
*/
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
  const HistAxis axNeutrinoQEEnergy("NeutrinoQEEnergy", binsEnergy, varNeutrinoQEEnergy);

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
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axNeutrinoQEEnergy, spillCut, cut));
  //====   Reco vs Truth
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axMuonRecoP, axMuonTruthP, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axMuonRecoCosineTheta, axMuonTruthCosineTheta, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axMuonRecoNuMICosineTheta, axMuonTruthNuMICosineTheta, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axNeutrinoCombinedEnergy, axNeutrinoTruthE, spillCut, cut));
  //====   Extra
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axMuonTruthP, axMuonBestmatchP, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axProtonTruthP, axProtonBestmatchP, spillCut, cut));
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

void HistoProducer::GetSmearMat(SpectrumLoader& loader, SpillCut spillCut, Cut cut){

  //==== Muon

  static vector<double> bins_ContainedMuonP = {0.0, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.70, 0.80, 0.90, 1.0, 2.0, 5.0};
  static vector<double> bins_ExitingMuonP = {0.0, 0.40, 0.60, 0.80, 0.90, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 5.0};
  static vector<double> bins_MuonLength = {0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 150.0, 160.0, 170.0, 180.0, 190.0, 200.0, 300.0, 500.0, 1000.0};

  const Binning MuonLBinning = Binning::Custom(bins_MuonLength);
  const Binning MuonPBinning = Binning::Simple(50., 0., 5.);
  const Binning ContainedMuonPBinning = Binning::Custom(bins_ContainedMuonP);
  const Binning ExitingMuonPBinning = Binning::Custom(bins_ExitingMuonP);
  const Binning MuonPResidualBinning = Binning::Simple(100, -5., 5.);
  const Binning MuonPResidualFractionBinning = Binning::Simple(1000, -2., 2.);
  const Binning binsXPosition = Binning::Simple(1000, -500., 500.);
  const Binning binsYPosition = Binning::Simple(400, -200., 200.);
  const Binning binsZPosition = Binning::Simple(2000, -1000., 1000.);
  const Binning binsCosineTheta = Binning::Simple(20, -1., 1.);
  const Binning binsChi2 = Binning::Simple(400, 0., 400);

  //==== Truth
  const HistAxis axMuonTruthP("MuonTruthP", MuonPBinning, varMuonTruthP);
  const HistAxis axContainedMuonTruthP("ContainedMuonTruthP", ContainedMuonPBinning, varMuonTruthP);
  const HistAxis axExitingMuonTruthP("ExitingMuonTruthP", ExitingMuonPBinning, varMuonTruthP);
  const HistAxis axMuonTruthNuMICosineTheta("MuonTruthNuMICosineTheta", binsCosineTheta, varMuonTruthNuMICosineTheta);
  //==== Matched Reco
  const HistAxis axTruthMuonMatchedTrackEndPositionX("TruthMuonMatchedTrackEndPositionX", binsXPosition, varTruthMuonMatchedTrackEndPositionX);
  const HistAxis axTruthMuonMatchedTrackEndPositionY("TruthMuonMatchedTrackEndPositionY", binsYPosition, varTruthMuonMatchedTrackEndPositionY);
  const HistAxis axTruthMuonMatchedTrackEndPositionZ("TruthMuonMatchedTrackEndPositionZ", binsZPosition, varTruthMuonMatchedTrackEndPositionZ);
  const HistAxis axTruthMuonMatchedTrackNuMICosineTheta("TruthMuonMatchedTrackNuMICosineTheta", binsCosineTheta, varTruthMuonMatchedTrackNuMICosineTheta);
  const HistAxis axTruthMuonMatchedTrackLength("TruthMuonMatchedTrackLength", MuonLBinning, varTruthMuonMatchedTrackLength);
  const HistAxis axTruthMuonMatchedTrackRangeP("TruthMuonMatchedTrackRangeP", ContainedMuonPBinning, varTruthMuonMatchedTrackRangeP);
  const HistAxis axTruthMuonMatchedTrackMCSP("TruthMuonMatchedTrackMCSP", ExitingMuonPBinning, varTruthMuonMatchedTrackMCSP);
  const HistAxis axTruthMuonMatchedTrackCombinedP("TruthMuonMatchedTrackCombinedP", ContainedMuonPBinning, varTruthMuonMatchedTrackCombinedP);
  const HistAxis axTruthMuonMatchedTrackChi2Muon("TruthMuonMatchedTrackChi2Muon", binsChi2, varTruthMuonMatchedTrackChi2Muon);
  const HistAxis axTruthMuonMatchedTrackChi2Proton("TruthMuonMatchedTrackChi2Proton", binsChi2, varTruthMuonMatchedTrackChi2Proton);
  const HistAxis axTruthMuonMatchedTrackChi2Pion("TruthMuonMatchedTrackChi2Pion", binsChi2, varTruthMuonMatchedTrackChi2Pion);
  //==== Residual
  const HistAxis axTruthMuonMatchedTrackRangePResidual("TruthMuonMatchedTrackRangePResidual", MuonPResidualBinning, varTruthMuonMatchedTrackRangePResidual);
  const HistAxis axTruthMuonMatchedTrackMCSPResidual("TruthMuonMatchedTrackMCSPResidual", MuonPResidualBinning, varTruthMuonMatchedTrackMCSPResidual);
  const HistAxis axTruthMuonMatchedTrackCombinedPResidual("TruthMuonMatchedTrackCombinedPResidual", MuonPResidualBinning, varTruthMuonMatchedTrackCombinedPResidual);
  //==== Residual frction
  const HistAxis axTruthMuonMatchedTrackRangePResidualFraction("TruthMuonMatchedTrackRangePResidualFraction", MuonPResidualFractionBinning, varTruthMuonMatchedTrackRangePResidualFraction);
  const HistAxis axTruthMuonMatchedTrackMCSPResidualFraction("TruthMuonMatchedTrackMCSPResidualFraction", MuonPResidualFractionBinning, varTruthMuonMatchedTrackMCSPResidualFraction);
  const HistAxis axTruthMuonMatchedTrackCombinedPResidualFraction("TruthMuonMatchedTrackCombinedPResidualFraction", MuonPResidualFractionBinning, varTruthMuonMatchedTrackCombinedPResidualFraction);

  //==== Truth
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axMuonTruthP, spillCut, cut));
  //==== Matched Reco
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axTruthMuonMatchedTrackEndPositionX, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axTruthMuonMatchedTrackEndPositionY, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axTruthMuonMatchedTrackEndPositionZ, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axTruthMuonMatchedTrackLength, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axTruthMuonMatchedTrackRangeP, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axTruthMuonMatchedTrackMCSP, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axTruthMuonMatchedTrackCombinedP, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axTruthMuonMatchedTrackChi2Muon, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axTruthMuonMatchedTrackChi2Proton, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axTruthMuonMatchedTrackChi2Pion, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum("TruthMuonMatchedTrackChi2Muon_Proton_Pion", loader, 
      Binning::Simple(20, 0., 100), varTruthMuonMatchedTrackChi2Muon,
      Binning::Simple(40, 0., 400), varTruthMuonMatchedTrackChi2Proton,
      Binning::Simple(20, 0., 100), varTruthMuonMatchedTrackChi2Pion,
    spillCut, cut)
  );
  //==== Residual
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axTruthMuonMatchedTrackRangePResidual, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axTruthMuonMatchedTrackMCSPResidual, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axTruthMuonMatchedTrackCombinedPResidual, spillCut, cut));
  //==== Residual frcation
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axTruthMuonMatchedTrackRangePResidualFraction, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axTruthMuonMatchedTrackMCSPResidualFraction, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axTruthMuonMatchedTrackCombinedPResidualFraction, spillCut, cut));
  //==== Truth vs Reco
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axMuonTruthNuMICosineTheta, axTruthMuonMatchedTrackNuMICosineTheta, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axMuonTruthP, axTruthMuonMatchedTrackLength, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axContainedMuonTruthP, axTruthMuonMatchedTrackRangeP, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axExitingMuonTruthP, axTruthMuonMatchedTrackMCSP, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axMuonTruthP, axTruthMuonMatchedTrackCombinedP, spillCut, cut));
  //==== 2D Truth vs ResidualFraction
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axContainedMuonTruthP, axTruthMuonMatchedTrackRangePResidualFraction, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axExitingMuonTruthP, axTruthMuonMatchedTrackMCSPResidualFraction, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axMuonTruthP, axTruthMuonMatchedTrackCombinedPResidualFraction, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axTruthMuonMatchedTrackLength, axTruthMuonMatchedTrackRangePResidualFraction, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axTruthMuonMatchedTrackLength, axTruthMuonMatchedTrackMCSPResidualFraction, spillCut, cut));

  //==== Proton

  const Binning ProtonLBinning = Binning::Simple(40., 0., 200.);
  const Binning ProtonPBinning = Binning::Simple(400, 0., 2.);
  const Binning ProtonPResidualBinning = Binning::Simple(100, -5., 5.);
  const Binning ProtonPResidualFractionBinning = Binning::Simple(1000, -5., 5.);

  //==== Truth
  const HistAxis axProtonTruthP("ProtonTruthP", ProtonPBinning, varProtonTruthP);
  //==== Matched Reco
  const HistAxis axTruthProtonMatchedTrackLength("TruthProtonMatchedTrackLength", ProtonLBinning, varTruthProtonMatchedTrackLength);
  const HistAxis axTruthProtonMatchedTrackRangeP("TruthProtonMatchedTrackRangeP", ProtonPBinning, varTruthProtonMatchedTrackRangeP);
  const HistAxis axTruthProtonMatchedTrackChi2Muon("TruthProtonMatchedTrackChi2Muon", binsChi2, varTruthProtonMatchedTrackChi2Muon);
  const HistAxis axTruthProtonMatchedTrackChi2Proton("TruthProtonMatchedTrackChi2Proton", binsChi2, varTruthProtonMatchedTrackChi2Proton);
  const HistAxis axTruthProtonMatchedTrackChi2Pion("TruthProtonMatchedTrackChi2Pion", binsChi2, varTruthProtonMatchedTrackChi2Pion);
  //==== Residual
  const HistAxis axTruthProtonMatchedTrackRangePResidual("TruthProtonMatchedTrackRangePResidual", ProtonPResidualBinning, varTruthProtonMatchedTrackRangePResidual);
  //==== Residual frction
  const HistAxis axTruthProtonMatchedTrackRangePResidualFraction("TruthProtonMatchedTrackRangePResidualFraction", ProtonPResidualFractionBinning, varTruthProtonMatchedTrackRangePResidualFraction);

  //==== Truth
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axProtonTruthP, spillCut, cut));
  //==== Matched Reco
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axTruthProtonMatchedTrackLength, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axTruthProtonMatchedTrackRangeP, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axTruthProtonMatchedTrackChi2Muon, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axTruthProtonMatchedTrackChi2Proton, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axTruthProtonMatchedTrackChi2Pion, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum("TruthProtonMatchedTrackChi2Muon_Proton_Pion", loader,
      Binning::Simple(20, 0., 100), varTruthProtonMatchedTrackChi2Muon,
      Binning::Simple(40, 0., 400), varTruthProtonMatchedTrackChi2Proton,
      Binning::Simple(20, 0., 100), varTruthProtonMatchedTrackChi2Pion,
    spillCut, cut)
  );

  //==== Residual
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axTruthProtonMatchedTrackRangePResidual, spillCut, cut));
  //==== Residual frcation
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axTruthProtonMatchedTrackRangePResidualFraction, spillCut, cut));
  //==== Truth vs Reco
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axProtonTruthP, axTruthProtonMatchedTrackLength, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axProtonTruthP, axTruthProtonMatchedTrackRangeP, spillCut, cut));
  //==== 2D Truth vs ResidualFraction
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axProtonTruthP, axTruthProtonMatchedTrackRangePResidualFraction, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axTruthProtonMatchedTrackLength, axTruthProtonMatchedTrackRangePResidualFraction, spillCut, cut));

  //==== Pion

  const Binning ChargedPionPBinning = Binning::Simple(400, 0., 2.);

  //==== Truth
  const HistAxis axChargedPionTruthP("ChargedPionTruthP", ChargedPionPBinning, varChargedPionTruthP);
  //==== Matched Reco
  const HistAxis axTruthChargedPionMatchedTrackChi2Muon("TruthChargedPionMatchedTrackChi2Muon", binsChi2, varTruthChargedPionMatchedTrackChi2Muon);
  const HistAxis axTruthChargedPionMatchedTrackChi2Proton("TruthChargedPionMatchedTrackChi2Proton", binsChi2, varTruthChargedPionMatchedTrackChi2Proton);
  const HistAxis axTruthChargedPionMatchedTrackChi2Pion("TruthChargedPionMatchedTrackChi2Pion", binsChi2, varTruthChargedPionMatchedTrackChi2Pion);
  //==== Truth
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axChargedPionTruthP, spillCut, cut));
  //==== Matched Reco
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axTruthChargedPionMatchedTrackChi2Muon, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axTruthChargedPionMatchedTrackChi2Proton, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axTruthChargedPionMatchedTrackChi2Pion, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum("TruthChargedPionMatchedTrackChi2Muon_Proton_Pion", loader,
      Binning::Simple(20, 0., 100), varTruthChargedPionMatchedTrackChi2Muon,
      Binning::Simple(40, 0., 400), varTruthChargedPionMatchedTrackChi2Proton,
      Binning::Simple(20, 0., 100), varTruthChargedPionMatchedTrackChi2Pion,
    spillCut, cut)
  );

}

void HistoProducer::GetChi2Eff(SpectrumLoader& loader, SpillCut spillCut, Cut cut){

  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum("TruthMuonMatchedTrackChi2Muon_Proton_Pion", loader, 
      Binning::Simple(100, 0., 100), varTruthMuonMatchedTrackChi2Muon,
      Binning::Simple(400, 0., 400), varTruthMuonMatchedTrackChi2Proton,
      Binning::Simple(100, 0., 100), varTruthMuonMatchedTrackChi2Pion,
    spillCut, cut)
  );

  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum("TruthProtonMatchedTrackChi2Muon_Proton_Pion", loader,
      Binning::Simple(100, 0., 100), varTruthProtonMatchedTrackChi2Muon,
      Binning::Simple(400, 0., 400), varTruthProtonMatchedTrackChi2Proton,
      Binning::Simple(100, 0., 100), varTruthProtonMatchedTrackChi2Pion,
    spillCut, cut)
  );

  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum("TruthChargedPionMatchedTrackChi2Muon_Proton_Pion", loader,
      Binning::Simple(100, 0., 100), varTruthChargedPionMatchedTrackChi2Muon,
      Binning::Simple(400, 0., 400), varTruthChargedPionMatchedTrackChi2Proton,
      Binning::Simple(100, 0., 100), varTruthChargedPionMatchedTrackChi2Pion,
    spillCut, cut)
  );

}

void HistoProducer::bookTEST(SpectrumLoader& loader, SpillCut spillCut, Cut cut){

  //====  Proton
  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("ProtonTruthT", Binning::Simple(500, 0., 1.0), loader, varProtonTruthT, spillCut, cut) );
  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("TruthProtonMatchedStubE", Binning::Simple(100, 0., 0.05), loader, varTruthProtonMatchedStubE, spillCut, cut) );
  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("TruthProtonMatchedStubLength", Binning::Simple(100, 0., 10.), loader, varTruthProtonMatchedStubLength, spillCut, cut) );
  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("ProtonTruthT_vs_TruthProtonMatchedStubLength", loader, Binning::Simple(50, 0., 0.2), varProtonTruthT, Binning::Simple(100, 0., 10.), varTruthProtonMatchedStubLength, spillCut, cut) );
  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("ProtonTruthT_vs_TruthProtonMatchedStubE", loader, Binning::Simple(50, 0., 0.2), varProtonTruthT, Binning::Simple(100, 0., 0.04), varTruthProtonMatchedStubE, spillCut, cut) );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum("TruthProtonMatchedObjectType", Binning::Simple(8, 0., 8.), loader, varTruthProtonMatchedObjectType, spillCut, cut));


/*
  //==== Test GENIE syst
  double xEnergyMin = 0.; // GeV
  double xEnergyMax = 10.;
  double dxEnergy = 0.1;
  const Binning binsEnergy = Binning::Simple( int( (xEnergyMax-xEnergyMin)/dxEnergy ), xEnergyMin, xEnergyMax );

  const HistAxis axNeutrinoTruthE("NeutrinoTruthE", binsEnergy, varNeutrinoTruthE);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axNeutrinoTruthE, spillCut, cut));
  //==== Systematic
  for(unsigned int i_syst=0; i_syst<IGENIESysts.size(); i_syst++){
    const ISyst* s = IGENIESysts.at(i_syst);
    addUpDownSystematic(loader, axNeutrinoTruthE, spillCut, cut, currentCutName, s);
  }

  if(vec_UniverseWeightsForEachGENIESource.size()==1){
    addEnsembleSpectrum(loader, axNeutrinoTruthE, spillCut, cut, currentCutName, vec_UniverseWeightsForEachGENIESource.at(0), "multisim");
  }
*/

/*
  Var thisWeight = wSterileOscProb;
  //thisWeight = kUnweighted;

  double xEnergyMin = 0.; // GeV
  double xEnergyMax = 10.;
  double dxEnergy = 0.1;
  const Binning binsEnergy = Binning::Simple( int( (xEnergyMax-xEnergyMin)/dxEnergy ), xEnergyMin, xEnergyMax );
  const HistAxis axNeutrinoTruthE("NeutrinoTruthE", binsEnergy, varNeutrinoTruthE);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axNeutrinoTruthE, spillCut, cut, kNoShift, thisWeight));

  double xMomentumMin = 0.; // GeV
  double xMomentumMax = 5.;
  double dxMomentum = 0.1;
  const Binning binsMomentum = Binning::Simple( int( (xMomentumMax-xMomentumMin)/dxMomentum ), xMomentumMin, xMomentumMax );
  const HistAxis axMuonRecoP("MuonRecoP", binsMomentum, varMuonRecoP);
  const HistAxis axMuonTruthP("MuonTruthP", binsMomentum, varMuonTruthP);

  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axMuonRecoP, spillCut, cut, kNoShift, thisWeight));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axMuonTruthP, spillCut, cut, kNoShift, thisWeight));

  double xCosineThetaMin = -1;
  double xCosineThetaMax = 1.;
  double dxCosineTheta = 0.1;
  const Binning binsCosineTheta = Binning::Simple( int( (xCosineThetaMax-xCosineThetaMin)/dxCosineTheta ), xCosineThetaMin, xCosineThetaMax );
  const HistAxis axMuonRecoNuMICosineTheta("MuonRecoNuMICosineTheta", binsCosineTheta, varMuonRecoNuMICosineTheta);
  const HistAxis axMuonTruthNuMICosineTheta("MuonTruthNuMICosineTheta", binsCosineTheta, varMuonTruthNuMICosineTheta);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axMuonRecoNuMICosineTheta, spillCut, cut, kNoShift, thisWeight));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axMuonTruthNuMICosineTheta, spillCut, cut, kNoShift, thisWeight));
*/

}

void HistoProducer::MichelStudy(SpectrumLoader& loader, SpillCut spillCut, Cut cut){

  const Binning binsChi2 = Binning::Simple(400, 0., 400);

  static vector<double> bins_LargedEdX = {0.01, 0.1, 0.5, 1., 10., 50., 100., 200., 300., 400, 500., 600, 700, 800, 900, 1000.};
  //const Binning LargedEdXBinning = Binning::Custom(bins_LargedEdX);
  const Binning LargedEdXBinning = Binning::Simple(300, 0., 300);

  //==== Muon

  const HistAxis axMuonTruthT(
    "MuonTruthT",
    Binning::Simple(500, 0., 5.),
    varMuonTruthT
  );
  const HistAxis axMuonTruthLength(
    "MuonTruthLength",
    Binning::Simple(1000, 0., 1000.),
    varMuonTruthLength
  );
  const HistAxis axTruthMuonMatchedTrackLength(
    "TruthMuonMatchedTrackLength",
    Binning::Simple(1000, 0., 1000.),
    varTruthMuonMatchedTrackLength
  );
  const HistAxis axTruthMuonMatchedTrackNDaughterTracks(
    "TruthMuonMatchedTrackNDaughterTracks",
    Binning::Simple(10, 0., 10.),
    varTruthMuonMatchedTrackNDaughterTracks
  );
  const HistAxis axTruthMuonMatchedTrackNDaughterShowers(
    "TruthMuonMatchedTrackNDaughterShowers",
    Binning::Simple(10, 0., 10.),
    varTruthMuonMatchedTrackNDaughterShowers
  );
  const HistAxis axTruthMuonMatchedTrackFrontLargedEdX(
    "TruthMuonMatchedTrackFrontLargedEdX",
    Binning::Simple(100, 0., 100.),
    varTruthMuonMatchedTrackFrontLargedEdX
  );

  const HistAxis axTruthMuonMatchedTrackStitchedTrackDistance(
    "TruthMuonMatchedTrackStitchedTrackDistance",
    Binning::Simple(300,0.,30.),
    varTruthMuonMatchedTrackStitchedTrackDistance
  );
  const HistAxis axTruthMuonMatchedTrackStitchedTrackPDG(
    "TruthMuonMatchedTrackStitchedTrackPDG",
    Binning::Simple(10000, -5000., 5000.),
    varTruthMuonMatchedTrackStitchedTrackPDG
  );
  const HistAxis axTruthMuonMatchedTrackStitchedTrackLength(
    "TruthMuonMatchedTrackStitchedTrackLength",
    Binning::Simple(100, 0., 100.),
    varTruthMuonMatchedTrackStitchedTrackLength
  );
  const HistAxis axTruthMuonMatchedTrackStitchedTrackCosine(
    "TruthMuonMatchedTrackStitchedTrackCosine",
    Binning::Simple(20, -1., 1.),
    varTruthMuonMatchedTrackStitchedTrackCosine
  );
  const HistAxis axTruthMuonMatchedTrackStitchedTrackChi2Muon(
    "TruthMuonMatchedTrackStitchedTrackChi2Muon",
    Binning::Simple(200, 0., 200.),
    varTruthMuonMatchedTrackStitchedTrackChi2Muon
  );
  const HistAxis axTruthMuonMatchedTrackStitchedTrackChi2Proton(
    "TruthMuonMatchedTrackStitchedTrackChi2Proton",
    Binning::Simple(200, 0., 200.),
    varTruthMuonMatchedTrackStitchedTrackChi2Proton
  );
  const HistAxis axTruthMuonMatchedTrackStitchedTrackChi2Pion(
    "TruthMuonMatchedTrackStitchedTrackChi2Pion",
    Binning::Simple(200, 0., 200.),
    varTruthMuonMatchedTrackStitchedTrackChi2Pion
  );

  const HistAxis axTruthMuonMatchedTrackStitchedShowerDistance(
    "TruthMuonMatchedTrackStitchedShowerDistance",
    Binning::Simple(300,0.,30.),
    varTruthMuonMatchedTrackStitchedShowerDistance
  );
  const HistAxis axTruthMuonMatchedTrackStitchedShowerPDG(
    "TruthMuonMatchedTrackStitchedShowerPDG",
    Binning::Simple(10000, -5000., 5000.),
    varTruthMuonMatchedTrackStitchedShowerPDG
  );
  const HistAxis axTruthMuonMatchedTrackStitchedShowerEnergy(
    "TruthMuonMatchedTrackStitchedShowerEnergy",
    Binning::Simple(200, 0., 2.),
    varTruthMuonMatchedTrackStitchedShowerEnergy
  );
  const HistAxis axTruthMuonMatchedTrackStitchedShowerCosine(
    "TruthMuonMatchedTrackStitchedShowerCosine",
    Binning::Simple(20, -1., 1.),
    varTruthMuonMatchedTrackStitchedShowerCosine
  );

  const HistAxis axTruthMuonMatchedTrackChi2Muon("TruthMuonMatchedTrackChi2Muon", binsChi2, varTruthMuonMatchedTrackChi2Muon);
  const HistAxis axTruthMuonMatchedTrackChi2Proton("TruthMuonMatchedTrackChi2Proton", binsChi2, varTruthMuonMatchedTrackChi2Proton);
  const HistAxis axTruthMuonMatchedTrackChi2Pion("TruthMuonMatchedTrackChi2Pion", binsChi2, varTruthMuonMatchedTrackChi2Pion);

  const HistAxis axTruthMuonMatchedTrackRangePResidualFraction(
    "TruthMuonMatchedTrackRangePResidualFraction", 
    Binning::Simple(1000, -2., 2.),
    varTruthMuonMatchedTrackRangePResidualFraction
  );

  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axMuonTruthT, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axMuonTruthLength, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axTruthMuonMatchedTrackLength, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axTruthMuonMatchedTrackNDaughterTracks, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axTruthMuonMatchedTrackNDaughterShowers, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axTruthMuonMatchedTrackFrontLargedEdX, spillCut, cut));

  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum("TruthMuonMatchedTrackFront_rr_vs_dedx", loader,
      Binning::Simple(30, 0., 30.), varTruthMuonMatchedTrackFrontrr,
      LargedEdXBinning, varTruthMuonMatchedTrackFrontdedx,
      spillCut, cut
    )
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum("TruthMuonMatchedTrackFront_rr_vs_dedxDiff", loader,
      Binning::Simple(30, 0., 30.), varTruthMuonMatchedTrackFrontrr,
      Binning::Simple(150, -5., 10.), varTruthMuonMatchedTrackFrontdedxDiff,
      spillCut, cut
    )
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum("TruthMuonMatchedTrackEnd_rr_vs_dedx", loader,
      Binning::Simple(26, 0., 26.), varTruthMuonMatchedTrackEndrr,
      Binning::Simple(100, 0., 10.), varTruthMuonMatchedTrackEnddedx,
      spillCut, cut
    )
  );

  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axTruthMuonMatchedTrackStitchedTrackDistance, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axTruthMuonMatchedTrackStitchedTrackPDG, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axTruthMuonMatchedTrackStitchedTrackLength, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axTruthMuonMatchedTrackStitchedTrackCosine, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axTruthMuonMatchedTrackStitchedTrackChi2Muon, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axTruthMuonMatchedTrackStitchedTrackChi2Proton, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axTruthMuonMatchedTrackStitchedTrackChi2Pion, spillCut, cut));

  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axTruthMuonMatchedTrackStitchedShowerDistance, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axTruthMuonMatchedTrackStitchedShowerPDG, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axTruthMuonMatchedTrackStitchedShowerEnergy, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axTruthMuonMatchedTrackStitchedShowerCosine, spillCut, cut));

  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axTruthMuonMatchedTrackChi2Muon, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axTruthMuonMatchedTrackChi2Proton, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axTruthMuonMatchedTrackChi2Pion, spillCut, cut));

  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axTruthMuonMatchedTrackRangePResidualFraction, spillCut, cut));


  //==== Pion

  const HistAxis axChargedPionTruthT(
    "ChargedPionTruthT",
    Binning::Simple(500, 0., 5.),
    varChargedPionTruthT
  );
  const HistAxis axChargedPionTruthLength(
    "ChargedPionTruthLength",
    Binning::Simple(1000, 0., 1000.),
    varChargedPionTruthLength
  );
  const HistAxis axTruthChargedPionMatchedTrackLength(
    "TruthChargedPionMatchedTrackLength",
    Binning::Simple(1000, 0., 1000.),
    varTruthChargedPionMatchedTrackLength
  );
  const HistAxis axTruthChargedPionMatchedTrackNDaughterTracks(
    "TruthChargedPionMatchedTrackNDaughterTracks",
    Binning::Simple(10, 0., 10.),
    varTruthChargedPionMatchedTrackNDaughterTracks
  );
  const HistAxis axTruthChargedPionMatchedTrackNDaughterShowers(
    "TruthChargedPionMatchedTrackNDaughterShowers",
    Binning::Simple(10, 0., 10.),
    varTruthChargedPionMatchedTrackNDaughterShowers
  );
  const HistAxis axTruthChargedPionMatchedTrackFrontLargedEdX(
    "TruthChargedPionMatchedTrackFrontLargedEdX",
    Binning::Simple(100, 0., 100.),
    varTruthChargedPionMatchedTrackFrontLargedEdX
  );


  const HistAxis axTruthChargedPionMatchedTrackStitchedTrackDistance(
    "TruthChargedPionMatchedTrackStitchedTrackDistance",
    Binning::Simple(300,0.,30.),
    varTruthChargedPionMatchedTrackStitchedTrackDistance
  );
  const HistAxis axTruthChargedPionMatchedTrackStitchedTrackPDG(
    "TruthChargedPionMatchedTrackStitchedTrackPDG",
    Binning::Simple(10000, -5000., 5000.),
    varTruthChargedPionMatchedTrackStitchedTrackPDG
  );
  const HistAxis axTruthChargedPionMatchedTrackStitchedTrackLength(
    "TruthChargedPionMatchedTrackStitchedTrackLength",
    Binning::Simple(100, 0., 100.),
    varTruthChargedPionMatchedTrackStitchedTrackLength
  );
  const HistAxis axTruthChargedPionMatchedTrackStitchedTrackCosine(
    "TruthChargedPionMatchedTrackStitchedTrackCosine",
    Binning::Simple(20, -1., 1.),
    varTruthChargedPionMatchedTrackStitchedTrackCosine
  );
  const HistAxis axTruthChargedPionMatchedTrackStitchedTrackChi2Muon(
    "TruthChargedPionMatchedTrackStitchedTrackChi2Muon",
    Binning::Simple(200, 0., 200.),
    varTruthChargedPionMatchedTrackStitchedTrackChi2Muon
  );
  const HistAxis axTruthChargedPionMatchedTrackStitchedTrackChi2Proton(
    "TruthChargedPionMatchedTrackStitchedTrackChi2Proton",
    Binning::Simple(200, 0., 200.),
    varTruthChargedPionMatchedTrackStitchedTrackChi2Proton
  );
  const HistAxis axTruthChargedPionMatchedTrackStitchedTrackChi2Pion(
    "TruthChargedPionMatchedTrackStitchedTrackChi2Pion",
    Binning::Simple(200, 0., 200.),
    varTruthChargedPionMatchedTrackStitchedTrackChi2Pion
  );

  const HistAxis axTruthChargedPionMatchedTrackStitchedShowerDistance(
    "TruthChargedPionMatchedTrackStitchedShowerDistance",
    Binning::Simple(300,0.,30.),
    varTruthChargedPionMatchedTrackStitchedShowerDistance
  );
  const HistAxis axTruthChargedPionMatchedTrackStitchedShowerPDG(
    "TruthChargedPionMatchedTrackStitchedShowerPDG",
    Binning::Simple(10000, -5000., 5000.),
    varTruthChargedPionMatchedTrackStitchedShowerPDG
  );
  const HistAxis axTruthChargedPionMatchedTrackStitchedShowerEnergy(
    "TruthChargedPionMatchedTrackStitchedShowerEnergy",
    Binning::Simple(200, 0., 2.),
    varTruthChargedPionMatchedTrackStitchedShowerEnergy
  );
  const HistAxis axTruthChargedPionMatchedTrackStitchedShowerCosine(
    "TruthChargedPionMatchedTrackStitchedShowerCosine",
    Binning::Simple(20, -1., 1.),
    varTruthChargedPionMatchedTrackStitchedShowerCosine
  );

  const HistAxis axTruthChargedPionMatchedTrackChi2Muon("TruthChargedPionMatchedTrackChi2Muon", binsChi2, varTruthChargedPionMatchedTrackChi2Muon);
  const HistAxis axTruthChargedPionMatchedTrackChi2Proton("TruthChargedPionMatchedTrackChi2Proton", binsChi2, varTruthChargedPionMatchedTrackChi2Proton);
  const HistAxis axTruthChargedPionMatchedTrackChi2Pion("TruthChargedPionMatchedTrackChi2Pion", binsChi2, varTruthChargedPionMatchedTrackChi2Pion);

  const HistAxis axTruthChargedPionMatchedTrackRangePResidualFraction(
    "TruthChargedPionMatchedTrackRangePResidualFraction",
    Binning::Simple(1000, -2., 2.), 
    varTruthChargedPionMatchedTrackRangePResidualFraction
  );

  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axChargedPionTruthT, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axChargedPionTruthLength, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axTruthChargedPionMatchedTrackLength, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axTruthChargedPionMatchedTrackNDaughterTracks, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axTruthChargedPionMatchedTrackNDaughterShowers, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axTruthChargedPionMatchedTrackFrontLargedEdX, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum("TruthChargedPionMatchedTrackFront_rr_vs_dedx", loader,
      Binning::Simple(30, 0., 30.), varTruthChargedPionMatchedTrackFrontrr,
      LargedEdXBinning, varTruthChargedPionMatchedTrackFrontdedx,
      spillCut, cut
    )
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum("TruthChargedPionMatchedTrackFront_rr_vs_dedxDiff", loader,
      Binning::Simple(30, 0., 30.), varTruthChargedPionMatchedTrackFrontrr,
      Binning::Simple(150, -5., 10.), varTruthChargedPionMatchedTrackFrontdedxDiff,
      spillCut, cut
    )
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum("TruthChargedPionMatchedTrackEnd_rr_vs_dedx", loader,
      Binning::Simple(26, 0., 26.), varTruthChargedPionMatchedTrackEndrr,
      Binning::Simple(100, 0., 10.), varTruthChargedPionMatchedTrackEnddedx,
      spillCut, cut
    )
  );

  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axTruthChargedPionMatchedTrackStitchedTrackDistance, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axTruthChargedPionMatchedTrackStitchedTrackPDG, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axTruthChargedPionMatchedTrackStitchedTrackLength, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axTruthChargedPionMatchedTrackStitchedTrackCosine, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axTruthChargedPionMatchedTrackStitchedTrackChi2Muon, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axTruthChargedPionMatchedTrackStitchedTrackChi2Proton, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axTruthChargedPionMatchedTrackStitchedTrackChi2Pion, spillCut, cut));

  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axTruthChargedPionMatchedTrackStitchedShowerDistance, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axTruthChargedPionMatchedTrackStitchedShowerPDG, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axTruthChargedPionMatchedTrackStitchedShowerEnergy, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axTruthChargedPionMatchedTrackStitchedShowerCosine, spillCut, cut));

  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axTruthChargedPionMatchedTrackChi2Muon, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axTruthChargedPionMatchedTrackChi2Proton, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axTruthChargedPionMatchedTrackChi2Pion, spillCut, cut));

  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axTruthChargedPionMatchedTrackRangePResidualFraction, spillCut, cut));

}

//==== 221014_PandoraCosmicTest
void HistoProducer::PandoraClearCosmicTest(SpectrumLoader& loader, SpillCut spillCut, Cut cut){

  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum("VertexRecoX", Binning::Simple(1000, -500., 500.), loader, varVertexRecoX, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum("VertexRecoY", Binning::Simple(400, -200., 200.), loader, varVertexRecoY, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum("VertexRecoZ", Binning::Simple(2000, -1000., 1000.), loader, varVertexRecoZ, spillCut, cut));

  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum("NeutrinoTruthE", Binning::Simple(100, 0., 10.), loader, varNeutrinoTruthE, spillCut, cut));

  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum("MuonTruthP", Binning::Simple(100, 0., 5.), loader, varMuonTruthP, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum("MuonTruthCosineTheta", Binning::Simple(100, -1., 1.), loader, varMuonTruthCosineTheta, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum("MuonTruthNuMICosineTheta", Binning::Simple(100, -1., 1.), loader, varMuonTruthNuMICosineTheta, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum("MuonTruthDirectionX", Binning::Simple(100, -1., 1.), loader, varMuonTruthDirectionX, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum("MuonTruthDirectionY", Binning::Simple(100, -1., 1.), loader, varMuonTruthDirectionY, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum("MuonTruthDirectionZ", Binning::Simple(100, -1., 1.), loader, varMuonTruthDirectionZ, spillCut, cut));

  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum("MuonBestmatchDirectionX", Binning::Simple(100, -1., 1.), loader, varMuonBestmatchDirectionX, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum("MuonBestmatchDirectionY", Binning::Simple(100, -1., 1.), loader, varMuonBestmatchDirectionY, spillCut, cut));
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum("MuonBestmatchDirectionZ", Binning::Simple(100, -1., 1.), loader, varMuonBestmatchDirectionZ, spillCut, cut));

  double xNuScoreMin = -0.2;
  double xNuScoreMax = 1.;
  double dxNuScore = 0.01;
  const Binning binsNuScore = Binning::Simple( int( (xNuScoreMax-xNuScoreMin)/dxNuScore ), xNuScoreMin, xNuScoreMax );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum("NuScore", binsNuScore, loader, varNuScore, spillCut, cut));

  const HistAxis axSliceNuNFinalStatePfos("SliceNuNFinalStatePfos", Binning::Simple(20,0.,20.), varSliceNuNFinalStatePfos);
  const HistAxis axSliceNuNHitsTotal("SliceNuNHitsTotal", Binning::Simple(500,0.,5000.), varSliceNuNHitsTotal);
  const HistAxis axSliceNuVertexY("SliceNuVertexY", Binning::Simple(400, -200., 200.), varSliceNuVertexY);
  const HistAxis axSliceNuWeightedDirZ("SliceNuWeightedDirZ", Binning::Simple(100,-1,1), varSliceNuWeightedDirZ);
  const HistAxis axSliceNuNSpacePointsInSphere("SliceNuNSpacePointsInSphere", Binning::Simple(300,0.,300.), varSliceNuNSpacePointsInSphere);
  const HistAxis axSliceNuEigenRatioInSphere("SliceNuEigenRatioInSphere", Binning::Simple(100,0.,1.), varSliceNuEigenRatioInSphere);
  const HistAxis axSliceCRLongestTrackDirY("SliceCRLongestTrackDirY", Binning::Simple(100,-1,1), varSliceCRLongestTrackDirY);
  const HistAxis axSliceCRLongestTrackDeflection("SliceCRLongestTrackDeflection", Binning::Simple(40,0.,2.), varSliceCRLongestTrackDeflection);
  const HistAxis axSliceCRFracHitsInLongestTrack("SliceCRFracHitsInLongestTrack", Binning::Simple(50,0.,1.), varSliceCRFracHitsInLongestTrack);
  const HistAxis axSliceCRNHitsMax("SliceCRNHitsMax", Binning::Simple(1000,0.,1000.), varSliceCRNHitsMax);
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


}

//==== 221027_PandoraCosmicVertexTest
void HistoProducer::PandoraCosmicVertexTest(SpectrumLoader& loader, SpillCut spillCut, Cut cut){
/*
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum("spillTEST", Binning::Simple(1, 0., 1.), loader, spillTEST, spillCut)
  );
*/

  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum("NuScore", Binning::Simple(100, -1., 1.), loader, varNuScore, spillCut, cut));


}

//==== 221121_CRTPMTMatching
void HistoProducer::CRTPMTMatchingStudy(SpectrumLoader& loader, SpillCut spillCut, Cut cut){

  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum("SpillCount", Binning::Simple(1, 0.,1.), loader, spillvarCountSpill, spillCut)
  );

  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum("OpFlashTime", Binning::Simple(2500, -5.,20.), loader, spillvarOpFlashTime, spillCut)
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum("InTimeOpFlashTime", Binning::Simple(2500, -5.,20.), loader, spillvarInTimeOpFlashTime, spillCut)
  );

  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum("TopCRTTime", Binning::Simple(2500, -5.,20.), loader, spillvarTopCRTTime, spillCut)
  );

  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum("CRTPMTTime", Binning::Simple(1500, -0.1, 0.05), loader, spillvarCRTPMTTime, spillCut)
  );

  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum("OpFlashTimeAllRange", Binning::Simple(6000, -3000.,3000.), loader, spillvarOpFlashTime, spillCut)
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum("TopCRTTimeAllRange", Binning::Simple(6000, -3000.,3000.), loader, spillvarTopCRTTime, spillCut)
  );

  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum("CRTPMTMatchingID", Binning::Simple(10, 0.,10.), loader, spillvarCRTPMTMatchingID, spillCut)
  );

}

//==== 221129_TrackBreakingTest
void HistoProducer::TrackBreakingTest(SpectrumLoader& loader, SpillCut spillCut, Cut cut){

  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum("spillLongestTrackStitchedTrackLength", Binning::Simple(1000, 0., 1000.), loader, spillLongestTrackStitchedTrackLength, spillCut)
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum("spillLongestTrackStitchedTrackDistance", Binning::Simple(1000, 0., 1000.), loader, spillLongestTrackStitchedTrackDistance, spillCut)
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum("spillLongestTrackStitchedTrackClosestMode", Binning::Simple(4, 0., 4.), loader, spillLongestTrackStitchedTrackClosestMode, spillCut)
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum("spillLongestTrackStitchedTrackDistanceSameCryo", Binning::Simple(1000, 0., 1000.), loader, spillLongestTrackStitchedTrackDistanceSameCryo, spillCut)
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum("spillLongestTrackStitchedTrackDistanceOtherCryo", Binning::Simple(1000, 0., 1000.), loader, spillLongestTrackStitchedTrackDistanceOtherCryo, spillCut)
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum("spillLongestTrackStitchedTrackDistanceSameTruthG4ID", Binning::Simple(1000, 0., 1000.), loader, spillLongestTrackStitchedTrackDistanceSameTruthG4ID, spillCut)
  );


}

//==== 221201_Cosmic_to_ClearCosmicTest
void HistoProducer::Cosmic_to_ClearCosmicTest(SpectrumLoader& loader, SpillCut spillCut, Cut cut){

  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum("IsTrueCosmic_vs_IsClearCosmic", loader,
      Binning::Simple(2, 0., 2.), varIsTrueCosmic,
      Binning::Simple(2, 0., 2.), varIsClearCosmic,
      spillCut, cut)
  );


}

void HistoProducer::ThreeTrackAnalysis(SpectrumLoader& loader, SpillCut spillCut, Cut cut){

  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("CountSlice", Binning::Simple(1, 0., 1.), loader, varCountSlice, spillCut, cut) );

  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("TruthE", Binning::Simple(20, 0., 5.), loader, varNeutrinoTruthE, spillCut, cut) );
  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("TruthQ2", Binning::Simple(20, 0., 2.), loader, varTruthQ2, spillCut, cut) );
  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("Truthq0_lab", Binning::Simple(20, 0., 2.), loader, varTruthq0_lab, spillCut, cut) );
  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("Truthmodq_lab", Binning::Simple(20, 0., 2.), loader, varTruthmodq_lab, spillCut, cut) );
  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("TruthW", Binning::Simple(30, 0., 3.), loader, varTruthW, spillCut, cut) );

  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("MuonTruthP", Binning::Simple(20, 0., 2.), loader, varMuonTruthP, spillCut, cut) );
  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("ProtonTruthP", Binning::Simple(20, 0., 2.), loader, varProtonTruthP, spillCut, cut) );

  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("NPrimaryTracks", Binning::Simple(10, 0., 10.), loader, TTAVAR_NPrimaryTracks, spillCut, cut) );

  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("MuonTrackLength", Binning::Simple(500, 0., 500), loader, TTAVAR_MuonTrackLength, spillCut, cut) );
  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("MuonTrackLengthMatchMuon", Binning::Simple(20, 0., 2.), loader, TTAVAR_MuonTrackLengthMatchMuon, spillCut, cut) );
  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("MuonTrackLengthMatchPionPlus", Binning::Simple(20, 0., 2.), loader, TTAVAR_MuonTrackLengthMatchPionPlus, spillCut, cut) );
  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("MuonTrackLengthMatchPionMinus", Binning::Simple(20, 0., 2.), loader, TTAVAR_MuonTrackLengthMatchPionMinus, spillCut, cut) );
  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("MuonTrackLengthMatchProton", Binning::Simple(20, 0., 2.), loader, TTAVAR_MuonTrackLengthMatchProton, spillCut, cut) );
  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("MuonTrackLengthMatchElse", Binning::Simple(20, 0., 2.), loader, TTAVAR_MuonTrackLengthMatchElse, spillCut, cut) );
  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("MuonTrackP", Binning::Simple(20, 0., 2.), loader, TTAVAR_MuonTrackP, spillCut, cut) );
  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("MuonTrackPMatchMuon", Binning::Simple(20, 0., 2.), loader, TTAVAR_MuonTrackPMatchMuon, spillCut, cut) );
  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("MuonTrackPMatchPionPlus", Binning::Simple(20, 0., 2.), loader, TTAVAR_MuonTrackPMatchPionPlus, spillCut, cut) );
  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("MuonTrackPMatchPionMinus", Binning::Simple(20, 0., 2.), loader, TTAVAR_MuonTrackPMatchPionMinus, spillCut, cut) );
  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("MuonTrackPMatchProton", Binning::Simple(20, 0., 2.), loader, TTAVAR_MuonTrackPMatchProton, spillCut, cut) );
  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("MuonTrackPMatchElse", Binning::Simple(20, 0., 2.), loader, TTAVAR_MuonTrackPMatchElse, spillCut, cut) );

  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("ProtonTrackP", Binning::Simple(20, 0., 2.), loader, TTAVAR_ProtonTrackP, spillCut, cut) );
  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("ProtonTrackPMatchMuon", Binning::Simple(20, 0., 2.), loader, TTAVAR_ProtonTrackPMatchMuon, spillCut, cut) );
  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("ProtonTrackPMatchPionPlus", Binning::Simple(20, 0., 2.), loader, TTAVAR_ProtonTrackPMatchPionPlus, spillCut, cut) );
  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("ProtonTrackPMatchPionMinus", Binning::Simple(20, 0., 2.), loader, TTAVAR_ProtonTrackPMatchPionMinus, spillCut, cut) );
  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("ProtonTrackPMatchProton", Binning::Simple(20, 0., 2.), loader, TTAVAR_ProtonTrackPMatchProton, spillCut, cut) );
  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("ProtonTrackPMatchElse", Binning::Simple(20, 0., 2.), loader, TTAVAR_ProtonTrackPMatchElse, spillCut, cut) );

  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("ThirdTrackLength", Binning::Simple(500, 0., 500), loader, TTAVAR_ThirdTrackLength, spillCut, cut) );
  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("ThirdTrackLengthMatchMuon", Binning::Simple(20, 0., 2.), loader, TTAVAR_ThirdTrackLengthMatchMuon, spillCut, cut) );
  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("ThirdTrackLengthMatchPionPlus", Binning::Simple(20, 0., 2.), loader, TTAVAR_ThirdTrackLengthMatchPionPlus, spillCut, cut) );
  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("ThirdTrackLengthMatchPionMinus", Binning::Simple(20, 0., 2.), loader, TTAVAR_ThirdTrackLengthMatchPionMinus, spillCut, cut) );
  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("ThirdTrackLengthMatchProton", Binning::Simple(20, 0., 2.), loader, TTAVAR_ThirdTrackLengthMatchProton, spillCut, cut) );
  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("ThirdTrackLengthMatchElse", Binning::Simple(20, 0., 2.), loader, TTAVAR_ThirdTrackLengthMatchElse, spillCut, cut) );
  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("ThirdTrackP", Binning::Simple(20, 0., 2.), loader, TTAVAR_ThirdTrackP, spillCut, cut) );
  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("ThirdTrackPMatchMuon", Binning::Simple(20, 0., 2.), loader, TTAVAR_ThirdTrackPMatchMuon, spillCut, cut) );
  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("ThirdTrackPMatchPionPlus", Binning::Simple(20, 0., 2.), loader, TTAVAR_ThirdTrackPMatchPionPlus, spillCut, cut) );
  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("ThirdTrackPMatchPionMinus", Binning::Simple(20, 0., 2.), loader, TTAVAR_ThirdTrackPMatchPionMinus, spillCut, cut) );
  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("ThirdTrackPMatchProton", Binning::Simple(20, 0., 2.), loader, TTAVAR_ThirdTrackPMatchProton, spillCut, cut) );
  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("ThirdTrackPMatchElse", Binning::Simple(20, 0., 2.), loader, TTAVAR_ThirdTrackPMatchElse, spillCut, cut) );

}

//==== 230207_TwoTrackAnalysis
void HistoProducer::TwoTrackAnalysis(SpectrumLoader& loader, SpillCut spillCut, Cut cut){

  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("CountSlice", Binning::Simple(1, 0., 1.), loader, varCountSlice, spillCut, cut) );
  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("SliceCRLongestTrackDirY", Binning::Simple(20,-1,1), loader, varSliceCRLongestTrackDirY, spillCut, cut) );

  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("TruthE", Binning::Simple(20, 0., 5.), loader, varNeutrinoTruthE, spillCut, cut) );
  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("TruthQ2", Binning::Simple(20, 0., 2.), loader, varTruthQ2, spillCut, cut) );
  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("Truthq0_lab", Binning::Simple(20, 0., 2.), loader, varTruthq0_lab, spillCut, cut) );
  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("Truthmodq_lab", Binning::Simple(20, 0., 2.), loader, varTruthmodq_lab, spillCut, cut) );
  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("TruthW", Binning::Simple(30, 0., 3.), loader, varTruthW, spillCut, cut) );

  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("TruthMuonMatchedTrackChi2Proton", Binning::Simple(15, 0., 150.), loader, varTruthMuonMatchedTrackChi2Proton, spillCut, cut) );
  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("TruthMuonMatchedTrackChi2Muon", Binning::Simple(15, 0., 150.), loader, varTruthMuonMatchedTrackChi2Muon, spillCut, cut) );

  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("TruthProtonMatchedTrackChi2Proton", Binning::Simple(15, 0., 150.), loader, varTruthProtonMatchedTrackChi2Proton, spillCut, cut) );
  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("TruthProtonMatchedTrackChi2Muon", Binning::Simple(15, 0., 150.), loader, varTruthProtonMatchedTrackChi2Muon, spillCut, cut) );

  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("MuonRecoP", Binning::Simple(50, 0., 5.), loader, varMuonRecoP, spillCut, cut) );
  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("MuonRecoNuMICosineTheta", Binning::Simple(20, -1., 1.), loader, varMuonRecoNuMICosineTheta, spillCut, cut) );
  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("MuonRecoLength", Binning::Simple(200, 0., 2000.), loader, varMuonLength, spillCut, cut) );
  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("MuonChi2Muon", Binning::Simple(15, 0., 150.), loader, varMuonChi2Muon, spillCut, cut) );
  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("MuonChi2Proton", Binning::Simple(15, 0., 150.), loader, varMuonChi2Proton, spillCut, cut) );

  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("ProtonRecoP", Binning::Simple(16, 0., 1.6), loader, varProtonRecoP, spillCut, cut) );
  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("ProtonRecoLength", Binning::Simple(100, 0., 100.), loader, varProtonLength, spillCut, cut) );

  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("NStub", Binning::Simple(5, 0., 5.), loader, varNStub, spillCut, cut) );

  map_cutName_to_vec_Spectrums[currentCutName].push_back( new Spectrum("MuonProtonCosineTheta",  Binning::Simple(20, -1., 1.), loader, varMuonProtonCosineTheta, spillCut, cut) );


}


void HistoProducer::bookTEMP(SpectrumLoader& loader, SpillCut spillCut, Cut cut){

/*
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "spillvarNeutrinoQEEnergyResidual",
      Binning::Simple(4000, -2, 2.),
      loader, spillvarNeutrinoQEEnergyResidual, spillCut)
  );

  const HistAxis axProtonTruthT("ProtonTruthT", Binning::Simple(1000., 0., 1.), varProtonTruthT);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axProtonTruthT, spillCut, cut));
*/

  //const HistAxis axMuonChi2Muon("MuonChi2Muon", Binning::Simple(400, 0., 400), varMuonChi2Muon);
  //map_cutName_to_vec_Spectrums[currentCutName].push_back(new Spectrum(loader, axMuonChi2Muon, spillCut, cut));

/*

  //==== 220831 dedx template comparision

  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum("TruthMuonMatchedTrackEnd_rr_vs_dedx", loader,
      Binning::Simple(26, 0., 26.), varTruthMuonMatchedTrackEndrr,
      Binning::Simple(300, 0., 30.), varTruthMuonMatchedTrackEnddedx,
      spillCut, cut
    )
  );

  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum("TruthProtonMatchedTrackEnd_rr_vs_dedx", loader,
      Binning::Simple(26, 0., 26.), varTruthProtonMatchedTrackEndrr,
      Binning::Simple(300, 0., 30.), varTruthProtonMatchedTrackEnddedx,
      spillCut, cut
    )
  );
*/

/*
  const Binning binsXPosition = Binning::Simple(1000, -500., 500.);
  const Binning binsYPosition = Binning::Simple(400, -200., 200.);
  const Binning binsZPosition = Binning::Simple(2000, -1000., 1000.);

  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "VertexRecoX",
      binsXPosition,
      loader,
      varVertexRecoX,
      spillCut, cut
    )
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "VertexRecoY",
      binsYPosition,
      loader,
      varVertexRecoY,
      spillCut, cut
    )
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "VertexRecoZ",
      binsZPosition,
      loader,
      varVertexRecoZ,
      spillCut, cut
    )
  );

  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "MuonRecoStartX",
      binsXPosition,
      loader,
      varMuonRecoStartX,
      spillCut, cut
    )
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "MuonRecoStartY",
      binsYPosition,
      loader,
      varMuonRecoStartY,
      spillCut, cut
    )
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "MuonRecoStartZ", 
      binsZPosition,
      loader,
      varMuonRecoStartZ,
      spillCut, cut
    )
  );

  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "MuonRecoEndX",
      binsXPosition,
      loader,
      varMuonRecoEndX,
      spillCut, cut
    )
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "MuonRecoEndY",
      binsYPosition,
      loader,
      varMuonRecoEndY,
      spillCut, cut
    )
  );
  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "MuonRecoEndZ",
      binsZPosition,
      loader,
      varMuonRecoEndZ,
      spillCut, cut
    )
  );


  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "MuonRecoDirectionY_vs_CRLongestTrackDirY", loader,
      Binning::Simple(40, -1, 1.), varMuonRecoDirectionY,
      Binning::Simple(40, -1, 1.), varSliceCRLongestTrackDirY,
      spillCut, cut
    )
  );
*/

  map_cutName_to_vec_Spectrums[currentCutName].push_back(
    new Spectrum(
      "NuMuCCEnergy_vs_NMatchedNuMuCCSlice", loader,
      Binning::Simple(50, 0., 5.), spillNuMuCCEnergy,
      Binning::Simple(5, 0., 5.), spillNMatchedNuMuCCSlice,
      spillCut
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
        else if(vec_Spectrums.at(i)->GetBinnings().size()==3){
          TH3 *h = vec_Spectrums.at(i)->ToTH3( vec_Spectrums.at(i)->POT() );
          cout << "[HistoProducer::saveHistograms]     Writing TH3, \"" << hName << "\"" << endl;
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
          TGraphAsymmErrors *grErrorBand = vec_SystEnsembleSpectrumPairs.at(i).second->ErrorBand(1., vec_SystEnsembleSpectrumPairs.at(i).second->POT() );
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

  std::vector<Var> weis;
  weis.reserve(1000);
  for(int i = 0; i < 1000; ++i) weis.push_back(GetUniverseWeight("multisim_Genie", i));
  vec_UniverseWeightsForEachGENIESource.push_back( weis );

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

