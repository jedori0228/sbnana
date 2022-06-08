#pragma once

#include "sbnana/SBNAna/Cuts/Cuts.h"
#include "sbnana/SBNAna/Vars/Vars.h"
#include "sbnana/SBNAna/Cuts/VolumeDefinitions.h"

#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_Contants.h"
#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_Utilities.h"
#include "sbnana/SBNAna/Cuts/NumuCutsIcarus202106.h"
#include "TVector3.h"
#include <iostream>

using namespace std;
using namespace ana;

namespace ICARUSNumuXsec{

  //==== Spill variable
  extern const SpillMultiVar spillvarCRTHitTime;
  extern const SpillMultiVar spillvarCRTHitT0;
  extern const SpillMultiVar spillvarCRTHitPosX;
  extern const SpillMultiVar spillvarCRTHitPosY;

  //==== Slice variables

  extern const Var varCountSlice;

  extern const Var varVertexRecoX;
  extern const Var varVertexRecoY;
  extern const Var varVertexRecoZ;

  extern const Var varFMScore;
  extern const Var varFMChargeQ;
  extern const Var varFMLightPE;

  extern const Var varFMChargeCenterX;
  extern const Var varFMChargeCenterY;
  extern const Var varFMChargeCenterZ;

  extern const Var varFMLightCenterX;
  extern const Var varFMLightCenterY;
  extern const Var varFMLightCenterZ;

  extern const Var varFMTime;
  extern const Var varTruthTime;

  extern const Var varSliceTrackNhitsPlane0;
  extern const Var varSliceTrackNhitsPlane1;
  extern const Var varSliceTrackNhitsPlane2;
  extern const Var varSliceShowerNhitsPlane0;
  extern const Var varSliceShowerNhitsPlane1;
  extern const Var varSliceShowerNhitsPlane2;
  extern const Var varSliceTrackChargePlane0;
  extern const Var varSliceTrackChargePlane1;
  extern const Var varSliceTrackChargePlane2;

  extern const MultiVar varAllTrackStartPositionX;
  extern const MultiVar varAllTrackStartPositionY;
  extern const MultiVar varAllTrackStartPositionZ;

  extern const Var varNuScore;
  extern const Var varSliceNuNFinalStatePfos;
  extern const Var varSliceNuNHitsTotal;
  extern const Var varSliceNuVertexY;
  extern const Var varSliceNuWeightedDirZ;
  extern const Var varSliceNuNSpacePointsInSphere;
  extern const Var varSliceNuEigenRatioInSphere;
  extern const Var varSliceCRLongestTrackDirY;
  extern const Var varSliceCRLongestTrackDeflection;
  extern const Var varSliceCRFracHitsInLongestTrack;
  extern const Var varSliceCRNHitsMax;

  //==== GENIE interaction code
  //==== https://internal.dunescience.org/doxygen/namespacesimb.html#a2cce734d1b71408bbc7d98d148ac4360
  extern const Var varGENIEIntCode;

  //==== Truth variables

  extern const Var varNeutrinoTruthE;
  extern const Var varTruthQ2;
  extern const Var varTruthq0_lab;
  extern const Var varTruthmodq_lab;

  extern const Var varTruthNNeutron;
  extern const Var varTruthNPiMinus;
  extern const Var varTruthNPiPlus;
  extern const Var varTruthNChargedPion;
  extern const Var varTruthNPiZero;
  extern const Var varTruthNProton;

  extern const Var varMuonTruthIndex;
  extern const Var varMuonTruthP;
  extern const Var varMuonTruthT;
  extern const Var varMuonTruthCosineTheta;
  extern const Var varMuonTruthNuMICosineTheta;

  extern const Var varProtonTruthIndex;
  extern const Var varProtonTruthP;
  extern const Var varProtonTruthT;
  extern const Var varProtonTruthCosineTheta;
  extern const Var varProtonTruthNuMICosineTheta;

  extern const Var varTruthMuonProtonCosineTheta;

  //==== For a given truth, find the matching Reco track
  //====   For a given true muon (truth_index), find a reco track whose best-matched is this muon

  extern const Var varTruthMuonMatchedTrackIndex;
  extern const Var varTruthMuonMatchedTrackContainedness;
  extern const Var varTruthMuonMatchedTrackChi2Proton;
  extern const Var varTruthMuonMatchedTrackReducedChi2Proton;
  extern const Var varTruthMuonMatchedTrackChi2Muon;
  extern const Var varTruthMuonMatchedTrackReducedChi2Muon;
  extern const Var varTruthMuonMatchedTrackLength;
  extern const Var varTruthMuonMatchedTrackRangeP;
  extern const Var varTruthMuonMatchedTrackMCSP;
  extern const Var varTruthMuonMatchedTrackCombinedP;
  extern const Var varTruthMuonMatchedTrackScore;
  extern const Var varTruthMuonMatchedTrackVertexDistance;

  extern const Var varTruthMuonMatchedShowerIndex;
  extern const Var varTruthMuonMatchedShowerScore;

  //====   For a given true proton (truth_index), find a reco track whose best-matched is this muon

  extern const Var varTruthProtonMatchedTrackIndex;
  extern const Var varTruthProtonMatchedTrackContainedness;
  extern const Var varTruthProtonMatchedTrackChi2Proton;
  extern const Var varTruthProtonMatchedTrackChi2ProtonCollection;
  extern const Var varTruthProtonMatchedTrackReducedChi2Proton;
  extern const Var varTruthProtonMatchedTrackChi2Muon;
  extern const Var varTruthProtonMatchedTrackReducedChi2Muon;
  extern const Var varTruthProtonMatchedTrackRangeP;
  extern const Var varTruthProtonMatchedTrackMCSP;
  extern const Var varTruthProtonMatchedTrackCombinedP;
  extern const Var varTruthProtonMatchedTrackLength;
  extern const Var varTruthProtonMatchedTrackScore;
  extern const Var varTruthProtonMatchedTrackVertexDistance;

  extern const Var varTruthProtonMatchedShowerIndex;
  extern const Var varTruthProtonMatchedShowerScore;

  //====   For a given true proton (truth_index), find a reco stub whose best-matched is this muon

  extern const Var varTruthProtonMatchedStubIndex;
  extern const Var varTruthProtonMatchedStubE;

  //==== Reco variables

  //====   Muon

  extern const MultiVar varAllMuonTrackIndices;
  extern const Var varMuonTrackInd;
  extern const Var varMuonRecoP;
  extern const Var varMuonCaloPlane0P;
  extern const Var varMuonCaloPlane1P;
  extern const Var varMuonCaloPlane2P;
  extern const Var varMuonLength;
  extern const Var varMuonChi2Muon;
  extern const Var varMuonChi2Proton;
  extern const Var varMuonReducedChi2Muon;
  extern const Var varMuonRecoStartX;
  extern const Var varMuonRecoStartY;
  extern const Var varMuonRecoStartZ;
  extern const Var varMuonTrackFromVertex;

  //==== Start from a reco Track, and look at its best match gen-particle.
  //==== This means that the get-particle may not be a true muon
  extern const Var varMuonBestmatchP;
  extern const Var varMuonBestmatchPDG;
  extern const Var varMuonPResidual;
  extern const Var varMuonPResidualFraction;
  extern const Var varMuonCaloPlane0PResidualFraction;
  extern const Var varMuonCaloPlane1PResidualFraction;
  extern const Var varMuonCaloPlane2PResidualFraction;
  extern const Var varMuonRecoCosineTheta;
  extern const Var varMuonRecoNuMICosineTheta;
  extern const Var varMuonBestmatchCosineTheta;
  extern const Var varMuonCosineThetaResidual;
  extern const Var varMuonCosineThetaResidualFraction;

  //==== Proton

  //====   Track-based

  extern const Var varNProtonCandTrack;
  extern const Var varNProtonCandMatched;
  extern const Var varProtonTrackInd;
  extern const Var varProtonCaloP;
  extern const Var varProtonRecoP;
  extern const Var varProtonLength;
  extern const Var varProtonChi2Proton;
  extern const Var varProtonReducedChi2Proton;
  extern const Var varProtonBestmatchP;
  extern const Var varProtonBestmatchPDG;
  extern const Var varProtonPResidual;
  extern const Var varProtonPResidualFraction;
  extern const Var varProtonRecoCosineTheta;
  extern const Var varProtonRecoNuMICosineTheta;
  extern const Var varProtonBestmatchCosineTheta;
  extern const Var varProtonCosineThetaResidual;
  extern const Var varProtonCosineThetaResidualFraction;

  //====   Stub-based

  extern const Var varNStub;

  //==== Neutrino

  extern const Var varNeutrinoCombinedEnergy;
  extern const Var varNeutrinoCombinedEnergyResidual;
  extern const Var varNeutrinoCombinedEnergyResidualFraction;

  //==== https://s3.cern.ch/inspire-prod-files-9/93642a13c46438d97680971700e2013c
  extern const Var varNeutrinoQE;
  extern const Var varNeutrinoQEResidual;
  extern const Var varNeutrinoQEResidualFraction;
  extern const Var varNeutrinoTestEnergy;


}

