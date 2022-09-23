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
  extern const SpillVar spillvarCountSpill;
  extern const SpillMultiVar spillvarCRTHitTime;
  extern const SpillMultiVar spillvarCRTHitT0;
  extern const SpillMultiVar spillvarCRTHitPosX;
  extern const SpillMultiVar spillvarCRTHitPosY;
  extern const SpillMultiVar spillvarMCNeutrinoE;
  extern const SpillVar spillvarNSlice;
  extern const SpillVar spillvarNTrack;
  extern const SpillMultiVar spillNuDirectionX;
  extern const SpillMultiVar spillNuDirectionY;
  extern const SpillMultiVar spillNuDirectionZ;
  //====   Spill-based Truth neutrino information
  extern const SpillMultiVar spillvarNeutrinoQEEnergyResidual;

  extern const SpillVar spillTEST;

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
  extern const MultiVar varAllTrackEndPositionX;
  extern const MultiVar varAllTrackEndPositionY;
  extern const MultiVar varAllTrackEndPositionZ;
  extern const MultiVar varAllTrackDirectionX;
  extern const MultiVar varAllTrackDirectionY;
  extern const MultiVar varAllTrackDirectionZ;

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
  //====   Scattering
  extern const Var varNeutrinoTruthE;
  extern const Var varTruthQ2;
  extern const Var varTruthq0_lab;
  extern const Var varTruthmodq_lab;
  //====   Vertex
  extern const Var varTruthVtxX;
  extern const Var varTruthVtxY;
  extern const Var varTruthVtxZ;
  //====     Vertex residual
  extern const Var varVtxResidualX;
  extern const Var varVtxResidualY;
  extern const Var varVtxResidualZ;
  //====   Number of particles
  extern const Var varTruthNNeutron;
  extern const Var varTruthNPiMinus;
  extern const Var varTruthNPiPlus;
  extern const Var varTruthNChargedPion;
  extern const Var varTruthNPiZero;
  extern const Var varTruthNProton;
  //====   Truth Muon
  extern const Var varMuonTruthIndex;
  extern const Var varMuonTruthP;
  extern const Var varMuonTruthT;
  extern const Var varMuonTruthCosineTheta;
  extern const Var varMuonTruthNuMICosineTheta;
  extern const Var varMuonTruthLength;
  //====   Truth Proton
  extern const Var varProtonTruthIndex;
  extern const Var varProtonTruthP;
  extern const Var varProtonTruthT;
  extern const Var varProtonTruthCosineTheta;
  extern const Var varProtonTruthNuMICosineTheta;
  //====   Truth Charged Pion
  extern const Var varChargedPionTruthIndex;
  extern const Var varChargedPionTruthP;
  extern const Var varChargedPionTruthT;
  extern const Var varChargedPionTruthLength;
  //====   Truth Muon+Proton
  extern const Var varTruthMuonProtonCosineTheta;

  //==== For a given true muon (truth_index), find a reco track whose best-matched is this muon
  extern const Var varTruthMuonMatchedTrackIndex;
  extern const Var varTruthMuonMatchedTrackContainedness;
  extern const Var varTruthMuonMatchedTrackChi2Proton;
  extern const Var varTruthMuonMatchedTrackReducedChi2Proton;
  extern const Var varTruthMuonMatchedTrackChi2Muon;
  extern const Var varTruthMuonMatchedTrackReducedChi2Muon;
  extern const Var varTruthMuonMatchedTrackChi2Pion;
  extern const Var varTruthMuonMatchedTrackReducedChi2Pion;
  //====   Matched reco track positions
  extern const Var varTruthMuonMatchedTrackEndPositionX;
  extern const Var varTruthMuonMatchedTrackEndPositionY;
  extern const Var varTruthMuonMatchedTrackEndPositionZ;
  //====   Matched reco angle
  extern const Var varTruthMuonMatchedTrackNuMICosineTheta;
  //====   Matched reco length
  extern const Var varTruthMuonMatchedTrackLength;
  //====   Matched reco momentum
  extern const Var varTruthMuonMatchedTrackRangeP;
  extern const Var varTruthMuonMatchedTrackMCSP;
  extern const Var varTruthMuonMatchedTrackCombinedP;
  //====   Matched reco momentum residual
  extern const Var varTruthMuonMatchedTrackRangePResidual;
  extern const Var varTruthMuonMatchedTrackMCSPResidual;
  extern const Var varTruthMuonMatchedTrackCombinedPResidual;
  //====   Matched reco momentum residual fraction
  extern const Var varTruthMuonMatchedTrackRangePResidualFraction;
  extern const Var varTruthMuonMatchedTrackMCSPResidualFraction;
  extern const Var varTruthMuonMatchedTrackCombinedPResidualFraction;
  //====   Matched reco track other properties
  extern const Var varTruthMuonMatchedTrackScore;
  extern const Var varTruthMuonMatchedTrackVertexDistance;
  //====   dEdX vs rr
  extern const MultiVar varTruthMuonMatchedTrackEnddedx;
  extern const MultiVar varTruthMuonMatchedTrackEndrr;
  extern const MultiVar varTruthMuonMatchedTrackFrontdedx;
  extern const MultiVar varTruthMuonMatchedTrackFrontdedxTemplate;
  extern const MultiVar varTruthMuonMatchedTrackFrontdedxDiff;
  extern const MultiVar varTruthMuonMatchedTrackFrontrr;
  extern const Var varTruthMuonMatchedTrackFrontLargedEdX;
  //====   Michel study (closest)
  extern const Var varTruthMuonMatchedTrackStitchedTrackIndex;
  extern const Var varTruthMuonMatchedTrackStitchedTrackDistance;
  extern const Var varTruthMuonMatchedTrackStitchedTrackPDG;
  extern const Var varTruthMuonMatchedTrackStitchedTrackLength;
  extern const Var varTruthMuonMatchedTrackStitchedTrackCosine;
  extern const Var varTruthMuonMatchedTrackStitchedTrackChi2Muon;
  extern const Var varTruthMuonMatchedTrackStitchedTrackChi2Proton;
  extern const Var varTruthMuonMatchedTrackStitchedTrackChi2Pion;
  extern const Var varTruthMuonMatchedTrackStitchedShowerIndex;
  extern const Var varTruthMuonMatchedTrackStitchedShowerDistance;
  extern const Var varTruthMuonMatchedTrackStitchedShowerPDG;
  extern const Var varTruthMuonMatchedTrackStitchedShowerEnergy;
  extern const Var varTruthMuonMatchedTrackStitchedShowerCosine;
  //====   Count all close-enought track/shower daughters
  extern const Var varTruthMuonMatchedTrackNDaughterTracks;
  extern const Var varTruthMuonMatchedTrackNDaughterShowers;
  //====   Matched reco shower
  extern const Var varTruthMuonMatchedShowerIndex;
  extern const Var varTruthMuonMatchedShowerScore;

  //==== For a given true proton (truth_index), find a reco track whose best-matched is this proton
  extern const Var varTruthProtonMatchedTrackIndex;
  extern const Var varTruthProtonMatchedTrackContainedness;
  extern const Var varTruthProtonMatchedTrackChi2Proton;
  extern const Var varTruthProtonMatchedTrackReducedChi2Proton;
  extern const Var varTruthProtonMatchedTrackChi2Muon;
  extern const Var varTruthProtonMatchedTrackReducedChi2Muon;
  extern const Var varTruthProtonMatchedTrackChi2Pion;
  extern const Var varTruthProtonMatchedTrackReducedChi2Pion;
  //====   Matched reco length
  extern const Var varTruthProtonMatchedTrackLength;
  //====   Matched reco momentum
  extern const Var varTruthProtonMatchedTrackRangeP;
  extern const Var varTruthProtonMatchedTrackMCSP;
  extern const Var varTruthProtonMatchedTrackCombinedP;
  //====   Matched reco momentum residual
  extern const Var varTruthProtonMatchedTrackRangePResidual;
  extern const Var varTruthProtonMatchedTrackMCSPResidual;
  extern const Var varTruthProtonMatchedTrackCombinedPResidual;
  //====   Matched reco momentum residual fraction
  extern const Var varTruthProtonMatchedTrackRangePResidualFraction;
  extern const Var varTruthProtonMatchedTrackMCSPResidualFraction;
  extern const Var varTruthProtonMatchedTrackCombinedPResidualFraction;
  //====   Matched reco track other properties
  extern const Var varTruthProtonMatchedTrackScore;
  extern const Var varTruthProtonMatchedTrackVertexDistance;
  //====   dEdX vs rr
  extern const MultiVar varTruthProtonMatchedTrackEnddedx;
  extern const MultiVar varTruthProtonMatchedTrackEndrr;
  //====   Matched reco shower
  extern const Var varTruthProtonMatchedShowerIndex;
  extern const Var varTruthProtonMatchedShowerScore;
  //====   Matched stub
  extern const Var varTruthProtonMatchedStubIndex;
  extern const Var varTruthProtonMatchedStubE;
  extern const Var varTruthProtonMatchedStubLength;

  //==== For a given true charged pion (truth_index), find a reco track whose best-matched is this charged pion
  extern const Var varTruthChargedPionMatchedTrackIndex;
  extern const Var varTruthChargedPionMatchedTrackContainedness;
  extern const Var varTruthChargedPionMatchedTrackChi2Proton;
  extern const Var varTruthChargedPionMatchedTrackReducedChi2Proton;
  extern const Var varTruthChargedPionMatchedTrackChi2Muon;
  extern const Var varTruthChargedPionMatchedTrackReducedChi2Muon;
  extern const Var varTruthChargedPionMatchedTrackChi2Pion;
  extern const Var varTruthChargedPionMatchedTrackReducedChi2Pion;
  //====   Matched reco length
  extern const Var varTruthChargedPionMatchedTrackLength;
  //====   Matched reco momentum
  extern const Var varTruthChargedPionMatchedTrackRangeP;
  //====   Matched reco momentum residual fraction
  extern const Var varTruthChargedPionMatchedTrackRangePResidualFraction;
  //====   dEdX vs rr
  extern const MultiVar varTruthChargedPionMatchedTrackEnddedx;
  extern const MultiVar varTruthChargedPionMatchedTrackEndrr;
  extern const MultiVar varTruthChargedPionMatchedTrackFrontdedx;
  extern const MultiVar varTruthChargedPionMatchedTrackFrontdedxTemplate;
  extern const MultiVar varTruthChargedPionMatchedTrackFrontdedxDiff;
  extern const MultiVar varTruthChargedPionMatchedTrackFrontrr;
  extern const Var varTruthChargedPionMatchedTrackFrontLargedEdX;
  //====   Decay study
  extern const Var varTruthChargedPionMatchedTrackStitchedTrackIndex;
  extern const Var varTruthChargedPionMatchedTrackStitchedTrackDistance;
  extern const Var varTruthChargedPionMatchedTrackStitchedTrackPDG;
  extern const Var varTruthChargedPionMatchedTrackStitchedTrackLength;
  extern const Var varTruthChargedPionMatchedTrackStitchedTrackCosine;
  extern const Var varTruthChargedPionMatchedTrackStitchedTrackChi2Muon;
  extern const Var varTruthChargedPionMatchedTrackStitchedTrackChi2Proton;
  extern const Var varTruthChargedPionMatchedTrackStitchedTrackChi2Pion;
  extern const Var varTruthChargedPionMatchedTrackStitchedShowerIndex;
  extern const Var varTruthChargedPionMatchedTrackStitchedShowerDistance;
  extern const Var varTruthChargedPionMatchedTrackStitchedShowerPDG;
  extern const Var varTruthChargedPionMatchedTrackStitchedShowerEnergy;
  extern const Var varTruthChargedPionMatchedTrackStitchedShowerCosine;
  //====   Count all close-enought track/shower daughters
  extern const Var varTruthChargedPionMatchedTrackNDaughterTracks;
  extern const Var varTruthChargedPionMatchedTrackNDaughterShowers;

  //==== Reco variables

  //====   Muon

  extern const MultiVar varAllMuonTrackIndices;
  extern const Var varNMuonCandTrack;
  extern const Var varMuonTrackInd;
  extern const Var varMuonRecoP;
  extern const Var varMuonCaloPlane0P;
  extern const Var varMuonCaloPlane1P;
  extern const Var varMuonCaloPlane2P;
  extern const Var varMuonLength;
  extern const Var varMuonChi2Muon;
  extern const Var varMuonChi2Proton;
  extern const Var varMuonChi2Pion;
  extern const Var varMuonChi2MuonReEval;
  extern const Var varMuonReducedChi2Muon;
  extern const Var varMuonRecoStartX;
  extern const Var varMuonRecoStartY;
  extern const Var varMuonRecoStartZ;
  extern const Var varMuonRecoEndX;
  extern const Var varMuonRecoEndY;
  extern const Var varMuonRecoEndZ;
  extern const Var varMuonRecoTrackFromVertex;
  extern const Var varMuonRecoDirX;
  extern const Var varMuonRecoDirY;
  extern const Var varMuonRecoDirZ;
  extern const Var varMuonRecoForceDownDirX;
  extern const Var varMuonRecoForceDownDirY;
  extern const Var varMuonRecoForceDownDirZ;
  //==== dedx 
  extern const MultiVar varMuonTrackCalodedx;
  extern const MultiVar varMuonTrackCalorr;
  //==== Start from a reco Track, and look at its best match gen-particle.
  //==== This means that the get-particle may not be a true muon
  extern const Var varMuonBestmatchP;
  extern const Var varMuonBestmatchPDG;
  extern const Var varMuonRecoCosineTheta;
  extern const Var varMuonRecoNuMICosineTheta;
  extern const Var varMuonBestmatchCosineTheta;

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
  extern const Var varProtonRecoCosineTheta;
  extern const Var varProtonRecoNuMICosineTheta;
  extern const Var varProtonBestmatchCosineTheta;

  //====   Stub-based

  extern const Var varNStub;

  //==== Muon+Proton

  extern const Var varMuonProtonCosineTheta;

  //==== Neutrino

  //====   Energy by Emu+Kp+Eb
  extern const Var varNeutrinoCombinedEnergy;
  //====   Matched residual
  extern const Var varNeutrinoCombinedEnergyResidual;
  extern const Var varNeutrinoCombinedEnergyResidualFraction;

  //==== https://s3.cern.ch/inspire-prod-files-9/93642a13c46438d97680971700e2013c
  extern const Var varNeutrinoQEEnergy;
  extern const Var varNeutrinoQEEnergyResidual;
  extern const Var varNeutrinoTestEnergy;


}

