#pragma once

#include "sbnana/SBNAna/Vars/NuMIXSecVars.h"
#include "sbnana/SBNAna/Cuts/NuMIXSecCuts.h"
#include "sbnana/SBNAna/Vars/NuMIFlux.h"

namespace ana{

  // True
  // - Interaction
  extern const Var kNuMITruePDG; //!< Neutrino pdg
  extern const Var kNuMITrueTarget; //!< Target pdg
  extern const Var kNuMITrueMode; //!< GENIE interaction code (https://internal.dunescience.org/doxygen/namespacesimb.html#a2cce734d1b71408bbc7d98d148ac4360)
  extern const Var kNuMITrueIsCC; //!< IsCC (0:NC, 1:CC, -1:Not neutrino)
  extern const Var kNuMITrueNProton; //!< Number of primary proton
  extern const Var kNuMITrueNNeutron; //!< Number of primary neutron
  extern const Var kNuMITrueNpip; //!< Number of primary pi+
  extern const Var kNuMITrueNpip_All; //!< Number of ALL pi+ (not just primary)
  extern const Var kNuMITrueNpim; //!< Number of primary pi-
  extern const Var kNuMITrueNpim_All; //!< Number of ALL pi- (not just primary)
  extern const Var kNuMITrueNpi0; //!< Number of primary pi0
  extern const Var kNuMITrueNpi0_All; //!< Number of ALL pi0 (not just primary)
  extern const Var kNuMITrueNuE; //!< Neutrino energy
  extern const Var kNuMITrueQ2; //!< Q2
  extern const Var kNuMITrueq0; //!< q0; energy transfer
  extern const Var kNuMITrueq3; //!< q3; momentum transfer
  extern const Var kNuMITruew; //!< w; hadronic mass
  extern const Var kNuMIIsFHC; //!< 0: RHC, 1: FHC
  extern const Var kNuMITrueProdVtxX;
  extern const Var kNuMITrueProdVtxY;
  extern const Var kNuMITrueProdVtxZ;
  // - Muon
  extern const Cut kNuMIHasTrueMuon;
  extern const Var kNuMITrueMuonKE; //!< True muon kinetic energy
  extern const Var kNuMITrueMuonNuCosineTheta; //!< True muon cosine angle w.r.t. neutrino
  extern const Var kNuMITrueMuonContained; //!< 1: contained, 0: not contained (-1: muon not found)
  // - Proton
  extern const Cut kNuMIHasTrueProton;
  extern const Var kNuMITrueProtonKE; //!< True proton kinetic energy
  extern const Var kNuMITrueProtonNuCosineTheta; //!< True muon cosine angle w.r.t. neutrino
  // - Charged pion
  extern const Cut kNuMIHasTrueChargedPion;
  extern const Var kNuMITrueChargedPionKE; //!< True pi+- kinetic energy
  extern const Var kNuMITrueChargedPionEndProcess;
  extern const Var kNuMITrueChargedPionLength;
  extern const Var kNuMITrueChargedPionHasMichel;

  // Reco
  // - Slice
  extern const Var kNuMIRecoVtxTPC;
  extern const Var kNuMIRecoVtxX;
  extern const Var kNuMIRecoVtxY;
  extern const Var kNuMIRecoVtxZ;
  // - Muon
  extern const Var kNuMIRecoMuonContained; //!< 0: Muon candidate track exiting, 1: Muon candidate track contained (-1: no muon candidate)
  extern const Var kNuMIRecoMuonTrackMatchType;
  extern const Cut kNuMIRecoMuonTrackMatchContainedNuMu;
  extern const Var kNuMISplitMuonCut;
  //   - Michel from muon (kTruth_MuonMichelIndex)
  extern const MultiVar kNuMIMuonMichelMatchedPfpIndices;
  extern const Cut kNuMIHasTrueMuonMichel;
  // - Charged pion
  //   - Track Var
  extern const MultiVar kNuMIChargedPionMatchedTrackIndices;
  extern const Var kNuMINChargedPionMatchedTracks;
  extern const MultiVar kNuMIChargedPionMatchedTrackScores;
  extern const MultiVar kNuMIChargedPionMatchedTrackLengths;
  extern const MultiVar kNuMIChargedPionMatchedTrackChi2Muons;
  extern const MultiVar kNuMIChargedPionMatchedTrackChi2Protons;
  extern const Var kNuMINChargedPionBestMatchedTrackByHitCompletenessIdx;
  extern const Var kNuMINChargedPionBestMatchedTrackByHitPurityIdx;
  extern const Var kNuMINChargedPionBestMatchedTrackHitCompleteness;
  extern const Var kNuMINChargedPionBestMatchedTrackHitPurity;
  extern const Var kNuMINChargedPionBestMatchedTrackByHitPurityChi2Muon;
  extern const Var kNuMINChargedPionBestMatchedTrackByHitPurityChi2Proton;
  extern const Var kNuMINChargedPionBestMatchedTrackByHitPurityChi2MIP;
  //   - Shower Var
  extern const MultiVar kNuMIChargedPionMatchedShowerIndices;
  extern const Var kNuMINChargedPionMatchedShowers;
  extern const MultiVar kNuMIChargedPionMatchedShowerScores;
  extern const MultiVar kNuMIChargedPionMatchedShowerGaps;
  extern const MultiVar kNuMIChargedPionMatchedShowerIsPrimaries;
  extern const MultiVar kNuMIChargedPionMatchedShowerEnergies;
  //   - Cut
  extern const Cut kNuMIOtherCCWithChargedPion;
  //   - Michel from pion (kTruth_ChargedPionMichelIndex)
  extern const Var kNuMITrueChargedPionMichelEnergy;
  extern const MultiVar kNuMIChargedPionMichelMatchedPfpIndices;
  extern const MultiVar kNuMIChargedPionMichelMatchedPfpScores;
  extern const MultiVar kNuMIChargedPionMichelMatchedPfpTrackLengths;
  extern const MultiVar kNuMIChargedPionMichelMatchedPfpTrackDistances;
  extern const MultiVar kNuMIChargedPionMichelMatchedPfpTrackIsPrimaries;
  extern const MultiVar kNuMIChargedPionMichelMatchedPfpTrackHitPurities;
  extern const MultiVar kNuMIChargedPionMichelMatchedPfpTrackHitCompletenesses;
  extern const MultiVar kNuMIChargedPionMichelMatchedPfpShowerEnergies;
  extern const MultiVar kNuMIChargedPionMichelMatchedPfpShowerLengths;
  extern const MultiVar kNuMIChargedPionMichelMatchedPfpShowerGaps;
  extern const MultiVar kNuMIChargedPionMichelMatchedPfpShowerOpeningAngles;
  extern const Var kNuMIChargedPionMichelMatchedPfpShowerEnergySum;
  // - Proton
  extern const Var kNuMIRecoProtonMatchedToTrueProton;

  // - Selection enum
  extern const Var kNuMIRecoSelectionFlag;

}
