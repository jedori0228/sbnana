#pragma once

#include "sbnana/CAFAna/Core/Cut.h"

namespace ana{

  // Cut
  extern const Cut kNuMICRLongestTrackDirYHardCut;

  // Muon
  // - Var
  extern const Var kNuMIRelaxedMuonTrackIdx;
  extern const Var kNuMIRelaxedMuonTrackLength;
  extern const Var kNuMIRelaxedMuonTrackChi2Muon;
  extern const Var kNuMIRelaxedMuonTrackChi2Proton;
  extern const Var kNuMIRelaxedMuonTrackScore;
  extern const Var kNuMIRelaxedMuonTrackMatchedTruthPDG;
  extern const Var kNuMIRelaxedMuonTrackMatchedTruthIntID;
  extern const Var kNuMIRelaxedMuonTrackMatchedTruthContained;
  extern const Var kNuMIRelaxedMuonTrackMatchType;
  // - Cut
  extern const Cut kNuMIRelaxedMuonIsContained;
  extern const Cut kNuMIRelaxedMuonLengthCut;
  extern const Cut kNuMIRelaxedMuonNotIsoChronous;
  extern const Cut kNuMIRelaxedMuonSelection;
  extern const Var kNuMIIsRelaxedMuonSelection;

  // Tag muon
  extern const Var kNuMITagMuonIdx;
  extern const Var kNuMITagMuonLength;

  // Proton
  // - Var
  extern const Var kNuMIRelaxedProtonTrackIdx;
  extern const Var kNuMIRelaxedProtonTrackP;
  extern const Var kNuMIRelaxedProtonTrackLength;
  extern const Var kNuMIRelaxedProtonTrackChi2Muon;
  extern const Var kNuMIRelaxedProtonTrackChi2Proton;
  extern const Var kNuMIRelaxedProtonTrackScore;
  extern const Var kNuMIRelaxedProtonTrackMatchedTruthPDG;
  extern const Var kNuMIRelaxedProtonTrackMatchedTruthIntID;
  extern const Var kNuMIRelaxedProtonTrackMatchedTruthContained;
  extern const Var kNuMIRelaxedProtonTrackMatchType;
  // - Cut
  extern const Cut kNuMIRelaxedProtonIsContained;
  extern const Cut kNuMIRelaxedProtonNotIsoChronous;
  extern const Cut kNuMIRelaxedProtonSelection;
  extern const Var kNuMIIsRelaxedProtonSelection;

  // ChargedPion
  // - Var
  extern const Var kNuMIRelaxedChargedPionTrackIdx;
  extern const Var kNuMIRelaxedChargedPionTrackLength;
  extern const Var kNuMIRelaxedChargedPionTrackChi2Muon;
  extern const Var kNuMIRelaxedChargedPionTrackChi2Proton;
  extern const Var kNuMIRelaxedChargedPionTrackScore;
  extern const Var kNuMIRelaxedChargedPionTrackMatchedTruthPDG;
  extern const Var kNuMIRelaxedChargedPionTrackMatchedTruthIntID;
  extern const Var kNuMIRelaxedChargedPionTrackMatchedTruthContained;
  extern const Var kNuMIRelaxedChargedPionTrackMatchType;
  // - Cut
  extern const Cut kNuMIRelaxedChargedPionIsContained;
  extern const Cut kNuMIRelaxedChargedPionNotIsoChronous;
  extern const Cut kNuMIRelaxedChargedPionSelection;
  extern const Var kNuMIIsRelaxedChargedPionSelection;

  // For skimming
  extern const Cut kNuMIRelaxedPIDCut;


}
