#pragma once

#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_TwoTrack_Variables.h"
#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_Cuts.h"

using namespace ana;
using namespace std;

namespace ICARUSNumuXsec{

namespace TwoTrack{

  // Test
  extern const SpillMultiVar TestVar;

  // Primray tracks
  extern const Cut HasTwoPrimaryTracks;
  extern const Cut HasOnlyTwoPrimaryTracks;

  // Reco muon track
  extern const Cut HasMuonTrack;
  extern const Cut MuonTrackContained;
  extern const Cut MuonTrackExiting;
  extern const Var MuonTrackType;
  extern const Cut MuonTrackOneMeter;
  extern const Cut MuonTrackBelowBlindP;
  // - Truth matching
  extern const Cut MuonTrackTruthContainedNuMuon;
  extern const Cut MuonTrackTruthExitingNuMuon;
  extern const Cut MuonTrackTruthCosmicMuon;
  extern const Cut MuonTrackTruthStoppingProton;
  extern const Cut MuonTrackTruthInelProton;
  extern const Cut MuonTrackTruthOtherProton;
  extern const Cut MuonTrackTruthStoppingChargedPion;
  extern const Cut MuonTrackTruthInelChargedPion;
  extern const Cut MuonTrackTruthOtherChargedPion;
  extern const Cut MuonTrackTruthOther;
  // - Comparing truth match to primary
  extern const Cut MuonTrackTruthMatchedPrimaryMuon;

  // Reco proton track
  extern const Cut HasProtonTrack;
  extern const Cut ProtonTrackPCut;
  // - Comparing truth match to primary
  extern const Cut ProtonTrackTruthMatchedPrimaryProton;

  // Hadron (non-muon) contrained
  extern const Cut HadronContained;

  // pion tagging
  // - charged
  extern const Cut HasChargedPionTrack;
  extern const Cut HasStoppedChargedPionTrack;
  extern const Cut StoppedChargedPionTrackLongEnough;
  extern const Cut HasInelasticChargedPionTrack;
  extern const Cut InelasticChargedPionTrackLongEnough;
  // - neutral
  extern const Cut HasNeutralPionPhotonShower;

  // Signal def
  extern const Cut SignalDef;
  extern const Cut SignalMuonContained;
  extern const Cut SignalMuonExiting;
  extern const Var IsSignal;

  // Sideband def
  extern const Cut ChargedPionSideBand;
  extern const Cut NeutralPionSideBand;

  // cut type
  extern const Var CutType;
  extern const Cut IsForTree;

  namespace Aux{
    extern const Cut HasRelaxedMuonTrack;
    extern const Cut HasRelaxedProtonTrack;
    extern const Cut HasRelaxedChargedPionTrack;
    extern const Cut RelaxedMuonTrackContained;
    extern const Cut RelaxedProtonTrackContained;
    extern const Cut RelaxedChargedPionTrackContained;

    extern const Cut RelaxedMuonTrackTruthContainedNuMuon;
    extern const Cut RelaxedMuonTrackTruthExitingNuMuon;
    extern const Cut RelaxedMuonTrackTruthCosmicMuon;
    extern const Cut RelaxedMuonTrackTruthStoppingProton;
    extern const Cut RelaxedMuonTrackTruthInelProton;
    extern const Cut RelaxedMuonTrackTruthOtherProton;
    extern const Cut RelaxedMuonTrackTruthStoppingChargedPion;
    extern const Cut RelaxedMuonTrackTruthInelChargedPion;
    extern const Cut RelaxedMuonTrackTruthOtherChargedPion;
    extern const Cut RelaxedMuonTrackTruthOther;

    extern const Cut RelaxedProtonTrackTruthContainedNuMuon;
    extern const Cut RelaxedProtonTrackTruthExitingNuMuon;
    extern const Cut RelaxedProtonTrackTruthCosmicMuon;
    extern const Cut RelaxedProtonTrackTruthStoppingProton;
    extern const Cut RelaxedProtonTrackTruthInelProton;
    extern const Cut RelaxedProtonTrackTruthOtherProton;
    extern const Cut RelaxedProtonTrackTruthStoppingChargedPion;
    extern const Cut RelaxedProtonTrackTruthInelChargedPion;
    extern const Cut RelaxedProtonTrackTruthOtherChargedPion;
    extern const Cut RelaxedProtonTrackTruthOther;

    extern const Cut RelaxedChargedPionTrackTruthContainedNuMuon;
    extern const Cut RelaxedChargedPionTrackTruthExitingNuMuon;
    extern const Cut RelaxedChargedPionTrackTruthCosmicMuon;
    extern const Cut RelaxedChargedPionTrackTruthStoppingProton;
    extern const Cut RelaxedChargedPionTrackTruthInelProton;
    extern const Cut RelaxedChargedPionTrackTruthOtherProton;
    extern const Cut RelaxedChargedPionTrackTruthStoppingChargedPion;
    extern const Cut RelaxedChargedPionTrackTruthInelChargedPion;
    extern const Cut RelaxedChargedPionTrackTruthOtherChargedPion;
    extern const Cut RelaxedChargedPionTrackTruthOther;

    // some presels
    // - muon
    extern const Cut RelaxedMuonTrackIsochronous;
    extern const Cut RelaxedMuonTrackDriftDirection;
    extern const Cut RelaxedMuonTrackLengthCut;
    extern const Cut RelaxedMuonTrackTruthPResFracTail;
    // - proton
    extern const Cut RelaxedProtonTrackIsochronous;
    extern const Cut RelaxedProtonTrackDriftDirection;
    extern const Cut RelaxedProtonTrackHasChi2MuonExcess;
    extern const Cut RelaxedProtonTrackHasChi2MuonDeficit;
    extern const Cut RelaxedProtonTrackShort;
    // - charged pion
    extern const Cut RelaxedChargedPionTrackIsochronous;

  }

} // end namespace TwoTrack

} // end namespace ICARUSNumuXsec
