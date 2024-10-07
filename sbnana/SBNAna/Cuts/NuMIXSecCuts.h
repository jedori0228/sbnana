#pragma once

#include "sbnana/CAFAna/Core/Cut.h"

namespace ana
{
  /// \ref SpillCut on valid trigger
  extern const SpillCut kNuMIValidTrigger;

  /// \ref Good run
  extern const SpillCut kNuMIGoodRun;

  /// \ref Cut on vertex reconstruced in FV
  extern const Cut kNuMIVertexInFV;

  /// \ref Cut on reco slice not tagged as clearly cosmic
  extern const Cut kNuMINotClearCosmic;

  /// \ref Cut on having muon candidate
  extern const Cut kNuMIHasMuonCandidate;

  extern const Var kNuMIMuonCandMatched;

  /// \ref Cut on having proton candidate
  extern const Cut kNuMIHasProtonCandidate;
  extern const Cut kNuMIProtonCandidateRecoPTreshold;

  /// \ref Cut on having contained primary hadrons
  extern const Cut kNuMIAllPrimaryHadronsContained;

  /// \ref Cut aimed at charged pion rejection
  extern const Cut kNuMINoSecondPrimaryMuonlikeTracks;

  /// \ref Cut aimed at pi0 rejection
  extern const Cut kNuMICutPhotons;

  /// Base selection of 1muNp from which the main selection and current sidebands branch
  extern const Cut kNuMISelection_1muNp_Base;

  /// Combined selection \ref Cut for 1muNp0pi with contained+exiting muons, without additional shower cut being used
  extern const Cut kNuMISelection_1muNp0pi_WithoutShowerCut;

  /// Combined selection \ref Cut for 1muNp0pi with contained+exiting muons
  extern const Cut kNuMISelection_1muNp0pi;

  /// \ref Cut aimed at muon candidate containment, if so desired
  extern const Cut kNuMIMuonCandidateContained;

  /// \ref Cut aimed at reconstruction quality (e.g. split tracks)
  extern const Cut kNuMIRejectSplitMuons;

  /// \ref Cut aimed at having TWO photons for better pi0 selection
  extern const Cut kNuMIHasTwoPhotons;

  /// \ref Cut pion sideband
  extern const Cut kNuMIChargedPionSideBand;
  extern const Cut kNuMINeutralPionSideBand;
  extern const Cut kNuMINeutralPion2phSideBand;

  /// \ref CutType; 1=Signal, 2=pi+- sideband, 3=pi0 sideband (0=other)
  extern const Var kNuMICutType;
  extern const Var kNuMICutTypeWithoutShowerCut;

  /// \ref Var for whether it would pass the RejectSplitMuons cut: 1 (passes, no split muon tagged), 0 (fails, split muon tagged)
  extern const Var kNuMIPassesSplitMuonCut;

  /// \ref Signal definitions: Neutrino Neutral Current
  extern const Cut kNuMI_IsSliceNuNC;
  /// \ref Not nu matched: i.e. cosmic, or noise, or not well-matched to an interaction
  extern const Cut kNuMI_IsSlcNotNu;

  /// \ref Check 1muNp0pi using vector of primaries; optionally apply phase space cut
  extern const Var kNuMITrueInteractionInFV;
  bool Is1muNp0pi(const caf::Proxy<caf::SRTrueInteraction>& true_int, bool ApplyPhaseSpcaeCut, bool printouts=false);
  inline bool Is1muNp0piWithPhaseSpaceCut(const caf::Proxy<caf::SRTrueInteraction>& true_int){ return Is1muNp0pi(true_int, true); }
  /// \ref Signal without phase space cut
  extern const Cut kNuMI_1muNp0piStudy_Signal_WithoutPhaseSpaceCut;
  /// \ref Signal with phase space cut = "Signal"
  extern const Cut kNuMI_1muNp0piStudy_Signal_WithPhaseSpaceCut;
  extern const Cut kNuMI_1muNp0piStudy_Signal_WithPhaseSpaceCutWithPrintouts;
  /// \ref TruthCut version of signal def
  extern const TruthCut kTruthCut_IsSignal;
  extern const TruthCut kTruthCut_IsSignalWithoutPhaseSpaceCut;
  extern const TruthVar kTruth_IsSignal;
  /// \ref Signal but fails phase space cut = "out of phase space" (OOPS)
  extern const Cut kNuMI_1muNp0piStudy_Signal_FailPhaseSpaceCut;
  /// \ref CC but NOT "Signal" or "OOPS"
  extern const Cut kNuMI_1muNp0piStudy_OtherNuCC;

  /// \ref Check if it's a neutrino in the fiducial volume
  bool IsNuInFV(const caf::Proxy<caf::SRTrueInteraction>& true_int);
  bool IsNuInAV(const caf::Proxy<caf::SRTrueInteraction>& true_int);

  /// \ref Var for slice type (CutType; 1=Signal, 2=OOPS, 3=OtherCC, 4=NuNC, 5=NotNu)
  extern const Var kNuMISliceSignalType;
  extern const Var kNuMISliceSignalTypeWithPrintouts;

  /// \ref Var for slice type (CutType; 1=Signal,         3=OtherCC, 4=NuNC, 5=NotNu)
  extern const Var kNuMISliceSignalTypeWithoutOOPS;

  // FSI cov
  extern const TruthCut kNuMITrueNuMuCCInFV;

  // mu-pi flip flag
  extern const Var kNuMIIsMuonPionFlipped;
  extern const Var kNuMIIsPionPionSelected;
  extern const Var kNuMIIsMuonPionCorrect;
  extern const Var kNuMI_AhtBY_m3sigma;
  extern const Var kNuMI_BhtBY_p3sigma;
  extern const Var kNuMI_CV1uBY_p3sigma;
  extern const Var kNuMI_CV2uBY_m3sigma;
  extern const Var kNuMI_NonRESBGvpCC2pi_m3sigma;
  extern const Var kNuMI_NonRESBGvpNC2pi_m3sigma;
  extern const Var kNuMI_NonRESBGvnCC2pi_m3sigma;
  extern const Var kNuMI_NonRESBGvnNC2pi_m3sigma;
  extern const Var kNuMI_NonRESBGvbarpCC2pi_m3sigma;
  extern const Var kNuMI_NonRESBGvbarpNC2pi_m3sigma;
  extern const Var kNuMI_NonRESBGvbarnCC2pi_m3sigma;
  extern const Var kNuMI_NonRESBGvbarnNC2pi_m3sigma;

  // NuSyst fake data studies
  // - 2p2h
  extern const Var kNuMI_Lars2p2h_XSecShape_CCMEC;
  extern const Var kNuMI_Lars2p2h_EnergyDependence_CCMEC;
  // - FSI
  extern const Var kNuMI_FSI_hN;
  extern const Var kNuMI_FSI_INCL;
  extern const Var kNuMI_FSI_G4BC;

}
