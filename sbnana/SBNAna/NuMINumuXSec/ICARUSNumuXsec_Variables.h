#pragma once

#include "sbnana/SBNAna/CRTPMTMatching/ICARUSCRTPMTMatching.h"
#include "sbnana/SBNAna/Cuts/NumuCutsIcarus202106.h"
#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_Contants.h"
#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_TruthMatch_Variables.h"
#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_TruthMatch_Cuts.h"
#include "sbnana/SBNAna/Vars/Vars.h"
#include "sbnana/SBNAna/Vars/NuMIXSecVars.h"
#include "sbnana/SBNAna/Cuts/NuMIXSecRecoEventVars.h"
#include "sbnana/SBNAna/Cuts/NuMIXSecTrueEventVars.h"
#include "sbnana/SBNAna/Cuts/NuMIXSecCuts.h"
#include "sbnana/CAFAna/Core/Cut.h"
#include "sbnana/CAFAna/Core/MultiVar.h"
#include <iostream>
#include <climits>

#include "sbnana/SBNAna/Cuts/NuMIXSecSysts.h"

using namespace std;
using namespace ana;

namespace ICARUSNumuXsec{

  // - PMT
  extern const SpillMultiVar OpFlashFirstTime;
  extern const SpillMultiVar OpFlashTime;
  extern const SpillMultiVar OpFlashTimeAfterSignalSelection;

  // - Recalc chi2

  extern const Var kNuMIRecoMuonChi2MuonPlusMichel5cmShift;
  extern const Var kNuMIRecoMuonChi2MuonPlusMichel5cmNoShift;
  extern const Var kNuMIRecoMuonChi2MuonPlusMichel10cmShift;
  extern const Var kNuMIRecoMuonChi2MuonPlusMichel10cmNoShift;

  extern const Var kNuMIRecoMuonFloatChi2MuonPlusMichel1cm;
  extern const Var kNuMIRecoMuonFloatChi2MuonPlusMichel2cm;
  extern const Var kNuMIRecoMuonFloatChi2MuonPlusMichel3cm;
  extern const Var kNuMIRecoMuonFloatChi2MuonPlusMichel4cm;
  extern const Var kNuMIRecoMuonFloatChi2MuonPlusMichel5cm;
  extern const Var kNuMIRecoMuonFloatChi2MuonPlusMichel10cm;
  extern const Var kNuMIRecoMuonFloatChi2MuonPlusMichel15cm;
  extern const Var kNuMIRecoMuonFloatChi2MuonPlusMichelDebug;

  extern const Var kNuMIRecoMuonFloatChi2MuonPlusMichel15cmDelta;
  extern const Var kNuMIRecoMuonFloatChi2MuonPlusMichelDebugDelta;

  extern const Var kNuMILeadingChargedPionCandidateChi2MuonRecalc0p5;
  extern const Var kNuMILeadingChargedPionCandidateChi2MuonRecalc1p0;

  // - Test
  extern const SpillMultiVar spillvarTest;

  // SpillVar
  extern const SpillVar spillvarCountSpill;
  // Trigger
  extern const SpillVar TriggerWithinGate;
  extern const SpillVar TriggerInfoTriggerType;
  extern const SpillVar TriggerInfoSourceType;
  // - Test
  extern const SpillVar spillvarNTrack;
  extern const SpillVar spillvarNShower;
  // - CRT Hit
  extern const SpillMultiVar spillvarSideCRTHitPe;
  extern const SpillMultiVar spillvarTopCRTHitPe;
  // - Pos
  extern const SpillMultiVar spillvarEastWestCRTHitPosX;
  extern const SpillMultiVar spillvarFlashPosX;
  // - OpFlash
  extern const SpillMultiVar spillvarOpFlashPeakToFirstTime;
  // - CRTHit
  extern const SpillMultiVar spillvarTopCRTHitTime;
  extern const SpillMultiVar spillvarSideCRTHitTime;
  // - PMT-CRT matching extra
  extern const SpillMultiVar spillvarTopCRTPMTTime;
  extern const SpillMultiVar spillvarSideCRTPMTTime;
  // Nu from spillvar
  extern const SpillVar TruthFirstNuEnergy;

  // Var
  extern const Var varNuDirectionX;
  extern const Var varNuDirectionY;
  extern const Var varNuDirectionZ;
  extern const Var varTruthQ2;
  extern const Var varTruthq0_lab;
  extern const Var varTruthmodq_lab;
  extern const Var varTruthW;
  // - Truth n particle
  extern const Var varTruthNProton;
  extern const Var varTruthNNeutron;
  extern const Var varTruthNPip;
  extern const Var varTruthNPim;
  extern const Var varTruthNPi0;
  // - Slice var
  extern const Var varCountSlice;
  extern const Var varIsClearCosmic;
  // - Flash matching
  extern const Var varFMScore;
  extern const Var varFMTime;
  // - NuID
  extern const Var varCRLongestTrackDirY;
  // - Long-enough tracks
  extern const MultiVar varLongTrackDirectionY;
  // - Primary tracks
  extern const MultiVar PrimaryTrackIndices;
  extern const Var NPrimaryTracks;
  // - Longest track
  //   - index
  extern const Var LongestTrackIndex;
  //   - length
  extern const Var LongestTrackLength;
  //   - direction
  extern const Var LongestTrackDirectionX;
  extern const Var LongestTrackDirectionY;
  extern const Var LongestTrackDirectionZ;
  extern const Var LongestTrackDirectionXZ;
  extern const Var LongestTrackForceDownDirectionX;
  extern const Var LongestTrackForceDownDirectionY;
  extern const Var LongestTrackForceDownDirectionZ;
  // - Reco muon
  // - Stub
  extern const Var NStubs;
  extern const MultiVar StubCollectionCharges;

  // - In case of neutrino overlapped with cosmic
  extern const Var NNuPFP;
  extern const Var NCosmicPFP;

  // - Test
  extern const Var SliceTestVar;

  // - For trigger eff study
  extern const Var SliceTruthTime;
  extern const SpillVar NuMuSliceLongestTrackLenForTriggerEff;
  extern const SpillVar InTimeCosmicSliceLongestTrackLenForTriggerEff;

  extern const Var Pass_VtxInFV;
  extern const Var Pass_NotClearCosmic;
  extern const Var Pass_HasMuon;
  extern const Var Pass_HasProton;
  extern const Var Pass_ProtonPCut;
  extern const Var Pass_PrimaryHadronContained;
  extern const Var Pass_NoChargedPionTrack;
  extern const Var Pass_NoNeutralPionShower;
  extern const Var Pass_MuonContained;

  extern const TruthVar kTruth_BNBDefaultWeight;

  extern const Var IsCosmicSlice;

  extern const TruthVar DummyTruthVar;

  extern const Var kNuMINUANCECode;

  extern const SpillVar kNuMIG3ChaseSpillWeightByClosesetNu;

  // LifetimeVariation
  extern const Var kNuMI_MuonMatchedTrack_TrackScore;
  extern const Var kNuMI_MuonMatchedTrack_Chi2Muon;
  extern const Var kNuMI_MuonMatchedTrack_Chi2Proton;
  extern const Var kNuMI_ProtonMatchedTrack_TrackScore;
  extern const Var kNuMI_ProtonMatchedTrack_Chi2Muon;
  extern const Var kNuMI_ProtonMatchedTrack_Chi2Proton;

  // Rock study
  extern const SpillCut kNuMI_HasSignalSelectionSlice;
  extern const SpillVar kNuMI_SignalSelectionSlice_CutType;
  extern const SpillMultiVar kNuMI_SignalSelectionSlices_Others_MuonTrackMatched_genT;
  extern const SpillMultiVar kNuMI_SignalSelectionSlices_Others_MuonTrackMatched_pdg;
  extern const SpillMultiVar kNuMI_SignalSelectionSlices_Others_MuonTrackMatched_interaction_id;

  extern const SpillMultiVar kNuMI_TrueNeutrino_PosX;
  extern const SpillMultiVar kNuMI_TrueNeutrino_PosY;
  extern const SpillMultiVar kNuMI_TrueNeutrino_PosZ;

  extern const SpillVar kNuMI_trigger_within_gate;
  extern const SpillVar kNuMI_NumberOfNeutrinos;
  extern const SpillVar kNuMI_TriggerNeutrino_Idx;
  extern const SpillVar kNuMI_TriggerNeutrino_time;
  extern const SpillVar kNuMI_TriggerNeutrino_PosX;
  extern const SpillVar kNuMI_TriggerNeutrino_PosY;
  extern const SpillVar kNuMI_TriggerNeutrino_PosZ;
  extern const SpillVar kNuMI_HasIntimeCosmic;
  extern const SpillVar kNuMI_NIntimeCosmic;
  extern const SpillVar kNuMI_IntimeCosmicClosestTime;
  extern const SpillVar kNuMI_TriggerNeutrino_PrimaryMuonIdx;
  extern const SpillVar kNuMI_TriggerNeutrino_PrimaryMuon_pdg;
  extern const SpillVar kNuMI_TriggerNeutrino_PrimaryMuon_genE;
  extern const SpillVar kNuMI_TriggerNeutrino_PrimaryMuon_GenX;
  extern const SpillVar kNuMI_TriggerNeutrino_PrimaryMuon_GenY;
  extern const SpillVar kNuMI_TriggerNeutrino_PrimaryMuon_GenZ;
  extern const SpillVar kNuMI_TriggerNeutrino_PrimaryMuon_StartX;
  extern const SpillVar kNuMI_TriggerNeutrino_PrimaryMuon_StartY;
  extern const SpillVar kNuMI_TriggerNeutrino_PrimaryMuon_StartZ;
  extern const SpillVar kNuMI_TriggerNeutrino_PrimaryMuon_EndX;
  extern const SpillVar kNuMI_TriggerNeutrino_PrimaryMuon_EndY;
  extern const SpillVar kNuMI_TriggerNeutrino_PrimaryMuon_EndZ;
  extern const SpillVar kNuMI_TriggerTrueParticle_Idx;
  extern const SpillVar kNuMI_TriggerTrueParticle_pdg;
  extern const SpillVar kNuMI_SignalSelectionSlice_Idx;
  extern const SpillVar kNuMI_SignalSelectionSlice_Truth_time;
  extern const SpillVar kNuMI_SignalSelectionSlice_Truth_vtx_x;
  extern const SpillVar kNuMI_SignalSelectionSlice_Truth_vtx_y;
  extern const SpillVar kNuMI_SignalSelectionSlice_Truth_vtx_z;
  extern const SpillVar kNuMI_MuonTrackMatchedTPIdx;
  extern const SpillVar kNuMI_MuonTrackMatchedTP_genT;
  extern const SpillVar kNuMI_MuonTrackMatchedTP_pdg;
  extern const SpillVar kNuMI_MuonTrackMatchedTP_genE;
  extern const SpillVar kNuMI_ProtonTrackMatchedTPIdx;
  extern const SpillVar kNuMI_ProtonTrackMatchedTP_genT;
  extern const SpillVar kNuMI_ProtonTrackMatchedTP_pdg;
  extern const SpillVar kNuMI_ProtonTrackMatchedTP_genE;



  extern const SpillMultiVar kNuMI_IntimeCosmics_genT;
  extern const SpillMultiVar kNuMI_IntimeCosmics_TimeFromTrig;

  extern const SpillCut kNuMI_Spill_NoIntimeCosmic;
  extern const SpillCut kNuMI_Spill_HasIntimeNu;

}
