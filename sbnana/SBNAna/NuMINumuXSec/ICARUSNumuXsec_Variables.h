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

using namespace std;
using namespace ana;

namespace ICARUSNumuXsec{

  // - PMT
  extern const SpillMultiVar OpFlashFirstTime;
  extern const SpillMultiVar OpFlashTime;

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

}
