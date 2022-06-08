#pragma once

#include "sbnana/SBNAna/Cuts/Cuts.h"
#include "sbnana/SBNAna/Vars/Vars.h"
#include "sbnana/SBNAna/Cuts/VolumeDefinitions.h"

#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_Variables.h"
#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_Utilities.h"

using namespace ana;
using namespace std;

namespace ICARUSNumuXsec{

  //==== Neutrino flavor
  extern const Cut cutIsNuMu;
  extern const Cut cutIsNuE;

  //==== GENIE Interaction code
  extern const Cut cutIsQE;
  extern const Cut cutIsRes;
  extern const Cut cutIsDIS;
  extern const Cut cutIsCoh;
  extern const Cut cutIsCohElastic;
  extern const Cut cutIsElectronScattering;
  extern const Cut cutIsIMDAnnihilation;
  extern const Cut cutIsInverseBetaDecay;
  extern const Cut cutIsGlashowResonance;
  extern const Cut cutIsAMNuGamma;
  extern const Cut cutIsMEC;
  extern const Cut cutIsDiffractive;
  extern const Cut cutIsEM;
  extern const Cut cutIsWeakMix;
  extern const Cut cutIsUnknownInteractionType1;
  extern const Cut cutIsUnknownInteractionType2;
  extern const Cut cutIsUnknownInteractionType3;

  //==== CC vs NC

  extern const Cut cutIsCC;
  extern const Cut cutIsNC;

  //==== NuMu-CC categories

  extern const Cut cutIsNuMuCC;
  extern const Cut cutIsNuMuCCQE;
  extern const Cut cutIsNuMuCCRes;
  extern const Cut cutIsNuMuCCMEC;
  extern const Cut cutIsNuMuCCDIS;

  //==== NuMu-NC categories

  extern const Cut cutIsNuMuNC;

  //==== NuE

  extern const Cut cutIsNuECC;

  //==== Number of truth pios

  extern const Cut cutTruthNoPiZero;
  extern const Cut cutTruthNoChargedPion;
  extern const Cut cutTruthNoNeutron;

  //==== Truth Fiducial cut for the signal definition

  extern const Cut cutTFiducial;

  //==== For CC, applying minimum kinetic energy cut on truth muon/proton
  extern const Cut cutTruthMuonTCut;
  extern const Cut cutTruthProtonTCut;

  //==== NuMu+P signal signal definition using truth
  extern const Cut cutIsNuMuCCSignalDef;

  //==== FV

  extern const Cut cutRFiducial;

  //==== Flash matching

  extern const Cut cutFMScore;
  extern const Cut cutFMTime;

  //==== NuScore

  extern const Cut cutNuScore;
  extern const Cut cutSliceCRLongestTrackDirY;
  extern const Cut cutSliceNuVertexYTop;

  //==== Muon related

  extern const Cut cutHasMuon;
  extern const Cut cutMuonContained;
  extern const Cut cutTruthMuonContained;
  extern const Cut cutRecoMuonTruthContained;
  //==== TEST
  extern const Cut cutZeroMomentum;

  extern const Cut cutMuonMatchedToMuon;
  extern const Cut cutMuonMatchedToProton;

  extern const Cut cutHasProton;
  extern const Cut cutProtonContained;
  extern const Cut cutProtonHighMomentum;
  extern const Cut cutProtonMatchedToMuon;
  extern const Cut cutProtonMatchedToProton;

  //==== CRT 

  extern const SpillCut spillcutFDTopCRTHitVeto;
  extern const SpillCut spillcutFDSideCRTHitVeto;

  extern const SpillCut spillcutFDTopCRTHitVetoTestMatching;
  extern const SpillCut spillcutFDSideCRTHitVetoTestMatching;

}

