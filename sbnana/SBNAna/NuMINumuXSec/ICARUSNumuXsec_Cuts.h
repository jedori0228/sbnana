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

  extern const Cut cutIsCC;
  extern const Cut cutIsNC;

  //==== Nue

  extern const Cut cutIsNuECC;

  //==== NC

  extern const Cut cutIsNuMuNC;

  extern const Cut cutTruthNoPiZero;
  extern const Cut cutTruthNoChargedPion;
  extern const Cut cutTruthNoNeutron;

  extern const Cut cutTFiducial;

  //==== NuMu+P signal signal definition using truth
  extern const Cut cutIsNuMuCCSignalDef;

  //==== FV

  extern const Cut cutRFiducial;

  //==== FMScore

  extern const Cut cutFMScore;

  //==== NuScore

  extern const Cut cutNuScore;

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
  extern const Cut cutProtonMatchedToMuon;
  extern const Cut cutProtonMatchedToProton;

  //==== CRT 

  extern const SpillCut spillcutSideCRTHitVetoFD;

  extern const SpillCut spillcutSideCRTHitVetoFDUpdatedT0;
  extern const SpillCut spillcutCRTHitVetoFDUpdatedT0;

}

