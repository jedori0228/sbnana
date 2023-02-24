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

  //==== Number of truth particles

  extern const Cut cutHasTruthMuon;
  extern const Cut cutHasTruthMuonHasRecoTrack;
  extern const Cut cutHasTruthMuonHasRecoTrackLongEnough;

  extern const Cut cutHasTruthProton;
  extern const Cut cutHasTruthProtonHasRecoTrack;
  extern const Cut cutHasTruthProtonHasRecoStub;

  extern const Cut cutHasTruthChargedPion;
  extern const Cut cutHasTruthChargedPionHasRecoTrack;
  extern const Cut cutHasTruthChargedPionHasRecoTrackLongEnough;

  extern const Cut cutTruthNoPiZero;
  extern const Cut cutTruthNoChargedPion;
  extern const Cut cutTruthNoNeutron;

  //==== CC vs NC

  extern const Cut cutIsCC;
  extern const Cut cutIsNC;

  //==== NuMu-CC categories

  extern const Cut cutIsNuMuCC;
  extern const Cut cutIsNuMuCCQE;
  extern const Cut cutIsNuMuCCRes;
  extern const Cut cutIsNuMuCCMEC;
  extern const Cut cutIsNuMuCCDIS;
  extern const Cut cutIsNuMuCCCoh;
  extern const Cut cutIsNuMuCCCohElastic;

  //==== NuMu-NC categories

  extern const Cut cutIsNuMuNC;

  //==== NuE

  extern const Cut cutIsNuECC;

  //==== Truth kinematic cuts

  extern const Cut cutTFiducial;
  extern const Cut cutTruthMuonTCut;
  extern const Cut cutTruthProtonTCut;

  //==== NuMu+P signal signal definition using truth
  extern const Cut cutIsNuMuCCSignalDef;

  //==== FV

  extern const Cut cutRFiducial;

  extern const Cut cutVertexCryoE;
  extern const Cut cutVertexCryoW;
  extern const Cut cutVertexTPCEE;
  extern const Cut cutVertexTPCEW;
  extern const Cut cutVertexTPCWE;
  extern const Cut cutVertexTPCWW;

  //==== Flash matching

  extern const Cut cutFMScore;
  extern const Cut cutFMTime;
  extern const Cut cutFMTimeDataTemp;

  //==== NuScore

  extern const Cut cutNuScore;
  extern const Cut cutSliceCRLongestTrackDirY;
  extern const Cut cutSliceNuVertexYTop;

  //==== 220811 Vertex test
  extern const Cut cutVertexYPos;
  extern const Cut cutVertexYNeg;

  //==== Longest track
  extern const Cut cutHasLongestTrack;
  extern const Cut cutSmallDirX;
  extern const Cut cutLargeDirX;
  extern const Cut cutLongestTrackShort;
  extern const Cut cutLongestTrackContained;
  extern const Cut cutLongestTrackExiting;

  //==== Muon related

  extern const Cut cutHasMuon;
  extern const Cut cutMuonContained;
  extern const Cut cutTruthMuonContained;
  extern const Cut cutRecoMuonTruthContained;
  //==== TEST
  extern const Cut cutTruthMuonMatchedTrackContained;
  extern const Cut cutMuonMatchedToMuon;
  extern const Cut cutMuonMatchedToProton;
  //==== Michel study
  extern const Cut cutTruthMuonMatchedTrackHasStitchedTrack;
  extern const Cut cutTruthMuonMatchedTrackStitchedTrackMatched;
  extern const Cut cutTruthMuonMatchedTrackHasStitchedShower;
  extern const Cut cutTruthMuonMatchedTrackStitchedShowerMatched;
  extern const Cut cutTruthMuonMatchedTrackIsMichelTagged;
  //==== Michel tagging for muon
  extern const Cut cutMuonIsMichelTagged;
  //==== Pion tagging for muon
  extern const Cut cutMuonIsPionTagged;

  //==== Proton related

  extern const Cut cutHasProton;
  extern const Cut cutProtonContained;
  extern const Cut cutTruthProtonContained;
  extern const Cut cutProtonHighMomentum;
  extern const Cut cutProtonMatchedToMuon;
  extern const Cut cutProtonMatchedToProton;

  extern const Cut cutTruthProtonLargePResidualFraction;

  //==== Charged pion related

  extern const Cut cutTruthChargedPionContained;
  extern const Cut cutTruthChargedPionMatchedTrackContained;
  //==== Michel study
  extern const Cut cutTruthChargedPionMatchedTrackHasStitchedTrack;
  extern const Cut cutTruthChargedPionMatchedTrackStitchedTrackMatched;
  extern const Cut cutTruthChargedPionMatchedTrackHasStitchedShower;
  extern const Cut cutTruthChargedPionMatchedTrackStitchedShowerMatched;

  //==== CRT 

  extern const SpillCut spillcutFDTopCRTHitVeto;
  extern const SpillCut spillcutFDSideCRTHitVeto;

  extern const SpillCut spillcutFDTopCRTHitVetoTestMatching;
  extern const SpillCut spillcutFDSideCRTHitVetoTestMatching;

  //==== CRTPMT matching

  extern const SpillCut spillcutHasValidFlash;
  extern const SpillCut spillcutHasInTimeFlash;
  extern const SpillCut spillcutCRTPMTCosmicByID;
  extern const SpillCut spillcutCRTPMTHasNegativeTOF;
  extern const SpillCut spillcutCRTPMTAllNegativeTOF;

  //==== combined

  extern const Cut cutNominal_ContainedMuon;
  extern const Cut cutNominal_ExitingMuon;

  extern const Cut cutMuonProtonCosineTheta;

  //==== Village
  extern const Cut kGraysProposedSampleCut;
  extern const Cut TTACUT_HasThreePrimaryTracks;
  extern const Cut TTACUT_HasMuonTrack;
  extern const Cut TTACUT_HasProtonTrack;
  extern const Cut TTACUT_HasMuonProtonPion;
  extern const Cut TTACUT_MuonContained;
  extern const Cut TTACUT_MuonExiting;

}

