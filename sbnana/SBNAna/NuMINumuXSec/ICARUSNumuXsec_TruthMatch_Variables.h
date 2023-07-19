#pragma once

#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_Contants.h"
#include "sbnana/SBNAna/Vars/Vars.h"
#include "sbnana/CAFAna/Core/Cut.h"
#include "sbnana/CAFAna/Core/MultiVar.h"
#include <iostream>

using namespace std;
using namespace ana;

namespace ICARUSNumuXsec{

namespace TruthMatch{

  // For a given true muon (truth_index), find a reco track whose best-matched is this particle
  extern const Var TruthMuonIndex;
  extern const Var TruthMuonLength;
  extern const Var TruthMuonNuCosineTheta;
  extern const Var TruthMuonKE;
  extern const Var TruthMuonP;
  extern const Var TruthMuonMatchedTrackIndex;
  extern const Var TruthMuonMatchedTrackContainedness;
  extern const Var TruthMuonMatchedTrackChi2Proton;
  extern const Var TruthMuonMatchedTrackChi2Muon;
  extern const Var TruthMuonMatchedTrackChi2Pion;
  extern const Var TruthMuonMatchedTrackCustomChi2InelasticPionCollection;
  // - Michel from muon
  extern const Var TruthMuonMichelIndex;
  extern const Var TruthMuonMichelStartProcess;
  extern const Var TruthMuonMichelKE;
  extern const Var TruthMuonMichelKEFromDecay;
  extern const Var TruthMuonMichelKEFromCapture;
  extern const Var TruthMuonMichelMatchedShowerIndex;
  extern const Var TruthMuonMichelMatchedShowerKE;
  extern const Var TruthMuonMichelMatchedShowerDistanceFromMuonEnd;

  // For a given true proton (truth_index), find a reco track whose best-matched is this particle
  extern const Var TruthProtonIndex;
  extern const Var TruthProtonLength;
  extern const Var TruthProtonNuCosineTheta;
  extern const Var TruthProtonKE;
  extern const Var TruthProtonP;
  extern const Var TruthProtonMatchedTrackIndex;
  extern const Var TruthProtonMatchedTrackContainedness;
  extern const Var TruthProtonMatchedTrackChi2Proton;
  extern const Var TruthProtonMatchedTrackChi2Muon;
  extern const Var TruthProtonMatchedTrackChi2Pion;
  extern const Var TruthProtonMatchedTrackCustomChi2InelasticPionCollection;

  // muon+proton
  extern const Var TruthMuonProtonCosineTheta;

  // For a given true charged pion (truth_index), find a reco track whose best-matched is this particle
  extern const Var TruthChargedPionIndex;
  extern const Var TruthChargedPionExistence;
  extern const Var TruthChargedPionLength;
  extern const Var TruthChargedPionKE;
  extern const Var TruthChargedPionContainedness;
  extern const Var TruthChargedPionMatchedTrackIndex;
  extern const Var TruthChargedPionMatchedTrackScore;
  extern const Var TruthChargedPionMatchedTrackContainedness;
  extern const Var TruthChargedPionMatchedTrackEndProcess;
  extern const Var TruthChargedPionMatchedTrackCustomChi2MuonCollection;
  extern const Var TruthChargedPionMatchedTrackCustomChi2ProtonCollection;
  extern const Var TruthChargedPionMatchedTrackCustomChi2InelasticPionCollection;
  extern const Var TruthChargedPionMatchedShowerIndex;
  extern const Var TruthChargedPionHasReco;
  extern const Var TruthChargedPionNDaughters;
  extern const MultiVar TruthChargedPionDaughterPDGs;
  // - Michel from pion
  extern const Var TruthChargedPionMichelIndex;
  extern const Var TruthChargedPionMichelExistence;
  extern const Var TruthChargedPionMichelStartProcess;
  extern const Var TruthChargedPionMichelKE;
  extern const Var TruthChargedPionMichelMatchedTrackIndex;
  extern const Var TruthChargedPionMichelMatchedShowerIndex;
  extern const Var TruthChargedPionMichelMatchedShowerKE;
  extern const Var TruthChargedPionMichelMatchedShowerDistanceFromMuonEnd;
  extern const Var TruthChargedPionMichelHasReco;

  // test
  extern const SpillMultiVar TruthChargedPionMichelMatchedSlice;

  namespace TKI{
    extern const Var deltaPT;
    extern const Var deltaPTx;
    extern const Var deltaPTy;
  }

} // end namespace TruthMatch

} // end namespace ICARUSNumuXsec
