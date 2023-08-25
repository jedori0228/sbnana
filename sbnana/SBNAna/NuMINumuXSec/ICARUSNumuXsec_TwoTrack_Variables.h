#pragma once

#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_Contants.h"
#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_Variables.h"
#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_TruthMatch_Variables.h"
#include "sbnana/SBNAna/Vars/Vars.h"
#include "sbnana/CAFAna/Core/Cut.h"
#include "sbnana/CAFAna/Core/MultiVar.h"
#include <iostream>

using namespace std;
using namespace ana;

namespace ICARUSNumuXsec{

namespace TwoTrack{

  // Muon
  extern const Var MuonTrackIndex;
  extern const Var MuonTrackLength;
  extern const Var MuonTrackDirX;
  extern const Var MuonTrackDirY;
  extern const Var MuonTrackDirZ;
  // truth match
  extern const Var MuonTrackTruthLength;
  extern const Var MuonTrackTruthP;
  extern const Var MuonTrackTruthNuMICosineTheta;
  extern const Var MuonTrackTruthPDG;
  extern const Var MuonTrackTruthMatchedPrimaryType;

  // Proton
  extern const Var ProtonTrackIndex;
  extern const Var ProtonTrackLength;
  extern const Var ProtonTrackDirX;
  extern const Var ProtonTrackDirY;
  extern const Var ProtonTrackDirZ;
  extern const Var ProtonTrackNHitsCollection;
  extern const Var ProtonTrackChi2MuonCollection;
  extern const Var ProtonTrackChi2ProtonCollection;
  // truth match
  extern const Var ProtonTrackTruthLength;
  extern const Var ProtonTrackTruthP;
  extern const Var ProtonTrackTruthNuMICosineTheta;
  extern const Var ProtonTrackTruthPDG;
  extern const Var ProtonTrackTruthMatchedPrimaryType;


  // All non-muon track indices, i.e. hadron candiates
  extern const MultiVar NonMuonTrackIndecies;

  // Charged pion, stopped
  extern const Var StoppedChargedPionTrackIndex;
  extern const Var StoppedChargedPionTrackLength;
  extern const Var StoppedChargedPionTrackP;
  extern const Var StoppedChargedPionTrackNuMICosineTheta;
  extern const Var StoppedChargedPionTrackNuMIToVtxCosineTheta;
  extern const Var StoppedChargedPionTrackChi2MuonCollection;
  extern const Var StoppedChargedPionTrackChi2ProtonCollection;
  // Charged pion, inelastic
  extern const Var InelasticChargedPionTrackIndex;
  extern const Var InelasticChargedPionTrackLength;
  extern const Var InelasticChargedPionTrackNuMICosineTheta;
  extern const Var InelasticChargedPionTrackNuMIToVtxCosineTheta;
  extern const Var InelasticChargedPionTrackChi2MIPCollection;

  // Neutral pion
  extern const MultiVar NeutralPionPhotonShowerIndices;
  extern const MultiVar NeutralPionPhotonShowerConvGaps;
  extern const MultiVar NeutralPionPhotonShowerLengths;
  extern const MultiVar NeutralPionPhotonShowerEnergies;
  extern const MultiVar NeutralPionPhotonShowerdEdxs;
  extern const Var NNeutralPionPhotonShower;
  extern const Var NeutralPionPhotonShowerSumEnergy;
  extern const Var NeutralPionPhotonShowerSumInvariantMass;

  namespace Aux{

    extern const SpillMultiVar TestSpillVar;
    extern const SpillMultiVar ForPrintingTruthInfos;

    extern const Var RelaxedMuonTrackIndex;
    extern const Var RelaxedMuonTrackLength;
    extern const Var RelaxedMuonTrackP;
    extern const Var RelaxedMuonTrackNuMICosineTheta;
    extern const Var RelaxedMuonTrackPosAbsX;
    extern const Var RelaxedMuonTrackDirX;
    extern const Var RelaxedMuonTrackDirY;
    extern const Var RelaxedMuonTrackDirZ;
    extern const Var RelaxedMuonTrackTrackScore;
    extern const Var RelaxedMuonTrackChi2Muon;
    extern const Var RelaxedMuonTrackChi2Proton;
    extern const Var RelaxedMuonTrackChi2MuonCollection;
    extern const Var RelaxedMuonTrackChi2ProtonCollection;
    extern const Var RelaxedMuonTrackCustomChi2MuonCollection;
    extern const Var RelaxedMuonTrackCustomChi2ProtonCollection;
    extern const MultiVar RelaxedMuonTrackCollectionRR;
    extern const MultiVar RelaxedMuonTrackCollectiondEdX;
    extern const MultiVar RelaxedMuonTrackCollectiondQdX;
    // truth match
    extern const Var RelaxedMuonTrackTruthLength;
    extern const Var RelaxedMuonTrackTruthP;
    extern const Var RelaxedMuonTrackTruthPResFrac;
    extern const Var RelaxedMuonTrackTruthOneOverP;
    extern const Var RelaxedMuonTrackTruthOneOverPResFrac;
    extern const Var RelaxedMuonTrackTruthNuMICosineTheta;
    extern const Var RelaxedMuonTrackTruthStartProcess;
    extern const Var RelaxedMuonTrackTruthEndProcess;

    extern const Var RelaxedProtonTrackIndex;
    extern const Var RelaxedProtonTrackLength;
    extern const Var RelaxedProtonTrackP;
    extern const Var RelaxedProtonTrackNuMICosineTheta;
    extern const Var RelaxedProtonTrackDirX;
    extern const Var RelaxedProtonTrackDirY;
    extern const Var RelaxedProtonTrackDirZ;
    extern const Var RelaxedProtonTrackTrackScore;
    extern const Var RelaxedProtonTrackChi2Muon;
    extern const Var RelaxedProtonTrackChi2Proton;
    extern const Var RelaxedProtonTrackChi2MuonCollection;
    extern const Var RelaxedProtonTrackChi2ProtonCollection;
    extern const Var RelaxedProtonTrackCustomChi2MuonCollection;
    extern const Var RelaxedProtonTrackCustomChi2ProtonCollection;
    extern const Var RelaxedProtonTrackNHitCollection;
    extern const MultiVar RelaxedProtonTrackCollectionRR;
    extern const MultiVar RelaxedProtonTrackCollectiondEdX;
    extern const MultiVar RelaxedProtonTrackCollectiondQdX;
    // truth match
    extern const Var RelaxedProtonTrackTruthLength;
    extern const Var RelaxedProtonTrackTruthP;
    extern const Var RelaxedProtonTrackTruthPResFrac;
    extern const Var RelaxedProtonTrackTruthNuMICosineTheta;
    extern const Var RelaxedProtonTrackTruthStartProcess;
    extern const Var RelaxedProtonTrackTruthEndProcess;

    // Pion
    extern const Var RelaxedChargedPionTrackIndex;
    extern const Var RelaxedChargedPionTrackLength;
    extern const Var RelaxedChargedPionTrackP;
    extern const Var RelaxedChargedPionTrackTrackScore;
    extern const Var RelaxedChargedPionTrackNuMICosineTheta;
    extern const Var RelaxedChargedPionTrackCustomChi2MuonCollection;
    extern const Var RelaxedChargedPionTrackCustomChi2ProtonCollection;
    extern const Var RelaxedChargedPionTrackCustomChi2InelasticPionCollection;
    extern const Var RelaxedChargedPionTrackFromVertex;
    extern const Var RelaxedChargedPionTrackNDaughter;
    extern const MultiVar RelaxedChargedPionTrackCollectionRR;
    extern const MultiVar RelaxedChargedPionTrackCollectiondEdX;
    extern const MultiVar RelaxedChargedPionTrackCollectiondQdX;
    extern const MultiVar RelaxedChargedPionTrackCollectionFrontRR;
    // truth match
    extern const Var RelaxedChargedPionTrackTruthLength;
    extern const Var RelaxedChargedPionTrackTruthP;
    extern const Var RelaxedChargedPionTrackTruthPDG;
    extern const Var RelaxedChargedPionTrackTruthStartProcess;
    extern const Var RelaxedChargedPionTrackTruthEndProcess;

  }

} // end namespace TwoTrack

} // end namespace ICARUSNumuXsec
