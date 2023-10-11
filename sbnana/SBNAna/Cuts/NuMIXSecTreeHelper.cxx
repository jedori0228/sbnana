#include "sbnana/SBNAna/Cuts/NuMIXSecTreeHelper.h"

using namespace ana;

namespace ana{

  std::vector<string> GetNuMITrueTreeLabels(){

    return {
      // Interaction
      "TruePDG/i", "TrueMode/i", "TrueTarget/i", "TrueIsCC/i",
      "TrueIsFHC/i",
      // Weight
      "FluxWeight",
      // Nu E
      "TrueE",
      // Muon
      "TrueMuonP",
      "TrueMuonPt",
      "TrueMuonCos",
      "TrueMuonCosBeam",
      "TrueMuonLength",
      "TrueMuonContained/i",
      // Proton
      "TrueProtonP",
      "TrueProtonPt",
      "TrueProtonCos",
      "TrueProtonLength",
      // Muon+Proton
      "TrueMuonProtonCos",
      // TKI
      "TruedeltaPT",
      "TruedeltaPTx",
      "TruedeltaPTy",
      "TruedeltaalphaT",
      "TruedeltaphiT",
    };

  }
  std::vector<TruthVar> GetNuMITrueTreeVars(){

    return {
      // Interaction
      kTruth_NeutrinoPDG, kTruth_NeutrinoMode, kTruth_Target, kTruth_IsCC,
      kTruth_IsFHC,
      // Weight
      kGetTruthNuMIFluxWeight,
      // Nu E
      kTruth_NeutrinoE,
      // Muon
      kTruth_MuonP,
      kTruth_MuonPt,
      kTruth_MuonNuCosineTheta,
      kTruth_MuonCosThBeam,
      kTruth_MuonLength,
      kTruth_MuonContained,
      // Proton
      kTruth_ProtonP,
      kTruth_ProtonPt,
      kTruth_ProtonNuCosineTheta,
      kTruth_ProtonLength,
      // Muon+Proton
      kTruth_CosThMuonProton,
      // TKI
      kTruth_deltaPT,
      kTruth_deltaPTx,
      kTruth_deltaPTy,
      kTruth_deltaalphaT,
      kTruth_deltaphiT,
    };

  }

  std::vector<std::string> GetNuMIRecoTreeLabels(){

    return {
      // CutType
      "CutType/i",
      // Selection flags for each step
      "SelectionFlag/i",
      // Intercation
      "TruePDG/i", "TrueMode/i", "TrueTarget/i", "TrueIsCC/i",
      "TrueIsFHC/i",
      // Number of primary particles
      "TrueNProton/i", "TrueNNeutron/i",
      "TrueNpip/i", "TrueNpim/i", "TrueNpi0/i",
      "TrueNpipAll/i", "TrueNpimAll/i", "TrueNpi0All/i",
      // True muon
      "TrueMuonContained/i",
      // True proton
      "RecoProtonMatchedToTrueProton/i",
      // Is Signal
      "IsSignal/i",
      // Weight
      "FluxWeight",
      // NuE
      "TrueE",
      // Vertex
      "RecoVtxTPC/i",
      "RecoVtxX", "RecoVtxY", "RecoVtxZ",
      "TrueProdVtxX", "TrueProdVtxY", "TrueProdVtxZ",
      // RecoMuon track notfound/contained/exiting
      "MuonTrackType/i",
      // Split muon cut
      "SplitMuonCut/i",
      // Muon momentum
      "RecoMuonP", "TrueMuonP",
      "RecoMuonRangeP", "RecoMuonMCSP",
      "RecoMuonPt", "TrueMuonPt",
      // Muon length
      "RecoMuonLength", "TrueMuonLength",
      // Muon angle
      "RecoMuonCos", "TrueMuonCos",
      "RecoMuonCosBeam", "TrueMuonCosBeam",
      // Proton momentum
      "RecoProtonP", "TrueProtonP",
      "RecoProtonPt", "TrueProtonPt",
      // Proton length
      "RecoProtonLength", "TrueProtonLength",
      // Proton angle
      "RecoProtonCos", "TrueProtonCos",
      // Muon,Proton angle
      "RecoMuonProtonCos", "TrueMuonProtonCos",
      // TKI
      "RecodeltaPT", "TruedeltaPT",
      "RecodeltaPTx", "TruedeltaPTx",
      "RecodeltaPTy", "TruedeltaPTy",
      "RecodeltaalphaT", "TruedeltaalphaT",
      "RecodeltaphiT", "TruedeltaphiT",
      // SideBand
      // pi+-
      // - True,
      "TrueChargedPionKE",
      "TrueChargedPionLength",
      "TrueChargedPionEndProcess/i",
      "TrueChargedPionHasMichel/i",
      "TrueChargedPionBestMatchedTrackHitCompleteness",
      "TrueChargedPionBestMatchedTrackHitPurity",
      "TrueChargedPionBestMatchedTrackHitPurityChi2Muon",
      "TrueChargedPionBestMatchedTrackHitPurityChi2Proton",
      "TrueChargedPionBestMatchedTrackHitPurityChi2MIP",
      // - Reco
      //   - Track
      "LeadingChargedPionCandidateLength",
      "LeadingChargedPionCandidateNDaughter/i",
      "LeadingChargedPionCandidateMatchedPDG/i",
      "LeadingChargedPionCandidateNCollectionHit/i",
      "LeadingChargedPionCandidateMIPChi2",
      // Michel
      "NMichelCandidates/i",
      "MichelCandidateTrackLength",
      "MichelCandidateTrackChi2Proton",
      "MichelCandidateTrackEnergyDensity",
      "MichelCandidateTrackEnergySum",
      "MichelCandidateTrackMatchedPDG/i",
      "MichelCandidateShowerOpeningAngle",
      "MichelCandidateShowerLength",
      "MichelCandidateShowerEnergySum",
      // pi0
      "LeadingPhotonCandidateE",
      "SecondaryPhotonCandidateE",
      "LeadingPhotonCandidateConvGap",
      "SecondaryPhotonCandidateConvGap",
      "PhotonCandidatesOpeningAngle",
      "DiPhotonMass",
    };

  }
  std::vector<Var> GetNuMIRecoTreeVars(){

    return {
      // CutType
      kNuMICutType,
      // Selection flags for each step
      kNuMIRecoSelectionFlag,
      // Intercation
      kNuMITruePDG, kNuMITrueMode, kNuMITrueTarget, kNuMITrueIsCC,
      kNuMIIsFHC,
      // Number of primary particles
      kNuMITrueNProton, kNuMITrueNNeutron,
      kNuMITrueNpip, kNuMITrueNpim, kNuMITrueNpi0,
      kNuMITrueNpip_All, kNuMITrueNpim_All, kNuMITrueNpi0_All,
      // True muon
      kNuMITrueMuonContained,
      // True proton
      kNuMIRecoProtonMatchedToTrueProton,
      // Is Signal
      kNuMISliceSignalType,
      // Weight
      kGetNuMIFluxWeight,
      // NuE
      kNuMITrueNuE,
      // Vertex
      kNuMIRecoVtxTPC,
      kNuMIRecoVtxX, kNuMIRecoVtxY, kNuMIRecoVtxZ,
      kNuMITrueProdVtxX, kNuMITrueProdVtxY, kNuMITrueProdVtxZ,
      // RecoMuon track notfound/contained/exiting
      kNuMIRecoMuonContained,
      // Split muon cut
      kNuMISplitMuonCut,
      // Muon momentum
      kNuMIMuonCandidateRecoP, kNuMIMuonTrueP,
      kNuMIMuonRecoRangeP, kNuMIMuonRecoMCSP,
      kNuMIRecoMuonPt, kNuMITrueMuonPt,
      // Muon length
      kNuMIRecoMuonLength, kNuMITrueMuonLength,
      // Muon angle
      kNuMIRecoCosThVtx, kNuMITrueCosThVtx,
      kNuMIRecoCosThBeam, kNuMITrueCosThBeam,
      // Proton momentum
      kNuMIProtonCandidateRecoP, kNuMIProtonTrueP,
      kNuMIRecoProtonPt, kNuMITrueProtonPt,
      // Proton legnth
      kNuMIRecoProtonLength, kNuMITrueProtonLength,
      // Proton angle
      kNuMIProtonRecoCosThVtx, kNuMIProtonTrueCosThVtx,
      // Muon,Proton angle
      kNuMIRecoCosThMuP, kNuMITrueCosThMuP,
      // TKI
      kNuMIRecodeltaPT, kNuMITruedeltaPT,
      kNuMIRecodeltaPTx, kNuMITruedeltaPTx,
      kNuMIRecodeltaPTy, kNuMITruedeltaPTy,
      kNuMIRecodeltaalphaT, kNuMITruedeltaalphaT,
      kNuMIRecodeltaphiT, kNuMITruedeltaphiT,
      // SideBand
      // pi+-
      // - True
      kNuMITrueChargedPionKE,
      kNuMITrueChargedPionLength,
      kNuMITrueChargedPionEndProcess,
      kNuMITrueChargedPionHasMichel,
      kNuMINChargedPionBestMatchedTrackHitCompleteness,
      kNuMINChargedPionBestMatchedTrackHitPurity,
      kNuMINChargedPionBestMatchedTrackByHitPurityChi2Muon,
      kNuMINChargedPionBestMatchedTrackByHitPurityChi2Proton,
      kNuMINChargedPionBestMatchedTrackByHitPurityChi2MIP,
      // - Reco
      //   - Track
      kNuMILeadingChargedPionCandidateLength,
      kNuMILeadingChargedPionCandidateNDaughter,
      kNuMILeadingChargedPionCandidateMatchedPDG,
      kNuMILeadingChargedPionCandidateNCollectionHit,
      kNuMILeadingChargedPionCandidateMIPChi2,
      // Michel
      kNuMINMichelCandidates,
      kNuMIMichelCandidateTrackLength,
      kNuMIMichelCandidateTrackChi2Proton,
      kNuMIMichelCandidateTrackEnergyDensity,
      kNuMIMichelCandidateTrackEnergySum,
      kNuMIMichelCandidateTrackMatchedPDG,
      kNuMIMichelCandidateShowerOpeningAngle,
      kNuMIMichelCandidateShowerLength,
      kNuMIMichelCandidateShowerEnergySum,
      // pi0
      kNuMILeadingPhotonCandidateE,
      kNuMISecondaryPhotonCandidateE,
      kNuMILeadingPhotonCandidateConvGap,
      kNuMISecondaryPhotonCandidateConvGap,
      kNuMIPhotonCandidatesOpeningAngle,
      kNuMIDiPhotonMass,
    };

  }

}
