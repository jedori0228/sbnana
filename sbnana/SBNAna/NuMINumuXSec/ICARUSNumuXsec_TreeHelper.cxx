#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_TreeHelper.h"

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
      "TrueMuonDirX",
      "TrueMuonDirY",
      "TrueMuonDirZ",
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
      kTruth_MuonDirX,
      kTruth_MuonDirY,
      kTruth_MuonDirZ,
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
      // Intercation
      "TruePDG/i", "TrueMode/i", "TrueTarget/i", "TrueIsCC/i",
      "TrueIsFHC/i",
      "TrueQ2",
      // Number of primary particles
      "TrueNProton/i", "TrueNNeutron/i",
      "TrueNProtonThreshold/i",
      "TrueNProtonAll/i",
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
      "TrueIntInFV/i",
      // Stub
      "NStub/i",
      "NExtraStub/i",
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
      // Muon track direction
      "RecoMuonDirX", "RecoMuonDirY", "RecoMuonDirZ",
      "RecoMuonStartX", "RecoMuonEndX",
      "RecoMuonThetaXW_Plane0",
      // Muon track chi2
      "RecoMuonChi2Muon", 
      "RecoMuonChi2Proton",
      // Proton momentum
      "RecoProtonP", "TrueProtonP",
      "RecoProtonPt", "TrueProtonPt",
      // G4Proton for G4 cov
      "TrueG4ProtonP",
      // Proton length
      "RecoProtonLength", "TrueProtonLength",
      // Proton angle
      "RecoProtonCos", "TrueProtonCos",
      // Proton track chi2
      "RecoProtonChi2Muon",
      "RecoProtonChi2Proton",
      // Proton angle
      "RecoProtonThetaXW_Plane0",
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
      "LeadingChargedPionCandidateMatchedEndProcess/i",
      "LeadingChargedPionCandidateMatchedLength",
      "LeadingChargedPionCandidateMatchedDirX",
      "LeadingChargedPionCandidateMatchedDirY",
      "LeadingChargedPionCandidateMatchedDirZ",
      "LeadingChargedPionCandidateNCollectionHit/i",
      "LeadingChargedPionCandidateMIPChi2",
      "LeadingChargedPionCandidateChi2Muon",
      "LeadingChargedPionCandidateChi2Proton",
      "LeadingChargedPionCandidateChi2MuonInd1",
      "LeadingChargedPionCandidateChi2ProtonInd1",
      "LeadingChargedPionCandidateChi2MuonInd2",
      "LeadingChargedPionCandidateChi2ProtonInd2",
      "LeadingChargedPionCandidateMindEdX",
      "LeadingChargedPionCandidateChi2MuonRecalc0p5",
      "LeadingChargedPionCandidateChi2MuonRecalc1p0",
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
      // Nu
      "ENu_Muon_LeadingProton",
      "ENu_Muon_AllProtons"
    };

  }
  std::vector<Var> GetNuMIRecoTreeVars(){

    return {
      // CutType
      kNuMICutType,
      // Intercation
      kNuMITruePDG, kNuMITrueMode, kNuMITrueTarget, kNuMITrueIsCC,
      kNuMIIsFHC,
      kNuMITrueQ2,
      // Number of primary particles
      kNuMITrueNProton, kNuMITrueNNeutron,
      kNuMITrueNProton_Threshold,
      kNuMITrueNProton_All,
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
      kNuMITrueInteractionInFV,
      // Stub
      kNuMINStub,
      kNuMINExtraStub,
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
      // Muon track direction
      kNuMIRecoMuonTrackDirX, kNuMIRecoMuonTrackDirY, kNuMIRecoMuonTrackDirZ,
      kNuMIRecoMuonTrackStartX, kNuMIRecoMuonTrackEndX,
      kNuMIRecoMuonTrackThetaXW_Plane0,
      // Muon track chi2
      kNuMIRecoMuonTrackChi2Muon,
      kNuMIRecoMuonTrackChi2Proton,
      // Proton momentum
      kNuMIProtonCandidateRecoP, kNuMIProtonTrueP,
      kNuMIRecoProtonPt, kNuMITrueProtonPt,
      // G4Proton for G4 cov
      kNuMITrueG4ProtonP,
      // Proton legnth
      kNuMIRecoProtonLength, kNuMITrueProtonLength,
      // Proton angle
      kNuMIProtonRecoCosThVtx, kNuMIProtonTrueCosThVtx,
      // Proton track chi2
      kNuMIRecoProtonTrackChi2Muon,
      kNuMIRecoProtonTrackChi2Proton,
      // Proton angle
      kNuMIRecoProtonTrackThetaXW_Plane0,
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
      kNuMILeadingChargedPionCandidateMatchedEndProcess,
      kNuMILeadingChargedPionCandidateMatchedLength,
      kNuMILeadingChargedPionCandidateMatchedDirX,
      kNuMILeadingChargedPionCandidateMatchedDirY,
      kNuMILeadingChargedPionCandidateMatchedDirZ,
      kNuMILeadingChargedPionCandidateNCollectionHit,
      kNuMILeadingChargedPionCandidateMIPChi2,
      kNuMILeadingChargedPionCandidateChi2Muon,
      kNuMILeadingChargedPionCandidateChi2Proton,
      kNuMILeadingChargedPionCandidateChi2MuonInd1,
      kNuMILeadingChargedPionCandidateChi2ProtonInd1,
      kNuMILeadingChargedPionCandidateChi2MuonInd2,
      kNuMILeadingChargedPionCandidateChi2ProtonInd2,
      kNuMILeadingChargedPionCandidateMindEdX,
      ICARUSNumuXsec::kNuMILeadingChargedPionCandidateChi2MuonRecalc0p5,
      ICARUSNumuXsec::kNuMILeadingChargedPionCandidateChi2MuonRecalc1p0,
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
      // Nu
      kNuMIRecoENu_Muon_LeadingProton,
      kNuMIRecoENu_Muon_AllProtons,
    };

  }

}
