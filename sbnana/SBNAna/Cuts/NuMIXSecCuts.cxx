#include "sbnana/SBNAna/Cuts/NuMIXSecCuts.h"
#include "sbnana/SBNAna/Vars/NuMIXSecVars.h"
#include <algorithm>
#include <iostream>
#include "sbnana/CAFAna/Systs/SBNWeightSysts.h"

namespace ana {

  // Cut on having valid trigger time in approx. beam window
  const SpillCut kNuMIValidTrigger ( [](const caf::SRSpillProxy *sr) {
    double spillTriggerTime = kNuMISpillTriggerTime(sr);
    //return spillTriggerTime > -0.1 && spillTriggerTime < 9.7;
    return spillTriggerTime > -5.;
  });

  /// \ref Good run
  const SpillCut kNuMIGoodRun( [](const caf::SRSpillProxy *sr) {

    if(sr->hdr.ismc) return true;

/*
    // used before 05/07/2024
    static std::vector<unsigned int> GoodRunList = {
8460, 8461, 8462, 8468, 8469, 8470, 8471,
8505, 8506, 8507, 8513, 8514, 8515, 8517, 8518, 8521, 8522, 8525, 8527, 8528, 8529, 8530, 8531, 8552, 8553,
9593, 9594, 9595, 9597, 9599,
9602, 9610, 9642, 9646, 9648, 9649, 9688, 9690, 9691, 9692, 9693, 9694, 9695, 9696, 9699,
9700, 9704, 9715, 9716, 9717, 9721, 9723, 9725, 9726, 9728, 9729, 9730, 9731, 9732, 9733,
9735, 9743, 9744, 9745, 9746, 9747, 9750, 9752, 9753, 9755, 9762, 9763, 9764, 9765, 9781,
9791, 9792, 9794, 9796, 9807, 9834, 9835, 9837, 9838, 9840, 9844, 9847, 9849, 9851, 9854,
9855, 9860, 9862, 9867, 9869, 9870, 9892, 9894, 9896, 9897, 9914, 9919, 9921, 9922, 9924,
9925, 9926, 9941, 9943, 9944, 9945, 9949, 9950, 9951, 9953, 9954, 9956, 9959, 9960, 9972,
9974, 9977, 9979, 9981, 10054, 10059, 10060, 10061, 10064, 10065, 10066, 10084, 10085, 10096, 10097
    };
*/
    // updated 05/07/2024
    static std::vector<unsigned int> GoodRunList = {
      8461, 8462, 8468, 8469, 8470, 8471,
8505, 8506, 8507, 8513, 8514, 8515,             8521, 8522,       8527, 8528, 8529, 8530, 8531, 8552, 8553,
9593, 9594, 9595, 9597, 9599,
9602, 9610, 9642, 9646, 9648, 9649, 9688, 9690, 9691, 9692, 9693, 9694, 9695, 9696, 9699,
9700, 9704, 9715, 9716, 9717, 9721, 9723, 9725, 9726, 9728, 9729, 9730, 9731, 9732, 9733,
9735, 9743, 9744, 9745, 9746, 9747, 9750, 9752, 9753, 9755, 9762, 9763, 9764, 9765, 9781,
9791, 9792, 9794, 9796, 9807, 9834, 9835, 9837, 9838, 9840, 9844, 9847, 9849, 9851, 9854,
9855, 9860, 9862, 9867, 9869, 9870, 9892, 9894, 9896, 9897, 9914, 9919, 9921, 9922, 9924,
9925, 9926,             9944, 9945, 9949, 9950, 9951, 9953, 9954, 9956, 9959, 9960, 9972,
9974, 9977, 9979, 9981, 10054, 10059, 10060, 10061, 10064, 10065, 10066, 10084, 10085, 10096, 10097
    };
    unsigned int RunNum = sr->hdr.run;
    auto it = std::find(GoodRunList.begin(), GoodRunList.end(), RunNum);
    if(it != GoodRunList.end()){
      //std::cout << "RunNum is in GoodRunList" << std::endl;
      return true;
    }
    else{
      //std::cout << "RunNum is not in GoodRunList" << std::endl;
      return false;
    }

  });

  // reco vertex fiducial volume cut
  const Cut kNuMIVertexInFV([](const caf::SRSliceProxy* slc) {
    if ( std::isnan(slc->vertex.x) || std::isnan(slc->vertex.y) || std::isnan(slc->vertex.z) ) return false;
    return isInFV(slc->vertex.x,slc->vertex.y,slc->vertex.z);
  });

  // reco tagged as clear cosmic cut
  const Cut kNuMINotClearCosmic([](const caf::SRSliceProxy* slc) {
    return !slc->is_clear_cosmic;
  });


  // Muon candidate
  const Cut kNuMIHasMuonCandidate([](const caf::SRSliceProxy* slc) {
    return ( kNuMIMuonCandidateIdx(slc) >= 0 );
  });

  const Var kNuMIMuonCandMatched([](const caf::SRSliceProxy* slc) -> int {
    int muonCandIdx = kNuMIMuonCandidateIdx(slc);
    int muonTrueIdx = kTruth_MuonIndex(&slc->truth);
    if(muonCandIdx<0) return -2;
    if(muonTrueIdx<0) return -1;
    if(slc->reco.pfp.at(muonCandIdx).trk.truth.p.G4ID == slc->truth.prim.at(muonTrueIdx).G4ID) return 1;
    else return 0;
  });

  // Proton candidate
  const Cut kNuMIHasProtonCandidate([](const caf::SRSliceProxy* slc) {
    return ( kNuMIProtonCandidateIdx(slc) >= 0 );
  });

  const Cut kNuMIProtonCandidateRecoPTreshold([](const caf::SRSliceProxy* slc) {
    float p(-5.f);

    if ( kNuMIProtonCandidateIdx(slc) >= 0 )
    {
      auto const& trk = slc->reco.pfp.at(kNuMIProtonCandidateIdx(slc)).trk;
      p = trk.rangeP.p_proton;
    }

    return (p > 0.4 && p < 1.0);
  });

  // Set of two cuts largely meant to cut charged pions
  const Cut kNuMIAllPrimaryHadronsContained([](const caf::SRSliceProxy* slc) {
    // Considers all track fits, not just track-like PFPs
    int primaryInd = kNuMIMuonCandidateIdx(slc);
    if ( primaryInd < 0 ) return false;

    unsigned int idxTrk = 0;
    while ( IsValidTrkIdx(slc, idxTrk) ) {
      int thisIdxInt = idxTrk;
      if ( thisIdxInt == primaryInd ) {
        idxTrk+=1;
        continue; // skip the particle which is the muon candidate!
      }
      auto const& trk = slc->reco.pfp.at(idxTrk).trk;
      unsigned int thisIdx = idxTrk;
      idxTrk+=1;

      if ( std::isnan(trk.start.x) || std::isnan(trk.len) || trk.len <= 0. ) continue;
      if ( std::isnan(slc->vertex.x) || std::isnan(slc->vertex.y) || std::isnan(slc->vertex.z) ) return false;
      const float Atslc = std::hypot(slc->vertex.x - trk.start.x,
                                     slc->vertex.y - trk.start.y,
                                     slc->vertex.z - trk.start.z);
      const bool isPrimCandidate = (Atslc < 10. && IsPrimaryPFP(slc,thisIdx));

      if ( !isPrimCandidate ) continue;
      const bool Contained = isContainedVol(trk.end.x,trk.end.y,trk.end.z);
      if ( !Contained ) return false;
    }

    return true;
  });

  const Cut kNuMINoSecondPrimaryMuonlikeTracks([](const caf::SRSliceProxy* slc) {
    int primaryInd = kNuMIMuonCandidateIdx(slc);
    if ( primaryInd < 0 ) return false;

    std::vector<double> chargedpion_indices = kNuMIChargedPionCandidateIdxs(slc);
    if(chargedpion_indices.size()==0) return true;
    else return false;

  });

  // Cut on showers aiming at rejecting pi0
  const Cut kNuMICutPhotons([](const caf::SRSliceProxy* slc) {

    int primaryInd = kNuMIMuonCandidateIdx(slc);
    if ( primaryInd < 0 ) return false;

    int primaryProtonInd = kNuMIProtonCandidateIdx(slc);
    if ( primaryProtonInd < 0 ) return false;

    std::vector<double> photon_indices = kNuMIPhotonCandidateIdxs(slc);
    if(photon_indices.size()==0) return true;
    else return false;

  });

  // Base selection common to side-bands (cuts back on number of entries one needs to carry):
  // For now the sidebands revolve around chg pi and pi0
  const Cut kNuMISelection_1muNp_Base = kNuMIVertexInFV && kNuMINotClearCosmic &&
                                        kNuMIHasMuonCandidate && kNuMIHasProtonCandidate &&
                                        kNuMIProtonCandidateRecoPTreshold &&
                                        kNuMIAllPrimaryHadronsContained;

  // Full selection: 1muNp0pi without caring about containment
  const Cut kNuMISelection_1muNp0pi_WithoutShowerCut = kNuMIVertexInFV && kNuMINotClearCosmic &&                                                /*Preselection*/
                                                       kNuMIHasMuonCandidate && kNuMIHasProtonCandidate && kNuMIProtonCandidateRecoPTreshold && /*Mu, P candidates*/
                                                       kNuMIAllPrimaryHadronsContained && kNuMINoSecondPrimaryMuonlikeTracks;                   /*Esp. chg pi rejection*/

  // Full selection: 1muNp0pi without caring about containment
  const Cut kNuMISelection_1muNp0pi = kNuMIVertexInFV && kNuMINotClearCosmic &&                                                /*Preselection*/
                                      kNuMIHasMuonCandidate && kNuMIHasProtonCandidate && kNuMIProtonCandidateRecoPTreshold && /*Mu, P candidates*/
                                      kNuMIAllPrimaryHadronsContained && kNuMINoSecondPrimaryMuonlikeTracks &&                 /*Esp. chg pi rejection*/
                                      kNuMICutPhotons;                                                                         /*Esp. pi0 rejection*/

  // Muon containment
  const Cut kNuMIMuonCandidateContained([](const caf::SRSliceProxy* slc) {
    int muonCandidate = kNuMIMuonCandidateIdx(slc);

    if ( muonCandidate < 0 ) return false;
    else {
      auto const& trk = slc->reco.pfp.at(kNuMIMuonCandidateIdx(slc)).trk;
      const bool Contained = isContainedVol(trk.end.x,trk.end.y,trk.end.z);
      return Contained;
    }

    return false;
  });

  // Reco quality (meant to help ID split tracks)
  const Cut kNuMIRejectSplitMuons([](const caf::SRSliceProxy* slc) {
    int primaryInd = kNuMIMuonCandidateIdx(slc);
    if ( primaryInd < 0 ) return false;

    unsigned int idxTrk = 0;
    while ( IsValidTrkIdx(slc, idxTrk) ) {
      int thisIdxInt = idxTrk;
      if ( thisIdxInt == primaryInd ) {
        idxTrk+=1;
        continue; // skip the particle which is the muon candidate!
      }
      auto const& trk = slc->reco.pfp.at(idxTrk).trk;
      idxTrk+=1;

      if ( std::isnan(trk.start.x) || std::isnan(trk.len) || trk.len <= 0. ) continue;
      if ( std::isnan(slc->vertex.x) || std::isnan(slc->vertex.y) || std::isnan(slc->vertex.z) ) return false;
      const float Atslc = std::hypot(slc->vertex.x - trk.start.x,
                                     slc->vertex.y - trk.start.y,
                                     slc->vertex.z - trk.start.z);
      const bool isSplitCandidate = (Atslc > 10.);

      if ( !isSplitCandidate || trk.calo[2].nhit < 5 ) continue;
      const float Chi2Proton = trk.chi2pid[2].chi2_proton;
      const float Chi2Muon = trk.chi2pid[2].chi2_muon;
      if ( !(Chi2Proton > 60. && Chi2Muon < 30.) ) continue;

      if ( trk.len > 10. ) return false;
    }

    return true;
  });

  const Cut kNuMIHasTwoPhotons([](const caf::SRSliceProxy* slc) {

    int primaryInd = kNuMIMuonCandidateIdx(slc);
    if ( primaryInd < 0 ) return false;

    int primaryProtonInd = kNuMIProtonCandidateIdx(slc);
    if ( primaryProtonInd < 0 ) return false;

    std::vector<double> photon_indices = kNuMIPhotonCandidateIdxs(slc);
    if(photon_indices.size()>=2) return true;
    else return false;

  });

  /// Pion sideband
  const Cut kNuMIChargedPionSideBand = kNuMIVertexInFV && kNuMINotClearCosmic &&                                                /*Preselection*/
                                       kNuMIHasMuonCandidate && kNuMIHasProtonCandidate && kNuMIProtonCandidateRecoPTreshold && /*Mu, P candidates*/
                                       kNuMIAllPrimaryHadronsContained && /*Hadron contained*/
                                       !kNuMINoSecondPrimaryMuonlikeTracks && /*Charged pion*/
                                       kNuMICutPhotons; /*Esp. pi0 rejection*/

  const Cut kNuMINeutralPionSideBand = kNuMIVertexInFV && kNuMINotClearCosmic &&                                                /*Preselection*/
                                       kNuMIHasMuonCandidate && kNuMIHasProtonCandidate && kNuMIProtonCandidateRecoPTreshold && /*Mu, P candidates*/
                                       kNuMIAllPrimaryHadronsContained && /*Hadron contained*/
                                       kNuMINoSecondPrimaryMuonlikeTracks && /*Esp. chg pi rejection*/
                                       !kNuMICutPhotons; /*Neutral pion*/

  /// Sideband with the 2ph cut --> this is a subset of the NeutralPionSideband
  const Cut kNuMINeutralPion2phSideBand = kNuMIVertexInFV && kNuMINotClearCosmic &&                                                /*Preselection*/
                                          kNuMIHasMuonCandidate && kNuMIHasProtonCandidate && kNuMIProtonCandidateRecoPTreshold && /*Mu, P candidates*/
                                          kNuMIAllPrimaryHadronsContained && /*Hadron contained*/
                                          kNuMINoSecondPrimaryMuonlikeTracks && /*Esp. chg pi rejection*/
                                          kNuMIHasTwoPhotons; /*Neutral pion*/

  /// CutType; 1=Signal, 2=pi+- sideband, 3=pi0 sideband, 0=other
  const Var kNuMICutType([](const caf::SRSliceProxy* slc) -> double {

    if( kNuMISelection_1muNp0pi(slc) ) return 1;
    else if( kNuMIChargedPionSideBand(slc) ) return 2;
    else if( kNuMINeutralPion2phSideBand(slc) ) return 3; // this should be a subset of 4! so if you want 4, use "3 || 4"
    else if( kNuMINeutralPionSideBand(slc) ) return 4;
    else if ( kNuMISelection_1muNp0pi_WithoutShowerCut(slc) ) return 5; // this is everything in "1" except it doesn't place a cut on showers...
    else return 0;

  });

  /// CutType without the showers cut in signal selection; 1=Signal, 2=pi+- sideband, 3=pi0 sideband, 0=other
  const Var kNuMICutTypeWithoutShowerCut([](const caf::SRSliceProxy* slc) -> double {

    if( kNuMISelection_1muNp0pi_WithoutShowerCut(slc) ) return 1;
    else if( kNuMIChargedPionSideBand(slc) ) return 2;
    else if( kNuMINeutralPion2phSideBand(slc) ) return 3; // this should be a subset of 4! so if you want 4, use "3 || 4"
    else if( kNuMINeutralPionSideBand(slc) ) return 4;
    else return 0;

  });

  /// Would pass the split muon cut or not?
  const Var kNuMIPassesSplitMuonCut([](const caf::SRSliceProxy* slc) -> double {
    if( kNuMIRejectSplitMuons(slc) ) return 1;
    else return 0;
  });

  // Signal definitions:
  // Neutrino NC
  const Cut kNuMI_IsSliceNuNC([](const caf::SRSliceProxy* slc) {
    if ( slc->truth.index < 0 ) return false;

    if ( slc->truth.iscc )
      return false; // IS CC

    return true;
  });

  // Not nu matched: i.e. cosmic, or noise, or not well-matched to an interaction
  const Cut kNuMI_IsSlcNotNu([](const caf::SRSliceProxy* slc) {
    return ( slc->truth.index < 0 );
  });
  const Var kNuMITrueInteractionInFV([](const caf::SRSliceProxy* slc) -> int {
    if ( slc->truth.index < 0 ) return -2;
    if ( std::isnan(slc->truth.position.x) || std::isnan(slc->truth.position.y) || std::isnan(slc->truth.position.z) ) return -1;

    if ( isInFV(slc->truth.position.x, slc->truth.position.y, slc->truth.position.z) ) return 1;
    else return 0;

  });
  /// \ref Check 1muNp0pi using vector of primaries
  bool Is1muNp0pi(const caf::Proxy<caf::SRTrueInteraction>& true_int, bool ApplyPhaseSpcaeCut, bool printouts){
    if ( printouts ) std::cout << "CHECKING Is1muNp0pi with printouts = true..." << std::endl;
    if ( true_int.index < 0 ) return false;

    if ( abs(true_int.pdg) != 14 ||
         !true_int.iscc ||
         std::isnan(true_int.position.x) || std::isnan(true_int.position.y) || std::isnan(true_int.position.z) ||
         !isInFV(true_int.position.x, true_int.position.y, true_int.position.z) )
      return false; // not signal

    unsigned int nMu(0), nP(0), nPi(0);
    unsigned int genieNPhotons(0), genieNMesons(0), genieNBaryonsAndPi0(0);
    double maxMomentumP = 0.;
    bool passProtonPCut = false;
    for ( auto const& prim : true_int.prim ) {
      if ( prim.start_process != 0 ) continue;

      double momentum = sqrt( (prim.genp.x*prim.genp.x) + (prim.genp.y*prim.genp.y) + (prim.genp.z*prim.genp.z) );

      bool PassMuonPCut = (momentum > 0.226);
      if ( abs(prim.pdg) == 13 ) {
        if (ApplyPhaseSpcaeCut ? PassMuonPCut : true) nMu+=1;
      }

      if ( abs(prim.pdg) == 2212 ) {
        nP+=1;
        if ( momentum > maxMomentumP ) {
          maxMomentumP = momentum;
          passProtonPCut = (momentum > 0.4 && momentum < 1.);
        }
      }

      if ( abs(prim.pdg) == 111 || abs(prim.pdg) == 211 ) nPi+=1;
      // CHECK A SIMILAR DEFINITION AS MINERVA FOR EXTRA REJECTION OF UNWANTED THINGS IN SIGNAL DEFN.
      if ( abs(prim.pdg) == 22 && prim.startE > 0.01 ) genieNPhotons+=1;
      else if ( abs(prim.pdg) == 211 || abs(prim.pdg) == 321 || abs(prim.pdg) == 323 ||
                prim.pdg == 111 || prim.pdg == 130 || prim.pdg == 310 || prim.pdg == 311 ||
                prim.pdg == 313 || abs(prim.pdg) == 221 || abs(prim.pdg) == 331 ) genieNMesons+=1;
      else if ( prim.pdg == 3112 || prim.pdg == 3122 || prim.pdg == 3212 || prim.pdg == 3222 ||
                prim.pdg == 4112 || prim.pdg == 4122 || prim.pdg == 4212 || prim.pdg == 4222 ||
                prim.pdg == 411 || prim.pdg == 421 || prim.pdg == 111 ) genieNBaryonsAndPi0+=1;
    }

    // if this would have passed signal cuts before ... without phase space consideration
    if ( printouts && !( nMu!=1 || nP==0 || nPi > 0 ) ) {
      std::cout << "---> THIS SLICE WOULD HAVE PASSED SIGNAL CRITERIA BEFORE!!" << std::endl;
    }

    if ( nMu!=1 || nP==0 || nPi > 0 || genieNPhotons > 0 || genieNMesons > 0 || genieNBaryonsAndPi0 > 0 ) {
      if ( printouts ) {
        std::cout << "---> !! This slice does not pass the current signal criteria." << std::endl;
        std::cout << "---> NPhotons > 10 MeV: " << genieNPhotons << ", NMesons: " << genieNMesons << ", NBaryonsPi0: " << genieNBaryonsAndPi0 << std::endl;
      }
      return false;
    }
    else if ( ApplyPhaseSpcaeCut ) return passProtonPCut; 

    return true;

  }
  /// \ref Signal without phase space cut
  const Cut kNuMI_1muNp0piStudy_Signal_WithoutPhaseSpaceCut([](const caf::SRSliceProxy* slc) {
    return Is1muNp0pi(slc->truth, false);
  });
  /// \ref Signal with phase space cut = "Signal"
  const Cut kNuMI_1muNp0piStudy_Signal_WithPhaseSpaceCut([](const caf::SRSliceProxy* slc) {
    return Is1muNp0pi(slc->truth, true);
  });
  const Cut kNuMI_1muNp0piStudy_Signal_WithPhaseSpaceCutWithPrintouts([](const caf::SRSliceProxy* slc) {
    return Is1muNp0pi(slc->truth, true, true);
  });
  /// \ref TruthCut version of signal def
  const TruthCut kTruthCut_IsSignal([](const caf::SRTrueInteractionProxy* nu) {
    return Is1muNp0pi(*nu, true);
  });
  const TruthCut kTruthCut_IsSignalWithoutPhaseSpaceCut([](const caf::SRTrueInteractionProxy* nu) {
    return Is1muNp0pi(*nu, false);
  });
  const TruthVar kTruth_IsSignal([](const caf::SRTrueInteractionProxy *nu) -> int {
    if( kTruthCut_IsSignal(nu) ) return 1;
    else return 0;
  });


  /// \ref Signal but fails phase space cut = "out of phase space" (OOPS)
  const Cut kNuMI_1muNp0piStudy_Signal_FailPhaseSpaceCut = kNuMI_1muNp0piStudy_Signal_WithoutPhaseSpaceCut && // signal,
                                                           !kNuMI_1muNp0piStudy_Signal_WithPhaseSpaceCut; // but fail phase cut
  /// \ref CC but NOT "Signal" or "OOPS"
  const Cut kNuMI_1muNp0piStudy_OtherNuCC([](const caf::SRSliceProxy* slc) {
    if ( slc->truth.index < 0 ) return false; // not Nu
    if ( !slc->truth.iscc ) return false; // not CC

    if ( kNuMI_1muNp0piStudy_Signal_WithoutPhaseSpaceCut(slc) ) return false; // covered by signal without phase space cut
    return true;
  });

  /// CutType; 1=Signal, 2=OOPS, 3=OtherCC, 4=NuNC, 5=NotNu
  /// can be expanded further
  const Var kNuMISliceSignalType([](const caf::SRSliceProxy* slc) -> int {
    if ( kNuMI_1muNp0piStudy_Signal_WithPhaseSpaceCut(slc) ) return 1; // Signal (with phase space cut)
    else if ( kNuMI_1muNp0piStudy_Signal_FailPhaseSpaceCut(slc) ) return 2; // Signal but out of phase space cut (OOPS)
    else if ( kNuMI_1muNp0piStudy_OtherNuCC(slc) ) return 3; // CC but not (signal without phase space cut)
    else if ( kNuMI_IsSliceNuNC(slc) ) return 4; // NC
    else if ( kNuMI_IsSlcNotNu(slc) ) return 5; // Not nu-slice Cosmic
    else return 0;
  });

  const Var kNuMISliceSignalTypeWithPrintouts([](const caf::SRSliceProxy* slc) -> int {
    if ( kNuMI_1muNp0piStudy_Signal_WithPhaseSpaceCutWithPrintouts(slc) ) return 1; // Signal (with phase space cut)
    else if ( kNuMI_1muNp0piStudy_Signal_FailPhaseSpaceCut(slc) ) return 2; // Signal but out of phase space cut (OOPS)
    else if ( kNuMI_1muNp0piStudy_OtherNuCC(slc) ) return 3; // CC but not (signal without phase space cut)
    else if ( kNuMI_IsSliceNuNC(slc) ) return 4; // NC
    else if ( kNuMI_IsSlcNotNu(slc) ) return 5; // Not nu-slice Cosmic
    else return 0;
  });

  /// \ref Var for slice type (CutType; 1=Signal,         3=OtherCC, 4=NuNC, 5=NotNu)
  const Var kNuMISliceSignalTypeWithoutOOPS([](const caf::SRSliceProxy* slc) -> int {
    if ( kNuMI_1muNp0piStudy_Signal_WithPhaseSpaceCut(slc) ) return 1; // Signal (with phase space cut)
    else if ( kNuMI_1muNp0piStudy_Signal_FailPhaseSpaceCut(slc) ) return 3; // Signal but out of phase space cut (OOPS)
    else if ( kNuMI_1muNp0piStudy_OtherNuCC(slc) ) return 3; // CC but not (signal without phase space cut)
    else if ( kNuMI_IsSliceNuNC(slc) ) return 4; // NC
    else if ( kNuMI_IsSlcNotNu(slc) ) return 5; // Not nu-slice Cosmic
    else return 0;
  });

  /// \ref Check if the nu interaction is in the fiducial volume only
  bool IsNuInFV(const caf::Proxy<caf::SRTrueInteraction>& true_int){
    if ( true_int.index < 0 ) return false;

    if ( std::isnan(true_int.position.x) || std::isnan(true_int.position.y) || std::isnan(true_int.position.z) ||
         !isInFV(true_int.position.x, true_int.position.y, true_int.position.z) )
      return false; // not signal

    return true;
  }

  /// \ref Check if the nu interaction is in the fiducial volume only
  bool IsNuInAV(const caf::Proxy<caf::SRTrueInteraction>& true_int){
    if ( true_int.index < 0 ) return false;

    if ( std::isnan(true_int.position.x) || std::isnan(true_int.position.y) || std::isnan(true_int.position.z) ||
         !isInAV(true_int.position.x, true_int.position.y, true_int.position.z) )
      return false; // not signal

    return true;
  }


  // FSI cov
  const TruthCut kNuMITrueNuMuCCInFV([](const caf::SRTrueInteractionProxy* nu) {

    if ( abs(nu->pdg) != 14 ||
         !nu->iscc ||
         std::isnan(nu->position.x) || std::isnan(nu->position.y) || std::isnan(nu->position.z) ||
         !isInFV(nu->position.x, nu->position.y, nu->position.z) ) return false;
    else return true;
  });

  // mu-pi flip flag
  const Var kNuMIIsMuonPionFlipped([](const caf::SRSliceProxy* slc) -> int {
    int muon_pdg = kNuMIRecoMuonCandidateTruePDG(slc);
    int pion_pdg = kNuMIRecoPionCandidateTruePDG(slc);

    if( abs(muon_pdg)==211 && abs(pion_pdg)==13 ) return 1;
    else return 0;
  });
  const Var kNuMIIsPionPionSelected([](const caf::SRSliceProxy* slc) -> int {
    int muon_pdg = kNuMIRecoMuonCandidateTruePDG(slc);
    int pion_pdg = kNuMIRecoPionCandidateTruePDG(slc);

    if( abs(muon_pdg)==211 && abs(pion_pdg)==211 ) return 1;
    else return 0;
  });
  const Var kNuMIIsMuonPionCorrect([](const caf::SRSliceProxy* slc) -> int {
    int muon_pdg = kNuMIRecoMuonCandidateTruePDG(slc);
    int pion_pdg = kNuMIRecoPionCandidateTruePDG(slc);

    if( abs(muon_pdg)==13 && abs(pion_pdg)==211 ) return 1;
    else return 0;
  });

  const Var kNuMI_AhtBY_m3sigma = GetUniverseWeight("GENIEReWeight_ICARUS_v2_multisigma_AhtBY", 4);
  const Var kNuMI_BhtBY_p3sigma = GetUniverseWeight("GENIEReWeight_ICARUS_v2_multisigma_BhtBY", 5);
  const Var kNuMI_CV1uBY_p3sigma = GetUniverseWeight("GENIEReWeight_ICARUS_v2_multisigma_CV1uBY", 5);
  const Var kNuMI_CV2uBY_m3sigma = GetUniverseWeight("GENIEReWeight_ICARUS_v2_multisigma_CV2uBY", 4);
  const Var kNuMI_NonRESBGvpCC2pi_m3sigma = GetUniverseWeight("GENIEReWeight_ICARUS_v2_multisigma_NonRESBGvpCC2pi", 4);
  const Var kNuMI_NonRESBGvpNC2pi_m3sigma = GetUniverseWeight("GENIEReWeight_ICARUS_v2_multisigma_NonRESBGvpNC2pi", 4);
  const Var kNuMI_NonRESBGvnCC2pi_m3sigma = GetUniverseWeight("GENIEReWeight_ICARUS_v2_multisigma_NonRESBGvnCC2pi", 4);
  const Var kNuMI_NonRESBGvnNC2pi_m3sigma = GetUniverseWeight("GENIEReWeight_ICARUS_v2_multisigma_NonRESBGvnNC2pi", 4);
  const Var kNuMI_NonRESBGvbarpCC2pi_m3sigma = GetUniverseWeight("GENIEReWeight_ICARUS_v2_multisigma_NonRESBGvbarpCC2pi", 4);
  const Var kNuMI_NonRESBGvbarpNC2pi_m3sigma = GetUniverseWeight("GENIEReWeight_ICARUS_v2_multisigma_NonRESBGvbarpNC2pi", 4);
  const Var kNuMI_NonRESBGvbarnCC2pi_m3sigma = GetUniverseWeight("GENIEReWeight_ICARUS_v2_multisigma_NonRESBGvbarnCC2pi", 4);
  const Var kNuMI_NonRESBGvbarnNC2pi_m3sigma = GetUniverseWeight("GENIEReWeight_ICARUS_v2_multisigma_NonRESBGvbarnNC2pi", 4);

  // NuSyst fake data studies
  // - 2p2h
  const Var kNuMI_Lars2p2h_XSecShape_CCMEC = GetUniverseWeight("GENIEReWeight_SBNNuSyst_GENIE_multisigma_XSecShape_CCMEC", 1);
  const Var kNuMI_Lars2p2h_EnergyDependence_CCMEC = GetUniverseWeight("GENIEReWeight_SBNNuSyst_GENIE_multisigma_EnergyDependence_CCMEC", 1);
  // - FSI
  const Var kNuMI_FSI_hN = GetUniverseWeight("FSIReweight_SBNNuSyst_FSI_hNReweight_multisigma_FSIReweight", 1);
  const Var kNuMI_FSI_INCL = GetUniverseWeight("FSIReweight_SBNNuSyst_FSI_INCLReweight_multisigma_FSIReweight", 1);
  const Var kNuMI_FSI_G4BC = GetUniverseWeight("FSIReweight_SBNNuSyst_FSI_G4BCReweight_multisigma_FSIReweight", 1);
  // - LQCD Zexp
  const Var kNuMI_LQCDZExpFit = GetUniverseWeight("GENIEReWeight_SBNNuSyst_LQCDZExpFit_multisim_ZExpAVariationResponse", 1);

}
