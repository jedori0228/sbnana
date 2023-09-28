#include "sbnana/SBNAna/Cuts/NuMIXSecRecoEventVars.h"

namespace ana{

  // Neutrino pdg
  const Var kNuMITruePDG([](const caf::SRSliceProxy* slc) -> int {
    if ( slc->truth.index < 0 ) return 0;
    return kTruth_NeutrinoPDG(&slc->truth);
  });
  // Target pdg
  const Var kNuMITrueTarget([](const caf::SRSliceProxy* slc) -> int {
    if ( slc->truth.index < 0 ) return 0;
    return kTruth_Target(&slc->truth);
  });
  // GENIE interaction code (https://internal.dunescience.org/doxygen/namespacesimb.html#a2cce734d1b71408bbc7d98d148ac4360)
  const Var kNuMITrueMode([](const caf::SRSliceProxy* slc) -> int {
    if ( slc->truth.index < 0 ) return -1;
    return kTruth_NeutrinoMode(&slc->truth);
  });
  // IsCC (0:NC, 1:CC, -1:Not neutrino)
  const Var kNuMITrueIsCC([](const caf::SRSliceProxy* slc) -> int {
    if ( slc->truth.index < 0 ) return -1;
    if (slc->truth.iscc) return 1;
    else return 0;
  });
  // Number of primary proton
  const Var kNuMITrueNProton([](const caf::SRSliceProxy* slc) -> int {
    if ( slc->truth.index < 0 ) return -1;
    return kTruth_NProton_Primary(&slc->truth);
  });
  // Number of primary neutron
  const Var kNuMITrueNNeutron([](const caf::SRSliceProxy* slc) -> int {
    if ( slc->truth.index < 0 ) return -1;
    return kTruth_NNeutron_Primary(&slc->truth);
  });
  // Number of primary pi+
  const Var kNuMITrueNpip([](const caf::SRSliceProxy* slc) -> int {
    if ( slc->truth.index < 0 ) return -1;
    return kTruth_Npip_Primary(&slc->truth);
  });
  // Number of ANY pi+ (not just primary)
  const Var kNuMITrueNpip_All([](const caf::SRSliceProxy* slc) -> int {
    if ( slc->truth.index < 0 ) return -1;
    return kTruth_Npip_All(&slc->truth);
  });
  // Number of primary pi-
  const Var kNuMITrueNpim([](const caf::SRSliceProxy* slc) -> int {
    if ( slc->truth.index < 0 ) return -1;
    return kTruth_Npim_Primary(&slc->truth);
  });
  // Number of ANY pi- (not just primary)
  const Var kNuMITrueNpim_All([](const caf::SRSliceProxy* slc) -> int {
    if ( slc->truth.index < 0 ) return -1;
    return kTruth_Npim_All(&slc->truth);
  });
  // Number of primary pi0
  const Var kNuMITrueNpi0([](const caf::SRSliceProxy* slc) -> int {
    if ( slc->truth.index < 0 ) return -1;
    return kTruth_Npi0_Primary(&slc->truth);
  });
  // Number of ANY pi0 (not just primary)
  const Var kNuMITrueNpi0_All([](const caf::SRSliceProxy* slc) -> int {
    if ( slc->truth.index < 0 ) return -1;
    return kTruth_Npi0_All(&slc->truth);
  });
  // Neutrino energy
  const Var kNuMITrueNuE([](const caf::SRSliceProxy* slc) -> double {
    if ( slc->truth.index < 0 ) return -5.; //TODO Define better dummy value
    return kTruth_NeutrinoE(&slc->truth);
  });
  // Q2
  const Var kNuMITrueQ2([](const caf::SRSliceProxy* slc) -> double {
    if ( slc->truth.index < 0 ) return -5.; //TODO Define better dummy value
    return kTruth_Q2(&slc->truth);
  });
  // q0; energy transfer
  const Var kNuMITrueq0([](const caf::SRSliceProxy* slc) -> double {
    if ( slc->truth.index < 0 ) return -5.; //TODO Define better dummy value
    return kTruth_q0(&slc->truth);
  });
  // q3; momentum transfer
  const Var kNuMITrueq3([](const caf::SRSliceProxy* slc) -> double {
    if ( slc->truth.index < 0 ) return -5.; //TODO Define better dummy value
    return kTruth_q3(&slc->truth);
  });
  // w; hadronic mass
  const Var kNuMITruew([](const caf::SRSliceProxy* slc) -> double {
    if ( slc->truth.index < 0 ) return -5.; //TODO Define better dummy value
    return kTruth_w(&slc->truth);
  });
  // 0: RHC, 1: FHC
  const Var kNuMIIsFHC([](const caf::SRSliceProxy* slc) -> int {
    if ( slc->truth.index < 0 ) return -5; //TODO Define better dummy value
    return kTruth_IsFHC(&slc->truth);
  });

  // Muon
  const Cut kNuMIHasTrueMuon([](const caf::SRSliceProxy* slc){
    if( kTruth_MuonIndex(&slc->truth)>=0 ) return true;
    else return false;
  });
  // True muon kinetic energy
  const Var kNuMITrueMuonKE([](const caf::SRSliceProxy* slc) -> double {
    if ( slc->truth.index < 0 ) return -5.; //TODO Define better dummy value
    return kTruth_MuonKE(&slc->truth);
  });
  // True muon cosine angle w.r.t. neutrino
  const Var kNuMITrueMuonNuCosineTheta([](const caf::SRSliceProxy* slc) -> double {
    if ( slc->truth.index < 0 ) return -5.; //TODO Define better dummy value
    return kTruth_MuonNuCosineTheta(&slc->truth);
  });
  // True muon contain?: 1: contained, 0: not contained (-1: muon not found)
  const Var kNuMITrueMuonContained([](const caf::SRSliceProxy* slc) -> int {
    if ( slc->truth.index < 0 ) return -1; //TODO Define better dummy value
    return kTruth_MuonContained(&slc->truth);
  });

  // Proton
  const Cut kNuMIHasTrueProton([](const caf::SRSliceProxy* slc){
    if( kTruth_ProtonIndex(&slc->truth)>=0 ) return true;
    else return false;
  });
  // True proton kinetic energy
  const Var kNuMITrueProtonKE([](const caf::SRSliceProxy* slc) -> double {
    if ( slc->truth.index < 0 ) return -5.; //TODO Define better dummy value
    return kTruth_ProtonKE(&slc->truth);
  });
  // True proton cosine angle w.r.t. neutrino
  const Var kNuMITrueProtonNuCosineTheta([](const caf::SRSliceProxy* slc) -> double {
    if ( slc->truth.index < 0 ) return -5.; //TODO Define better dummy value
    return kTruth_ProtonNuCosineTheta(&slc->truth);
  });

  // Chargedpion
  const Cut kNuMIHasTrueChargedPion([](const caf::SRSliceProxy* slc){
    if( kTruth_ChargedPionIndex(&slc->truth)>=0 ) return true;
    else return false;
  });
  // True Charged pion kinetic energy
  const Var kNuMITrueChargedPionKE([](const caf::SRSliceProxy* slc) -> double {
    if ( slc->truth.index < 0 ) return -5.; //TODO Define better dummy value
    return kTruth_ChargedPionKE(&slc->truth);
  });
  const Var kNuMITrueChargedPionEndProcess([](const caf::SRSliceProxy* slc) -> int {
    int truthChargedPionIndex = kTruth_ChargedPionIndex(&slc->truth);
    if( truthChargedPionIndex >= 0 ){
      const auto& prim_ChargedPion = slc->truth.prim.at(truthChargedPionIndex);
      return prim_ChargedPion.end_process;
    }
    else{
      return -5;
    }
  });
  const Var kNuMITrueChargedPionLength([](const caf::SRSliceProxy* slc) -> int {
    int truthChargedPionIndex = kTruth_ChargedPionIndex(&slc->truth);
    if( truthChargedPionIndex >= 0 ){
      const auto& prim_ChargedPion = slc->truth.prim.at(truthChargedPionIndex);
      return prim_ChargedPion.length;
    }
    else{
      return -5;
    }
  });
  const Var kNuMITrueChargedPionHasMichel([](const caf::SRSliceProxy* slc) -> int {
    int michelidx = kTruth_ChargedPionMichelIndex(&slc->truth);
    if(michelidx>=0) return 1;
    else return 0;
  });

  // 0: Muon candidate track exiting, 1: Muon candidate track contained (-1: no muon candidate)
  const Var kNuMIRecoMuonContained([](const caf::SRSliceProxy* slc) -> int {
    if( kNuMIMuonCandidateIdx(slc) >= 0 ){
      auto const& trk = slc->reco.pfp.at(kNuMIMuonCandidateIdx(slc)).trk;
      const bool Contained = isContainedVol(trk.end.x,trk.end.y,trk.end.z);
      if(Contained) return 1;
      else return 0;
    }
    else{
      return -1;
    }
  });

  // - Charged pion
  const MultiVar kNuMIChargedPionMatchedTrackIndices([](const caf::SRSliceProxy* slc) -> std::vector<double> {
    std::vector<double> rets;
    int truthChargedPionIndex = kTruth_ChargedPionIndex(&slc->truth);
    if( truthChargedPionIndex >= 0 ){
      const auto& prim_ChargedPion = slc->truth.prim.at(truthChargedPionIndex);
      for(unsigned int i_pfp=0; i_pfp<slc->reco.pfp.size(); i_pfp++){
        auto const& pfp = slc->reco.pfp.at(i_pfp);
        if(pfp.trk.truth.p.G4ID==prim_ChargedPion.G4ID){
          rets.push_back(i_pfp);
        }
      }
    }
    return rets;
  });
  const Var kNuMINChargedPionMatchedTracks([](const caf::SRSliceProxy* slc) -> int {
    return kNuMIChargedPionMatchedTrackIndices(slc).size();
  });

  const MultiVar kNuMIChargedPionMatchedTrackScores([](const caf::SRSliceProxy* slc) -> std::vector<double> {
    std::vector<double> matchedTrackIndices = kNuMIChargedPionMatchedTrackIndices(slc);
    std::vector<double> rets;
    for(const auto& idx: matchedTrackIndices){
      auto const& pfp = slc->reco.pfp.at(round(idx));
      rets.push_back( pfp.trackScore );
    }
    return rets;
  });
  const MultiVar kNuMIChargedPionMatchedTrackLengths([](const caf::SRSliceProxy* slc) -> std::vector<double> {
    std::vector<double> matchedTrackIndices = kNuMIChargedPionMatchedTrackIndices(slc);
    std::vector<double> rets;
    for(const auto& idx: matchedTrackIndices){
      auto const& pfp = slc->reco.pfp.at(round(idx));
      rets.push_back( pfp.trk.len );
    }
    return rets;
  });
  const MultiVar kNuMIChargedPionMatchedTrackChi2Muons([](const caf::SRSliceProxy* slc) -> std::vector<double> {
    std::vector<double> matchedTrackIndices = kNuMIChargedPionMatchedTrackIndices(slc);
    std::vector<double> rets;
    for(const auto& idx: matchedTrackIndices){
      auto const& pfp = slc->reco.pfp.at(round(idx));
      rets.push_back( pfp.trk.chi2pid[2].chi2_muon );
    }
    return rets;
  });
  const MultiVar kNuMIChargedPionMatchedTrackChi2Protons([](const caf::SRSliceProxy* slc) -> std::vector<double> {
    std::vector<double> matchedTrackIndices = kNuMIChargedPionMatchedTrackIndices(slc);
    std::vector<double> rets;
    for(const auto& idx: matchedTrackIndices){
      auto const& pfp = slc->reco.pfp.at(round(idx));
      rets.push_back( pfp.trk.chi2pid[2].chi2_proton );
    }
    return rets;
  });
  const Var kNuMINChargedPionBestMatchedTrackByHitCompletenessIdx([](const caf::SRSliceProxy* slc) -> int {
    int idx = -1;
    double max_match_val = -1.;
    int truthChargedPionIndex = kTruth_ChargedPionIndex(&slc->truth);
    if( truthChargedPionIndex >= 0 ){
      const auto& prim_ChargedPion = slc->truth.prim.at(truthChargedPionIndex);
      for(unsigned int i_pfp=0; i_pfp<slc->reco.pfp.size(); i_pfp++){
        auto const& pfp = slc->reco.pfp.at(i_pfp);
        auto const& trk = pfp.trk;
        double this_pfp_match_val = -1.;
        for(const auto& match: trk.truth.matches){
          if(match.G4ID==prim_ChargedPion.G4ID){
            this_pfp_match_val = match.hit_completeness;
            break;
          }
        }

        if(this_pfp_match_val > max_match_val){
          max_match_val = this_pfp_match_val;
          idx = i_pfp;
        }

      }
    }
    return idx;
  });
  const Var kNuMINChargedPionBestMatchedTrackByHitPurityIdx([](const caf::SRSliceProxy* slc) -> int {
    int idx = -1;
    double max_match_val = -1.;
    int truthChargedPionIndex = kTruth_ChargedPionIndex(&slc->truth);
    if( truthChargedPionIndex >= 0 ){
      const auto& prim_ChargedPion = slc->truth.prim.at(truthChargedPionIndex);
      for(unsigned int i_pfp=0; i_pfp<slc->reco.pfp.size(); i_pfp++){
        auto const& pfp = slc->reco.pfp.at(i_pfp);
        auto const& trk = pfp.trk;
        double this_pfp_match_val = -1.;
        for(const auto& match: trk.truth.matches){
          if(match.G4ID==prim_ChargedPion.G4ID){
            this_pfp_match_val = match.hit_purity;
            break;
          }
        }

        if(this_pfp_match_val > max_match_val){
          max_match_val = this_pfp_match_val;
          idx = i_pfp;
        }

      }
    }
    return idx;
  });
  const Var kNuMINChargedPionBestMatchedTrackHitCompleteness([](const caf::SRSliceProxy* slc) -> double {
    double max_match_val = -1.;
    int truthChargedPionIndex = kTruth_ChargedPionIndex(&slc->truth);
    if( truthChargedPionIndex >= 0 ){
      const auto& prim_ChargedPion = slc->truth.prim.at(truthChargedPionIndex);
      for(unsigned int i_pfp=0; i_pfp<slc->reco.pfp.size(); i_pfp++){
        auto const& pfp = slc->reco.pfp.at(i_pfp);
        auto const& trk = pfp.trk;
        double this_pfp_match_val = -1.;
        for(const auto& match: trk.truth.matches){
          if(match.G4ID==prim_ChargedPion.G4ID){
            this_pfp_match_val = match.hit_completeness;
            break;
          }
        }

        if(this_pfp_match_val > max_match_val){
          max_match_val = this_pfp_match_val;
        }

      }
    }
    return max_match_val;
  });
  const Var kNuMINChargedPionBestMatchedTrackHitPurity([](const caf::SRSliceProxy* slc) -> double {
    double max_match_val = -1.;
    int truthChargedPionIndex = kTruth_ChargedPionIndex(&slc->truth);
    if( truthChargedPionIndex >= 0 ){
      const auto& prim_ChargedPion = slc->truth.prim.at(truthChargedPionIndex);
      for(unsigned int i_pfp=0; i_pfp<slc->reco.pfp.size(); i_pfp++){
        auto const& pfp = slc->reco.pfp.at(i_pfp);
        auto const& trk = pfp.trk;
        double this_pfp_match_val = -1.;
        for(const auto& match: trk.truth.matches){
          if(match.G4ID==prim_ChargedPion.G4ID){
            this_pfp_match_val = match.hit_purity;
            break;
          }
        }

        if(this_pfp_match_val > max_match_val){
          max_match_val = this_pfp_match_val;
        }

      }
    }
    return max_match_val;
  });
  const Var kNuMINChargedPionBestMatchedTrackByHitPurityChi2Muon([](const caf::SRSliceProxy* slc) -> double {
    int trackIdx = kNuMINChargedPionBestMatchedTrackByHitPurityIdx(slc);
    if(trackIdx>=0){
      auto const& pfp = slc->reco.pfp.at(trackIdx);
      const auto& trk = pfp.trk;
      return trk.chi2pid[2].chi2_muon;
    }
    else{
      return -5.;
    }
  });
  const Var kNuMINChargedPionBestMatchedTrackByHitPurityChi2Proton([](const caf::SRSliceProxy* slc) -> double {
    int trackIdx = kNuMINChargedPionBestMatchedTrackByHitPurityIdx(slc);
    if(trackIdx>=0){
      auto const& pfp = slc->reco.pfp.at(trackIdx);
      const auto& trk = pfp.trk;
      return trk.chi2pid[2].chi2_proton;
    }
    else{
      return -5.;
    }
  });
  const Var kNuMINChargedPionBestMatchedTrackByHitPurityChi2MIP([](const caf::SRSliceProxy* slc) -> double {
    int trackIdx = kNuMINChargedPionBestMatchedTrackByHitPurityIdx(slc);
    if(trackIdx>=0){
      return GetChi2MIP(slc->reco.pfp[trackIdx].trk.calo[2]);
    }
    else{
      return -5.;
    }
  });

  //   - Shower Var
  const MultiVar kNuMIChargedPionMatchedShowerIndices([](const caf::SRSliceProxy* slc) -> std::vector<double> {
    std::vector<double> rets;
    int truthChargedPionIndex = kTruth_ChargedPionIndex(&slc->truth);
    if( truthChargedPionIndex >= 0 ){
      const auto& prim_ChargedPion = slc->truth.prim.at(truthChargedPionIndex);
      for(unsigned int i_pfp=0; i_pfp<slc->reco.pfp.size(); i_pfp++){
        auto const& pfp = slc->reco.pfp.at(i_pfp);
        if(pfp.shw.truth.p.G4ID==prim_ChargedPion.G4ID){
          rets.push_back(i_pfp);
        }
      }
    }
    return rets;
  });
  const Var kNuMINChargedPionMatchedShowers([](const caf::SRSliceProxy* slc) -> int {
    return kNuMIChargedPionMatchedShowerIndices(slc).size();
  });
  const MultiVar kNuMIChargedPionMatchedShowerScores([](const caf::SRSliceProxy* slc) -> std::vector<double> {
    std::vector<double> matchedShowerIndices = kNuMIChargedPionMatchedShowerIndices(slc);
    std::vector<double> rets;
    for(const auto& idx: matchedShowerIndices){
      auto const& pfp = slc->reco.pfp.at(round(idx));
      rets.push_back( pfp.trackScore );
    }
    return rets;
  });
  const MultiVar kNuMIChargedPionMatchedShowerGaps([](const caf::SRSliceProxy* slc) -> std::vector<double> {
    std::vector<double> matchedShowerIndices = kNuMIChargedPionMatchedShowerIndices(slc);
    std::vector<double> rets;
    for(const auto& idx: matchedShowerIndices){
      auto const& pfp = slc->reco.pfp.at(round(idx));
      rets.push_back( pfp.shw.conversion_gap );
    }
    return rets;
  });
  const MultiVar kNuMIChargedPionMatchedShowerIsPrimaries([](const caf::SRSliceProxy* slc) -> std::vector<double> {
    std::vector<double> matchedShowerIndices = kNuMIChargedPionMatchedShowerIndices(slc);
    std::vector<double> rets;
    for(const auto& idx: matchedShowerIndices){
      auto const& pfp = slc->reco.pfp.at(round(idx));
      if( pfp.parent_is_primary ) rets.push_back( 1 );
      else rets.push_back( 0 );
    }
    return rets;
  });
  const MultiVar kNuMIChargedPionMatchedShowerEnergies([](const caf::SRSliceProxy* slc) -> std::vector<double> {
    std::vector<double> matchedShowerIndices = kNuMIChargedPionMatchedShowerIndices(slc);
    std::vector<double> rets;
    for(const auto& idx: matchedShowerIndices){
      auto const& pfp = slc->reco.pfp.at(round(idx));
      rets.push_back( pfp.shw.plane[2].energy );
    }
    return rets;
  });
  //   - Cut
  const Cut kNuMIOtherCCWithChargedPion([](const caf::SRSliceProxy* slc) {
    if(!kNuMI_1muNp0piStudy_OtherNuCC(slc)) return false;
    int n_chargedpion = kNuMITrueNpip(slc)+kNuMITrueNpim(slc);
    if( n_chargedpion== 0 ) return false;
    return true;
  });
  const Var kNuMITrueChargedPionMichelEnergy([](const caf::SRSliceProxy* slc) -> double {
    double ret = -5.;
    int truthChargedPionMichelIndex = kTruth_ChargedPionMichelIndex(&slc->truth);
    if( truthChargedPionMichelIndex >= 0 ){
      const auto& prim_ChargedPionMichel = slc->truth.prim.at(truthChargedPionMichelIndex);
      ret = prim_ChargedPionMichel.genE;
    }
    return ret;
  });
  //   - Michel from pion (kTruth_ChargedPionMichelIndex)
  const MultiVar kNuMIChargedPionMichelMatchedPfpIndices([](const caf::SRSliceProxy* slc) -> std::vector<double> {
    std::vector<double> rets;
    int truthChargedPionMichelIndex = kTruth_ChargedPionMichelIndex(&slc->truth);
    if( truthChargedPionMichelIndex >= 0 ){
      const auto& prim_ChargedPionMichel = slc->truth.prim.at(truthChargedPionMichelIndex);
      for(unsigned int i_pfp=0; i_pfp<slc->reco.pfp.size(); i_pfp++){
        auto const& pfp = slc->reco.pfp.at(i_pfp);
        if(pfp.trk.truth.p.G4ID==prim_ChargedPionMichel.G4ID){
          rets.push_back(i_pfp);
        }
      }
    }
    return rets;
  });
  const MultiVar kNuMIChargedPionMichelMatchedPfpScores([](const caf::SRSliceProxy* slc) -> std::vector<double> {
    std::vector<double> matchedPfpIndices = kNuMIChargedPionMichelMatchedPfpIndices(slc);
    std::vector<double> rets;
    for(const auto& idx: matchedPfpIndices){
      auto const& pfp = slc->reco.pfp.at(round(idx));
      rets.push_back( pfp.trackScore );
    }
    return rets;
  });
  const MultiVar kNuMIChargedPionMichelMatchedPfpTrackLengths([](const caf::SRSliceProxy* slc) -> std::vector<double> {
    std::vector<double> matchedPfpIndices = kNuMIChargedPionMichelMatchedPfpIndices(slc);
    std::vector<double> rets;
    for(const auto& idx: matchedPfpIndices){
      auto const& pfp = slc->reco.pfp.at(round(idx));
      rets.push_back( pfp.trk.len );
    }
    return rets;
  });
  const MultiVar kNuMIChargedPionMichelMatchedPfpTrackDistances([](const caf::SRSliceProxy* slc) -> std::vector<double> {
    std::vector<double> matchedPfpIndices = kNuMIChargedPionMichelMatchedPfpIndices(slc);
    std::vector<double> rets;
    for(const auto& idx: matchedPfpIndices){
      auto const& pfp = slc->reco.pfp.at(round(idx));
      const float Atslc = std::hypot(slc->vertex.x - pfp.trk.start.x,
                                    slc->vertex.y - pfp.trk.start.y,
                                    slc->vertex.z - pfp.trk.start.z);
      rets.push_back( Atslc );
    }
    return rets;
  });
  const MultiVar kNuMIChargedPionMichelMatchedPfpTrackIsPrimaries([](const caf::SRSliceProxy* slc) -> std::vector<double> {
    std::vector<double> matchedPfpIndices = kNuMIChargedPionMichelMatchedPfpIndices(slc);
    std::vector<double> rets;
    for(const auto& idx: matchedPfpIndices){
      auto const& pfp = slc->reco.pfp.at(round(idx));
      if(pfp.parent_is_primary) rets.push_back(1);
      else rets.push_back(0);
    }
    return rets;
  });
  const MultiVar kNuMIChargedPionMichelMatchedPfpTrackHitPurities([](const caf::SRSliceProxy* slc) -> std::vector<double> {
    std::vector<double> matchedPfpIndices = kNuMIChargedPionMichelMatchedPfpIndices(slc);
    std::vector<double> rets;
    for(const auto& idx: matchedPfpIndices){
      auto const& pfp = slc->reco.pfp.at(round(idx));

      auto const& bestmatch = pfp.trk.truth.p;
      auto const& matches = pfp.trk.truth.matches;
      double this_match_val = -1.;
      for(const auto& match: matches){
        if(match.G4ID == bestmatch.G4ID){
          this_match_val = match.hit_purity;
        }
      }
      rets.push_back(this_match_val);
    }
    return rets;
  });
  const MultiVar kNuMIChargedPionMichelMatchedPfpTrackHitCompletenesses([](const caf::SRSliceProxy* slc) -> std::vector<double> {
    std::vector<double> matchedPfpIndices = kNuMIChargedPionMichelMatchedPfpIndices(slc);
    std::vector<double> rets;
    for(const auto& idx: matchedPfpIndices){
      auto const& pfp = slc->reco.pfp.at(round(idx));

      auto const& bestmatch = pfp.trk.truth.p;
      auto const& matches = pfp.trk.truth.matches;
      double this_match_val = -1.;
      for(const auto& match: matches){
        if(match.G4ID == bestmatch.G4ID){
          this_match_val = match.hit_completeness;
        }
      }
      rets.push_back(this_match_val);
    }
    return rets;
  });
  const MultiVar kNuMIChargedPionMichelMatchedPfpShowerEnergies([](const caf::SRSliceProxy* slc) -> std::vector<double> {
    std::vector<double> matchedPfpIndices = kNuMIChargedPionMichelMatchedPfpIndices(slc);
    std::vector<double> rets;
    for(const auto& idx: matchedPfpIndices){
      auto const& pfp = slc->reco.pfp.at(round(idx));
      rets.push_back( pfp.shw.plane[2].energy );
    }
    return rets;
  });
  const MultiVar kNuMIChargedPionMichelMatchedPfpShowerLengths([](const caf::SRSliceProxy* slc) -> std::vector<double> {
    std::vector<double> matchedPfpIndices = kNuMIChargedPionMichelMatchedPfpIndices(slc);
    std::vector<double> rets;
    for(const auto& idx: matchedPfpIndices){
      auto const& pfp = slc->reco.pfp.at(round(idx));
      rets.push_back( pfp.shw.len );
    }
    return rets;
  });
  const MultiVar kNuMIChargedPionMichelMatchedPfpShowerGaps([](const caf::SRSliceProxy* slc) -> std::vector<double> {
    std::vector<double> matchedPfpIndices = kNuMIChargedPionMichelMatchedPfpIndices(slc);
    std::vector<double> rets;
    for(const auto& idx: matchedPfpIndices){
      auto const& pfp = slc->reco.pfp.at(round(idx));
      rets.push_back( pfp.shw.conversion_gap );
    }
    return rets;
  });
  const MultiVar kNuMIChargedPionMichelMatchedPfpShowerOpeningAngles([](const caf::SRSliceProxy* slc) -> std::vector<double> {
    std::vector<double> matchedPfpIndices = kNuMIChargedPionMichelMatchedPfpIndices(slc);
    std::vector<double> rets;
    for(const auto& idx: matchedPfpIndices){
      auto const& pfp = slc->reco.pfp.at(round(idx));
      rets.push_back( pfp.shw.open_angle );
    }
    return rets;
  });
  const Var kNuMIChargedPionMichelMatchedPfpShowerEnergySum([](const caf::SRSliceProxy* slc) -> double {
    std::vector<double> matchedPfpIndices = kNuMIChargedPionMichelMatchedPfpIndices(slc);
    if(matchedPfpIndices.size()==0) return -1.;

    double sumE = 0.;
    for(const auto& idx: matchedPfpIndices){
      auto const& pfp = slc->reco.pfp.at(round(idx));
      sumE += pfp.shw.plane[2].energy;
    }

    return sumE;

  });
  // - Proton
  const Var kNuMIRecoProtonMatchedToTrueProton([](const caf::SRSliceProxy* slc) -> int {
    int true_leading_proton_index = kTruth_ProtonIndex(&slc->truth);
    if(true_leading_proton_index<0) return -1;
    const auto& prim_leading_proton = slc->truth.prim.at(true_leading_proton_index);

    int reco_proton_index = kNuMIProtonCandidateIdx(slc);
    if(reco_proton_index<0) return -2;
    const auto& trk = slc->reco.pfp.at(reco_proton_index).trk;
    const auto& bestmatched = trk.truth.p;

    if( prim_leading_proton.G4ID == bestmatched.G4ID ) return 1;
    else return 0;
  });

}
