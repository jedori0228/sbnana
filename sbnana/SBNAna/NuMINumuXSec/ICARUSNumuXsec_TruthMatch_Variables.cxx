#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_TruthMatch_Variables.h"

using namespace std;
using namespace ana;
using namespace ana::PrimaryUtil;

namespace ICARUSNumuXsec{

namespace TruthMatch{

  // For a given true muon (truth_index), find a reco track whose best-matched is this particle
  const Var TruthMuonIndex([](const caf::SRSliceProxy* slc) -> double {
    return ana::PrimaryUtil::MuonIndex_True(slc->truth);
  });
  const Var TruthMuonLength([](const caf::SRSliceProxy* slc) -> double {
    int truth_idx = TruthMuonIndex(slc);
    if(truth_idx>=0){
      if(isnan(slc->truth.prim.at(truth_idx).length)) return -999.;
      else return slc->truth.prim.at(truth_idx).length;
    }
    else{
      return -999.;
    }
  });
  const Var TruthMuonMatchedTrackIndex([](const caf::SRSliceProxy* slc) -> double {
    int truth_idx = TruthMuonIndex(slc);
    return GetMatchedRecoTrackIndex(slc, truth_idx);
  });
  const Var TruthMuonMatchedTrackContainedness([](const caf::SRSliceProxy* slc) -> double {
    int muonTrackIndex = TruthMuonMatchedTrackIndex(slc);
    // -1 : No track found
    // 0 : Exiting
    // 1 : Contained
    if(muonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
      bool isContained = fv_track.isContained(trk.end.x, trk.end.y, trk.end.z);
      if(isContained) return 1;
      else return 0;
    }
    else{
      return -1;
    }
  });
  const Var TruthMuonMatchedTrackChi2Proton([](const caf::SRSliceProxy* slc) -> double {
    int muonTrackIndex = TruthMuonMatchedTrackIndex(slc);
    if(muonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
      int bp = trk.bestplane;
      float Chi2Proton = trk.chi2pid[bp].chi2_proton;
      if(isnan(Chi2Proton)) return -999.;
      else return Chi2Proton;
    }
    else{
      return 999999;
    }
  });

  const Var TruthMuonMatchedTrackChi2Muon([](const caf::SRSliceProxy* slc) -> double {
    int muonTrackIndex = TruthMuonMatchedTrackIndex(slc);
    if(muonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
      int bp = trk.bestplane;
      float Chi2Muon = trk.chi2pid[bp].chi2_muon;
      if(isnan(Chi2Muon)) return -999.;
      else return Chi2Muon;
    }
    else{
      return 999999;
    }
  });
  const Var TruthMuonMatchedTrackChi2Pion([](const caf::SRSliceProxy* slc) -> double {
    int muonTrackIndex = TruthMuonMatchedTrackIndex(slc);
    if(muonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
      int bp = trk.bestplane;
      float Chi2Pion = trk.chi2pid[bp].chi2_pion;
      if(isnan(Chi2Pion)) return -999.;
      else return Chi2Pion;
    }
    else{
      return 999999;
    }
  });
  const Var TruthMuonMatchedTrackCustomChi2InelasticPionCollection([](const caf::SRSliceProxy* slc) -> double {
    int muonTrackIndex = TruthMuonMatchedTrackIndex(slc);
    if(muonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
      if(trk.calo[2].nhit==0) return -999.;
      double new_chi2 = dedxtempt.CalculateInelasticPionChi2(trk.calo[2], 0);
      return new_chi2;
    }
    else{
      return -999.;
    }
  });
  // - Michel from muon
  const Var TruthMuonMichelIndex([](const caf::SRSliceProxy* slc) -> double {
    int truth_idx = TruthMuonIndex(slc);
    if(truth_idx>=0){
      const auto& mu_prim = slc->truth.prim[truth_idx];
      int truth_el_idx(-1);
      for(std::size_t i(0); i < slc->truth.prim.size(); ++i){
        if( abs(slc->truth.prim[i].pdg)==11 && (int)slc->truth.prim[i].parent==(int)mu_prim.G4ID ){
          truth_el_idx = i;
/*
          if(slc->truth.prim[i].start_process==41){
          std::cout << "[JSKIMDEBUG][TruthMuonMichelIndex] Michel electron from MuonCapture:" << std::endl;
          PrintPrimaries(slc);
          }
*/
        }
      }
      return truth_el_idx;
    }
    else{
      return -999.;
    }
  });
  const Var TruthMuonMichelStartProcess([](const caf::SRSliceProxy* slc) -> double {
    int truth_idx = TruthMuonMichelIndex(slc);
    if(truth_idx>=0){
      return slc->truth.prim[truth_idx].start_process;
    }
    else{
      return -999.;
    }
  });
  const Var TruthMuonMichelKE([](const caf::SRSliceProxy* slc) -> double {
    int truth_idx = TruthMuonMichelIndex(slc);
    if(truth_idx>=0){
      return slc->truth.prim[truth_idx].genE - M_ELECTRON;
    }
    else{
      return -999.;
    }
  });
  const Var TruthMuonMichelKEFromDecay([](const caf::SRSliceProxy* slc) -> double {
    int truth_idx = TruthMuonMichelIndex(slc);
    if(truth_idx>=0){
      if( slc->truth.prim[truth_idx].start_process == 3 ) return slc->truth.prim[truth_idx].genE - M_ELECTRON;
      else return -999.;
    }
    else{
      return -999.;
    }
  });
  const Var TruthMuonMichelKEFromCapture([](const caf::SRSliceProxy* slc) -> double {
    int truth_idx = TruthMuonMichelIndex(slc);
    if(truth_idx>=0){
      if( slc->truth.prim[truth_idx].start_process == 41 ) return slc->truth.prim[truth_idx].genE - M_ELECTRON;
      else return -999.;
    }
    else{
      return -999.;
    }
  });
  const Var TruthMuonMichelMatchedShowerIndex([](const caf::SRSliceProxy* slc) -> double {
    int truth_idx = TruthMuonMichelIndex(slc);
    return GetMatchedRecoShowerIndex(slc, truth_idx);
  });
  const Var TruthMuonMichelMatchedShowerKE([](const caf::SRSliceProxy* slc) -> double {
    int muonMichelShowerIndex = TruthMuonMichelMatchedShowerIndex(slc);
    if(muonMichelShowerIndex>=0){
      const auto& shw = slc->reco.pfp.at(muonMichelShowerIndex).shw;
      return shw.plane[2].energy;
    }
    else{
      return -1;
    }
  });
  const Var TruthMuonMichelMatchedShowerDistanceFromMuonEnd([](const caf::SRSliceProxy* slc) -> double {
    int muonTrackIndex = TruthMuonMatchedTrackIndex(slc);
    int muonMichelShowerIndex = TruthMuonMichelMatchedShowerIndex(slc);
    if(muonTrackIndex>=0 && muonMichelShowerIndex>=0){
      const auto& trk_mu = slc->reco.pfp.at(muonTrackIndex).trk;
      const auto& shw_el = slc->reco.pfp.at(muonMichelShowerIndex).shw;

      TVector3 vec_trk_mu_end(trk_mu.end.x, trk_mu.end.y, trk_mu.end.z);
      TVector3 vec_shw_el_start(shw_el.start.x, shw_el.start.y, shw_el.start.z);

      TVector3 vec_mu_to_el = vec_shw_el_start-vec_trk_mu_end;
      return vec_mu_to_el.Mag();

    }
    else{
      return -1;
    }
  });

  // For a given true proton (truth_index), find a reco track whose best-matched is this particle
  const Var TruthProtonIndex([](const caf::SRSliceProxy* slc) -> double {
    return ana::PrimaryUtil::ProtonIndex_True(slc->truth);
  });
  const Var TruthProtonLength([](const caf::SRSliceProxy* slc) -> double {
    int truth_idx = TruthProtonIndex(slc);
    if(truth_idx>=0){
      if(isnan(slc->truth.prim.at(truth_idx).length)) return -999.;
      else return slc->truth.prim.at(truth_idx).length;
    }
    else{
      return -999.;
    }
  });
  const Var TruthProtonMatchedTrackIndex([](const caf::SRSliceProxy* slc) -> double {
    int truth_idx = TruthProtonIndex(slc);
    return GetMatchedRecoTrackIndex(slc, truth_idx);
  });
  const Var TruthProtonMatchedTrackContainedness([](const caf::SRSliceProxy* slc) -> double {
    int protonTrackIndex = TruthProtonMatchedTrackIndex(slc);
    // -1 : No track found
    // 0 : Exiting
    // 1 : Contained
    if(protonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
      bool isContained = fv_track.isContained(trk.end.x, trk.end.y, trk.end.z);
      if(isContained) return 1;
      else return 0;
    }
    else{
      return -1;
    }
  });
  const Var TruthProtonMatchedTrackChi2Proton([](const caf::SRSliceProxy* slc) -> double {
    int protonTrackIndex = TruthProtonMatchedTrackIndex(slc);
    if(protonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
      int bp = trk.bestplane;
      float Chi2Proton = trk.chi2pid[bp].chi2_proton;
      if(isnan(Chi2Proton)) return -999.;
      else return Chi2Proton;
    }
    else{
      return 999999;
    }
  });

  const Var TruthProtonMatchedTrackChi2Muon([](const caf::SRSliceProxy* slc) -> double {
    int protonTrackIndex = TruthProtonMatchedTrackIndex(slc);
    if(protonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
      int bp = trk.bestplane;
      float Chi2Muon = trk.chi2pid[bp].chi2_muon;
      if(isnan(Chi2Muon)) return -999.;
      else return Chi2Muon;
    }
    else{
      return 999999;
    }
  });
  const Var TruthProtonMatchedTrackChi2Pion([](const caf::SRSliceProxy* slc) -> double {
    int protonTrackIndex = TruthProtonMatchedTrackIndex(slc);
    if(protonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
      int bp = trk.bestplane;
      float Chi2Pion = trk.chi2pid[bp].chi2_pion;
      if(isnan(Chi2Pion)) return -999.;
      else return Chi2Pion;
    }
    else{
      return 999999;
    }
  });
  const Var TruthProtonMatchedTrackCustomChi2InelasticPionCollection([](const caf::SRSliceProxy* slc) -> double {
    int protonTrackIndex = TruthProtonMatchedTrackIndex(slc);
    if(protonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
      if(trk.calo[2].nhit==0) return -999.;
      double new_chi2 = dedxtempt.CalculateInelasticPionChi2(trk.calo[2], 0);
      return new_chi2;
    }
    else{ 
      return -999.;
    }
  });

  // For a given true charged pion (truth_index), find a reco track whose best-matched is this particle
  const Var TruthChargedPionIndex([](const caf::SRSliceProxy* slc) -> double {
    double max_E(-999);
    int truth_idx(-1);
    for(std::size_t i(0); i < slc->truth.prim.size(); ++i){
      if( abs(slc->truth.prim.at(i).pdg)==211 && slc->truth.prim.at(i).start_process==0 ){
        if(isnan(slc->truth.prim.at(i).genE)) continue;
        double this_E = slc->truth.prim.at(i).genE;
        if(this_E>max_E){
          max_E = this_E;
          truth_idx = i;
        }
      }
    }
    return truth_idx;
  });
  const Var TruthChargedPionExistence([](const caf::SRSliceProxy* slc) -> double {
    int truth_idx = TruthChargedPionIndex(slc);
    if(truth_idx>=0) return 1.;
    else return 0.;
  });
  const Var TruthChargedPionLength([](const caf::SRSliceProxy* slc) -> double {
    int truth_idx = TruthChargedPionIndex(slc);
    if(truth_idx>=0){
      if(isnan(slc->truth.prim.at(truth_idx).length)) return -999.;
      else return slc->truth.prim.at(truth_idx).length;
    }
    else{
      return -999.;
    }
  });
  const Var TruthChargedPionKE([](const caf::SRSliceProxy* slc) -> double {
    int truth_idx = TruthChargedPionIndex(slc);
    if(truth_idx>=0){
      return slc->truth.prim.at(truth_idx).genE - M_CHARGEDPION;
    }
    else{
      return -999.;
    }
  });
  const Var TruthChargedPionContainedness([](const caf::SRSliceProxy* slc) -> double {
    int truth_idx = TruthChargedPionIndex(slc);
    if(truth_idx>=0){
      return slc->truth.prim.at(truth_idx).contained;
    }
    else{
      return -999.;
    }
  });
  const Var TruthChargedPionMatchedTrackIndex([](const caf::SRSliceProxy* slc) -> double {
    int truth_idx = TruthChargedPionIndex(slc);
    if(truth_idx<0) return -999.;
    return GetMatchedRecoTrackIndex(slc, truth_idx, 0.);
  });
  const Var TruthChargedPionMatchedTrackScore([](const caf::SRSliceProxy* slc) -> double {
    int cpionTrackIndex = TruthChargedPionMatchedTrackIndex(slc);
    if(cpionTrackIndex>=0){
      return slc->reco.pfp[cpionTrackIndex].trackScore;
    }
    else{
      return -999.;
    }
  });
  const Var TruthChargedPionMatchedTrackContainedness([](const caf::SRSliceProxy* slc) -> double {
    int cpionTrackIndex = TruthChargedPionMatchedTrackIndex(slc);
    // -1 : No track found
    // 0 : Exiting
    // 1 : Contained
    if(cpionTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(cpionTrackIndex).trk;
      bool isContained = fv_track.isContained(trk.end.x, trk.end.y, trk.end.z);
      if(isContained) return 1;
      else return 0;
    }
    else{
      return -1;
    }
  });
  const Var TruthChargedPionMatchedTrackEndProcess([](const caf::SRSliceProxy* slc) -> double {
    int cpionTrackIndex = TruthChargedPionMatchedTrackIndex(slc);
    if(cpionTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(cpionTrackIndex).trk;
      return trk.truth.p.end_process;
    }
    else{
      return -999.;
    }
  });
  const Var TruthChargedPionMatchedTrackCustomChi2MuonCollection([](const caf::SRSliceProxy* slc) -> double {
    int cpionTrackIndex = TruthChargedPionMatchedTrackIndex(slc);
    if(cpionTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(cpionTrackIndex).trk;
      if(trk.calo[2].nhit==0) return -999.;
      double new_chi2 = dedxtempt.CalculateChi2(trk.calo[2], 3);
      return new_chi2;
    }
    else{
      return -999.;
    }
  });
  const Var TruthChargedPionMatchedTrackCustomChi2ProtonCollection([](const caf::SRSliceProxy* slc) -> double {
    int cpionTrackIndex = TruthChargedPionMatchedTrackIndex(slc);
    if(cpionTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(cpionTrackIndex).trk;
      if(trk.calo[2].nhit==0) return -999.;
      double new_chi2 = dedxtempt.CalculateChi2(trk.calo[2], 0);
      return new_chi2;
    }
    else{
      return -999.;
    }
  });
  const Var TruthChargedPionMatchedTrackCustomChi2InelasticPionCollection([](const caf::SRSliceProxy* slc) -> double {
    int cpionTrackIndex = TruthChargedPionMatchedTrackIndex(slc);
    if(cpionTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(cpionTrackIndex).trk;
      if(trk.calo[2].nhit==0) return -999.;
      double new_chi2 = dedxtempt.CalculateInelasticPionChi2(trk.calo[2], 0);
      return new_chi2;
    } 
    else{
      return -999.;
    } 
  });
  const Var TruthChargedPionMatchedShowerIndex([](const caf::SRSliceProxy* slc) -> double {
    int truth_idx = TruthChargedPionIndex(slc);
    if(truth_idx<0) return -999.;
    return GetMatchedRecoShowerIndex(slc, truth_idx);
  });
  const Var TruthChargedPionHasReco([](const caf::SRSliceProxy* slc) -> double {
    bool TruePionFound = TruthChargedPionIndex(slc)>=0;
    if(!TruePionFound) return -1;
    // -1    : true pion not found
    // 00 = 0: no match
    // 01 = 1: track match, no shower
    // 10 = 2: no track, shower match
    // 11 = 3: all match
    //      4: should not happen
    bool trkMatch = TruthChargedPionMatchedTrackIndex(slc)>=0;
    bool shwMatch = TruthChargedPionMatchedShowerIndex(slc)>=0;
    if(!trkMatch && !shwMatch) return 0.;
    if(trkMatch && !shwMatch) return 1.;
    if(!trkMatch && shwMatch) return 2.;
    if(trkMatch && shwMatch) return 3.;
    return 4.;

  });


  const Var TruthChargedPionNDaughters([](const caf::SRSliceProxy* slc) -> double {
    int cpionTrackIndex = TruthChargedPionMatchedTrackIndex(slc);
    if(cpionTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(cpionTrackIndex).trk;
      return trk.truth.p.daughters.size();
    }
    else{
      return -999.;
    }
  });
  const MultiVar TruthChargedPionDaughterPDGs([](const caf::SRSliceProxy* slc) -> vector<double> {
    vector<double> rets;
    int cpionTrackIndex = TruthChargedPionMatchedTrackIndex(slc);
    if(cpionTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(cpionTrackIndex).trk;
      for(const auto& d: trk.truth.p.daughters){
        for(const auto& prim: slc->truth.prim){
          if(d == (unsigned int)prim.G4ID){
            rets.push_back( prim.pdg );
          }
        }
      }
    }
    return rets;
  });
  // - Michel from pion
  const Var TruthChargedPionMichelIndex([](const caf::SRSliceProxy* slc) -> double {

    int truth_idx = TruthChargedPionIndex(slc);
    if(truth_idx>=0){
      const auto& pi_prim = slc->truth.prim[truth_idx];

      int truth_el_idx(-1);
      for(std::size_t i(0); i < slc->truth.prim.size(); ++i){
        if( abs(slc->truth.prim[i].pdg)==11 && (int)slc->truth.prim[i].parent==(int)pi_prim.G4ID ){
          truth_el_idx = i;
          //std::cout << "[PionMichel] Electron from pion found" << std::endl;
        }
      }

      int truth_mu_idx(-1);
      for(std::size_t i(0); i < slc->truth.prim.size(); ++i){
        if( abs(slc->truth.prim[i].pdg)==13 && (int)slc->truth.prim[i].parent==(int)pi_prim.G4ID ){
          truth_mu_idx = i;
          //std::cout << "[PionMichel] Muon from pion found" << std::endl;
        }
      }

      if(truth_mu_idx>=0){

        const auto& mu_prim = slc->truth.prim[truth_mu_idx];
        for(std::size_t i(0); i < slc->truth.prim.size(); ++i){
          if( abs(slc->truth.prim[i].pdg)==11 && (int)slc->truth.prim[i].parent==(int)mu_prim.G4ID ){
            truth_el_idx = i;
            //std::cout << "[PionMichel] Electron from muon which was from pion found" << std::endl;
          }
        }

     }

      return truth_el_idx;
    }
    else{
      return -999.;
    }
  });
  const Var TruthChargedPionMichelExistence([](const caf::SRSliceProxy* slc) -> double {
    int truth_pi_idx = TruthChargedPionIndex(slc);
    int truth_el_idx = TruthChargedPionMichelIndex(slc);
    if(!truth_pi_idx) return -1.;
    else{
      if(truth_el_idx) return 1.;
      else return 0.;
    }
  });
  const Var TruthChargedPionMichelStartProcess([](const caf::SRSliceProxy* slc) -> double {
    int truth_idx = TruthChargedPionMichelIndex(slc);
    if(truth_idx>=0){
      return slc->truth.prim[truth_idx].start_process;
    }
    else{
      return -999.;
    }
  });
  const Var TruthChargedPionMichelKE([](const caf::SRSliceProxy* slc) -> double {
    int truth_idx = TruthChargedPionMichelIndex(slc);
    if(truth_idx>=0){
      return slc->truth.prim[truth_idx].genE - M_ELECTRON;
    }
    else{
      return -999.;
    }
  });
  const Var TruthChargedPionMichelMatchedTrackIndex([](const caf::SRSliceProxy* slc) -> double {
    int truth_idx = TruthChargedPionMichelIndex(slc);
    return GetMatchedRecoTrackIndex(slc, truth_idx, 0.);
  });
  const Var TruthChargedPionMichelMatchedShowerIndex([](const caf::SRSliceProxy* slc) -> double {
    int truth_idx = TruthChargedPionMichelIndex(slc);
    return GetMatchedRecoShowerIndex(slc, truth_idx, 1.0);
  });
  const Var TruthChargedPionMichelMatchedShowerKE([](const caf::SRSliceProxy* slc) -> double {
    int pionMichelShowerIndex = TruthChargedPionMichelMatchedShowerIndex(slc);
    if(pionMichelShowerIndex>=0){
      const auto& shw = slc->reco.pfp.at(pionMichelShowerIndex).shw;
      return shw.plane[2].energy;
    }
    else{
      return -1;
    }
  });
  const Var TruthChargedPionMichelMatchedShowerDistanceFromMuonEnd([](const caf::SRSliceProxy* slc) -> double {
    int pionTrackIndex = TruthChargedPionMatchedTrackIndex(slc);
    int pionMichelShowerIndex = TruthChargedPionMichelMatchedShowerIndex(slc);
    if(pionTrackIndex>=0 && pionMichelShowerIndex>=0){
      const auto& trk_pi = slc->reco.pfp.at(pionTrackIndex).trk;
      const auto& shw_el = slc->reco.pfp.at(pionMichelShowerIndex).shw;

      TVector3 vec_trk_pi_end(trk_pi.end.x, trk_pi.end.y, trk_pi.end.z);
      TVector3 vec_shw_el_start(shw_el.start.x, shw_el.start.y, shw_el.start.z);

      TVector3 vec_pi_to_el = vec_shw_el_start-vec_trk_pi_end;
      return vec_pi_to_el.Mag();

    }
    else{
      return -1;
    }
  });
  const Var TruthChargedPionMichelHasReco([](const caf::SRSliceProxy* slc) -> double {
    bool TrueMichelFound = TruthChargedPionMichelIndex(slc)>=0;
    if(!TrueMichelFound) return -1;
    // -1    : true michel (from pion) not found
    // 00 = 0: no match
    // 01 = 1: track match, no shower
    // 10 = 2: no track, shower match
    // 11 = 3: all match
    //      4: should not happen
    bool trkMatch = TruthChargedPionMichelMatchedTrackIndex(slc)>=0;
    bool shwMatch = TruthChargedPionMichelMatchedShowerIndex(slc)>=0;
    if(!trkMatch && !shwMatch) return 0.;
    if(trkMatch && !shwMatch) return 1.;
    if(!trkMatch && shwMatch) return 2.;
    if(trkMatch && shwMatch) return 3.;
    return 4.;

  });

  // - Neutral pion
  // For a given true neutral pion (truth_index), find a reco track whose best-matched is this particle
  const Var TruthNeutralPionIndex([](const caf::SRSliceProxy* slc) -> double {
    double max_E(-999);
    int truth_idx(-1);
    for(std::size_t i(0); i < slc->truth.prim.size(); ++i){
      if( abs(slc->truth.prim.at(i).pdg)==111 && slc->truth.prim.at(i).start_process==0 ){
        if(isnan(slc->truth.prim.at(i).genE)) continue;
        double this_E = slc->truth.prim.at(i).genE;
        if(this_E>max_E){
          max_E = this_E;
          truth_idx = i;
        }
      }
    }
    return truth_idx;
  });
  const Var TruthNeutralPionKE([](const caf::SRSliceProxy* slc) -> double {
    int truth_idx = TruthNeutralPionIndex(slc);
    if(truth_idx>=0){
      return slc->truth.prim.at(truth_idx).genE - M_NEUTRALPION;
    }
    else{
      return -999.;
    }

  });
  const MultiVar TruthNeutralPionMatchedShowerIndicies([](const caf::SRSliceProxy* slc) -> vector<double> {
    int truth_idx = TruthNeutralPionIndex(slc);

    vector<double> rets;

    if(truth_idx>=0){

      //return GetMatchedRecoShowerIndices(slc, truth_idx, 1.0);
      const auto& prim = slc->truth.prim.at(truth_idx);

      for(std::size_t i(0); i < slc->reco.pfp.size(); ++i){
        const auto& pfp = slc->reco.pfp.at(i);
        const auto& shw = pfp.shw;

        const auto& shw_pdg = shw.truth.p.pdg;
        // photon
        if(shw_pdg!=22) continue;

        const auto& shw_ParentG4ID = shw.truth.p.parent;

        if((int)shw_ParentG4ID==prim.G4ID){
          rets.push_back(i);
        }
      }

    }

    return rets;

  });
  const MultiVar TruthNeutralPionMatchedShowerTrackScores([](const caf::SRSliceProxy* slc) -> vector<double> {
    vector<double> shw_indices = TruthNeutralPionMatchedShowerIndicies(slc);
    vector<double> rets;
    for(const auto& shw_index: shw_indices){
      rets.push_back( slc->reco.pfp[shw_index].trackScore );
    }
    return rets;
  });
  const MultiVar TruthNeutralPionMatchedShowerEnergies([](const caf::SRSliceProxy* slc) -> vector<double> {
    vector<double> shw_indices = TruthNeutralPionMatchedShowerIndicies(slc);
    vector<double> rets;
    for(const auto& shw_index: shw_indices){
      rets.push_back( slc->reco.pfp[shw_index].shw.plane[2].energy );
    }
    return rets;
  });
  const MultiVar TruthNeutralPionMatchedShowerdEdxs([](const caf::SRSliceProxy* slc) -> vector<double> {
    vector<double> shw_indices = TruthNeutralPionMatchedShowerIndicies(slc);
    vector<double> rets;
    for(const auto& shw_index: shw_indices){
      rets.push_back( slc->reco.pfp[shw_index].shw.plane[2].dEdx );
    }
    return rets;
  });
  const MultiVar TruthNeutralPionMatchedShowerConvGaps([](const caf::SRSliceProxy* slc) -> vector<double> {
    vector<double> shw_indices = TruthNeutralPionMatchedShowerIndicies(slc);
    vector<double> rets;
    for(const auto& shw_index: shw_indices){
      rets.push_back( slc->reco.pfp[shw_index].shw.conversion_gap );
    }
    return rets;
  });
  const MultiVar TruthNeutralPionMatchedShowerLengths([](const caf::SRSliceProxy* slc) -> vector<double> {
    vector<double> shw_indices = TruthNeutralPionMatchedShowerIndicies(slc);
    vector<double> rets;
    for(const auto& shw_index: shw_indices){
      rets.push_back( slc->reco.pfp[shw_index].shw.len );
    }
    return rets;
  });
  const Var TruthNeutralPionNMatchedShower([](const caf::SRSliceProxy* slc) -> double {
    return TruthNeutralPionMatchedShowerIndicies(slc).size();
  });
  const Var TruthNeutralPionNMatchedShowerWithCut([](const caf::SRSliceProxy* slc) -> double {
    vector<double> shw_indices = TruthNeutralPionMatchedShowerIndicies(slc);
    int Nshw = 0;
    for(const auto& shw_index: shw_indices){
      const auto& pfp = slc->reco.pfp.at(shw_index);
      const auto& shw = pfp.shw;

      if(pfp.trackScore>0.5) continue;
      if(shw.plane[2].dEdx<3) continue;

      Nshw++;
    }
    return Nshw;
  });

  const Var TruthNeutralPionMatchedShowerSumEnergy([](const caf::SRSliceProxy* slc) -> double {
    vector<double> shw_indices = TruthNeutralPionMatchedShowerIndicies(slc);
    double esum = 0.;
    for(const auto& shw_index: shw_indices){
      const auto& shw = slc->reco.pfp.at(shw_index).shw;
      esum += shw.plane[2].energy;
    }

    return esum;

  });
  const Var TruthNeutralPionMatchedShowerSumInvariantMass([](const caf::SRSliceProxy* slc) -> double {
    vector<double> shw_indices = TruthNeutralPionMatchedShowerIndicies(slc);
    double esum = 0.;
    TVector3 psum(0., 0., 0.);
    for(const auto& shw_index: shw_indices){
      const auto& shw = slc->reco.pfp.at(shw_index).shw;

      TVector3 this_p(shw.dir.x.GetValue(), shw.dir.y.GetValue(), shw.dir.z.GetValue());
      this_p *= shw.plane[2].energy;

      esum += shw.plane[2].energy;
      psum += this_p;
    }

    double m2 = esum*esum - psum.Mag2();

    return m2>0 ? sqrt(m2) : 0.;

  });


  // test
  const SpillMultiVar TruthChargedPionMichelMatchedSlice([](const caf::SRSpillProxy *sr) -> vector<double> {

    vector<double> rets;

    for(std::size_t i(0); i < sr->slc.size(); ++i){
      const auto& slc = sr->slc.at(i);
      int pi_idx = TruthChargedPionIndex(&slc);
      int el_idx = TruthChargedPionMichelIndex(&slc);
      if(pi_idx<0 || el_idx<0){
        continue;
      }

      const auto& prim_el = slc.truth.prim[el_idx];

      bool InSameSlice = false;
      bool InOtherSlice = false;

      for(std::size_t j(0); j < sr->slc.size(); ++j){
        const auto& slc_2 = sr->slc.at(j);
        for(const auto& pfp: slc_2.reco.pfp){
          const auto& trk = pfp.trk;
          const auto& shw = pfp.shw;
          bool trkMatched = (trk.truth.p.G4ID == prim_el.G4ID);
          bool shwMatched = (shw.truth.p.G4ID == prim_el.G4ID);
          if(trkMatched||shwMatched){
            if(i==j) InSameSlice = true;
            if(i!=j) InOtherSlice = true;
          }
        }
      }

      int ret = -1;
      if(!InSameSlice && !InOtherSlice){
        ret = 0;
        //std::cout << "[JSKIMDEBUG] Michel electron exists, but not matched to any pfps in the spill" << std::endl;
      }
      else if(InSameSlice && !InOtherSlice){
        ret = 1;
        //std::cout << "[JSKIMDEBUG] Michel electron exists and found in the same slice of pion" << std::endl;
      }
      else if(!InSameSlice && InOtherSlice){
        ret = 2;
        //std::cout << "[JSKIMDEBUG] Michel electron exists but found in other slice to pion" << std::endl;
      }
      else{
        ret = 3;
        //std::cout << "[JSKIMDEBUG] Michel electron exists and found in multiple slices" << std::endl;
      }
      rets.push_back(ret);

    }

    return rets;

  });


} // end namespace TruthMatch

} // end namespace ICARUSNumuXsec
