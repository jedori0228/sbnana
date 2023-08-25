#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_TwoTrack_Cuts.h"

using namespace ana;
using namespace std;

namespace ICARUSNumuXsec{

namespace TwoTrack{

  // Test

  const SpillMultiVar TestVar([](const caf::SRSpillProxy *sr) -> vector<double> {

    std::vector<double> rets;
/*
    for(std::size_t i(0); i < sr->slc.size(); ++i){
      const auto& slc = sr->slc.at(i);
    }
*/
    return rets;


  });

  const Cut MuonTrackContained([](const caf::SRSliceProxy* slc) {
    int muonTrackIndex = kNuMIMuonCandidateIdx(slc);
    if(muonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
      bool isContained = fv_track.isContained(trk.end.x, trk.end.y, trk.end.z);
      return isContained;
    }
    else{
      return false;
    }
  });
  const Cut MuonTrackExiting([](const caf::SRSliceProxy* slc) {
    int muonTrackIndex = kNuMIMuonCandidateIdx(slc);
    if(muonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
      bool isContained = fv_track.isContained(trk.end.x, trk.end.y, trk.end.z);
      return !isContained;
    }
    else{
      return false;
    }
  });
  // - Truth matching
  const Cut MuonTrackTruthContainedNuMuon([](const caf::SRSliceProxy* slc) {
    int muonTrackIndex = kNuMIMuonCandidateIdx(slc);
    if(muonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
      const auto& trk_truth = trk.truth.p;
      int intID = trk_truth.interaction_id;
      if(intID==-1) return false;

      int trk_truth_pdg = trk_truth.pdg;
      if(isnan(trk_truth.end.x) || isnan(trk_truth.end.y) || isnan(trk_truth.end.z)) return false;
      bool isTrkTruthContained = fv_track.isContained(trk_truth.end.x, trk_truth.end.y, trk_truth.end.z);
      return ( abs(trk_truth_pdg)==13 && isTrkTruthContained );
    }
    else{
      return false;
    }
  });
  const Cut MuonTrackTruthExitingNuMuon([](const caf::SRSliceProxy* slc) {
    int muonTrackIndex = kNuMIMuonCandidateIdx(slc);
    if(muonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
      const auto& trk_truth = trk.truth.p;
      int intID = trk_truth.interaction_id;
      if(intID==-1) return false;
      
      int trk_truth_pdg = trk_truth.pdg;
      if(isnan(trk_truth.end.x) || isnan(trk_truth.end.y) || isnan(trk_truth.end.z)) return false;
      bool isTrkTruthContained = fv_track.isContained(trk_truth.end.x, trk_truth.end.y, trk_truth.end.z);
      return ( abs(trk_truth_pdg)==13 && !isTrkTruthContained );
    }
    else{
      return false;
    }
  });
  const Cut MuonTrackTruthCosmicMuon([](const caf::SRSliceProxy* slc) {
    int muonTrackIndex = kNuMIMuonCandidateIdx(slc);
    if(muonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
      const auto& trk_truth = trk.truth.p;
      int intID = trk_truth.interaction_id;
      return (intID==-1);
    }
    else{
      return false;
    }
  });
  const Cut MuonTrackTruthStoppingProton([](const caf::SRSliceProxy* slc) {
    int muonTrackIndex = kNuMIMuonCandidateIdx(slc);
    if(muonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
      const auto& trk_truth = trk.truth.p;
      const int& endProc = trk_truth.end_process;

      int trk_truth_pdg = trk_truth.pdg;
      if(isnan(trk_truth.end.x) || isnan(trk_truth.end.y) || isnan(trk_truth.end.z)) return false;
      bool isTrkTruthContained = fv_track.isContained(trk_truth.end.x, trk_truth.end.y, trk_truth.end.z);
      bool isStopped = (endProc==2);
      return ( trk_truth_pdg==2212 && isTrkTruthContained && isStopped );
    }
    else{
      return false;
    }
  });
  const Cut MuonTrackTruthInelProton([](const caf::SRSliceProxy* slc) {
    int muonTrackIndex = kNuMIMuonCandidateIdx(slc);
    if(muonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
      const auto& trk_truth = trk.truth.p;
      const int& endProc = trk_truth.end_process;

      int trk_truth_pdg = trk_truth.pdg;
      if(isnan(trk_truth.end.x) || isnan(trk_truth.end.y) || isnan(trk_truth.end.z)) return false;
      bool isTrkTruthContained = fv_track.isContained(trk_truth.end.x, trk_truth.end.y, trk_truth.end.z);
      bool isInel = (endProc==7);
      return ( trk_truth_pdg==2212 && isTrkTruthContained && isInel );
    }
    else{
      return false;
    }
  });
  const Cut MuonTrackTruthOtherProton([](const caf::SRSliceProxy* slc) {
    int muonTrackIndex = kNuMIMuonCandidateIdx(slc);
    if(muonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
      const auto& trk_truth = trk.truth.p;
      const int& endProc = trk_truth.end_process;

      int trk_truth_pdg = trk_truth.pdg;
      if(isnan(trk_truth.end.x) || isnan(trk_truth.end.y) || isnan(trk_truth.end.z)) return false;
      bool isTrkTruthContained = fv_track.isContained(trk_truth.end.x, trk_truth.end.y, trk_truth.end.z);
      bool isStopped = (endProc==2);
      bool isInel = (endProc==7);

      return ( trk_truth_pdg==2212 && 
               ( !isTrkTruthContained ||
                 (isTrkTruthContained && !(isStopped||isInel))
               )
             );
    }
    else{
      return false;
    }
  });
  const Cut MuonTrackTruthStoppingChargedPion([](const caf::SRSliceProxy* slc) {
    int muonTrackIndex = kNuMIMuonCandidateIdx(slc);
    if(muonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
      const auto& trk_truth = trk.truth.p;
      const int& endProc = trk_truth.end_process;

      int trk_truth_pdg = trk_truth.pdg;
      if(isnan(trk_truth.end.x) || isnan(trk_truth.end.y) || isnan(trk_truth.end.z)) return false;
      bool isTrkTruthContained = fv_track.isContained(trk_truth.end.x, trk_truth.end.y, trk_truth.end.z);
      bool isStopped = (endProc==2);
      return ( abs(trk_truth_pdg)==211 && isTrkTruthContained && isStopped );
    }
    else{
      return false;
    }
  });
  const Cut MuonTrackTruthInelChargedPion([](const caf::SRSliceProxy* slc) {
    int muonTrackIndex = kNuMIMuonCandidateIdx(slc);
    if(muonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
      const auto& trk_truth = trk.truth.p;
      const int& endProc = trk_truth.end_process;

      int trk_truth_pdg = trk_truth.pdg;
      if(isnan(trk_truth.end.x) || isnan(trk_truth.end.y) || isnan(trk_truth.end.z)) return false;
      bool isTrkTruthContained = fv_track.isContained(trk_truth.end.x, trk_truth.end.y, trk_truth.end.z);
      bool isInel = (endProc==9) || (endProc==10);
      return ( abs(trk_truth_pdg)==211 && isTrkTruthContained && isInel );
    }
    else{
      return false;
    }
  });
  const Cut MuonTrackTruthOtherChargedPion([](const caf::SRSliceProxy* slc) {
    int muonTrackIndex = kNuMIMuonCandidateIdx(slc);
    if(muonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
      const auto& trk_truth = trk.truth.p;
      const int& endProc = trk_truth.end_process;

      int trk_truth_pdg = trk_truth.pdg;
      if(isnan(trk_truth.end.x) || isnan(trk_truth.end.y) || isnan(trk_truth.end.z)) return false;
      bool isTrkTruthContained = fv_track.isContained(trk_truth.end.x, trk_truth.end.y, trk_truth.end.z);
      bool isStopped = (endProc==2);
      bool isInel = (endProc==9) || (endProc==10);

      return ( abs(trk_truth_pdg)==211 &&
               ( !isTrkTruthContained ||
                 (isTrkTruthContained && !(isStopped||isInel))
               )
             );
    }
    else{
      return false;
    }
  });
  const Cut MuonTrackTruthOther([](const caf::SRSliceProxy* slc) {
    bool IsCategorized = MuonTrackTruthContainedNuMuon(slc);
    IsCategorized |= MuonTrackTruthExitingNuMuon(slc);
    IsCategorized |= MuonTrackTruthCosmicMuon(slc);
    IsCategorized |= MuonTrackTruthStoppingProton(slc);
    IsCategorized |= MuonTrackTruthInelProton(slc);
    IsCategorized |= MuonTrackTruthOtherProton(slc);
    IsCategorized |= MuonTrackTruthStoppingChargedPion(slc);
    IsCategorized |= MuonTrackTruthInelChargedPion(slc);
    IsCategorized |= MuonTrackTruthOtherChargedPion(slc);
    return !IsCategorized;
  });
  // - Comparing truth match to primary
  const Cut MuonTrackTruthMatchedPrimary([](const caf::SRSliceProxy* slc) {
    int muonTrackIndex = kNuMIMuonCandidateIdx(slc); // index of slc->reco.pfp
    int muonPrimrayIndex = ana::PrimaryUtil::MuonIndex_True(slc->truth); // index of slc->truth.prim
    if(muonTrackIndex>=0 && muonPrimrayIndex>=0){
      // reco track
      const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
      const auto& trk_truth = trk.truth.p;
      const auto& trk_truth_G4ID = trk_truth.G4ID;
      // primray
      const auto& prim = slc->truth.prim.at(muonPrimrayIndex);
      const auto& prim_G4ID = prim.G4ID;

      if(trk_truth_G4ID==prim_G4ID) return true;
      else return false;
    }
    else{
      return false;
    }
  });

  // - Comparing truth match to primary
  const Cut ProtonTrackTruthMatchedPrimaryProton([](const caf::SRSliceProxy* slc) {

    int protonTrackIndex = kNuMIProtonCandidateIdx(slc); // index of slc->reco.pfp
    int protonPrimrayIndex = ana::PrimaryUtil::ProtonIndex_True(slc->truth); // index of slc->truth.prim
    if(protonTrackIndex>=0 && protonPrimrayIndex>=0){
      // reco track
      const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
      const auto& trk_truth = trk.truth.p;
      const auto& trk_truth_G4ID = trk_truth.G4ID;
      // primray
      const auto& prim = slc->truth.prim.at(protonPrimrayIndex);
      const auto& prim_G4ID = prim.G4ID;

      if(trk_truth_G4ID==prim_G4ID) return true;
      else return false;
    }
    else{
      return false;
    }
  });

  const Cut IsForTree([](const caf::SRSliceProxy* slc) {
    int cuttpye = kNuMICutType(slc);
    if(cuttpye==0) return false;
    else return true;
  });

  namespace Aux{

    const Cut HasRelaxedMuonTrack([](const caf::SRSliceProxy* slc) {
      return RelaxedMuonTrackIndex(slc)>=0;
    });
    const Cut HasRelaxedProtonTrack([](const caf::SRSliceProxy* slc) {
      return RelaxedProtonTrackIndex(slc)>=0;
    });
    const Cut HasRelaxedChargedPionTrack([](const caf::SRSliceProxy* slc) {
      return RelaxedChargedPionTrackIndex(slc)>=0;
    });

    const Cut RelaxedMuonTrackContained([](const caf::SRSliceProxy* slc) {
      int muonTrackIndex = RelaxedMuonTrackIndex(slc);
      if(muonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
        if(isnan(trk.end.x) || isnan(trk.end.y) || isnan(trk.end.z)) return false;
        bool isTrkContained = fv_track.isContained(trk.end.x, trk.end.y, trk.end.z);
        return isTrkContained;
      }
      else{
        return false;
      }
    });
    const Cut RelaxedProtonTrackContained([](const caf::SRSliceProxy* slc) {
      int protonTrackIndex = RelaxedProtonTrackIndex(slc);
      if(protonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
        if(isnan(trk.end.x) || isnan(trk.end.y) || isnan(trk.end.z)) return false;
        bool isTrkContained = fv_track.isContained(trk.end.x, trk.end.y, trk.end.z);
        return isTrkContained;
      }
      else{
        return false;
      }
    });
    const Cut RelaxedChargedPionTrackContained([](const caf::SRSliceProxy* slc) {
      int cpionTrackIndex = RelaxedChargedPionTrackIndex(slc);
      if(cpionTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(cpionTrackIndex).trk;
        if(isnan(trk.end.x) || isnan(trk.end.y) || isnan(trk.end.z)) return false;
        bool isTrkContained = fv_track.isContained(trk.end.x, trk.end.y, trk.end.z);
        return isTrkContained;
      }
      else{
        return false;
      }
    });

    const Cut RelaxedMuonTrackTruthContainedNuMuon([](const caf::SRSliceProxy* slc) {
      int muonTrackIndex = RelaxedMuonTrackIndex(slc);
      if(muonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
        const auto& trk_truth = trk.truth.p;
        int intID = trk_truth.interaction_id;
        if(intID==-1) return false;

        int trk_truth_pdg = trk_truth.pdg;
        if(isnan(trk_truth.end.x) || isnan(trk_truth.end.y) || isnan(trk_truth.end.z)) return false;
        bool isTrkTruthContained = fv_track.isContained(trk_truth.end.x, trk_truth.end.y, trk_truth.end.z);
        return ( abs(trk_truth_pdg)==13 && isTrkTruthContained );
      }
      else{
        return false;
      }
    });
    const Cut RelaxedMuonTrackTruthExitingNuMuon([](const caf::SRSliceProxy* slc) {
      int muonTrackIndex = RelaxedMuonTrackIndex(slc);
      if(muonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
        const auto& trk_truth = trk.truth.p;
        int intID = trk_truth.interaction_id;
        if(intID==-1) return false;
        
        int trk_truth_pdg = trk_truth.pdg;
        if(isnan(trk_truth.end.x) || isnan(trk_truth.end.y) || isnan(trk_truth.end.z)) return false;
        bool isTrkTruthContained = fv_track.isContained(trk_truth.end.x, trk_truth.end.y, trk_truth.end.z);
        return ( abs(trk_truth_pdg)==13 && !isTrkTruthContained );
      }
      else{
        return false;
      }
    });
    const Cut RelaxedMuonTrackTruthCosmicMuon([](const caf::SRSliceProxy* slc) {
      int muonTrackIndex = RelaxedMuonTrackIndex(slc);
      if(muonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
        const auto& trk_truth = trk.truth.p;
        int intID = trk_truth.interaction_id;
        return (intID==-1);
      }
      else{
        return false;
      }
    });
    const Cut RelaxedMuonTrackTruthStoppingProton([](const caf::SRSliceProxy* slc) {
      int muonTrackIndex = RelaxedMuonTrackIndex(slc);
      if(muonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
        const auto& trk_truth = trk.truth.p;
        const int& endProc = trk_truth.end_process;

        int trk_truth_pdg = trk_truth.pdg;
        if(isnan(trk_truth.end.x) || isnan(trk_truth.end.y) || isnan(trk_truth.end.z)) return false;
        bool isTrkTruthContained = fv_track.isContained(trk_truth.end.x, trk_truth.end.y, trk_truth.end.z);
        bool isStopped = (endProc==2);
        return ( trk_truth_pdg==2212 && isTrkTruthContained && isStopped );
      }
      else{
        return false;
      }
    });
    const Cut RelaxedMuonTrackTruthInelProton([](const caf::SRSliceProxy* slc) {
      int muonTrackIndex = RelaxedMuonTrackIndex(slc);
      if(muonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
        const auto& trk_truth = trk.truth.p;
        const int& endProc = trk_truth.end_process;

        int trk_truth_pdg = trk_truth.pdg;
        if(isnan(trk_truth.end.x) || isnan(trk_truth.end.y) || isnan(trk_truth.end.z)) return false;
        bool isTrkTruthContained = fv_track.isContained(trk_truth.end.x, trk_truth.end.y, trk_truth.end.z);
        bool isInel = (endProc==7);
        return ( trk_truth_pdg==2212 && isTrkTruthContained && isInel );
      }
      else{
        return false;
      }
    });
    const Cut RelaxedMuonTrackTruthOtherProton([](const caf::SRSliceProxy* slc) {
      int muonTrackIndex = RelaxedMuonTrackIndex(slc);
      if(muonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
        const auto& trk_truth = trk.truth.p;
        const int& endProc = trk_truth.end_process;

        int trk_truth_pdg = trk_truth.pdg;
        if(isnan(trk_truth.end.x) || isnan(trk_truth.end.y) || isnan(trk_truth.end.z)) return false;
        bool isTrkTruthContained = fv_track.isContained(trk_truth.end.x, trk_truth.end.y, trk_truth.end.z);
        bool isStopped = (endProc==2);
        bool isInel = (endProc==7);

        return ( trk_truth_pdg==2212 && 
                 ( !isTrkTruthContained ||
                   (isTrkTruthContained && !(isStopped||isInel))
                 )
               );
      }
      else{
        return false;
      }
    });
    const Cut RelaxedMuonTrackTruthStoppingChargedPion([](const caf::SRSliceProxy* slc) {
      int muonTrackIndex = RelaxedMuonTrackIndex(slc);
      if(muonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
        const auto& trk_truth = trk.truth.p;
        const int& endProc = trk_truth.end_process;

        int trk_truth_pdg = trk_truth.pdg;
        if(isnan(trk_truth.end.x) || isnan(trk_truth.end.y) || isnan(trk_truth.end.z)) return false;
        bool isTrkTruthContained = fv_track.isContained(trk_truth.end.x, trk_truth.end.y, trk_truth.end.z);
        bool isStopped = (endProc==2);
        return ( abs(trk_truth_pdg)==211 && isTrkTruthContained && isStopped );
      }
      else{
        return false;
      }
    });
    const Cut RelaxedMuonTrackTruthInelChargedPion([](const caf::SRSliceProxy* slc) {
      int muonTrackIndex = RelaxedMuonTrackIndex(slc);
      if(muonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
        const auto& trk_truth = trk.truth.p;
        const int& endProc = trk_truth.end_process;

        int trk_truth_pdg = trk_truth.pdg;
        if(isnan(trk_truth.end.x) || isnan(trk_truth.end.y) || isnan(trk_truth.end.z)) return false;
        bool isTrkTruthContained = fv_track.isContained(trk_truth.end.x, trk_truth.end.y, trk_truth.end.z);
        bool isInel = (endProc==9) || (endProc==10);
        return ( abs(trk_truth_pdg)==211 && isTrkTruthContained && isInel );
      }
      else{
        return false;
      }
    });
    const Cut RelaxedMuonTrackTruthOtherChargedPion([](const caf::SRSliceProxy* slc) {
      int muonTrackIndex = RelaxedMuonTrackIndex(slc);
      if(muonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
        const auto& trk_truth = trk.truth.p;
        const int& endProc = trk_truth.end_process;

        int trk_truth_pdg = trk_truth.pdg;
        if(isnan(trk_truth.end.x) || isnan(trk_truth.end.y) || isnan(trk_truth.end.z)) return false;
        bool isTrkTruthContained = fv_track.isContained(trk_truth.end.x, trk_truth.end.y, trk_truth.end.z);
        bool isStopped = (endProc==2);
        bool isInel = (endProc==9) || (endProc==10);

        return ( abs(trk_truth_pdg)==211 &&
                 ( !isTrkTruthContained ||
                   (isTrkTruthContained && !(isStopped||isInel))
                 )
               );
      }
      else{
        return false;
      }
    });


    const Cut RelaxedMuonTrackTruthOther([](const caf::SRSliceProxy* slc) {
      bool IsCategorized = RelaxedMuonTrackTruthContainedNuMuon(slc);
      IsCategorized |= RelaxedMuonTrackTruthExitingNuMuon(slc);
      IsCategorized |= RelaxedMuonTrackTruthCosmicMuon(slc);
      IsCategorized |= RelaxedMuonTrackTruthStoppingProton(slc);
      IsCategorized |= RelaxedMuonTrackTruthInelProton(slc);
      IsCategorized |= RelaxedMuonTrackTruthOtherProton(slc);
      IsCategorized |= RelaxedMuonTrackTruthStoppingChargedPion(slc);
      IsCategorized |= RelaxedMuonTrackTruthInelChargedPion(slc);
      IsCategorized |= RelaxedMuonTrackTruthOtherChargedPion(slc);
      return !IsCategorized;
    });


    const Cut RelaxedProtonTrackTruthContainedNuMuon([](const caf::SRSliceProxy* slc) {
      int protonTrackIndex = RelaxedProtonTrackIndex(slc);
      if(protonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
        const auto& trk_truth = trk.truth.p;
        int intID = trk_truth.interaction_id;
        if(intID==-1) return false;

        int trk_truth_pdg = trk_truth.pdg;
        if(isnan(trk_truth.end.x) || isnan(trk_truth.end.y) || isnan(trk_truth.end.z)) return false;
        bool isTrkTruthContained = fv_track.isContained(trk_truth.end.x, trk_truth.end.y, trk_truth.end.z);
        return ( abs(trk_truth_pdg)==13 && isTrkTruthContained );
      }
      else{
        return false;
      }
    });
    const Cut RelaxedProtonTrackTruthExitingNuMuon([](const caf::SRSliceProxy* slc) {
      int protonTrackIndex = RelaxedProtonTrackIndex(slc);
      if(protonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
        const auto& trk_truth = trk.truth.p;
        int intID = trk_truth.interaction_id;
        if(intID==-1) return false;
        
        int trk_truth_pdg = trk_truth.pdg;
        if(isnan(trk_truth.end.x) || isnan(trk_truth.end.y) || isnan(trk_truth.end.z)) return false;
        bool isTrkTruthContained = fv_track.isContained(trk_truth.end.x, trk_truth.end.y, trk_truth.end.z);
        return ( abs(trk_truth_pdg)==13 && !isTrkTruthContained );
      }
      else{
        return false;
      }
    });
    const Cut RelaxedProtonTrackTruthCosmicMuon([](const caf::SRSliceProxy* slc) {
      int protonTrackIndex = RelaxedProtonTrackIndex(slc);
      if(protonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
        const auto& trk_truth = trk.truth.p;
        int intID = trk_truth.interaction_id;
        return (intID==-1);
      }
      else{
        return false;
      }
    });
    const Cut RelaxedProtonTrackTruthStoppingProton([](const caf::SRSliceProxy* slc) {
      int protonTrackIndex = RelaxedProtonTrackIndex(slc);
      if(protonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
        const auto& trk_truth = trk.truth.p;
        const int& endProc = trk_truth.end_process;

        int trk_truth_pdg = trk_truth.pdg;
        if(isnan(trk_truth.end.x) || isnan(trk_truth.end.y) || isnan(trk_truth.end.z)) return false;
        bool isTrkTruthContained = fv_track.isContained(trk_truth.end.x, trk_truth.end.y, trk_truth.end.z);
        bool isStopped = (endProc==2);
        return ( trk_truth_pdg==2212 && isTrkTruthContained && isStopped );
      }
      else{
        return false;
      }
    });
    const Cut RelaxedProtonTrackTruthInelProton([](const caf::SRSliceProxy* slc) {
      int protonTrackIndex = RelaxedProtonTrackIndex(slc);
      if(protonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
        const auto& trk_truth = trk.truth.p;
        const int& endProc = trk_truth.end_process;
        
        int trk_truth_pdg = trk_truth.pdg;
        if(isnan(trk_truth.end.x) || isnan(trk_truth.end.y) || isnan(trk_truth.end.z)) return false;
        bool isTrkTruthContained = fv_track.isContained(trk_truth.end.x, trk_truth.end.y, trk_truth.end.z);
        bool isInel = (endProc==7);
        return ( trk_truth_pdg==2212 && isTrkTruthContained && isInel );
      }
      else{
        return false;
      }
    });
    const Cut RelaxedProtonTrackTruthOtherProton([](const caf::SRSliceProxy* slc) {
      int protonTrackIndex = RelaxedProtonTrackIndex(slc);
      if(protonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
        const auto& trk_truth = trk.truth.p;
        const int& endProc = trk_truth.end_process;

        int trk_truth_pdg = trk_truth.pdg;
        if(isnan(trk_truth.end.x) || isnan(trk_truth.end.y) || isnan(trk_truth.end.z)) return false;
        bool isTrkTruthContained = fv_track.isContained(trk_truth.end.x, trk_truth.end.y, trk_truth.end.z);
        bool isStopped = (endProc==2);
        bool isInel = (endProc==7);

        return ( trk_truth_pdg==2212 &&
                 ( !isTrkTruthContained ||
                   (isTrkTruthContained && !(isStopped||isInel))
                 )
               );
      }
      else{
        return false;
      }
    });
    const Cut RelaxedProtonTrackTruthStoppingChargedPion([](const caf::SRSliceProxy* slc) {
      int protonTrackIndex = RelaxedProtonTrackIndex(slc);
      if(protonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
        const auto& trk_truth = trk.truth.p;
        const int& endProc = trk_truth.end_process;

        int trk_truth_pdg = trk_truth.pdg;
        if(isnan(trk_truth.end.x) || isnan(trk_truth.end.y) || isnan(trk_truth.end.z)) return false;
        bool isTrkTruthContained = fv_track.isContained(trk_truth.end.x, trk_truth.end.y, trk_truth.end.z);
        bool isStopped = (endProc==2);
        return ( abs(trk_truth_pdg)==211 && isTrkTruthContained && isStopped );
      }
      else{
        return false;
      }
    });
    const Cut RelaxedProtonTrackTruthInelChargedPion([](const caf::SRSliceProxy* slc) {
      int protonTrackIndex = RelaxedProtonTrackIndex(slc);
      if(protonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
        const auto& trk_truth = trk.truth.p;
        const int& endProc = trk_truth.end_process;

        int trk_truth_pdg = trk_truth.pdg;
        if(isnan(trk_truth.end.x) || isnan(trk_truth.end.y) || isnan(trk_truth.end.z)) return false;
        bool isTrkTruthContained = fv_track.isContained(trk_truth.end.x, trk_truth.end.y, trk_truth.end.z);
        bool isInel = (endProc==9) || (endProc==10);
        return ( abs(trk_truth_pdg)==211 && isTrkTruthContained && isInel );
      }
      else{
        return false;
      }
    });
    const Cut RelaxedProtonTrackTruthOtherChargedPion([](const caf::SRSliceProxy* slc) {
      int protonTrackIndex = RelaxedProtonTrackIndex(slc);
      if(protonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
        const auto& trk_truth = trk.truth.p;
        const int& endProc = trk_truth.end_process;

        int trk_truth_pdg = trk_truth.pdg;
        if(isnan(trk_truth.end.x) || isnan(trk_truth.end.y) || isnan(trk_truth.end.z)) return false;
        bool isTrkTruthContained = fv_track.isContained(trk_truth.end.x, trk_truth.end.y, trk_truth.end.z);
        bool isStopped = (endProc==2);
        bool isInel = (endProc==9) || (endProc==10);

        return ( abs(trk_truth_pdg)==211 &&
                 ( !isTrkTruthContained ||
                   (isTrkTruthContained && !(isStopped||isInel))
                 )
               );
      }
      else{
        return false;
      }
    });
    const Cut RelaxedProtonTrackTruthOther([](const caf::SRSliceProxy* slc) {
      bool IsCategorized = RelaxedProtonTrackTruthContainedNuMuon(slc);
      IsCategorized |= RelaxedProtonTrackTruthExitingNuMuon(slc);
      IsCategorized |= RelaxedProtonTrackTruthCosmicMuon(slc);
      IsCategorized |= RelaxedProtonTrackTruthStoppingProton(slc);
      IsCategorized |= RelaxedProtonTrackTruthInelProton(slc);
      IsCategorized |= RelaxedProtonTrackTruthOtherProton(slc);
      IsCategorized |= RelaxedProtonTrackTruthStoppingChargedPion(slc);
      IsCategorized |= RelaxedProtonTrackTruthInelChargedPion(slc);
      IsCategorized |= RelaxedProtonTrackTruthOtherChargedPion(slc);
      return !IsCategorized;
    });

    const Cut RelaxedChargedPionTrackTruthContainedNuMuon([](const caf::SRSliceProxy* slc) {
      int cpionTrackIndex = RelaxedChargedPionTrackIndex(slc);
      if(cpionTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(cpionTrackIndex).trk;
        const auto& trk_truth = trk.truth.p;
        int intID = trk_truth.interaction_id;
        if(intID==-1) return false;

        int trk_truth_pdg = trk_truth.pdg;
        if(isnan(trk_truth.end.x) || isnan(trk_truth.end.y) || isnan(trk_truth.end.z)) return false;
        bool isTrkTruthContained = fv_track.isContained(trk_truth.end.x, trk_truth.end.y, trk_truth.end.z);
        return ( abs(trk_truth_pdg)==13 && isTrkTruthContained );
      }
      else{
        return false;
      }
    });
    const Cut RelaxedChargedPionTrackTruthExitingNuMuon([](const caf::SRSliceProxy* slc) {
      int cpionTrackIndex = RelaxedChargedPionTrackIndex(slc);
      if(cpionTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(cpionTrackIndex).trk;
        const auto& trk_truth = trk.truth.p;
        int intID = trk_truth.interaction_id;
        if(intID==-1) return false;
        
        int trk_truth_pdg = trk_truth.pdg;
        if(isnan(trk_truth.end.x) || isnan(trk_truth.end.y) || isnan(trk_truth.end.z)) return false;
        bool isTrkTruthContained = fv_track.isContained(trk_truth.end.x, trk_truth.end.y, trk_truth.end.z);
        return ( abs(trk_truth_pdg)==13 && !isTrkTruthContained );
      }
      else{
        return false;
      }
    });
    const Cut RelaxedChargedPionTrackTruthCosmicMuon([](const caf::SRSliceProxy* slc) {
      int cpionTrackIndex = RelaxedChargedPionTrackIndex(slc);
      if(cpionTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(cpionTrackIndex).trk;
        const auto& trk_truth = trk.truth.p;
        int intID = trk_truth.interaction_id;
        return (intID==-1);
      }
      else{
        return false;
      }
    });
    const Cut RelaxedChargedPionTrackTruthStoppingProton([](const caf::SRSliceProxy* slc) {
      int cpionTrackIndex = RelaxedChargedPionTrackIndex(slc);
      if(cpionTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(cpionTrackIndex).trk;
        const auto& trk_truth = trk.truth.p;
        const int& endProc = trk_truth.end_process;

        int trk_truth_pdg = trk_truth.pdg;
        if(isnan(trk_truth.end.x) || isnan(trk_truth.end.y) || isnan(trk_truth.end.z)) return false;
        bool isTrkTruthContained = fv_track.isContained(trk_truth.end.x, trk_truth.end.y, trk_truth.end.z);
        bool isStopped = (endProc==2);
        return ( trk_truth_pdg==2212 && isTrkTruthContained && isStopped );
      }
      else{
        return false;
      }
    });
    const Cut RelaxedChargedPionTrackTruthInelProton([](const caf::SRSliceProxy* slc) {
      int cpionTrackIndex = RelaxedChargedPionTrackIndex(slc);
      if(cpionTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(cpionTrackIndex).trk;
        const auto& trk_truth = trk.truth.p;
        const int& endProc = trk_truth.end_process;

        int trk_truth_pdg = trk_truth.pdg;
        if(isnan(trk_truth.end.x) || isnan(trk_truth.end.y) || isnan(trk_truth.end.z)) return false;
        bool isTrkTruthContained = fv_track.isContained(trk_truth.end.x, trk_truth.end.y, trk_truth.end.z);
        bool isInel = (endProc==7);
        return ( trk_truth_pdg==2212 && isTrkTruthContained && isInel );
      }
      else{
        return false;
      }
    });
    const Cut RelaxedChargedPionTrackTruthOtherProton([](const caf::SRSliceProxy* slc) {
      int cpionTrackIndex = RelaxedChargedPionTrackIndex(slc);
      if(cpionTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(cpionTrackIndex).trk;
        const auto& trk_truth = trk.truth.p;
        const int& endProc = trk_truth.end_process;

        int trk_truth_pdg = trk_truth.pdg;
        if(isnan(trk_truth.end.x) || isnan(trk_truth.end.y) || isnan(trk_truth.end.z)) return false;
        bool isTrkTruthContained = fv_track.isContained(trk_truth.end.x, trk_truth.end.y, trk_truth.end.z);
        bool isStopped = (endProc==2);
        bool isInel = (endProc==7);

        return ( trk_truth_pdg==2212 && 
                 ( !isTrkTruthContained ||
                   (isTrkTruthContained && !(isStopped||isInel))
                 )
               );
      }
      else{
        return false;
      }
    });
    const Cut RelaxedChargedPionTrackTruthStoppingChargedPion([](const caf::SRSliceProxy* slc) {
      int cpionTrackIndex = RelaxedChargedPionTrackIndex(slc);
      if(cpionTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(cpionTrackIndex).trk;
        const auto& trk_truth = trk.truth.p;
        const int& endProc = trk_truth.end_process;

        int trk_truth_pdg = trk_truth.pdg;
        if(isnan(trk_truth.end.x) || isnan(trk_truth.end.y) || isnan(trk_truth.end.z)) return false;
        bool isTrkTruthContained = fv_track.isContained(trk_truth.end.x, trk_truth.end.y, trk_truth.end.z);
        bool isStopped = (endProc==2);
        return ( abs(trk_truth_pdg)==211 && isTrkTruthContained && isStopped );
      }
      else{
        return false;
      }
    });
    const Cut RelaxedChargedPionTrackTruthInelChargedPion([](const caf::SRSliceProxy* slc) {
      int cpionTrackIndex = RelaxedChargedPionTrackIndex(slc);
      if(cpionTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(cpionTrackIndex).trk;
        const auto& trk_truth = trk.truth.p;
        const int& endProc = trk_truth.end_process;

        int trk_truth_pdg = trk_truth.pdg;
        if(isnan(trk_truth.end.x) || isnan(trk_truth.end.y) || isnan(trk_truth.end.z)) return false;
        bool isTrkTruthContained = fv_track.isContained(trk_truth.end.x, trk_truth.end.y, trk_truth.end.z);
        bool isInel = (endProc==9) || (endProc==10);
        return ( abs(trk_truth_pdg)==211 && isTrkTruthContained && isInel );
      }
      else{
        return false;
      }
    });
    const Cut RelaxedChargedPionTrackTruthOtherChargedPion([](const caf::SRSliceProxy* slc) {
      int cpionTrackIndex = RelaxedChargedPionTrackIndex(slc);
      if(cpionTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(cpionTrackIndex).trk;
        const auto& trk_truth = trk.truth.p;
        const int& endProc = trk_truth.end_process;

        int trk_truth_pdg = trk_truth.pdg;
        if(isnan(trk_truth.end.x) || isnan(trk_truth.end.y) || isnan(trk_truth.end.z)) return false;
        bool isTrkTruthContained = fv_track.isContained(trk_truth.end.x, trk_truth.end.y, trk_truth.end.z);
        bool isStopped = (endProc==2);
        bool isInel = (endProc==9) || (endProc==10);

        return ( abs(trk_truth_pdg)==211 &&
                 ( !isTrkTruthContained ||
                   (isTrkTruthContained && !(isStopped||isInel))
                 )
               );
      }
      else{
        return false;
      }
    });

    const Cut RelaxedChargedPionTrackTruthOther([](const caf::SRSliceProxy* slc) {
      bool IsCategorized = RelaxedChargedPionTrackTruthContainedNuMuon(slc);
      IsCategorized |= RelaxedChargedPionTrackTruthExitingNuMuon(slc);
      IsCategorized |= RelaxedChargedPionTrackTruthCosmicMuon(slc);
      IsCategorized |= RelaxedChargedPionTrackTruthStoppingProton(slc);
      IsCategorized |= RelaxedChargedPionTrackTruthInelProton(slc);
      IsCategorized |= RelaxedChargedPionTrackTruthOtherProton(slc);
      IsCategorized |= RelaxedChargedPionTrackTruthStoppingChargedPion(slc);
      IsCategorized |= RelaxedChargedPionTrackTruthInelChargedPion(slc);
      IsCategorized |= RelaxedChargedPionTrackTruthOtherChargedPion(slc);

      return !IsCategorized;
    });

    // some presels
    // - muon
    const Cut RelaxedMuonTrackIsochronous([](const caf::SRSliceProxy* slc) { 
      int muonTrackIndex = RelaxedMuonTrackIndex(slc);
      if(muonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
        return ( fabs(trk.dir.x)<0.1 );

      }
      else{
        return false;
      }
    });

    const Cut RelaxedMuonTrackDriftDirection([](const caf::SRSliceProxy* slc) { 
      int muonTrackIndex = RelaxedMuonTrackIndex(slc);
      if(muonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
        return ( fabs(trk.dir.x)>0.8 );

      }
      else{
        return false;
      }
    });
    const Cut RelaxedMuonTrackLengthCut([](const caf::SRSliceProxy* slc) {
      int muonTrackIndex = RelaxedMuonTrackIndex(slc);
      if(muonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
        if(isnan(trk.len)) return false;
        else return (trk.len>50.0);
        //else return (trk.len>0.);
      }
      else{
        return false;
      }
    });
    const Cut RelaxedMuonTrackTruthPResFracTail([](const caf::SRSliceProxy* slc) {
      double prf = RelaxedMuonTrackTruthPResFrac(slc);
      return (prf>-999 && prf<-0.2);
    });
    // - proton
    const Cut RelaxedProtonTrackIsochronous([](const caf::SRSliceProxy* slc) {
      int protonTrackIndex = RelaxedProtonTrackIndex(slc);
      if(protonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
        return ( fabs(trk.dir.x)<0.1 );
       
      }
      else{
        return false;
      }
    });

    const Cut RelaxedProtonTrackDriftDirection([](const caf::SRSliceProxy* slc) {
      int protonTrackIndex = RelaxedProtonTrackIndex(slc);
      if(protonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
        return ( fabs(trk.dir.x)>0.8 );

      }
      else{
        return false;
      }
    });

    const Cut RelaxedProtonTrackHasChi2MuonExcess([](const caf::SRSliceProxy* slc) {
      int protonTrackIndex = RelaxedProtonTrackIndex(slc);
      if(protonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
        if(trk.bestplane == -1) return false;
        const float Chi2Muon = trk.chi2pid[2].chi2_muon;
        return (Chi2Muon>40. && Chi2Muon<50.);
      }
      else{
        return false;
      }
    });

    const Cut RelaxedProtonTrackHasChi2MuonDeficit([](const caf::SRSliceProxy* slc) {
      int protonTrackIndex = RelaxedProtonTrackIndex(slc);
      if(protonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
        if(trk.bestplane == -1) return false;
        const float Chi2Muon = trk.chi2pid[2].chi2_muon;
        return (Chi2Muon>50.);
      }
      else{
        return false;
      }
    });
    const Cut RelaxedProtonTrackShort([](const caf::SRSliceProxy* slc) {
      int protonTrackIndex = RelaxedProtonTrackIndex(slc);
      if(protonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
        return (trk.len<10.);
      }
      else{
        return false;
      }
    });
    // - charged pion
    const Cut RelaxedChargedPionTrackIsochronous([](const caf::SRSliceProxy* slc) {
      int cpionTrackIndex = RelaxedChargedPionTrackIndex(slc);
      if(cpionTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(cpionTrackIndex).trk;
        return ( fabs(trk.dir.x)<0.1 );

      }
      else{
        return false;
      }
    });

  } // end namespace Aux


} // end namespace TwoTrack

} // end namespace ICARUSNumuXsec
