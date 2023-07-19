#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_TwoTrack_Cuts.h"

using namespace ana;
using namespace std;

namespace ICARUSNumuXsec{

namespace TwoTrack{

  const Cut HasTwoPrimaryTracks([](const caf::SRSliceProxy* slc) {
    return ICARUSNumuXsec::PrimaryTrackIndices(slc).size()>=2;
  });
  const Cut HasOnlyTwoPrimaryTracks([](const caf::SRSliceProxy* slc) {
    return ICARUSNumuXsec::PrimaryTrackIndices(slc).size()==2;
  });
  const Cut HasMuonTrack([](const caf::SRSliceProxy* slc) {
    return MuonTrackIndex(slc)>=0.;
  });

  const Cut HasProtonTrack([](const caf::SRSliceProxy* slc) {
    return ProtonTrackIndex(slc)>=0;
  });

  const Cut MuonTrackContained([](const caf::SRSliceProxy* slc) {
    int muonTrackIndex = MuonTrackIndex(slc);
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
    int muonTrackIndex = MuonTrackIndex(slc);
    if(muonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
      bool isContained = fv_track.isContained(trk.end.x, trk.end.y, trk.end.z);
      return !isContained;
    }
    else{
      return false;
    }
  });
  const Cut MuonTrackOneMeter([](const caf::SRSliceProxy* slc) {
    int muonTrackIndex = MuonTrackIndex(slc);
    if(muonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
      return trk.len>100.;
    }
    else{
      return false;
    }
  });
  const Cut MuonTrackBelowBlindP([](const caf::SRSliceProxy* slc) {
    double muonTrackRecoP = MuonTrackP(slc);
    return muonTrackRecoP<0.6;
  });

  const Cut MuonTrackTruthContainedNuMuon([](const caf::SRSliceProxy* slc) {
    int muonTrackIndex = MuonTrackIndex(slc);
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
    int muonTrackIndex = MuonTrackIndex(slc);
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
    int muonTrackIndex = MuonTrackIndex(slc);
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
    int muonTrackIndex = MuonTrackIndex(slc);
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
    int muonTrackIndex = MuonTrackIndex(slc);
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
    int muonTrackIndex = MuonTrackIndex(slc);
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
    int muonTrackIndex = MuonTrackIndex(slc);
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
    int muonTrackIndex = MuonTrackIndex(slc);
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
    int muonTrackIndex = MuonTrackIndex(slc);
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


  const Cut ProtonPCut([](const caf::SRSliceProxy* slc) {
    double protonP = ProtonTrackP(slc);
    if(protonP<0){
      return false;
    }
    else{
      return (protonP>0.6);
    }
  });

  // Hadron (non-muon) contrained
  const Cut HadronContained([](const caf::SRSliceProxy* slc) {
/*
    vector<double> nonMuonTrackIndcies = NonMuonTrackIndecies(slc);
    if(nonMuonTrackIndcies.size()==0){
      return true;
    }
    else{
      for(const auto& trkIdx: nonMuonTrackIndcies){
        const auto& pfp = slc->reco.pfp.at(trkIdx);
        const auto& trk = pfp.trk;

        if(isnan(trk.start.x)) continue;

        const float Atslc = std::hypot(slc->vertex.x - trk.start.x,
                                       slc->vertex.y - trk.start.y,
                                       slc->vertex.z - trk.start.z);
        const bool AtSlice = ( Atslc < 10.0 && pfp.parent_is_primary);
        const bool Contained = fv_track.isContained(trk.end.x, trk.end.y, trk.end.z);

        if(AtSlice){
          if(!Contained) return false;
        }

      }
      return true;
    }
*/

    int muonTrackIndex = MuonTrackIndex(slc);
    if(muonTrackIndex>=0){

      for(unsigned int i_pfp=0; i_pfp<slc->reco.pfp.size(); i_pfp++){

        const auto& pfp = slc->reco.pfp.at(i_pfp);
        const auto& trk = pfp.trk;

        if(i_pfp==(unsigned int )muonTrackIndex) continue;
        if(isnan(trk.start.x)) continue;

        const float Atslc = std::hypot(slc->vertex.x - trk.start.x,
                                       slc->vertex.y - trk.start.y,
                                       slc->vertex.z - trk.start.z);
        const bool AtSlice = ( Atslc < 10.0 && pfp.parent_is_primary);
        const bool Contained = fv_track.isContained(trk.end.x, trk.end.y, trk.end.z);

        if(AtSlice){
          if(!Contained) return false;
        }

      }
      return true;

    }
    else{
      return false;
    }

  });

  // pion tagging
  const Cut HasChargedPionTrack([](const caf::SRSliceProxy* slc) {
    return ChargedPionTrackIndex(slc)>=0;
  });
  const Cut HasStoppedChargedPionTrack([](const caf::SRSliceProxy* slc) {
    return StoppedChargedPionTrackIndex(slc)>=0;
  });
  const Cut StoppedChargedPionTrackLongEnough([](const caf::SRSliceProxy* slc) {
    int stoppedCPionIndex = StoppedChargedPionTrackIndex(slc);
    if(stoppedCPionIndex>=0){
      const auto& trk = slc->reco.pfp.at(stoppedCPionIndex).trk;
      if(trk.len>20.) return true;
      else return false;
    }
    else{
      return false;
    }
  });

  const Cut HasInelasticChargedPionTrack([](const caf::SRSliceProxy* slc) {
    return InelasticChargedPionTrackIndex(slc)>=0;
  });
  const Cut InelasticChargedPionTrackLongEnough([](const caf::SRSliceProxy* slc) {
    int inelCPionIndex = InelasticChargedPionTrackIndex(slc);
    if(inelCPionIndex>=0){
      const auto& trk = slc->reco.pfp.at(inelCPionIndex).trk;
      if(trk.len>20.) return true;
      else return false;
    }
    else{
      return false;
    }
  });

  // Signal def

  const Cut SignalDef([](const caf::SRSliceProxy* slc) {

    bool bk_UseGHepRecord = intt.UseGHepRecord;
    intt.UseGHepRecord = false;
    ICARUSNumuXsec::InteractionTool::NParticles nptls = intt.GetNParticles(slc);
    intt.UseGHepRecord = bk_UseGHepRecord;

    bool isCCQELike = (nptls.NMuon==1 && nptls.NProton>=1 && nptls.NPip+nptls.NPim+nptls.NPi0==0);
    if(!isCCQELike) return false;

    bool VertexContained = false;
    if( !isnan(slc->truth.position.x) ){
      VertexContained = fv.isContained(slc->truth.position.x, slc->truth.position.y, slc->truth.position.z);
    }
    if(!VertexContained) return false;

    const auto& mu_ghep = slc->truth.prim[ intt.MuonIndices[0] ];

    double MuonP = sqrt(mu_ghep.genp.x*mu_ghep.genp.x + mu_ghep.genp.y*mu_ghep.genp.y + mu_ghep.genp.z*mu_ghep.genp.z);

    double minMuonP = mu_ghep.contained ? 0.226 : 0.340;

    if(MuonP<minMuonP) return false;

    double LargestProtonP = -999.;
    bool HasExitingProton = false;
    for(const auto& protonInx: intt.ProtonIndices){
      const auto& p_ghep = slc->truth.prim[ protonInx ];
      double ProtonP = sqrt(p_ghep.genp.x*p_ghep.genp.x + p_ghep.genp.y*p_ghep.genp.y + p_ghep.genp.z*p_ghep.genp.z);
      LargestProtonP = max(LargestProtonP, ProtonP);

      bool IsProtonContained = fv_track.isContained(p_ghep.end.x, p_ghep.end.y, p_ghep.end.z);
      if(!IsProtonContained) HasExitingProton = true;

    }
    if(LargestProtonP<0.400) return false;
    if(HasExitingProton) return false;


    return true;

  });
  const Cut SignalMuonContained([](const caf::SRSliceProxy* slc) {

    bool bk_UseGHepRecord = intt.UseGHepRecord;
    intt.UseGHepRecord = false;
    ICARUSNumuXsec::InteractionTool::NParticles nptls = intt.GetNParticles(slc);
    intt.UseGHepRecord = bk_UseGHepRecord;

    if(nptls.NMuon!=1) return false;

    const auto& mu_ghep = slc->truth.prim[ intt.MuonIndices[0] ];
/*
    bool MuonStartContained = fv_track.isContained(mu_ghep.start.x, mu_ghep.start.y, mu_ghep.start.z);
    bool MuonEndContained = fv_track.isContained(mu_ghep.end.x, mu_ghep.end.y, mu_ghep.end.z);

    return MuonStartContained&&MuonEndContained;
*/
    if(mu_ghep.contained) return true;
    else return false;
  });

  const Cut SignalMuonExiting([](const caf::SRSliceProxy* slc) {

    bool bk_UseGHepRecord = intt.UseGHepRecord;
    intt.UseGHepRecord = false;
    ICARUSNumuXsec::InteractionTool::NParticles nptls = intt.GetNParticles(slc);
    intt.UseGHepRecord = bk_UseGHepRecord;
    
    if(nptls.NMuon!=1) return false;
    
    const auto& mu_ghep = slc->truth.prim[ intt.MuonIndices[0] ];

    //double MuonP = sqrt(mu_ghep.genp.x*mu_ghep.genp.x + mu_ghep.genp.y*mu_ghep.genp.y + mu_ghep.genp.z*mu_ghep.genp.z);
    //if(MuonP<0.340) return false;

/*
    bool MuonStartContained = fv_track.isContained(mu_ghep.start.x, mu_ghep.start.y, mu_ghep.start.z);
    bool MuonEndContained = fv_track.isContained(mu_ghep.end.x, mu_ghep.end.y, mu_ghep.end.z);
    
    return !(MuonStartContained&&MuonEndContained);
*/

    if(mu_ghep.contained) return false;
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
