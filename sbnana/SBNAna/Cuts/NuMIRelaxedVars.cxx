#include "sbnana/SBNAna/Cuts/NuMIRelaxedVars.h"
#include "sbnana/SBNAna/Cuts/NuMIXSecCuts.h"
#include "sbnana/SBNAna/Vars/NuMIXSecVars.h"

#include <iostream>

namespace ana {

  // Cut
  const Cut kNuMICRLongestTrackDirYHardCut([](const caf::SRSliceProxy* slc) {
    return (slc->nuid.crlongtrkdiry>-0.4);
  });

  // Muon
  // - Var
  const Var kNuMIRelaxedMuonTrackIdx([](const caf::SRSliceProxy* slc) -> int {
    if(slc->reco.pfp.size()==0){
      return -999.;
    }
    else{
      int muonIdx = -1;
      double muonMaxLength = -1;
      for(unsigned int i_pfp=0; i_pfp<slc->reco.pfp.size(); i_pfp++){
        const auto& pfp = slc->reco.pfp.at(i_pfp);
        const auto& trk = pfp.trk;

        if(trk.bestplane == -1) continue;

        // First we calculate the distance of each track to the slice vertex.
        if(isnan(trk.start.x)) continue;
        const float Atslc = std::hypot(slc->vertex.x - trk.start.x,
                                       slc->vertex.y - trk.start.y,
                                       slc->vertex.z - trk.start.z);

        // We require that the distance of the track from the slice is less than
        // 10 cm and that the parent of the track has been marked as the primary.
        const bool AtSlice = ( Atslc < 10.0 && pfp.parent_is_primary);

        if ( AtSlice && trk.len > muonMaxLength ) {
          muonIdx = i_pfp;
          muonMaxLength = trk.len;
        }
      }
      return muonIdx;

    }
  });
  const Var kNuMIRelaxedMuonTrackLength([](const caf::SRSliceProxy* slc) -> double {
    int muonTrackIndex = kNuMIRelaxedMuonTrackIdx(slc);
    if(muonTrackIndex>=0){
      return slc->reco.pfp.at(muonTrackIndex).trk.len;
    }
    else{
      return -999.;
    }
  });
  const Var kNuMIRelaxedMuonTrackChi2Muon([](const caf::SRSliceProxy* slc) -> double {
    int muonTrackIndex = kNuMIRelaxedMuonTrackIdx(slc);
    if(muonTrackIndex>=0){
      return slc->reco.pfp.at(muonTrackIndex).trk.chi2pid[2].chi2_muon;
    }
    else{
      return -999.;
    }
  });
  const Var kNuMIRelaxedMuonTrackChi2Proton([](const caf::SRSliceProxy* slc) -> double {
    int muonTrackIndex = kNuMIRelaxedMuonTrackIdx(slc);
    if(muonTrackIndex>=0){
      return slc->reco.pfp.at(muonTrackIndex).trk.chi2pid[2].chi2_proton;
    }
    else{
      return -999.;
    }
  });
  const Var kNuMIRelaxedMuonTrackMatchedTruthPDG([](const caf::SRSliceProxy* slc) -> int {
    int muonTrackIndex = kNuMIRelaxedMuonTrackIdx(slc);
    if(muonTrackIndex>=0){
      return slc->reco.pfp.at(muonTrackIndex).trk.truth.p.pdg;
    }
    else{
      return -99999999;
    }
  });
  const Var kNuMIRelaxedMuonTrackMatchedTruthIntID([](const caf::SRSliceProxy* slc) -> int {
    int muonTrackIndex = kNuMIRelaxedMuonTrackIdx(slc);
    if(muonTrackIndex>=0){
      return slc->reco.pfp.at(muonTrackIndex).trk.truth.p.interaction_id;
    }
    else{
      return -3;
    }
  });
  const Var kNuMIRelaxedMuonTrackMatchedTruthContained([](const caf::SRSliceProxy* slc) -> int {
    int muonTrackIndex = kNuMIRelaxedMuonTrackIdx(slc);
    if(muonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
      const auto& trk_truth = trk.truth.p;
      if(isnan(trk_truth.end.x) || isnan(trk_truth.end.y) || isnan(trk_truth.end.z)) return false;
      bool isTrkTruthContained = isContainedVol(trk_truth.end.x, trk_truth.end.y, trk_truth.end.z);
      if(isTrkTruthContained) return 1;
      else return 0;
    }
    else{
      return -1;
    }
  });
  const Var kNuMIRelaxedMuonTrackMatchType([](const caf::SRSliceProxy* slc) -> int {
    int muonTrackIndex = kNuMIRelaxedMuonTrackIdx(slc);
    if(muonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
      const auto& trk_truth = trk.truth.p;
      int intID = trk_truth.interaction_id;
      int trk_truth_pdg = trk_truth.pdg;
      bool isTrkTruthContained = ( isnan(trk_truth.end.x) || isnan(trk_truth.end.y) || isnan(trk_truth.end.z) ) ? false : isContainedVol(trk_truth.end.x, trk_truth.end.y, trk_truth.end.z);
      const int& endProc = trk_truth.end_process;
      bool isStopped = (endProc==2);
      bool isProtonInel = (endProc==7);
      bool isChargedPionInel = (endProc==9) || (endProc==10);

      // 0: Contained nu-mu
      if( intID!= -1 && abs(trk_truth_pdg)==13 && isTrkTruthContained ) return 0;
      // 1: Exiting nu-mu
      else if( intID!= -1 && abs(trk_truth_pdg)==13 && !isTrkTruthContained ) return 1;
      // 2: Cosmic muon
      else if( intID==-1 && abs(trk_truth_pdg)==13 ) return 2;
      // 3: stoppting proton
      else if( trk_truth_pdg==2212 && isTrkTruthContained && isStopped ) return 3;
      // 4: Inel. proton
      else if( trk_truth_pdg==2212 && isTrkTruthContained && isProtonInel ) return 4;
      // 5: Other proton
      else if( trk_truth_pdg==2212 && ( !(isTrkTruthContained && isStopped) && !(isTrkTruthContained && isProtonInel) ) ) return 5;
      // 6: Stopping charged pion
      else if( abs(trk_truth_pdg)==211 && isTrkTruthContained && isStopped ) return 6;
      // 7: Inel. charged pion
      else if( abs(trk_truth_pdg)==211 && isTrkTruthContained && isChargedPionInel ) return 7;
      // 8: Other charged pion
      else if( abs(trk_truth_pdg)==211 && ( !(isTrkTruthContained && isStopped) && !(isTrkTruthContained && isChargedPionInel) ) ) return 8;
      // 9: Other
      else return 9;
    }
    else{
      return -1;
    }
  });

  // - Cut
  const Cut kNuMIRelaxedMuonIsContained([](const caf::SRSliceProxy* slc) {
    int muonTrackIndex = kNuMIRelaxedMuonTrackIdx(slc);
    if(muonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
      const bool Contained = isContainedVol(trk.end.x,trk.end.y,trk.end.z);
      if( Contained ) return true;
      else return false;
    }
    else{
      return false;
    }
  });
  const Cut kNuMIRelaxedMuonNotIsoChronous([](const caf::SRSliceProxy* slc) {
    int muonTrackIndex = kNuMIRelaxedMuonTrackIdx(slc);
    if(muonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
      return ( fabs(trk.dir.x)>=0.1 );
    }
    else{
      return false;
    }
  });
  const Cut kNuMIRelaxedMuonLengthCut([](const caf::SRSliceProxy* slc) {
    if( kNuMIRelaxedMuonTrackLength(slc)>50. ) return true;
    else return false;
  });
  const Cut kNuMIRelaxedMuonSelection([](const caf::SRSliceProxy* slc) {
    if(!kNuMIVertexInFV(slc)) return false;
    if(!kNuMINotClearCosmic(slc)) return false;
    if(!kNuMICRLongestTrackDirYHardCut(slc)) return false;
    if(!kNuMIRelaxedMuonIsContained(slc)) return false;
    if(!kNuMIRelaxedMuonNotIsoChronous(slc)) return false;
    if(!kNuMIRelaxedMuonLengthCut(slc)) return false;
    return true;
  });
  const Var kNuMIIsRelaxedMuonSelection([](const caf::SRSliceProxy* slc) -> int {
    if( kNuMIRelaxedMuonSelection(slc) ) return 1;
    else return 0;
  });

  // Tag muon
  const Var kNuMITagMuonIdx([](const caf::SRSliceProxy* slc) -> int {
    double Longest = 0.;
    int PTrackInd = -1;
    for(unsigned int i_pfp=0; i_pfp<slc->reco.pfp.size(); i_pfp++){
      const auto& pfp = slc->reco.pfp.at(i_pfp);
      if( pfp.trackScore<0.5 ) continue;
      const auto& trk = pfp.trk;
      if(trk.bestplane == -1) continue;
      if(isnan(trk.start.x)) continue;
      const float Atslc = std::hypot(slc->vertex.x - trk.start.x,
                                     slc->vertex.y - trk.start.y,
                                     slc->vertex.z - trk.start.z);

      const bool isPrimCandidate = (Atslc < 10. && pfp.parent_is_primary);

      if(!isPrimCandidate) continue;
      const bool Contained = isContainedVol(trk.end.x,trk.end.y,trk.end.z);
      const float Chi2Proton = trk.chi2pid[2].chi2_proton;
      const float Chi2Muon = trk.chi2pid[2].chi2_muon;
      if ( (!Contained && trk.len > 100.) || (Contained && trk.calo[2].nhit>=5 && trk.len > 50. && Chi2Proton > 60. && Chi2Muon < 30.) ) {
        if ( trk.len <= Longest ) continue;
        Longest = trk.len;
        PTrackInd = i_pfp;
      }
    }
    return PTrackInd;
  });
  const Var kNuMITagMuonLength([](const caf::SRSliceProxy* slc) -> double {
    int muonTrackIndex = kNuMITagMuonIdx(slc);
    if(muonTrackIndex>=0){
      return slc->reco.pfp.at(muonTrackIndex).trk.len;
    }
    else{
      return -999.;
    }
  });

  // Proton
  // - Var
  const Var kNuMIRelaxedProtonTrackIdx([](const caf::SRSliceProxy* slc) -> int {
    if(slc->reco.pfp.size()==0){
      return -999.;
    }
    else{
      // Requiring good muon
      int muonTrackIndex = kNuMITagMuonIdx(slc);
      int protonIdx = -1;
      if(muonTrackIndex>=0){
        double protonMaxLength = -1.;
        // find the longest track
        for(unsigned int i_pfp=0; i_pfp<slc->reco.pfp.size(); i_pfp++){
          const auto& pfp = slc->reco.pfp.at(i_pfp);
          const auto& trk = pfp.trk;
          if(i_pfp==(unsigned int)muonTrackIndex) continue;
          if(trk.bestplane == -1) continue;
          if(isnan(trk.start.x)) continue;
          // First we calculate the distance of each track to the slice vertex.
          const float Atslc = std::hypot(slc->vertex.x - trk.start.x,
                                         slc->vertex.y - trk.start.y,
                                         slc->vertex.z - trk.start.z);
          // We require that the distance of the track from the slice is less than
          // 10 cm and that the parent of the track has been marked as the primary.
          const bool AtSlice = ( Atslc < 10.0 && pfp.parent_is_primary);
          const unsigned int idxPrim = (unsigned int)muonTrackIndex;
          const TVector3 muDir( slc->reco.pfp[idxPrim].trk.dir.x, slc->reco.pfp[idxPrim].trk.dir.y, slc->reco.pfp[idxPrim].trk.dir.z );
          const TVector3 pDir( trk.dir.x, trk.dir.y, trk.dir.z );
          const float angle = TMath::Cos(muDir.Angle(pDir));
          if( AtSlice && angle>=-0.9 && trk.len > protonMaxLength ){
            protonMaxLength = trk.len;
            protonIdx = i_pfp;
          }
        } // end loop pfp
      } // end if tag muon exist
      return protonIdx;
    }
  });
  const Var kNuMIRelaxedProtonTrackP([](const caf::SRSliceProxy* slc) -> double {
    int protonTrackIndex = kNuMIRelaxedProtonTrackIdx(slc);
    if(protonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
      bool Contained = isContainedVol(trk.end.x, trk.end.y, trk.end.z);
      if(Contained){
        if(isnan(trk.rangeP.p_proton)) return -999.;
        return trk.rangeP.p_proton;
      }
      else{
        if(isnan(trk.mcsP.fwdP_proton)) return -999.;
        return trk.mcsP.fwdP_proton;
      }

    }
    else{
      return -999.;
    }
  });
  const Var kNuMIRelaxedProtonTrackMatchedTruthPDG([](const caf::SRSliceProxy* slc) -> int {
    int protonTrackIndex = kNuMIRelaxedProtonTrackIdx(slc);
    if(protonTrackIndex>=0){
      return slc->reco.pfp.at(protonTrackIndex).trk.truth.p.pdg;
    }
    else{
      return -99999999;
    }
  });
  const Var kNuMIRelaxedProtonTrackMatchedTruthIntID([](const caf::SRSliceProxy* slc) -> int {
    int protonTrackIndex = kNuMIRelaxedProtonTrackIdx(slc);
    if(protonTrackIndex>=0){
      return slc->reco.pfp.at(protonTrackIndex).trk.truth.p.interaction_id;
    }
    else{
      return -3;
    }
  });
  const Var kNuMIRelaxedProtonTrackMatchedTruthContained([](const caf::SRSliceProxy* slc) -> int {
    int protonTrackIndex = kNuMIRelaxedProtonTrackIdx(slc);
    if(protonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
      const auto& trk_truth = trk.truth.p;
      if(isnan(trk_truth.end.x) || isnan(trk_truth.end.y) || isnan(trk_truth.end.z)) return false;
      bool isTrkTruthContained = isContainedVol(trk_truth.end.x, trk_truth.end.y, trk_truth.end.z);
      if(isTrkTruthContained) return 1;
      else return 0;
    }
    else{
      return -1;
    }
  });
  const Var kNuMIRelaxedProtonTrackMatchType([](const caf::SRSliceProxy* slc) -> int {
    int protonTrackIndex = kNuMIRelaxedProtonTrackIdx(slc);
    if(protonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
      const auto& trk_truth = trk.truth.p;
      int intID = trk_truth.interaction_id;
      int trk_truth_pdg = trk_truth.pdg;
      bool isTrkTruthContained = ( isnan(trk_truth.end.x) || isnan(trk_truth.end.y) || isnan(trk_truth.end.z) ) ? false : isContainedVol(trk_truth.end.x, trk_truth.end.y, trk_truth.end.z);
      const int& endProc = trk_truth.end_process;
      bool isStopped = (endProc==2);
      bool isProtonInel = (endProc==7);
      bool isChargedPionInel = (endProc==9) || (endProc==10);
      
      // 0: Contained nu-mu
      if( intID!= -1 && abs(trk_truth_pdg)==13 && isTrkTruthContained ) return 0;
      // 1: Exiting nu-mu 
      else if( intID!= -1 && abs(trk_truth_pdg)==13 && !isTrkTruthContained ) return 1;
      // 2: Cosmic muon
      else if( intID==-1 && abs(trk_truth_pdg)==13 ) return 2;
      // 3: stoppting proton
      else if( trk_truth_pdg==2212 && isTrkTruthContained && isStopped ) return 3;
      // 4: Inel. proton
      else if( trk_truth_pdg==2212 && isTrkTruthContained && isProtonInel ) return 4;
      // 5: Other proton
      else if( trk_truth_pdg==2212 && ( !(isTrkTruthContained && isStopped) && !(isTrkTruthContained && isProtonInel) ) ) return 5;
      // 6: Stopping charged pion
      else if( abs(trk_truth_pdg)==211 && isTrkTruthContained && isStopped ) return 6;
      // 7: Inel. charged pion
      else if( abs(trk_truth_pdg)==211 && isTrkTruthContained && isChargedPionInel ) return 7;
      // 8: Other charged pion
      else if( abs(trk_truth_pdg)==211 && ( !(isTrkTruthContained && isStopped) && !(isTrkTruthContained && isChargedPionInel) ) ) return 8;
      // 9: Other 
      else return 9;
    }
    else{
      return -1;
    }
  });
  // - Cut
  const Cut kNuMIRelaxedProtonIsContained([](const caf::SRSliceProxy* slc) {
    int protonTrackIndex = kNuMIRelaxedProtonTrackIdx(slc);
    if(protonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
      const bool Contained = isContainedVol(trk.end.x,trk.end.y,trk.end.z);
      if( Contained ) return true;
      else return false;
    }
    else{
      return false;
    }
  });
  const Cut kNuMIRelaxedProtonNotIsoChronous([](const caf::SRSliceProxy* slc) {
    int protonTrackIndex = kNuMIRelaxedProtonTrackIdx(slc);
    if(protonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
      return ( fabs(trk.dir.x)>=0.1 );
    }
    else{
      return false;
    }
  });
  const Cut kNuMIRelaxedProtonSelection([](const caf::SRSliceProxy* slc) {
    if(!kNuMIVertexInFV(slc)) return false;
    if(!kNuMINotClearCosmic(slc)) return false;
    if(!kNuMICRLongestTrackDirYHardCut(slc)) return false;
    if(!kNuMIRelaxedProtonIsContained(slc)) return false;
    if(!kNuMIRelaxedProtonNotIsoChronous(slc)) return false;
    return true;
  });
  const Var kNuMIIsRelaxedProtonSelection([](const caf::SRSliceProxy* slc) -> int {
    if( kNuMIRelaxedProtonSelection(slc) ) return 1;
    else return 0;
  });

  // ChargedPion
  // - Var
  const Var kNuMIRelaxedChargedPionTrackIdx([](const caf::SRSliceProxy* slc) -> int {
    if(slc->reco.pfp.size()==0){
     return -999.;
    }
    else{
      // Requiring good muon
      int muonTrackIndex = kNuMITagMuonIdx(slc);
      int PTrackInd(-1);
      if(muonTrackIndex>=0){
        double Longest(0);
        for(unsigned int i_pfp=0; i_pfp<slc->reco.pfp.size(); i_pfp++){
          const auto& pfp = slc->reco.pfp.at(i_pfp);
          const auto& trk = pfp.trk;
          if(i_pfp==(unsigned int)muonTrackIndex) continue;
          if(trk.bestplane==-1) continue;
          if(isnan(trk.start.x)) continue;
          // First we calculate the distance of each track to the slice vertex.
          const float Atslc = std::hypot(slc->vertex.x - trk.start.x,
                                         slc->vertex.y - trk.start.y,
                                         slc->vertex.z - trk.start.z);
          // We require that the distance of the track from the slice is less than
          // 10 cm and that the parent of the track has been marked as the primary.
          const bool AtSlice = ( Atslc < 10.0 && pfp.parent_is_primary);
          // pid from collection
          const float Chi2Proton = trk.chi2pid[2].chi2_proton;
          const float Chi2Muon = trk.chi2pid[2].chi2_muon;
          const bool Contained = isContainedVol(trk.end.x, trk.end.y, trk.end.z);
          const bool MaybeExiting = ( !Contained && trk.len > 100);
          const bool MaybeContained = ( Contained && trk.calo[2].nhit!=0 && Chi2Proton > 60 && Chi2Muon < 30 && trk.len > 20 );
          if ( AtSlice && ( MaybeExiting || MaybeContained ) && trk.len > Longest )
          {
            Longest = trk.len;
            PTrackInd = i_pfp;
          }
        } // END pfp loop
      } // END if muon track exist
      return PTrackInd;
    } // END if slice has primary tracks
  });
  const Var kNuMIRelaxedChargedPionTrackLength([](const caf::SRSliceProxy* slc) -> double {
    int cpionTrackIndex = kNuMIRelaxedChargedPionTrackIdx(slc);
    if(cpionTrackIndex>=0){
      return slc->reco.pfp.at(cpionTrackIndex).trk.len;
    }
    else{
      return -999.;
    }
  });
  const Var kNuMIRelaxedChargedPionTrackMatchedTruthPDG([](const caf::SRSliceProxy* slc) -> int {
    int cpionTrackIndex = kNuMIRelaxedChargedPionTrackIdx(slc);
    if(cpionTrackIndex>=0){
      return slc->reco.pfp.at(cpionTrackIndex).trk.truth.p.pdg;
    }
    else{
      return -99999999;
    }
  });
  const Var kNuMIRelaxedChargedPionTrackMatchedTruthIntID([](const caf::SRSliceProxy* slc) -> int {
    int cpionTrackIndex = kNuMIRelaxedChargedPionTrackIdx(slc);
    if(cpionTrackIndex>=0){
      return slc->reco.pfp.at(cpionTrackIndex).trk.truth.p.interaction_id;
    }
    else{
      return -3;
    }
  });
  const Var kNuMIRelaxedChargedPionTrackMatchedTruthContained([](const caf::SRSliceProxy* slc) -> int {
    int cpionTrackIndex = kNuMIRelaxedChargedPionTrackIdx(slc);
    if(cpionTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(cpionTrackIndex).trk;
      const auto& trk_truth = trk.truth.p;
      if(isnan(trk_truth.end.x) || isnan(trk_truth.end.y) || isnan(trk_truth.end.z)) return false;
      bool isTrkTruthContained = isContainedVol(trk_truth.end.x, trk_truth.end.y, trk_truth.end.z);
      if(isTrkTruthContained) return 1;
      else return 0;
    }
    else{
      return -1;
    }
  });
  const Var kNuMIRelaxedChargedPionTrackMatchType([](const caf::SRSliceProxy* slc) -> int {
    int cpionTrackIndex = kNuMIRelaxedChargedPionTrackIdx(slc);
    if(cpionTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(cpionTrackIndex).trk;
      const auto& trk_truth = trk.truth.p;
      int intID = trk_truth.interaction_id;
      int trk_truth_pdg = trk_truth.pdg;
      bool isTrkTruthContained = ( isnan(trk_truth.end.x) || isnan(trk_truth.end.y) || isnan(trk_truth.end.z) ) ? false : isContainedVol(trk_truth.end.x, trk_truth.end.y, trk_truth.end.z);
      const int& endProc = trk_truth.end_process;
      bool isStopped = (endProc==2);
      bool isProtonInel = (endProc==7);
      bool isChargedPionInel = (endProc==9) || (endProc==10);
      
      // 0: Contained nu-mu
      if( intID!= -1 && abs(trk_truth_pdg)==13 && isTrkTruthContained ) return 0;
      // 1: Exiting nu-mu 
      else if( intID!= -1 && abs(trk_truth_pdg)==13 && !isTrkTruthContained ) return 1;
      // 2: Cosmic muon
      else if( intID==-1 && abs(trk_truth_pdg)==13 ) return 2;
      // 3: stoppting proton
      else if( trk_truth_pdg==2212 && isTrkTruthContained && isStopped ) return 3;
      // 4: Inel. proton
      else if( trk_truth_pdg==2212 && isTrkTruthContained && isProtonInel ) return 4;
      // 5: Other proton
      else if( trk_truth_pdg==2212 && ( !(isTrkTruthContained && isStopped) && !(isTrkTruthContained && isProtonInel) ) ) return 5;
      // 6: Stopping charged pion
      else if( abs(trk_truth_pdg)==211 && isTrkTruthContained && isStopped ) return 6;
      // 7: Inel. charged pion
      else if( abs(trk_truth_pdg)==211 && isTrkTruthContained && isChargedPionInel ) return 7;
      // 8: Other charged pion
      else if( abs(trk_truth_pdg)==211 && ( !(isTrkTruthContained && isStopped) && !(isTrkTruthContained && isChargedPionInel) ) ) return 8;
      // 9: Other 
      else return 9;
    }
    else{
      return -1;
    }
  });
  // - Cut
  const Cut kNuMIRelaxedChargedPionIsContained([](const caf::SRSliceProxy* slc) {
    int cpionTrackIndex = kNuMIRelaxedChargedPionTrackIdx(slc);
    if(cpionTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(cpionTrackIndex).trk;
      const bool Contained = isContainedVol(trk.end.x,trk.end.y,trk.end.z);
      if( Contained ) return true;
      else return false;
    }
    else{
      return false;
    }
  });
  const Cut kNuMIRelaxedChargedPionNotIsoChronous([](const caf::SRSliceProxy* slc) {
    int cpionTrackIndex = kNuMIRelaxedChargedPionTrackIdx(slc);
    if(cpionTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(cpionTrackIndex).trk;
      return ( fabs(trk.dir.x)>=0.1 );
    }
    else{
      return false;
    }
  });
  const Cut kNuMIRelaxedChargedPionSelection([](const caf::SRSliceProxy* slc) {
    if(!kNuMIVertexInFV(slc)) return false;
    if(!kNuMINotClearCosmic(slc)) return false;
    if(!kNuMICRLongestTrackDirYHardCut(slc)) return false;
    if(!kNuMIRelaxedChargedPionIsContained(slc)) return false;
    if(!kNuMIRelaxedChargedPionNotIsoChronous(slc)) return false;
    return true;
  });
  const Var kNuMIIsRelaxedChargedPionSelection([](const caf::SRSliceProxy* slc) -> int {
    if( kNuMIRelaxedChargedPionSelection(slc) ) return 1;
    else return 0;
  });

  // For skimming
  const Cut kNuMIRelaxedPIDCut([](const caf::SRSliceProxy* slc) {
    bool ms = kNuMIRelaxedMuonSelection(slc);
    bool ps = kNuMIRelaxedProtonSelection(slc);
    bool cps = kNuMIRelaxedChargedPionSelection(slc);
    if( ms || ps || cps ) return true;
    else return false;
  });



} // end namespace ana
