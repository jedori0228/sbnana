#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_TwoTrack_Variables.h"

using namespace std;
using namespace ana;

namespace ICARUSNumuXsec{

namespace TwoTrack{

  // Muon
  const Var MuonTrackIndex([](const caf::SRSliceProxy* slc) -> double {
    vector<double> primTrackIndices = ICARUSNumuXsec::PrimaryTrackIndices(slc);
    if(primTrackIndices.size()==0){
      return -999.;
    }
    else{
      float Longest(0);
      int PTrackInd(-1);
      for(const auto& trkIdx: primTrackIndices){
        const auto& pfp = slc->reco.pfp.at(trkIdx);
        const auto& trk = pfp.trk;

        if(trk.bestplane == -1) continue;
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

        const bool Contained = fv_track.isContained(trk.end.x, trk.end.y, trk.end.z);

        const bool MaybeMuonExiting = ( !Contained && trk.len > 100);
        const bool MaybeMuonContained = ( Contained && trk.calo[2].nhit>=5 && Chi2Proton > 60 && Chi2Muon < 30 && trk.len > 50 );
        if ( AtSlice && ( MaybeMuonExiting || MaybeMuonContained ) && trk.len > Longest )
        {
          Longest = trk.len;
          PTrackInd = trkIdx;
        }

      }

      return PTrackInd;

    }
  });
  const Var MuonTrackLength([](const caf::SRSliceProxy* slc) -> double {
    int muonTrackIndex = MuonTrackIndex(slc);
    if(muonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
      return trk.len;
    }
    else{
      return -999.;
    }

  });
  const Var MuonTrackP([](const caf::SRSliceProxy* slc) -> double {
    int muonTrackIndex = MuonTrackIndex(slc);
    if(muonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
      bool Contained = fv_track.isContained(trk.end.x, trk.end.y, trk.end.z);
      if(Contained){
        if(isnan(trk.rangeP.p_muon)) return -999.;
        return trk.rangeP.p_muon;
      }
      else{
        if(isnan(trk.mcsP.fwdP_muon)) return -999.;
        return trk.mcsP.fwdP_muon;
      }
    }
    else{
      return -999.;
    }

  });
  const Var MuonTrackDirX([](const caf::SRSliceProxy* slc) -> double {
    int muonTrackIndex = MuonTrackIndex(slc);
    if(muonTrackIndex>=0){
      return slc->reco.pfp.at(muonTrackIndex).trk.dir.x;
    }
    else{
      return -999.;
    }
  });
  const Var MuonTrackDirY([](const caf::SRSliceProxy* slc) -> double {
    int muonTrackIndex = MuonTrackIndex(slc);
    if(muonTrackIndex>=0){
      return slc->reco.pfp.at(muonTrackIndex).trk.dir.y;
    }
    else{
      return -999.;
    }
  });
  const Var MuonTrackDirZ([](const caf::SRSliceProxy* slc) -> double {
    int muonTrackIndex = MuonTrackIndex(slc);
    if(muonTrackIndex>=0){
      return slc->reco.pfp.at(muonTrackIndex).trk.dir.z;
    }
    else{
      return -999.;
    }
  });
  const Var MuonTrackNuMICosineTheta([](const caf::SRSliceProxy* slc) -> double {
    int muonTrackIndex = MuonTrackIndex(slc);
    if(muonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
      TVector3 v3_trk(trk.dir.x, trk.dir.y, trk.dir.z);
      double angleNuMI = v3_trk.Angle(NuDirection_NuMI);
      return TMath::Cos(angleNuMI);
    }
    else{
      return -999.;
    }
  });
  // truth match
  const Var MuonTrackNuMIToVtxCosineTheta([](const caf::SRSliceProxy* slc) -> double {
    int muonTrackIndex = MuonTrackIndex(slc);
    if(muonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
      const auto& vtx = slc->vertex;

      TVector3 v3_trk(trk.dir.x, trk.dir.y, trk.dir.z);
      TVector3 v3_vtx(vtx.x, vtx.y, vtx.z);
      static const TVector3 dFromNuMI(315.120380, 33.644912, 733.632532);

      TVector3 v3_numi_to_vtx = dFromNuMI+v3_vtx;

      double angleNuMI = v3_trk.Angle(v3_numi_to_vtx);
      return TMath::Cos(angleNuMI);
    }
    else{
      return -999.;
    }
  });
  // truth match
  const Var MuonTrackTruthLength([](const caf::SRSliceProxy* slc) -> double {
    int muonTrackIndex = MuonTrackIndex(slc);
    if(muonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
      if( isnan(trk.truth.p.length) ) return -999.;
      return trk.truth.p.length;
    }
    else{
      return -999.;
    }
  });
  const Var MuonTrackTruthP([](const caf::SRSliceProxy* slc) -> double {
    int muonTrackIndex = MuonTrackIndex(slc);
    if(muonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
      const auto& truth_p = trk.truth.p;
      if( isnan(truth_p.genp.x) ) return -999.;
      TVector3 v3_gen_p(truth_p.genp.x, truth_p.genp.y, truth_p.genp.z);
      return v3_gen_p.Mag();
    }
    else{
      return -999.;
    }
  });
  const Var MuonTrackTruthNuMICosineTheta([](const caf::SRSliceProxy* slc) -> double {
    int muonTrackIndex = MuonTrackIndex(slc);
    if(muonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
      const auto& truth_p = trk.truth.p;
      if( isnan(truth_p.genp.x) ) return -999.;
      TVector3 v3_gen_p(truth_p.genp.x, truth_p.genp.y, truth_p.genp.z);
      double angleNuMI = v3_gen_p.Angle(NuDirection_NuMI);
      return TMath::Cos(angleNuMI);
    }
    else{
      return -999.;
    }
  });

  // Proton
  const Var ProtonTrackIndex([](const caf::SRSliceProxy* slc) -> double {
    vector<double> primTrackIndices = ICARUSNumuXsec::PrimaryTrackIndices(slc);
    if(primTrackIndices.size()==0){
      return -999.;
    }
    else{
      float Longest(0);
      int PTrackInd(-1);
      int muonTrackIndex = MuonTrackIndex(slc);
      for(const auto& trkIdx: primTrackIndices){
        const auto& pfp = slc->reco.pfp.at(trkIdx);
        const auto& trk = pfp.trk;

        if(trkIdx==muonTrackIndex) continue;
        if(trk.bestplane == -1) continue;

        // First we calculate the distance of each track to the slice vertex.
        if(isnan(trk.start.x)) continue;
        const float Atslc = std::hypot(slc->vertex.x - trk.start.x,
                                       slc->vertex.y - trk.start.y,
                                       slc->vertex.z - trk.start.z);

        // We require that the distance of the track from the slice is less than
        // 10 cm and that the parent of the track has been marked as the primary.
        const bool AtSlice = ( Atslc < 10.0 && pfp.parent_is_primary);
        // pid from collection only
        const float Chi2Proton = trk.chi2pid[2].chi2_proton;
        const float Chi2Muon = trk.chi2pid[2].chi2_muon;

        const bool Contained = fv_track.isContained(trk.end.x, trk.end.y, trk.end.z);

        float angle = -5.0;
        if ( muonTrackIndex >= 0 ) {
          const unsigned int idxPrim = (unsigned int)muonTrackIndex;
          TVector3 muDir( slc->reco.pfp[idxPrim].trk.dir.x, slc->reco.pfp[idxPrim].trk.dir.y, slc->reco.pfp[idxPrim].trk.dir.z );
          TVector3 pDir( slc->reco.pfp[trkIdx].trk.dir.x, slc->reco.pfp[trkIdx].trk.dir.y, slc->reco.pfp[trkIdx].trk.dir.z );
          angle = TMath::Cos(muDir.Angle(pDir));
        }

        bool PassPCut = false;
        if( !isnan(trk.rangeP.p_proton) ){
          PassPCut = trk.rangeP.p_proton>0.4;
        }

        if ( AtSlice && Contained && trk.calo[2].nhit>=5 && Chi2Proton <= 90 && Chi2Muon >= 30 && angle >= -0.9 && trk.len > Longest && PassPCut ) {
        //if ( AtSlice && Contained && trk.calo[2].nhit!=0 && Chi2Proton <= 100 && angle >= -0.9 && trk.len > Longest ) {
          Longest = trk.len;
          PTrackInd = trkIdx;
        }

      }

      return PTrackInd;

    }
  });

  const Var ProtonTrackLength([](const caf::SRSliceProxy* slc) -> double {
    int protonTrackIndex = ProtonTrackIndex(slc);
    if(protonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
      return trk.len;
    }
    else{
      return -999.;
    }

  });
  const Var ProtonTrackP([](const caf::SRSliceProxy* slc) -> double {
    int protonTrackIndex = ProtonTrackIndex(slc);
    if(protonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
      bool Contained = fv_track.isContained(trk.end.x, trk.end.y, trk.end.z);
      if(Contained){
        if(isnan(trk.rangeP.p_proton)) return -999.;
        return trk.rangeP.p_proton;
      }
      else return -999.;
    }
    else{
      return -999.;
    }

  });
  const Var ProtonTrackDirX([](const caf::SRSliceProxy* slc) -> double {
    int protonTrackIndex = ProtonTrackIndex(slc);
    if(protonTrackIndex>=0){
        return slc->reco.pfp.at(protonTrackIndex).trk.dir.x;
    }
    else{
      return -999.;
    }
  });
  const Var ProtonTrackDirY([](const caf::SRSliceProxy* slc) -> double {
    int protonTrackIndex = ProtonTrackIndex(slc);
    if(protonTrackIndex>=0){
        return slc->reco.pfp.at(protonTrackIndex).trk.dir.y;
    }
    else{
      return -999.;
    }
  });
  const Var ProtonTrackDirZ([](const caf::SRSliceProxy* slc) -> double {
    int protonTrackIndex = ProtonTrackIndex(slc);
    if(protonTrackIndex>=0){
        return slc->reco.pfp.at(protonTrackIndex).trk.dir.z;
    }
    else{
      return -999.;
    }
  });
  const Var ProtonTrackNuMICosineTheta([](const caf::SRSliceProxy* slc) -> double {
    int protonTrackIndex = ProtonTrackIndex(slc);
    if(protonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
      TVector3 v3_trk(trk.dir.x, trk.dir.y, trk.dir.z);
      double angleNuMI = v3_trk.Angle(NuDirection_NuMI);
      return TMath::Cos(angleNuMI);
    }
    else{
      return -999.;
    }
  });
  const Var ProtonTrackNuMIToVtxCosineTheta([](const caf::SRSliceProxy* slc) -> double {
    int protonTrackIndex = ProtonTrackIndex(slc);
    if(protonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
      const auto& vtx = slc->vertex;
      
      TVector3 v3_trk(trk.dir.x, trk.dir.y, trk.dir.z);
      TVector3 v3_vtx(vtx.x, vtx.y, vtx.z);
      static const TVector3 dFromNuMI(315.120380, 33.644912, 733.632532);
      
      TVector3 v3_numi_to_vtx = dFromNuMI+v3_vtx;
      
      double angleNuMI = v3_trk.Angle(v3_numi_to_vtx);
      return TMath::Cos(angleNuMI);
    }
    else{
      return -999.;
    }
  });
  const Var ProtonTrackNHitsCollection([](const caf::SRSliceProxy* slc) -> double {
    int protonTrackIndex = ProtonTrackIndex(slc);
    if(protonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
      return trk.calo[2].nhit;
    }
    else{
      return -999.;
    }
  });


  const Var ProtonTrackChi2MuonCollection([](const caf::SRSliceProxy* slc) -> double {
    int protonTrackIndex = ProtonTrackIndex(slc);
    if(protonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
      if(trk.calo[2].nhit==0) return -999.;
      else return trk.chi2pid[2].chi2_muon;
    }
    else{
      return -999.;
    }
  });
  const Var ProtonTrackChi2ProtonCollection([](const caf::SRSliceProxy* slc) -> double {
    int protonTrackIndex = ProtonTrackIndex(slc);
    if(protonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
      if(trk.calo[2].nhit==0) return -999.;
      else return trk.chi2pid[2].chi2_proton;
    }
    else{
      return -999.;
    }
  });
  // truth match
  const Var ProtonTrackTruthLength([](const caf::SRSliceProxy* slc) -> double {
    int protonTrackIndex = ProtonTrackIndex(slc);
    if(protonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
      if( isnan(trk.truth.p.length) ) return -999.;
      return trk.truth.p.length;
    }
    else{
      return -999.;
    }
  }); 
  const Var ProtonTrackTruthP([](const caf::SRSliceProxy* slc) -> double {
    int protonTrackIndex = ProtonTrackIndex(slc);
    if(protonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
      const auto& truth_p = trk.truth.p; 
      if( isnan(truth_p.genp.x) ) return -999.;
      TVector3 v3_gen_p(truth_p.genp.x, truth_p.genp.y, truth_p.genp.z);
      return v3_gen_p.Mag();
    }
    else{
      return -999.;
    }
  });
  const Var ProtonTrackTruthNuMICosineTheta([](const caf::SRSliceProxy* slc) -> double {
    int protonTrackIndex = ProtonTrackIndex(slc);
    if(protonTrackIndex>=0){
      const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
      const auto& truth_p = trk.truth.p;
      if( isnan(truth_p.genp.x) ) return -999.;
      TVector3 v3_gen_p(truth_p.genp.x, truth_p.genp.y, truth_p.genp.z);
      double angleNuMI = v3_gen_p.Angle(NuDirection_NuMI);
      return TMath::Cos(angleNuMI);
    }
    else{
      return -999.;
    }
  });

  // Muon+Proton
  const Var MuonProtonCosineTheta([](const caf::SRSliceProxy* slc) -> double {
    int muonTrackIndex = MuonTrackIndex(slc);
    int protonTrackIndex = ProtonTrackIndex(slc);
    if(muonTrackIndex>=0 && protonTrackIndex>=0){
      const auto& trk_mu = slc->reco.pfp.at(muonTrackIndex).trk;
      const auto& trk_pro = slc->reco.pfp.at(protonTrackIndex).trk;

      TVector3 v3_trk_mu(trk_mu.dir.x, trk_mu.dir.y, trk_mu.dir.z);
      TVector3 v3_trk_pro(trk_pro.dir.x, trk_pro.dir.y, trk_pro.dir.z);

      double angle = v3_trk_mu.Angle(v3_trk_pro);
      return TMath::Cos(angle);

    }
    else{
      return -999.;
    }
  });

  namespace TKI{
    const Var deltaPT([](const caf::SRSliceProxy* slc) -> double {
      int muonTrackIndex = MuonTrackIndex(slc);
      int protonTrackIndex = ProtonTrackIndex(slc);
      if(muonTrackIndex>=0 && protonTrackIndex>=0){
        const auto& trk_mu = slc->reco.pfp.at(muonTrackIndex).trk;
        const auto& trk_pro = slc->reco.pfp.at(protonTrackIndex).trk;
        const auto& vtx = slc->vertex;

        TVector3 vec_mu(trk_mu.dir.x, trk_mu.dir.y, trk_mu.dir.z);
        TVector3 vec_pro(trk_pro.dir.x, trk_pro.dir.y, trk_pro.dir.z);

        double p_mu = MuonTrackP(slc);
        double p_pro = ProtonTrackP(slc);

        vec_mu *= p_mu;
        vec_pro *= p_pro;

        TVector3 vec_vtx_icarus(vtx.x, vtx.y, vtx.z);
        static const TVector3 vec_NuMI_to_ICARUS(315.120380, 33.644912, 733.632532);

        TVector3 unit_numi_to_vtx = (vec_NuMI_to_ICARUS+vec_vtx_icarus).Unit();

        TVector3 pt_mu = vec_mu - (vec_mu.Dot(unit_numi_to_vtx))*unit_numi_to_vtx ;
        TVector3 pt_pro = vec_pro - (vec_pro.Dot(unit_numi_to_vtx))*unit_numi_to_vtx;

        TVector3 vec_deltaPT = pt_mu+pt_pro;

        return vec_deltaPT.Mag();
      }
      else{
        return -9999999.;
      }
    });
    const Var deltaPTx([](const caf::SRSliceProxy* slc) -> double {
      int muonTrackIndex = MuonTrackIndex(slc);
      int protonTrackIndex = ProtonTrackIndex(slc);
      if(muonTrackIndex>=0 && protonTrackIndex>=0){
        const auto& trk_mu = slc->reco.pfp.at(muonTrackIndex).trk;
        const auto& trk_pro = slc->reco.pfp.at(protonTrackIndex).trk;
        const auto& vtx = slc->vertex;

        TVector3 vec_mu(trk_mu.dir.x, trk_mu.dir.y, trk_mu.dir.z);
        TVector3 vec_pro(trk_pro.dir.x, trk_pro.dir.y, trk_pro.dir.z);

        double p_mu = MuonTrackP(slc);
        double p_pro = ProtonTrackP(slc);

        vec_mu *= p_mu;
        vec_pro *= p_pro;

        TVector3 vec_vtx_icarus(vtx.x, vtx.y, vtx.z);
        static const TVector3 vec_NuMI_to_ICARUS(315.120380, 33.644912, 733.632532);

        TVector3 unit_numi_to_vtx = (vec_NuMI_to_ICARUS+vec_vtx_icarus).Unit();

        TVector3 pt_mu = vec_mu - (vec_mu.Dot(unit_numi_to_vtx))*unit_numi_to_vtx ;
        TVector3 pt_pro = vec_pro - (vec_pro.Dot(unit_numi_to_vtx))*unit_numi_to_vtx;

        TVector3 vec_deltaPT = pt_mu+pt_pro;

        double deltaPT_x = ( unit_numi_to_vtx.Cross(pt_mu.Unit()) ).Dot(vec_deltaPT);

        return deltaPT_x;
      }
      else{
        return -9999999.;
      }
    });
    const Var deltaPTy([](const caf::SRSliceProxy* slc) -> double {
      int muonTrackIndex = MuonTrackIndex(slc);
      int protonTrackIndex = ProtonTrackIndex(slc);
      if(muonTrackIndex>=0 && protonTrackIndex>=0){
        const auto& trk_mu = slc->reco.pfp.at(muonTrackIndex).trk;
        const auto& trk_pro = slc->reco.pfp.at(protonTrackIndex).trk;
        const auto& vtx = slc->vertex;

        TVector3 vec_mu(trk_mu.dir.x, trk_mu.dir.y, trk_mu.dir.z);
        TVector3 vec_pro(trk_pro.dir.x, trk_pro.dir.y, trk_pro.dir.z);

        double p_mu = MuonTrackP(slc);
        double p_pro = ProtonTrackP(slc);

        vec_mu *= p_mu;
        vec_pro *= p_pro;

        TVector3 vec_vtx_icarus(vtx.x, vtx.y, vtx.z);
        static const TVector3 vec_NuMI_to_ICARUS(315.120380, 33.644912, 733.632532);

        TVector3 unit_numi_to_vtx = (vec_NuMI_to_ICARUS+vec_vtx_icarus).Unit();

        TVector3 pt_mu = vec_mu - (vec_mu.Dot(unit_numi_to_vtx))*unit_numi_to_vtx ;
        TVector3 pt_pro = vec_pro - (vec_pro.Dot(unit_numi_to_vtx))*unit_numi_to_vtx;

        TVector3 vec_deltaPT = pt_mu+pt_pro;

        double deltaPT_y = -1.*(pt_mu.Unit().Dot(vec_deltaPT));

        return deltaPT_y;
      }
      else{
        return -9999999.;
      }
    });
  } // END namespace TKI

  // All non-muon track indices, i.e. hadron candiates
  const MultiVar NonMuonTrackIndecies([](const caf::SRSliceProxy* slc) -> vector<double> {

    vector<double> rets;

    vector<double> primTrackIndices = ICARUSNumuXsec::PrimaryTrackIndices(slc);
    if(primTrackIndices.size()==0){
      return rets;
    }
    else{

      // Requiring good muon
      int muonTrackIndex = MuonTrackIndex(slc);
      if(muonTrackIndex>=0){

        for(const auto& trkIdx: primTrackIndices){
          const auto& pfp = slc->reco.pfp.at(trkIdx);
          const auto& trk = pfp.trk;

          if(trkIdx==(unsigned int)muonTrackIndex) continue;
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

          const bool Contained = fv_track.isContained(trk.end.x, trk.end.y, trk.end.z);

          const bool MaybeExiting = ( !Contained && trk.len > 50);
          const bool MaybeContained = ( Contained && trk.calo[2].nhit>=5 && Chi2Proton > 60 && Chi2Muon < 30 && trk.len > 20 );

          bool isMuonCand = AtSlice && ( MaybeExiting || MaybeContained );

          if(!isMuonCand) rets.push_back(trkIdx);

        } // END pfp loop

      } // END If muon exist

      return rets;

    }

  });

  // Charged pion
  const Var ChargedPionTrackIndex([](const caf::SRSliceProxy* slc) -> double {
    vector<double> primTrackIndices = ICARUSNumuXsec::PrimaryTrackIndices(slc);
    if(primTrackIndices.size()==0){
      return -999.;
    }
    else{

      // Requiring good muon
      int muonTrackIndex = MuonTrackIndex(slc);
      int PTrackInd(-1);
      if(muonTrackIndex>=0){

        int protonTrackIndex = ProtonTrackIndex(slc);

        float Longest(0);
        for(const auto& trkIdx: primTrackIndices){
          const auto& pfp = slc->reco.pfp.at(trkIdx);
          const auto& trk = pfp.trk;

          if(trkIdx==(unsigned int)muonTrackIndex) continue;
          if(trkIdx==(unsigned int)protonTrackIndex) continue;
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

          const bool Contained = fv_track.isContained(trk.end.x, trk.end.y, trk.end.z);

          const bool MaybeExiting = ( !Contained && trk.len > 50);
          const bool MaybeContained = ( Contained && trk.calo[2].nhit>=5 && Chi2Proton > 60 && Chi2Muon < 30 && trk.len > 20 );

          if ( AtSlice && ( MaybeExiting || MaybeContained ) && trk.len > Longest )
          {
            Longest = trk.len;
            PTrackInd = trkIdx;
          }
        } // END pfp loop

      } // END if muon track exist

      return PTrackInd;

    } // END if slice has primary tracks
  });
  // Charged pion, stopped
  const Var StoppedChargedPionTrackIndex([](const caf::SRSliceProxy* slc) -> double {
    if(slc->reco.pfp.size()==0){
      return -999.;
    }
    else{

      // Requiring good muon
      int muonTrackIndex = MuonTrackIndex(slc);
      int PTrackInd(-1);
      if(muonTrackIndex>=0){

        int protonTrackIndex = ProtonTrackIndex(slc);

        float Longest(0);
        for(unsigned int i_pfp=0; i_pfp<slc->reco.pfp.size(); i_pfp++){
          const auto& pfp = slc->reco.pfp.at(i_pfp);
          const auto& trk = pfp.trk;

          if(i_pfp==(unsigned int)muonTrackIndex) continue;
          if(i_pfp==(unsigned int)protonTrackIndex) continue;
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

          const bool Contained = fv_track.isContained(trk.end.x, trk.end.y, trk.end.z);

          const bool MaybeExiting = ( !Contained && trk.len > 50);
          //TODO
          const bool MaybeContained = ( Contained && trk.calo[2].nhit>=5 && Chi2Proton > 60 && Chi2Muon < 30 && trk.len > 0. );

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
  const Var StoppedChargedPionTrackLength([](const caf::SRSliceProxy* slc) -> double {
    int stoppedCPionIndex = StoppedChargedPionTrackIndex(slc);
    if(stoppedCPionIndex>=0){
      const auto& trk = slc->reco.pfp.at(stoppedCPionIndex).trk;
      return trk.len;
    }
    else{
      return -999.;
    }
  });
  const Var StoppedChargedPionTrackP([](const caf::SRSliceProxy* slc) -> double {
    int stoppedCPionIndex = StoppedChargedPionTrackIndex(slc);
    if(stoppedCPionIndex>=0){
      const auto& trk = slc->reco.pfp.at(stoppedCPionIndex).trk;
      bool Contained = fv_track.isContained(trk.end.x, trk.end.y, trk.end.z);
      if(Contained){
        if(isnan(trk.rangeP.p_pion)) return -999.;
        return trk.rangeP.p_pion;
      }
      else{
        if(isnan(trk.mcsP.fwdP_pion)) return -999.;
        return trk.mcsP.fwdP_pion;
      }
    }
    else{
      return -999.;
    }
  });
  const Var StoppedChargedPionTrackNuMICosineTheta([](const caf::SRSliceProxy* slc) -> double {
    int stoppedCPionIndex = StoppedChargedPionTrackIndex(slc);
    if(stoppedCPionIndex>=0){
      const auto& trk = slc->reco.pfp.at(stoppedCPionIndex).trk;
      TVector3 v3_trk(trk.dir.x, trk.dir.y, trk.dir.z);
      double angleNuMI = v3_trk.Angle(NuDirection_NuMI);
      return TMath::Cos(angleNuMI);
    }
    else{
      return -999.;
    }
  });
  const Var StoppedChargedPionTrackNuMIToVtxCosineTheta([](const caf::SRSliceProxy* slc) -> double {
    int stoppedCPionIndex = StoppedChargedPionTrackIndex(slc);
    if(stoppedCPionIndex>=0){
      const auto& trk = slc->reco.pfp.at(stoppedCPionIndex).trk;
      const auto& vtx = slc->vertex;

      TVector3 v3_trk(trk.dir.x, trk.dir.y, trk.dir.z);
      TVector3 v3_vtx(vtx.x, vtx.y, vtx.z);
      static const TVector3 dFromNuMI(315.120380, 33.644912, 733.632532);

      TVector3 v3_numi_to_vtx = dFromNuMI+v3_vtx;

      double angleNuMI = v3_trk.Angle(v3_numi_to_vtx);
      return TMath::Cos(angleNuMI);
    }
    else{
      return -999.;
    }
  });
  const Var StoppedChargedPionTrackChi2MuonCollection([](const caf::SRSliceProxy* slc) -> double {
    int stoppedCPionIndex = StoppedChargedPionTrackIndex(slc);
    if(stoppedCPionIndex>=0){
      const auto& trk = slc->reco.pfp.at(stoppedCPionIndex).trk;
      if(trk.calo[2].nhit==0) return -999.;
      else return trk.chi2pid[2].chi2_muon;
    }
    else{
      return -999.;
    }
  });
  const Var StoppedChargedPionTrackChi2ProtonCollection([](const caf::SRSliceProxy* slc) -> double {
    int stoppedCPionIndex = StoppedChargedPionTrackIndex(slc);
    if(stoppedCPionIndex>=0){
      const auto& trk = slc->reco.pfp.at(stoppedCPionIndex).trk;
      if(trk.calo[2].nhit==0) return -999.;
      else return trk.chi2pid[2].chi2_proton;
    }
    else{
      return -999.;
    }
  }); 
  // charged pion, inelastic
  const Var InelasticChargedPionTrackIndex([](const caf::SRSliceProxy* slc) -> double {
    if(slc->reco.pfp.size()==0){
      return -999.;
    }
    else{

      // Requiring good muon
      int muonTrackIndex = MuonTrackIndex(slc);
      int PTrackInd(-1);
      if(muonTrackIndex>=0){

        int protonTrackIndex = ProtonTrackIndex(slc);

        float Longest(0);
        for(unsigned int i_pfp=0; i_pfp<slc->reco.pfp.size(); i_pfp++){
          const auto& pfp = slc->reco.pfp.at(i_pfp);
          const auto& trk = pfp.trk;

          if(i_pfp==(unsigned int)muonTrackIndex) continue;
          if(i_pfp==(unsigned int)protonTrackIndex) continue;
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

          const bool Contained = fv_track.isContained(trk.end.x, trk.end.y, trk.end.z);

          bool isPionInelastic = false;
          if(trk.calo[2].nhit!=0){
            double new_chi2 = dedxtempt.CalculateInelasticPionChi2(trk.calo[2], 3);
            isPionInelastic = (Contained && trk.calo[2].nhit>=5 && new_chi2<10.);
          }

          if ( AtSlice && ( isPionInelastic ) && trk.len > Longest )
          {
            Longest = trk.len;
            PTrackInd = i_pfp;
          }
        } // END pfp loop

      } // END if muon track exist

      return PTrackInd;

    } // END if slice has primary tracks
  });
  const Var InelasticChargedPionTrackLength([](const caf::SRSliceProxy* slc) -> double {
    int inelCPionIndex = InelasticChargedPionTrackIndex(slc);
    if(inelCPionIndex>=0){
      const auto& trk = slc->reco.pfp.at(inelCPionIndex).trk;
      return trk.len;
    }
    else{
      return -999.;
    }
  });
  const Var InelasticChargedPionTrackNuMICosineTheta([](const caf::SRSliceProxy* slc) -> double {
    int inelCPionIndex = InelasticChargedPionTrackIndex(slc);
    if(inelCPionIndex>=0){
      const auto& trk = slc->reco.pfp.at(inelCPionIndex).trk;
      TVector3 v3_trk(trk.dir.x, trk.dir.y, trk.dir.z);
      double angleNuMI = v3_trk.Angle(NuDirection_NuMI);
      return TMath::Cos(angleNuMI);
    }
    else{
      return -999.;
    }
  });
  const Var InelasticChargedPionTrackNuMIToVtxCosineTheta([](const caf::SRSliceProxy* slc) -> double {
    int inelCPionIndex = InelasticChargedPionTrackIndex(slc);
    if(inelCPionIndex>=0){
      const auto& trk = slc->reco.pfp.at(inelCPionIndex).trk;
      const auto& vtx = slc->vertex;

      TVector3 v3_trk(trk.dir.x, trk.dir.y, trk.dir.z);
      TVector3 v3_vtx(vtx.x, vtx.y, vtx.z);
      static const TVector3 dFromNuMI(315.120380, 33.644912, 733.632532);

      TVector3 v3_numi_to_vtx = dFromNuMI+v3_vtx;

      double angleNuMI = v3_trk.Angle(v3_numi_to_vtx);
      return TMath::Cos(angleNuMI);
    }
    else{
      return -999.;
    }
  });
  const Var InelasticChargedPionTrackChi2MIPCollection([](const caf::SRSliceProxy* slc) -> double {
    int inelCPionIndex = InelasticChargedPionTrackIndex(slc);
    if(inelCPionIndex>=0){
      const auto& trk = slc->reco.pfp.at(inelCPionIndex).trk;
      if(trk.calo[2].nhit==0) return -999.;
      else return dedxtempt.CalculateInelasticPionChi2(trk.calo[2], 3);
    }
    else{
      return -999.;
    }
  });

  // Neutral pion
  const MultiVar NeutralPionPhotonShowerIndices([](const caf::SRSliceProxy* slc) -> vector<double> {
    vector<double> rets;
    for(unsigned int i_pfp=0; i_pfp<slc->reco.pfp.size(); i_pfp++){
      const auto& pfp = slc->reco.pfp.at(i_pfp);

      if( !IsPFPShower(pfp) ) continue; // TODO

      const auto& shw = pfp.shw;

      // Check if shower fit even seems kind-of valid:
      if ( std::isnan(shw.start.x) || (shw.start.x > -5.5 && shw.start.x < -4.5) ||
           std::isnan(shw.len) || shw.len <= 0. ) continue;

      // if it meets this then we're not going to cut on it...
      if ( std::isnan(shw.plane[2].energy) || std::isinf(shw.plane[2].energy) || shw.plane[2].energy <= 0.04 ) continue;

      // and... if it meets then then we're not going to cut on it...
      if ( std::isnan(shw.conversion_gap) || std::isinf(shw.conversion_gap) || shw.conversion_gap <= 5. ) continue;

      rets.push_back(i_pfp);
    }

    return rets;

  });
  const MultiVar NeutralPionPhotonShowerConvGaps([](const caf::SRSliceProxy* slc) -> vector<double> {
    vector<double> shw_indices = NeutralPionPhotonShowerIndices(slc);
    vector<double> rets;
    for(const auto& shw_index: shw_indices){
      const auto& shw = slc->reco.pfp.at(shw_index).shw;
      rets.push_back( shw.conversion_gap );
    }
    return rets;
  });
  const MultiVar NeutralPionPhotonShowerLengths([](const caf::SRSliceProxy* slc) -> vector<double> {
    vector<double> shw_indices = NeutralPionPhotonShowerIndices(slc);
    vector<double> rets;
    for(const auto& shw_index: shw_indices){
      const auto& shw = slc->reco.pfp.at(shw_index).shw;
      rets.push_back( shw.len );
    }
    return rets;
  });
  const MultiVar NeutralPionPhotonShowerEnergies([](const caf::SRSliceProxy* slc) -> vector<double> {
    vector<double> shw_indices = NeutralPionPhotonShowerIndices(slc);
    vector<double> rets;
    for(const auto& shw_index: shw_indices){
      const auto& shw = slc->reco.pfp.at(shw_index).shw;
      if( isnan(shw.plane[2].energy) ) rets.push_back(0.);
      else rets.push_back( shw.plane[2].energy );
    }
    return rets;
  });
  const MultiVar NeutralPionPhotonShowerdEdxs([](const caf::SRSliceProxy* slc) -> vector<double> {
    vector<double> shw_indices = NeutralPionPhotonShowerIndices(slc);
    vector<double> rets;
    for(const auto& shw_index: shw_indices){
      const auto& shw = slc->reco.pfp.at(shw_index).shw;
      if( isnan(shw.plane[2].dEdx) ) rets.push_back(0.);
      else rets.push_back( shw.plane[2].dEdx );
    }
    return rets;
  });
  const Var NNeutralPionPhotonShower([](const caf::SRSliceProxy* slc) -> double {
    return NeutralPionPhotonShowerIndices(slc).size();
  });
  const Var NeutralPionPhotonShowerSumEnergy([](const caf::SRSliceProxy* slc) -> double {
    vector<double> shw_indices = NeutralPionPhotonShowerIndices(slc);
    double esum = 0.;
    for(const auto& shw_index: shw_indices){
      const auto& shw = slc->reco.pfp.at(shw_index).shw;
      if( !isnan(shw.plane[2].energy) ) esum += shw.plane[2].energy;
    }
    return esum;
  });
  const Var NeutralPionPhotonShowerSumInvariantMass([](const caf::SRSliceProxy* slc) -> double {
    vector<double> shw_indices = NeutralPionPhotonShowerIndices(slc);
    double esum = 0.;
    TVector3 psum(0., 0., 0.);
    for(const auto& shw_index: shw_indices){
      const auto& shw = slc->reco.pfp.at(shw_index).shw;

      TVector3 this_p(shw.dir.x.GetValue(), shw.dir.y.GetValue(), shw.dir.z.GetValue());
      if( isnan(shw.plane[2].energy) ) continue;
      this_p *= shw.plane[2].energy;

      esum += shw.plane[2].energy;
      psum += this_p;
    }

    double m2 = esum*esum - psum.Mag2();
    return m2>0 ? sqrt(m2) : 0.;
  });

  namespace Aux{

    const SpillMultiVar TestSpillVar([](const caf::SRSpillProxy *sr) -> vector<double> {

      vector<double> rets;

      for(std::size_t i(0); i < sr->slc.size(); ++i){
        const auto& slc = sr->slc.at(i);
        int truth_index = ICARUSNumuXsec::TruthMatch::TruthNeutralPionIndex(&slc);
        if(truth_index>=0){
          std::cout << "[JSKIMDEBUG][Neutral pion event found]" << std::endl;
          PrintPrimaries(&slc);
        }
      }

      return rets;

    });
    const SpillMultiVar ForPrintingTruthInfos([](const caf::SRSpillProxy *sr) -> vector<double> {

      vector<double> rets;

      printf("\n(run, subrun, event) = (%d, %d, %d), nSlice = %ld\n", sr->hdr.run.GetValue(), sr->hdr.subrun.GetValue(), sr->hdr.evt.GetValue(),sr->slc.size());
      printf("- Number of truth nu = %ld\n",sr->mc.nu.size());
      for(std::size_t i(0); i < sr->mc.nu.size(); ++i){
        const auto& nu = sr->mc.nu[i];
        printf("- %ld-th nu\n",i);
        if(nu.iscc) printf("  - GENIE mode = %d, CC\n", nu.genie_mode.GetValue());
        else if(nu.isnc) printf("  - GENIE mode = %d, NC\n", nu.genie_mode.GetValue());
        else printf("  - GENIE mode = %d, not CC nor NC?\n", nu.genie_mode.GetValue());
        printf("  - time = %1.3f\n",nu.time.GetValue());

        rets.push_back( nu.time );

        printf("  - vertex = (%1.2f, %1.2f, %1.2f)\n", nu.vtx.x.GetValue(), nu.vtx.y.GetValue(), nu.vtx.z.GetValue());
        printf("  - nu pos = (%1.2f, %1.2f, %1.2f)\n", nu.position.x.GetValue(), nu.position.y.GetValue(), nu.position.z.GetValue());
        //for(std::size_t j(0); j < nu.prim.size(); ++j){
        //  const auto& prim = nu.prim[j];
        for(std::size_t j(0); j < nu.ghepptl.size(); ++j){
          const auto& prim = nu.ghepptl[j];
          printf("  - %ld-th prim\n",j);
          printf("    - pdg = %d\n", prim.pdg.GetValue());
          const double this_mass = ptlt.GetMass(prim.pdg.GetValue());
          printf("    - Mass = %1.3f\n", this_mass);
          printf("    - GStatus = %d\n", prim.gstatus.GetValue());
          printf("    - Parent ID = %d\n", prim.parent.GetValue());
          printf("    - ndaughters = %ld\n", prim.daughters.size());
          printf("    - start = (%1.2f, %1.2f, %1.2f)\n", prim.start.x.GetValue(), prim.start.y.GetValue(), prim.start.z.GetValue());
          printf("    - end = (%1.2f, %1.2f, %1.2f)\n", prim.end.x.GetValue(), prim.end.y.GetValue(), prim.end.z.GetValue());
          const float dist = std::hypot(prim.end.x - prim.start.x, prim.end.y - prim.start.y, prim.end.z - prim.start.z);
          printf("    - |end-start| = %1.2f\n", dist);
          printf("    - genp = (%1.2f, %1.2f, %1.2f)\n", prim.genp.x.GetValue(), prim.genp.y.GetValue(), prim.genp.z.GetValue());
          printf("    - startE = %1.3f\n", prim.startE.GetValue());
          printf("    - startE-Mass = %1.3f\n", prim.startE.GetValue()-this_mass);
          printf("    - startE-endE = %1.3f\n", prim.startE.GetValue()-prim.endE.GetValue());
          printf("    - end_process = %d\n", prim.end_process.GetValue());

        }
      }

      for(std::size_t i(0); i < sr->slc.size(); ++i){
        const auto& slc = sr->slc.at(i);
/*
        printf("- %ld-th slice\n",i);
        const auto& relaxedMuonIndex = RelaxedMuonTrackIndex(&slc);
        const auto& relaxedMuonTrack =  slc.reco.pfp[relaxedMuonIndex].trk;
        const auto& relaxedMuonTrackTruth = relaxedMuonTrack.truth;
        printf("  - Truth muon start = (%1.2f, %1.2f, %1.2f)\n", relaxedMuonTrackTruth.p.start.x.GetValue(), relaxedMuonTrackTruth.p.start.y.GetValue(), relaxedMuonTrackTruth.p.start.z.GetValue());
        printf("  - Truth muon end = (%1.2f, %1.2f, %1.2f)\n", relaxedMuonTrackTruth.p.end.x.GetValue(), relaxedMuonTrackTruth.p.end.y.GetValue(), relaxedMuonTrackTruth.p.end.z.GetValue());
        printf("  - Truth muon P = %1.2f GeV\n", RelaxedMuonTrackTruthP(&slc) );
        printf("  - Truth muon length = %1.2f\n", RelaxedMuonTrackTruthLength(&slc) );
        printf("  - Muon candidate index = %d:\n",int(relaxedMuonIndex));
        printf("  - Reco vertex = (%1.2f, %1.2f, %1.2f)\n", slc.vertex.x.GetValue(), slc.vertex.y.GetValue(), slc.vertex.z.GetValue());
        for(std::size_t ip(0); ip < slc.reco.pfp.size(); ++ip){
          printf("  - %ld-th pfp:\n",ip);
          printf("    - pfp id = %d\n", slc.reco.pfp[ip].id.GetValue());
          printf("    - trackScore = %1.2f\n", slc.reco.pfp[ip].trackScore.GetValue());
          printf("    - Track start = (%1.2f, %1.2f, %1.2f)\n", slc.reco.pfp[ip].trk.start.x.GetValue(), slc.reco.pfp[ip].trk.start.y.GetValue(), slc.reco.pfp[ip].trk.start.z.GetValue());
          printf("    - Track end = (%1.2f, %1.2f, %1.2f)\n", slc.reco.pfp[ip].trk.end.x.GetValue(), slc.reco.pfp[ip].trk.end.y.GetValue(), slc.reco.pfp[ip].trk.end.z.GetValue());

          const float dist = std::hypot(slc.reco.pfp[ip].trk.end.x - slc.reco.pfp[ip].trk.start.x, slc.reco.pfp[ip].trk.end.y - slc.reco.pfp[ip].trk.start.y, slc.reco.pfp[ip].trk.end.z - slc.reco.pfp[ip].trk.start.z);
          printf("    - Track |end-start| = %1.2f\n", dist);

          printf("    - Track length =  %1.2f\n", slc.reco.pfp[ip].trk.len.GetValue());
          printf("    - Track (chi2_muon, chi2_proton) = (%1.2f, %1.2f)\n", slc.reco.pfp[ip].trk.chi2pid[2].chi2_muon.GetValue(), slc.reco.pfp[ip].trk.chi2pid[2].chi2_proton.GetValue());
          printf("    - bestmatch pdg = %d\n", slc.reco.pfp[ip].trk.truth.p.pdg.GetValue());
          printf("    - bsetmatch start = (%1.2f, %1.2f, %1.2f)\n", slc.reco.pfp[ip].trk.truth.p.start.x.GetValue(), slc.reco.pfp[ip].trk.truth.p.start.y.GetValue(), slc.reco.pfp[ip].trk.truth.p.start.z.GetValue());
          printf("    - bestmatch end = (%1.2f, %1.2f, %1.2f)\n", slc.reco.pfp[ip].trk.truth.p.end.x.GetValue(), slc.reco.pfp[ip].trk.truth.p.end.y.GetValue(), slc.reco.pfp[ip].trk.truth.p.end.z.GetValue());
        }
*/

        const auto& relaxedChargedPionIndex = RelaxedChargedPionTrackIndex(&slc);
        if(relaxedChargedPionIndex<0) continue;
        const auto& relaxedChargedPionTrack =  slc.reco.pfp[relaxedChargedPionIndex].trk;
        const auto& relaxedChargedPionTrackTruth = relaxedChargedPionTrack.truth;
        bool isConvPion = (relaxedChargedPionTrackTruth.p.end_process==36);
        if(!isConvPion) continue;
        printf("- %ld-th slice\n",i);
        printf("  - Truth pi+- start = (%1.2f, %1.2f, %1.2f)\n", relaxedChargedPionTrackTruth.p.start.x.GetValue(), relaxedChargedPionTrackTruth.p.start.y.GetValue(), relaxedChargedPionTrackTruth.p.start.z.GetValue());
        printf("  - Truth pi+- end = (%1.2f, %1.2f, %1.2f)\n", relaxedChargedPionTrackTruth.p.end.x.GetValue(), relaxedChargedPionTrackTruth.p.end.y.GetValue(), relaxedChargedPionTrackTruth.p.end.z.GetValue());
        printf("  - Truth pi+- end_process = %d\n", relaxedChargedPionTrackTruth.p.end_process.GetValue());

      }

      return rets;

    });

    const Var RelaxedMuonTrackIndex([](const caf::SRSliceProxy* slc) -> double {
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

    const Var RelaxedMuonTrackLength([](const caf::SRSliceProxy* slc) -> double {
      int muonTrackIndex = RelaxedMuonTrackIndex(slc);
      if(muonTrackIndex>=0){
        return slc->reco.pfp.at(muonTrackIndex).trk.len;
      }
      else{
        return -999.;
      }
    });
    const Var RelaxedMuonTrackP([](const caf::SRSliceProxy* slc) -> double {
      int muonTrackIndex = RelaxedMuonTrackIndex(slc);
      if(muonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
        bool Contained = fv_track.isContained(trk.end.x, trk.end.y, trk.end.z);
        if(Contained){
          if(isnan(trk.rangeP.p_muon)) return -999.;
          return trk.rangeP.p_muon;
        }
        else{
          if(isnan(trk.mcsP.fwdP_muon)) return -999.;
          return trk.mcsP.fwdP_muon;
        }

      }
      else{
        return -999.;
      }
    });
    const Var RelaxedMuonTrackNuMICosineTheta([](const caf::SRSliceProxy* slc) -> double {
      int muonTrackIndex = RelaxedMuonTrackIndex(slc);
      if(muonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
        TVector3 v3_trk(trk.dir.x, trk.dir.y, trk.dir.z);
        double angleNuMI = v3_trk.Angle(NuDirection_NuMI);
        return TMath::Cos(angleNuMI);
      }
      else{
        return -999.;
      }
    });
    const Var RelaxedMuonTrackPosAbsX([](const caf::SRSliceProxy* slc) -> double {
      int muonTrackIndex = RelaxedMuonTrackIndex(slc);
      if(muonTrackIndex>=0){
        return abs(slc->reco.pfp.at(muonTrackIndex).trk.end.x);
      }
      else{
        return -999.;
      }
    });
    const Var RelaxedMuonTrackDirX([](const caf::SRSliceProxy* slc) -> double {
      int muonTrackIndex = RelaxedMuonTrackIndex(slc);
      if(muonTrackIndex>=0){
        return slc->reco.pfp.at(muonTrackIndex).trk.dir.x;
      }
      else{
        return -999.;
      }
    });
    const Var RelaxedMuonTrackDirY([](const caf::SRSliceProxy* slc) -> double {
      int muonTrackIndex = RelaxedMuonTrackIndex(slc);
      if(muonTrackIndex>=0){
        return slc->reco.pfp.at(muonTrackIndex).trk.dir.y;
      }
      else{
        return -999.;
      }
    });
    const Var RelaxedMuonTrackDirZ([](const caf::SRSliceProxy* slc) -> double {
      int muonTrackIndex = RelaxedMuonTrackIndex(slc);
      if(muonTrackIndex>=0){
        return slc->reco.pfp.at(muonTrackIndex).trk.dir.z;
      }
      else{
        return -999.;
      }
    });

    const Var RelaxedMuonTrackTrackScore([](const caf::SRSliceProxy* slc) -> double {
      int muonTrackIndex = RelaxedMuonTrackIndex(slc);
      if(muonTrackIndex>=0){
        return slc->reco.pfp.at(muonTrackIndex).trackScore;
      }
      else{
        return -999.;
      }
    });
    const Var RelaxedMuonTrackChi2Muon([](const caf::SRSliceProxy* slc) -> double {
      int muonTrackIndex = RelaxedMuonTrackIndex(slc);
      if(muonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
        if(trk.bestplane == -1) return -999.;

        return trk.chi2pid[trk.bestplane].chi2_muon;
      }
      else{
        return -999.;
      }
    });
    const Var RelaxedMuonTrackChi2Proton([](const caf::SRSliceProxy* slc) -> double {
      int muonTrackIndex = RelaxedMuonTrackIndex(slc);
      if(muonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
        if(trk.bestplane == -1) return -999.;

        return trk.chi2pid[trk.bestplane].chi2_proton;
      }
      else{
        return -999.;
      }
    });
    const Var RelaxedMuonTrackChi2MuonCollection([](const caf::SRSliceProxy* slc) -> double {
      int muonTrackIndex = RelaxedMuonTrackIndex(slc);
      if(muonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
        if(trk.bestplane == -1) return -999.;
        if(trk.calo[2].nhit==0) return -999.;
        
        return trk.chi2pid[2].chi2_muon;
      }
      else{
        return -999.;
      }
    });
    const Var RelaxedMuonTrackChi2ProtonCollection([](const caf::SRSliceProxy* slc) -> double {
      int muonTrackIndex = RelaxedMuonTrackIndex(slc);
      if(muonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
        if(trk.bestplane == -1) return -999.;
        if(trk.calo[2].nhit==0) return -999.;
        
        return trk.chi2pid[2].chi2_proton;
      }
      else{
        return -999.;
      }
    });
    const Var RelaxedMuonTrackCustomChi2MuonCollection([](const caf::SRSliceProxy* slc) -> double {
      int muonTrackIndex = RelaxedMuonTrackIndex(slc);
      if(muonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
        if(trk.calo[2].nhit==0) return -999.;
        double new_chi2 = dedxtempt.CalculateChi2(trk.calo[2], 3);
        return new_chi2;
      }
      else{
        return -999.;
      }
    });
    const Var RelaxedMuonTrackCustomChi2ProtonCollection([](const caf::SRSliceProxy* slc) -> double {
      int muonTrackIndex = RelaxedMuonTrackIndex(slc);
      if(muonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
        if(trk.calo[2].nhit==0) return -999.;
        double new_chi2 = dedxtempt.CalculateChi2(trk.calo[2], 0);
        return new_chi2;
      }
      else{
        return -999.;
      }
    });
    const MultiVar RelaxedMuonTrackCollectionRR([](const caf::SRSliceProxy* slc) -> vector<double> {
      int muonTrackIndex = RelaxedMuonTrackIndex(slc);
      vector<double> rets;
      if(muonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
        for(const auto& pt: trk.calo[2].points){
          rets.push_back( pt.rr );
        }
      }
      return rets;
    });
    const MultiVar RelaxedMuonTrackCollectiondEdX([](const caf::SRSliceProxy* slc) -> vector<double> {
      int muonTrackIndex = RelaxedMuonTrackIndex(slc);
      vector<double> rets;
      if(muonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
        for(const auto& pt: trk.calo[2].points){
          rets.push_back( pt.dedx );
        }
      }
      return rets;
    });
    const MultiVar RelaxedMuonTrackCollectiondQdX([](const caf::SRSliceProxy* slc) -> vector<double> {
      int muonTrackIndex = RelaxedMuonTrackIndex(slc);
      vector<double> rets;
      if(muonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
        for(const auto& pt: trk.calo[2].points){
          rets.push_back( pt.dqdx );
        }
      }
      return rets;
    });
    // truth match
    const Var RelaxedMuonTrackTruthLength([](const caf::SRSliceProxy* slc) -> double {
      int muonTrackIndex = RelaxedMuonTrackIndex(slc);
      if(muonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
        if( isnan(trk.truth.p.length) ) return -999.;
        return trk.truth.p.length;
      }
      else{
        return -999.;
      }
    });

    const Var RelaxedMuonTrackTruthP([](const caf::SRSliceProxy* slc) -> double {
      int muonTrackIndex = RelaxedMuonTrackIndex(slc);
      if(muonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
        if( isnan(trk.truth.p.genE) ) return -999.;
        double P2 = trk.truth.p.genE*trk.truth.p.genE - M_MUON*M_MUON;
        if(P2>0) return std::sqrt(P2);
        else return -999.;
      }
      else{
        return -999.;
      }
    });
    const Var RelaxedMuonTrackTruthPResFrac([](const caf::SRSliceProxy* slc) -> double {
      double Reco = RelaxedMuonTrackP(slc);
      double True = RelaxedMuonTrackTruthP(slc);
      if(Reco>0. && True>0.){
        return (Reco-True)/True;
      }
      else{
        return -999.;
      }
    });
    const Var RelaxedMuonTrackTruthOneOverP([](const caf::SRSliceProxy* slc) -> double {
      double P = RelaxedMuonTrackTruthP(slc);
      if(P>0.) return 1./P;
      else return -999.;
    });
    const Var RelaxedMuonTrackTruthOneOverPResFrac([](const caf::SRSliceProxy* slc) -> double {
      double RecoP = RelaxedMuonTrackP(slc);
      double True = RelaxedMuonTrackTruthOneOverP(slc);
      if(RecoP>0. && True>0.){
        double Reco = 1./RecoP;
        return (Reco-True)/True;
      }
      else{
        return -999.;
      }
    });

    const Var RelaxedMuonTrackTruthNuMICosineTheta([](const caf::SRSliceProxy* slc) -> double {
      int muonTrackIndex = RelaxedMuonTrackIndex(slc);
      if(muonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
        const auto& truth_p = trk.truth.p;
        if( isnan(truth_p.genp.x) ) return -999.;
        TVector3 v3_gen_p(truth_p.genp.x, truth_p.genp.y, truth_p.genp.z);
        double angleNuMI = v3_gen_p.Angle(NuDirection_NuMI);
        return TMath::Cos(angleNuMI);
      }
      else{
        return -999.;
      }
    });
    const Var RelaxedMuonTrackTruthStartProcess([](const caf::SRSliceProxy* slc) -> double {
      int muonTrackIndex = RelaxedMuonTrackIndex(slc);
      if(muonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
        return trk.truth.p.start_process;
      }
      else{
        return -999.;
      }
    });
    const Var RelaxedMuonTrackTruthEndProcess([](const caf::SRSliceProxy* slc) -> double {
      int muonTrackIndex = RelaxedMuonTrackIndex(slc);
      if(muonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(muonTrackIndex).trk;
        return trk.truth.p.end_process;
      }
      else{
        return -999.;
      }
    });

    const Var RelaxedProtonTrackIndex([](const caf::SRSliceProxy* slc) -> double {
      if(slc->reco.pfp.size()==0){
        return -999.;
      }
      else{

        // Requiring good muon
        int muonTrackIndex = MuonTrackIndex(slc);
        int protonIdx = -1;

        //const int chargedpionTrackIndex = ChargedPionTrackIndex(slc);
        //const bool HasPionCand = (chargedpionTrackIndex>=0);
        const bool HasPionCand = false;

        if(muonTrackIndex>=0 && !HasPionCand){

          double protonMaxLength = -1.;

          // find the longest track

          for(unsigned int i_pfp=0; i_pfp<slc->reco.pfp.size(); i_pfp++){
            const auto& pfp = slc->reco.pfp.at(i_pfp);
            const auto& trk = pfp.trk;

            if(i_pfp==(unsigned int )muonTrackIndex) continue;
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


          }

        }

        if(HasPionCand) return -1;
        return protonIdx;

      }
    });
    const Var RelaxedProtonTrackLength([](const caf::SRSliceProxy* slc) -> double {
      int protonTrackIndex = RelaxedProtonTrackIndex(slc);
      if(protonTrackIndex>=0){
        return slc->reco.pfp.at(protonTrackIndex).trk.len;
      }
      else{
        return -999.;
      }
    });
    const Var RelaxedProtonTrackP([](const caf::SRSliceProxy* slc) -> double {
      int protonTrackIndex = RelaxedProtonTrackIndex(slc);
      if(protonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
        bool Contained = fv_track.isContained(trk.end.x, trk.end.y, trk.end.z);
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
    const Var RelaxedProtonTrackNuMICosineTheta([](const caf::SRSliceProxy* slc) -> double {
      int protonTrackIndex = RelaxedProtonTrackIndex(slc);
      if(protonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
        TVector3 v3_trk(trk.dir.x, trk.dir.y, trk.dir.z);
        double angleNuMI = v3_trk.Angle(NuDirection_NuMI);
        return TMath::Cos(angleNuMI);
      }
      else{
        return -999.;
      }
    });
    const Var RelaxedProtonTrackDirX([](const caf::SRSliceProxy* slc) -> double {
      int protonTrackIndex = RelaxedProtonTrackIndex(slc);
      if(protonTrackIndex>=0){
        return slc->reco.pfp.at(protonTrackIndex).trk.dir.x;
      }
      else{
        return -999.;
      }
    });
    const Var RelaxedProtonTrackDirY([](const caf::SRSliceProxy* slc) -> double {
      int protonTrackIndex = RelaxedProtonTrackIndex(slc);
      if(protonTrackIndex>=0){
        return slc->reco.pfp.at(protonTrackIndex).trk.dir.y;
      }
      else{
        return -999.;
      }
    });
    const Var RelaxedProtonTrackDirZ([](const caf::SRSliceProxy* slc) -> double {
      int protonTrackIndex = RelaxedProtonTrackIndex(slc);
      if(protonTrackIndex>=0){
        return slc->reco.pfp.at(protonTrackIndex).trk.dir.z;
      }
      else{
        return -999.;
      }
    });
    const Var RelaxedProtonTrackTrackScore([](const caf::SRSliceProxy* slc) -> double {
      int protonTrackIndex = RelaxedProtonTrackIndex(slc);
      if(protonTrackIndex>=0){
        return slc->reco.pfp.at(protonTrackIndex).trackScore;
      }
      else{
        return -999.;
      }
    });
    const Var RelaxedProtonTrackChi2Muon([](const caf::SRSliceProxy* slc) -> double {
      int protonTrackIndex = RelaxedProtonTrackIndex(slc);
      if(protonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
        if(trk.bestplane == -1) return -999.;

        return trk.chi2pid[trk.bestplane].chi2_muon;
      }
      else{
        return -999.;
      }
    });
    const Var RelaxedProtonTrackChi2Proton([](const caf::SRSliceProxy* slc) -> double {
      int protonTrackIndex = RelaxedProtonTrackIndex(slc);
      if(protonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
        if(trk.bestplane == -1) return -999.;

        return trk.chi2pid[trk.bestplane].chi2_proton;
      }
      else{
        return -999.;
      }
    });
    const Var RelaxedProtonTrackChi2MuonCollection([](const caf::SRSliceProxy* slc) -> double {
      int protonTrackIndex = RelaxedProtonTrackIndex(slc);
      if(protonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
        if(trk.bestplane == -1) return -999.;
        if(trk.calo[2].nhit==0) return -999.;
        return trk.chi2pid[2].chi2_muon;
      }
      else{
        return -999.;
      }
    });
    const Var RelaxedProtonTrackChi2ProtonCollection([](const caf::SRSliceProxy* slc) -> double {
      int protonTrackIndex = RelaxedProtonTrackIndex(slc);
      if(protonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
        if(trk.bestplane == -1) return -999.;
        if(trk.calo[2].nhit==0) return -999.;

        return trk.chi2pid[2].chi2_proton;
      }
      else{
        return -999.;
      }
    });
    const Var RelaxedProtonTrackCustomChi2MuonCollection([](const caf::SRSliceProxy* slc) -> double {
      int protonTrackIndex = RelaxedProtonTrackIndex(slc);
      if(protonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
        if(trk.calo[2].nhit==0) return -999.;
        double new_chi2 = dedxtempt.CalculateChi2(trk.calo[2], 3);
        return new_chi2;
      }
      else{
        return -999.;
      }
    });
    const Var RelaxedProtonTrackCustomChi2ProtonCollection([](const caf::SRSliceProxy* slc) -> double {
      int protonTrackIndex = RelaxedProtonTrackIndex(slc);
      if(protonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
        if(trk.calo[2].nhit==0) return -999.;
        double new_chi2 = dedxtempt.CalculateChi2(trk.calo[2], 0);
        return new_chi2;
      }
      else{
        return -999.;
      }
    });
    const Var RelaxedProtonTrackNHitCollection([](const caf::SRSliceProxy* slc) -> double {
      int protonTrackIndex = RelaxedProtonTrackIndex(slc);
      if(protonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
        return trk.calo[2].nhit;
      }
      else{
        return -999.;
      }
    });
    const MultiVar RelaxedProtonTrackCollectionRR([](const caf::SRSliceProxy* slc) -> vector<double> {
      int protonTrackIndex = RelaxedProtonTrackIndex(slc);
      vector<double> rets;
      if(protonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
        for(const auto& pt: trk.calo[2].points){
          rets.push_back( pt.rr );
        }
      }
      return rets;
    });
    const MultiVar RelaxedProtonTrackCollectiondEdX([](const caf::SRSliceProxy* slc) -> vector<double> {
      int protonTrackIndex = RelaxedProtonTrackIndex(slc);
      vector<double> rets;
      if(protonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
        for(const auto& pt: trk.calo[2].points){
          rets.push_back( pt.dedx );
        }
      }
      return rets;
    });
    const MultiVar RelaxedProtonTrackCollectiondQdX([](const caf::SRSliceProxy* slc) -> vector<double> {
      int protonTrackIndex = RelaxedProtonTrackIndex(slc);
      vector<double> rets;
      if(protonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
        for(const auto& pt: trk.calo[2].points){
          rets.push_back( pt.dqdx );
        }
      }
      return rets;
    });
    // truth match
    const Var RelaxedProtonTrackTruthLength([](const caf::SRSliceProxy* slc) -> double {
      int protonTrackIndex = RelaxedProtonTrackIndex(slc);
      if(protonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
        if( isnan(trk.truth.p.length) ) return -999.;
        return trk.truth.p.length;
      }
      else{
        return -999.;
      }
    });
    const Var RelaxedProtonTrackTruthP([](const caf::SRSliceProxy* slc) -> double {
      int protonTrackIndex = RelaxedProtonTrackIndex(slc);
      if(protonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
        if( isnan(trk.truth.p.genE) ) return -999.;
        double P2 = trk.truth.p.genE*trk.truth.p.genE - M_PROTON*M_PROTON;
        if(P2>0) return std::sqrt(P2);
        else return -999.;
      }
      else{
        return -999.;
      }
    });
    const Var RelaxedProtonTrackTruthPResFrac([](const caf::SRSliceProxy* slc) -> double {
      double Reco = RelaxedProtonTrackP(slc);
      double True = RelaxedProtonTrackTruthP(slc);
      if(Reco>0. && True>0.){
        return (Reco-True)/True;
      }
      else{
        return -999.;
      }
    });
    const Var RelaxedProtonTrackTruthNuMICosineTheta([](const caf::SRSliceProxy* slc) -> double {
      int protonTrackIndex = RelaxedProtonTrackIndex(slc);
      if(protonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
        const auto& truth_p = trk.truth.p;
        if( isnan(truth_p.genp.x) ) return -999.;
        TVector3 v3_gen_p(truth_p.genp.x, truth_p.genp.y, truth_p.genp.z);
        double angleNuMI = v3_gen_p.Angle(NuDirection_NuMI);
        return TMath::Cos(angleNuMI);
      }
      else{
        return -999.;
      }
    });
    const Var RelaxedProtonTrackTruthStartProcess([](const caf::SRSliceProxy* slc) -> double {
      int protonTrackIndex = RelaxedProtonTrackIndex(slc);
      if(protonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
        return trk.truth.p.start_process;
      }
      else{
        return -999.;
      }
    });
    const Var RelaxedProtonTrackTruthEndProcess([](const caf::SRSliceProxy* slc) -> double {
      int protonTrackIndex = RelaxedProtonTrackIndex(slc);
      if(protonTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(protonTrackIndex).trk;
        return trk.truth.p.end_process;
      }
      else{
        return -999.;
      }
    });

    const Var RelaxedChargedPionTrackIndex([](const caf::SRSliceProxy* slc) -> double {
      if(slc->reco.pfp.size()==0){
        return -999.;
      }
      else{

        // Requiring good muon
        int muonTrackIndex = MuonTrackIndex(slc);
        int PTrackInd(-1);
        if(muonTrackIndex>=0){

          float Longest(0);
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

            const bool Contained = fv_track.isContained(trk.end.x, trk.end.y, trk.end.z);

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
    const Var RelaxedChargedPionTrackLength([](const caf::SRSliceProxy* slc) -> double {
      int cpionTrackIndex = RelaxedChargedPionTrackIndex(slc);
      if(cpionTrackIndex>=0){
        return slc->reco.pfp.at(cpionTrackIndex).trk.len;
      }
      else{
        return -999.;
      }
    });
    const Var RelaxedChargedPionTrackP([](const caf::SRSliceProxy* slc) -> double {
      int cpionTrackIndex = RelaxedChargedPionTrackIndex(slc);
      if(cpionTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(cpionTrackIndex).trk;
        bool Contained = fv_track.isContained(trk.end.x, trk.end.y, trk.end.z);
        if(Contained){
          if(isnan(trk.rangeP.p_pion)) return -999.;
          return trk.rangeP.p_pion;
        }
        else{
          if(isnan(trk.mcsP.fwdP_pion)) return -999.;
          return trk.mcsP.fwdP_pion;
        }

      }
      else{
        return -999.;
      }
    });
    const Var RelaxedChargedPionTrackTrackScore([](const caf::SRSliceProxy* slc) -> double {
      int cpionTrackIndex = RelaxedChargedPionTrackIndex(slc);
      if(cpionTrackIndex>=0){
        return slc->reco.pfp.at(cpionTrackIndex).trackScore;
      }
      else{
        return -999.;
      }
    });
    const Var RelaxedChargedPionTrackNuMICosineTheta([](const caf::SRSliceProxy* slc) -> double {
      int cpionTrackIndex = RelaxedChargedPionTrackIndex(slc);
      if(cpionTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(cpionTrackIndex).trk;
        TVector3 v3_trk(trk.dir.x, trk.dir.y, trk.dir.z);
        double angleNuMI = v3_trk.Angle(NuDirection_NuMI);
        return TMath::Cos(angleNuMI);
      }
      else{
        return -999.;
      }
    });
    const Var RelaxedChargedPionTrackCustomChi2MuonCollection([](const caf::SRSliceProxy* slc) -> double {
      int cpionTrackIndex = RelaxedChargedPionTrackIndex(slc);
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
    const Var RelaxedChargedPionTrackCustomChi2ProtonCollection([](const caf::SRSliceProxy* slc) -> double {
      int cpionTrackIndex = RelaxedChargedPionTrackIndex(slc);
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
    const Var RelaxedChargedPionTrackCustomChi2InelasticPionCollection([](const caf::SRSliceProxy* slc) -> double {
      int cpionTrackIndex = RelaxedChargedPionTrackIndex(slc);
      if(cpionTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(cpionTrackIndex).trk;
        if(trk.calo[2].nhit==0) return -999.;
        double new_chi2 = dedxtempt.CalculateInelasticPionChi2(trk.calo[2], 3);
        return new_chi2;
      }
      else{
        return -999.;
      }
    });
    const Var RelaxedChargedPionTrackFromVertex([](const caf::SRSliceProxy* slc) -> double {
      int cpionTrackIndex = RelaxedChargedPionTrackIndex(slc);
      if(cpionTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(cpionTrackIndex).trk;
        const float Atslc = std::hypot(slc->vertex.x - trk.start.x,
                                       slc->vertex.y - trk.start.y,
                                       slc->vertex.z - trk.start.z);
        return Atslc;

      }
      else{
        return -999.;
      }
    });
    const Var RelaxedChargedPionTrackNDaughter([](const caf::SRSliceProxy* slc) -> double {
      int cpionTrackIndex = RelaxedChargedPionTrackIndex(slc);
      if(cpionTrackIndex>=0){
        return slc->reco.pfp.at(cpionTrackIndex).ndaughters;
      }
      else{
        return -999.;
      }
    });
    const MultiVar RelaxedChargedPionTrackCollectionRR([](const caf::SRSliceProxy* slc) -> vector<double> {
      int cpionTrackIndex = RelaxedChargedPionTrackIndex(slc);
      vector<double> rets;
      if(cpionTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(cpionTrackIndex).trk;
        for(const auto& pt: trk.calo[2].points){
          rets.push_back( pt.rr );
        }
      }
      return rets;
    });
    const MultiVar RelaxedChargedPionTrackCollectiondEdX([](const caf::SRSliceProxy* slc) -> vector<double> {
      int cpionTrackIndex = RelaxedChargedPionTrackIndex(slc);
      vector<double> rets;
      if(cpionTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(cpionTrackIndex).trk;
        for(const auto& pt: trk.calo[2].points){
          rets.push_back( pt.dedx );
        }
      }
      return rets;
    });
    const MultiVar RelaxedChargedPionTrackCollectiondQdX([](const caf::SRSliceProxy* slc) -> vector<double> {
      int cpionTrackIndex = RelaxedChargedPionTrackIndex(slc);
      vector<double> rets;
      if(cpionTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(cpionTrackIndex).trk;
        for(const auto& pt: trk.calo[2].points){
          rets.push_back( pt.dqdx );
        }
      }
      return rets;
    });
    const MultiVar RelaxedChargedPionTrackCollectionFrontRR([](const caf::SRSliceProxy* slc) -> std::vector<double> {
      int cpionTrackIndex = RelaxedChargedPionTrackIndex(slc);
      vector<double> rets;
      if(cpionTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(cpionTrackIndex).trk;
        vector<float> rrs;
        for(const auto& pt: trk.calo[2].points){
          rrs.push_back( pt.rr );
        }
        float rrmax = !rrs.empty() ? *std::max_element(rrs.begin(), rrs.end()) : 0.;
        for( auto const& rr : rrs ){
          rets.push_back(rrmax-rr);
        }
      }
      return rets;
    });
    // truth match
    const Var RelaxedChargedPionTrackTruthLength([](const caf::SRSliceProxy* slc) -> double {
      int cpionTrackIndex = RelaxedChargedPionTrackIndex(slc);
      if(cpionTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(cpionTrackIndex).trk;
        if( isnan(trk.truth.p.length) ) return -999.;
        return trk.truth.p.length;
      }
      else{
        return -999.;
      }
    });
    const Var RelaxedChargedPionTrackTruthP([](const caf::SRSliceProxy* slc) -> double {
      int cpionTrackIndex = RelaxedChargedPionTrackIndex(slc);
      if(cpionTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(cpionTrackIndex).trk;
        if( isnan(trk.truth.p.genE) ) return -999.;
        double P2 = trk.truth.p.genE*trk.truth.p.genE - M_CHARGEDPION*M_CHARGEDPION;
        if(P2>0) return std::sqrt(P2);
        else return -999.;
      }
      else{
        return -999.;
      }
    });
    const Var RelaxedChargedPionTrackTruthPDG([](const caf::SRSliceProxy* slc) -> double {
      int cpionTrackIndex = RelaxedChargedPionTrackIndex(slc);
      if(cpionTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(cpionTrackIndex).trk;
        return trk.truth.p.pdg;
      }
      else{
        return -9999999.;
      }
    });
    const Var RelaxedChargedPionTrackTruthStartProcess([](const caf::SRSliceProxy* slc) -> double {
      int cpionTrackIndex = RelaxedChargedPionTrackIndex(slc);
      if(cpionTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(cpionTrackIndex).trk;
        return trk.truth.p.start_process;
      }
      else{
        return -999.;
      }
    });
    const Var RelaxedChargedPionTrackTruthEndProcess([](const caf::SRSliceProxy* slc) -> double {
      int cpionTrackIndex = RelaxedChargedPionTrackIndex(slc);
      if(cpionTrackIndex>=0){
        const auto& trk = slc->reco.pfp.at(cpionTrackIndex).trk;
        return trk.truth.p.end_process;
      }
      else{
        return -999.;
      }
    });


  } // end namespace Aux


} // end namespace TwoTrack

} // end namespace ICARUSNumuXsec
