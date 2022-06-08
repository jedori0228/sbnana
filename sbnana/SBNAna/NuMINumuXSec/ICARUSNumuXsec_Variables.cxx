#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_Variables.h"

using namespace std;
using namespace ana;

namespace ICARUSNumuXsec{

  //==== Spill variable
  const SpillMultiVar spillvarCRTHitTime([](const caf::SRSpillProxy *sr)
  {
    std::vector<double> times;
    for(const auto& hit : sr->crt_hits){
      times.push_back(hit.time);
    }
    return times;
  });

  const SpillMultiVar spillvarCRTHitT0([](const caf::SRSpillProxy *sr)
  {
    std::vector<double> times;
    for(const auto& hit : sr->crt_hits){
      times.push_back(hit.t0);
    }
    return times;
  });

  //==== only for south
  const SpillMultiVar spillvarCRTHitPosX([](const caf::SRSpillProxy *sr)
  {
    std::vector<double> poses;
    for(const auto& hit : sr->crt_hits){
      if(hit.plane==46){
        poses.push_back(hit.position.x);
      }
    }
    return poses;
  });
  const SpillMultiVar spillvarCRTHitPosY([](const caf::SRSpillProxy *sr)
  {
    std::vector<double> poses;
    for(const auto& hit : sr->crt_hits){
      if(hit.plane==46){
        poses.push_back(hit.position.y);
      }
    }
    return poses;
  });

  //==== Slice variables

  const Var varCountSlice([](const caf::SRSliceProxy* slc) ->int {
    return 0.;
  });

  const Var varVertexRecoX([](const caf::SRSliceProxy* slc) -> double {
    if(isnan(slc->vertex.x)) return -9999.;
    else return slc->vertex.x;
  });

  const Var varVertexRecoY([](const caf::SRSliceProxy* slc) -> double {
    if(isnan(slc->vertex.y)) return -9999.;
    else return slc->vertex.y;
  });

  const Var varVertexRecoZ([](const caf::SRSliceProxy* slc) -> double {
    if(isnan(slc->vertex.z)) return -9999.;
    else return slc->vertex.z;
  });

  const Var varFMScore([](const caf::SRSliceProxy* slc) -> double {
    if(isnan(slc->fmatch.score)) return -2.;
    else if(slc->fmatch.score<0) return -1.;
    else return slc->fmatch.score;
  });

  const Var varFMChargeQ([](const caf::SRSliceProxy* slc) -> double {
    if(isnan(slc->fmatch.chargeQ)) return -2.;
    else if(slc->fmatch.chargeQ<0) return -1.;
    else return slc->fmatch.chargeQ/1000.;
  });

  const Var varFMLightPE([](const caf::SRSliceProxy* slc) -> double {
    if(isnan(slc->fmatch.lightPE)) return -2.;
    else if(slc->fmatch.lightPE<0) return -1.;
    else return slc->fmatch.lightPE;
  });

  const Var varFMChargeCenterX([](const caf::SRSliceProxy* slc) -> double {
    if(isnan(slc->fmatch.chargeCenter.x)) return -9999.;
    else return slc->fmatch.chargeCenter.x;
  });

  const Var varFMChargeCenterY([](const caf::SRSliceProxy* slc) -> double {
    if(isnan(slc->fmatch.chargeCenter.y)) return -9999.;
    else return slc->fmatch.chargeCenter.y;
  });

  const Var varFMChargeCenterZ([](const caf::SRSliceProxy* slc) -> double {
    if(isnan(slc->fmatch.chargeCenter.z)) return -9999.;
    else return slc->fmatch.chargeCenter.z;
  });

  const Var varFMLightCenterX([](const caf::SRSliceProxy* slc) -> double {
    if(isnan(slc->fmatch.lightCenter.x)) return -9999.;
    else return slc->fmatch.lightCenter.x;
  });

  const Var varFMLightCenterY([](const caf::SRSliceProxy* slc) -> double {
    if(isnan(slc->fmatch.lightCenter.y)) return -9999.;
    else return slc->fmatch.lightCenter.y;
  });

  const Var varFMLightCenterZ([](const caf::SRSliceProxy* slc) -> double {
    if(isnan(slc->fmatch.lightCenter.z)) return -9999.;
    else return slc->fmatch.lightCenter.z;
  });

  const Var varFMTime([](const caf::SRSliceProxy* slc) -> double {
    if(isnan(slc->fmatch.time)) return -62.;
    else if(slc->fmatch.time<-50) return -61.;
    else return slc->fmatch.time;
  });

  const Var varTruthTime([](const caf::SRSliceProxy* slc) -> double {
    if(isnan(slc->truth.time)) return -999.;
    else return slc->truth.time;
  });

  const Var varSliceTrackNhitsPlane0([](const caf::SRSliceProxy* slc) -> int {
    int nHits(0);
    for(std::size_t i(0); i < slc->reco.trk.size(); ++i){
      auto const& trk = slc->reco.trk.at(i);
      nHits += trk.calo[0].nhit;
    }
    return nHits;
  });
  const Var varSliceTrackNhitsPlane1([](const caf::SRSliceProxy* slc) -> int {
    int nHits(0);
    for(std::size_t i(0); i < slc->reco.trk.size(); ++i){
      auto const& trk = slc->reco.trk.at(i);
      nHits += trk.calo[1].nhit;
    }
    return nHits;
  });
  const Var varSliceTrackNhitsPlane2([](const caf::SRSliceProxy* slc) -> int {
    int nHits(0);
    for(std::size_t i(0); i < slc->reco.trk.size(); ++i){
      auto const& trk = slc->reco.trk.at(i);
      nHits += trk.calo[2].nhit;
    }
    return nHits;
  });
  const Var varSliceShowerNhitsPlane0([](const caf::SRSliceProxy* slc) -> int {
    int nHits(0);
    for(std::size_t i(0); i < slc->reco.shw.size(); ++i){
      auto const& shw = slc->reco.shw.at(i);
      nHits += shw.plane[0].nHits;
    }
    return nHits;
  });
  const Var varSliceShowerNhitsPlane1([](const caf::SRSliceProxy* slc) -> int {
    int nHits(0);
    for(std::size_t i(0); i < slc->reco.shw.size(); ++i){
      auto const& shw = slc->reco.shw.at(i);
      nHits += shw.plane[1].nHits;
    }
    return nHits;
  });
  const Var varSliceShowerNhitsPlane2([](const caf::SRSliceProxy* slc) -> int {
    int nHits(0);
    for(std::size_t i(0); i < slc->reco.shw.size(); ++i){
      auto const& shw = slc->reco.shw.at(i);
      nHits += shw.plane[2].nHits;
    }
    return nHits;
  });
  const Var varSliceTrackChargePlane0([](const caf::SRSliceProxy* slc) -> double {
    double sumCharge(0);
    for(std::size_t i(0); i < slc->reco.trk.size(); ++i){
      auto const& trk = slc->reco.trk.at(i);
      sumCharge += trk.calo[0].charge;
      //std::cout << "[varSliceTrackChargePlane0] trk.calo[0].charge = " << trk.calo[0].charge << std::endl;
    }
    //std::cout << "[varSliceTrackChargePlane0] -> sumCharge = " << sumCharge << std::endl;
    return sumCharge/1000.;
  });
  const Var varSliceTrackChargePlane1([](const caf::SRSliceProxy* slc) -> double {
    double sumCharge(0);
    for(std::size_t i(0); i < slc->reco.trk.size(); ++i){
      auto const& trk = slc->reco.trk.at(i);
      sumCharge += trk.calo[1].charge;
    }
    return sumCharge/1000.;
  });
  const Var varSliceTrackChargePlane2([](const caf::SRSliceProxy* slc) -> double {
    double sumCharge(0);
    for(std::size_t i(0); i < slc->reco.trk.size(); ++i){
      auto const& trk = slc->reco.trk.at(i);
      sumCharge += trk.calo[2].charge;
    }
    return sumCharge/1000.;
  });

  const MultiVar varAllTrackStartPositionX([](const caf::SRSliceProxy* slc) -> std::vector<double> {
    std::vector<double> out;
    for (auto const& trk : slc->reco.trk) {
      out.push_back(trk.start.x);
    }
    return out;
  });

  const MultiVar varAllTrackStartPositionY([](const caf::SRSliceProxy* slc) -> std::vector<double> {
    std::vector<double> out;
    for (auto const& trk : slc->reco.trk) {
      out.push_back(trk.start.y);
    }
    return out;
  });

  const MultiVar varAllTrackStartPositionZ([](const caf::SRSliceProxy* slc) -> std::vector<double> {
    std::vector<double> out;
    for (auto const& trk : slc->reco.trk) {
      out.push_back(trk.start.z);
    }
    return out;
  });

  const Var varNuScore([](const caf::SRSliceProxy* slc) -> double {
    if(isnan(slc->nu_score)) return -0.1;
    else if(slc->nu_score<0) return -0.2;
    else return slc->nu_score;
  });
  const Var varSliceNuNFinalStatePfos([](const caf::SRSliceProxy* slc) -> double {
    return slc->nuid.nufspfos;
  });
  const Var varSliceNuNHitsTotal([](const caf::SRSliceProxy* slc) -> double {
    return slc->nuid.nutothits;
  });
  const Var varSliceNuVertexY([](const caf::SRSliceProxy* slc) -> double {
    return slc->nuid.nuvtxy;
  });
  const Var varSliceNuWeightedDirZ([](const caf::SRSliceProxy* slc) -> double {
    return slc->nuid.nuwgtdirz;
  });
  const Var varSliceNuNSpacePointsInSphere([](const caf::SRSliceProxy* slc) -> double {
    return slc->nuid.nusps;
  });
  const Var varSliceNuEigenRatioInSphere([](const caf::SRSliceProxy* slc) -> double {
    return slc->nuid.nueigen;
  });
  const Var varSliceCRLongestTrackDirY([](const caf::SRSliceProxy* slc) -> double {
    return slc->nuid.crlongtrkdiry;
  });
  const Var varSliceCRLongestTrackDeflection([](const caf::SRSliceProxy* slc) -> double {
    return slc->nuid.crlongtrkdef;
  });
  const Var varSliceCRFracHitsInLongestTrack([](const caf::SRSliceProxy* slc) -> double {
    return slc->nuid.crlongtrkhitfrac;
  });
  const Var varSliceCRNHitsMax([](const caf::SRSliceProxy* slc) -> double {
    return slc->nuid.crmaxhits;
  });

  //==== GENIE interaction code
  //==== https://internal.dunescience.org/doxygen/namespacesimb.html#a2cce734d1b71408bbc7d98d148ac4360
  const Var varGENIEIntCode([](const caf::SRSliceProxy* slc) -> int {
    //if(!isnan(slc->truth.genie_mode)){
    if(slc->truth.index >= 0.f){
      if(slc->truth.genie_mode<0) return -1;
      else if(slc->truth.genie_mode>13) return 14;
      else return slc->truth.genie_mode;
    }
    else{
      return -2;
    }
  });

  //==== Truth variables

  const Var varNeutrinoTruthE([](const caf::SRSliceProxy* slc) -> double {
    if(isnan(slc->truth.E)) return -999.;
    else return slc->truth.E;
  });
  const Var varTruthQ2([](const caf::SRSliceProxy* slc) -> double {
    if(isnan(slc->truth.Q2)) return -999.;
    else return slc->truth.Q2;
  });
  const Var varTruthq0_lab([](const caf::SRSliceProxy* slc) -> double {
    if(isnan(slc->truth.q0_lab)) return -999.;
    else return slc->truth.q0_lab;
  });
  const Var varTruthmodq_lab([](const caf::SRSliceProxy* slc) -> double {
    double Q2 = varTruthQ2(slc);
    double q0_lab = varTruthq0_lab(slc);
    if(isnan(Q2)||isnan(q0_lab)) return -999.;
    else{
      return sqrt(Q2*Q2+q0_lab*q0_lab);
    }
  });

  const Var varTruthNNeutron = SIMPLEVAR(truth.nneutron);
  const Var varTruthNPiMinus = SIMPLEVAR(truth.npiminus);
  const Var varTruthNPiPlus = SIMPLEVAR(truth.npiplus);
  const Var varTruthNChargedPion = varTruthNPiMinus+varTruthNPiPlus;
  const Var varTruthNPiZero = SIMPLEVAR(truth.npizero);
  const Var varTruthNProton = SIMPLEVAR(truth.nproton);

  const Var varMuonTruthIndex([](const caf::SRSliceProxy* slc) -> double {

    double max_E(-999);
    int truth_idx(-1);
    for(std::size_t i(0); i < slc->truth.prim.size(); ++i){
      if( abs(slc->truth.prim.at(i).pdg)== 13 ){
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

  const Var varMuonTruthP([](const caf::SRSliceProxy* slc) -> double {

    int truth_idx = varMuonTruthIndex(slc);
    if(truth_idx>=0){
      double genE = slc->truth.prim.at(truth_idx).genE;
      double genP2 = genE*genE-M_MUON*M_MUON;
      if(genP2>=0) return sqrt(genP2);
      else return 0;
    }
    else{
      return -999.;
    }

  });

  const Var varMuonTruthT([](const caf::SRSliceProxy* slc) -> double {

    int truth_idx = varMuonTruthIndex(slc);
    if(truth_idx>=0){
      double genE = slc->truth.prim.at(truth_idx).genE;
      return (genE-M_MUON);
    }
    else{
      return -999.;
    }

  });

  const Var varMuonTruthCosineTheta([](const caf::SRSliceProxy* slc) -> double {

    int truth_idx = varMuonTruthIndex(slc);
    if(truth_idx>=0){
      TVector3 v3(slc->truth.prim.at(truth_idx).genp.x, slc->truth.prim.at(truth_idx).genp.y, slc->truth.prim.at(truth_idx).genp.z);
      return v3.CosTheta();
    }
    else{
      return -999.;
    }

  });

  const Var varMuonTruthNuMICosineTheta([](const caf::SRSliceProxy* slc) -> double {

    int truth_idx = varMuonTruthIndex(slc);
    if(truth_idx>=0){
      TVector3 v3(slc->truth.prim.at(truth_idx).genp.x, slc->truth.prim.at(truth_idx).genp.y, slc->truth.prim.at(truth_idx).genp.z);
      double angleNuMI = v3.Angle(NuDirection_NuMI);
      return TMath::Cos(angleNuMI);
    }
    else{
      return -999.;
    }

  });

  const Var varProtonTruthIndex([](const caf::SRSliceProxy* slc) -> double {

    double max_E(-999);
    int truth_idx(-1);
    for(std::size_t i(0); i < slc->truth.prim.size(); ++i){
      if( abs(slc->truth.prim.at(i).pdg)== 2212 ){
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


  const Var varProtonTruthP([](const caf::SRSliceProxy* slc) -> double {

    int truth_idx = varProtonTruthIndex(slc);
    if(truth_idx>=0){
      double genE = slc->truth.prim.at(truth_idx).genE;
      double genP2 = genE*genE-M_PROTON*M_PROTON;
      if(genP2>=0) return sqrt(genP2);
      else return 0;
    }
    else{
      return -999.;
    }

  });

  const Var varProtonTruthT([](const caf::SRSliceProxy* slc) -> double {

    int truth_idx = varProtonTruthIndex(slc);
    if(truth_idx>=0){
      double genE = slc->truth.prim.at(truth_idx).genE;
      return (genE-M_PROTON);
    }
    else{
      return -999.;
    }

  });

  const Var varProtonTruthCosineTheta([](const caf::SRSliceProxy* slc) -> double {

    int truth_idx = varProtonTruthIndex(slc);
    if(truth_idx>=0){
      TVector3 v3(slc->truth.prim.at(truth_idx).genp.x, slc->truth.prim.at(truth_idx).genp.y, slc->truth.prim.at(truth_idx).genp.z);
      return v3.CosTheta();
    }
    else{
      return -999.;
    }

  });

  const Var varProtonTruthNuMICosineTheta([](const caf::SRSliceProxy* slc) -> double {

    int truth_idx = varProtonTruthIndex(slc);
    if(truth_idx>=0){
      TVector3 v3(slc->truth.prim.at(truth_idx).genp.x, slc->truth.prim.at(truth_idx).genp.y, slc->truth.prim.at(truth_idx).genp.z);
      double angleNuMI = v3.Angle(NuDirection_NuMI);
      return TMath::Cos(angleNuMI);

    }
    else{
      return -999.;
    }

  });

  const Var varTruthMuonProtonCosineTheta([](const caf::SRSliceProxy* slc) -> double {

    int truth_muon_idx = varMuonTruthIndex(slc);
    int truth_proton_idx = varProtonTruthIndex(slc);

    if(truth_muon_idx>=0 && truth_proton_idx>=0){
      TVector3 v3_muon(slc->truth.prim.at(truth_muon_idx).genp.x, slc->truth.prim.at(truth_muon_idx).genp.y, slc->truth.prim.at(truth_muon_idx).genp.z);
      TVector3 v3_proton(slc->truth.prim.at(truth_proton_idx).genp.x, slc->truth.prim.at(truth_proton_idx).genp.y, slc->truth.prim.at(truth_proton_idx).genp.z);
      return TMath::Cos( v3_muon.Angle(v3_proton) );
    }
    else{
      return -999.;
    }

  });

  //==== For a given truth, find the matching Reco track

  //====   For a given true muon (truth_index), find a reco track whose best-matched is this muon

  const Var varTruthMuonMatchedTrackIndex([](const caf::SRSliceProxy* slc) -> double {

    int truth_idx = varMuonTruthIndex(slc);

    return GetMatchedRecoTrackIndex(slc, truth_idx);

  });

  const Var varTruthMuonMatchedTrackContainedness([](const caf::SRSliceProxy* slc) -> int {
    int muonTrackIndex = varTruthMuonMatchedTrackIndex(slc);
    //==== -1 : No track found
    //==== 0 : Exiting
    //==== 1 : Contained
    if(muonTrackIndex>=0){

      auto const& trk = slc->reco.trk.at(muonTrackIndex);
      bool isContained = fv.isContained(trk.end.x, trk.end.y, trk.end.z);
      if(isContained) return 1;
      else return 0;
    }
    else{
      return -1;
    }

  });

  const Var varTruthMuonMatchedTrackChi2Proton([](const caf::SRSliceProxy* slc) -> double {
    int muonTrackIndex = varTruthMuonMatchedTrackIndex(slc);
    if(muonTrackIndex>=0){

      auto const& trk = slc->reco.trk.at(muonTrackIndex);
      int bp = trk.bestplane;
      float Chi2Proton = trk.chi2pid[bp].chi2_proton;
      if(isnan(Chi2Proton)) return -999.;
      else return Chi2Proton;

    }
    else{
      return 999999;
    }
  });

  const Var varTruthMuonMatchedTrackReducedChi2Proton([](const caf::SRSliceProxy* slc) -> double {
    int muonTrackIndex = varTruthMuonMatchedTrackIndex(slc);
    if(muonTrackIndex>=0){

      auto const& trk = slc->reco.trk.at(muonTrackIndex);
      int bp = trk.bestplane;
      float Chi2Proton = trk.chi2pid[bp].chi2_proton;
      int Chi2Ndof = trk.chi2pid[bp].pid_ndof;
      if(isnan(Chi2Proton)) return -999.;
      else if(Chi2Ndof==0) return -999.;
      else return Chi2Proton/double(Chi2Ndof);

    }
    else{
      return 999999;
    }
  });

  const Var varTruthMuonMatchedTrackChi2Muon([](const caf::SRSliceProxy* slc) -> double {
    int muonTrackIndex = varTruthMuonMatchedTrackIndex(slc);
    if(muonTrackIndex>=0){

      auto const& trk = slc->reco.trk.at(muonTrackIndex);
      int bp = trk.bestplane;
      float Chi2Muon = trk.chi2pid[bp].chi2_muon;
      if(isnan(Chi2Muon)) return -999.;
      else return Chi2Muon;

    }
    else{
      return 999999;
    }

  });

  const Var varTruthMuonMatchedTrackReducedChi2Muon([](const caf::SRSliceProxy* slc) -> double {
    int muonTrackIndex = varTruthMuonMatchedTrackIndex(slc);
    if(muonTrackIndex>=0){

      auto const& trk = slc->reco.trk.at(muonTrackIndex);
      int bp = trk.bestplane;
      float Chi2Muon = trk.chi2pid[bp].chi2_muon;
      int Chi2Ndof = trk.chi2pid[bp].pid_ndof;
      if(isnan(Chi2Muon)) return -999.;
      else if(Chi2Ndof==0) return -999.;
      else return Chi2Muon/double(Chi2Ndof);

    }
    else{
      return 999999;
    }

  });

  const Var varTruthMuonMatchedTrackLength([](const caf::SRSliceProxy* slc) -> double {
    int muonTrackIndex = varTruthMuonMatchedTrackIndex(slc);
    if(muonTrackIndex>=0){
      auto const& trk = slc->reco.trk.at(muonTrackIndex);
      return trk.len;
    }
    else{
      return -999;
    }

  });

  const Var varTruthMuonMatchedTrackRangeP([](const caf::SRSliceProxy* slc) -> double {
    int muonTrackIndex = varTruthMuonMatchedTrackIndex(slc);
    if(muonTrackIndex>=0){
      auto const& trk = slc->reco.trk.at(muonTrackIndex);
      return trk.rangeP.p_muon;
    }
    else{
      return -999;
    }

  });
  const Var varTruthMuonMatchedTrackMCSP([](const caf::SRSliceProxy* slc) -> double {
    int muonTrackIndex = varTruthMuonMatchedTrackIndex(slc);
    if(muonTrackIndex>=0){
      auto const& trk = slc->reco.trk.at(muonTrackIndex);
      double P(-999.);
      if(!isnan(trk.mcsP.fwdP_muon)){
        P = trk.mcsP.fwdP_muon;
      }
      return P;
    }
    else{
      return -999;
    }

  });
  const Var varTruthMuonMatchedTrackCombinedP([](const caf::SRSliceProxy* slc) -> double {
    int muonTrackIndex = varTruthMuonMatchedTrackIndex(slc);
    if(muonTrackIndex>=0){
      int mtc = varTruthMuonMatchedTrackContainedness(slc); // 0 : exiting, 1 : contained
      if(mtc==0) return varTruthMuonMatchedTrackMCSP(slc);
      else if(mtc==1) return varTruthMuonMatchedTrackRangeP(slc);
      else{
        cout << "[varTruthMuonMatchedTrackCombinedP] wtf?" << endl;
        exit(EXIT_FAILURE);
        return varTruthMuonMatchedTrackMCSP(slc); //
      }
    }
    else{
      return -999;
    }

  });

  const Var varTruthMuonMatchedTrackScore([](const caf::SRSliceProxy* slc) -> double {
    int muonTrackIndex = varTruthMuonMatchedTrackIndex(slc);
    if(muonTrackIndex>=0){
      return slc->reco.trk.at(muonTrackIndex).pfp.trackScore;
    }
    else{
      return -999;
    }

  });

  const Var varTruthMuonMatchedTrackVertexDistance([](const caf::SRSliceProxy* slc) -> double {
    int muonTrackIndex = varTruthMuonMatchedTrackIndex(slc);
    if(muonTrackIndex>=0){
      auto const& trk = slc->reco.trk.at(muonTrackIndex);
      return sqrt( pow( slc->vertex.x - trk.start.x, 2.0 ) +
              pow( slc->vertex.y - trk.start.y, 2.0 ) +
              pow( slc->vertex.z - trk.start.z, 2.0 ) );

    }
    else{
      return -999;
    }

  });

  const Var varTruthMuonMatchedShowerIndex([](const caf::SRSliceProxy* slc) -> double {

    int truth_idx = varMuonTruthIndex(slc);

    return GetMatchedRecoShowerIndex(slc, truth_idx);

  });

  const Var varTruthMuonMatchedShowerScore([](const caf::SRSliceProxy* slc) -> double {
    int muonShowerIndex = varTruthMuonMatchedShowerIndex(slc);
    if(muonShowerIndex>=0){
      return slc->reco.shw.at(muonShowerIndex).pfp.trackScore;
    }
    else{
      return -999;
    }

  });

  const Var varTruthProtonMatchedTrackScore([](const caf::SRSliceProxy* slc) -> double {
    int protonTrackIndex = varTruthProtonMatchedTrackIndex(slc);
    if(protonTrackIndex>=0){
      return slc->reco.trk.at(protonTrackIndex).pfp.trackScore;
    }
    else{
      return -999;
    }

  });

  const Var varTruthProtonMatchedTrackVertexDistance([](const caf::SRSliceProxy* slc) -> double {
    int protonTrackIndex = varTruthProtonMatchedTrackIndex(slc);
    if(protonTrackIndex>=0){
      auto const& trk = slc->reco.trk.at(protonTrackIndex);
      return sqrt( pow( slc->vertex.x - trk.start.x, 2.0 ) +
              pow( slc->vertex.y - trk.start.y, 2.0 ) +
              pow( slc->vertex.z - trk.start.z, 2.0 ) );

    }
    else{
      return -999;
    }

  });

  const Var varTruthProtonMatchedShowerIndex([](const caf::SRSliceProxy* slc) -> double {

    int truth_idx = varProtonTruthIndex(slc);

    return GetMatchedRecoShowerIndex(slc, truth_idx);

  });

  const Var varTruthProtonMatchedShowerScore([](const caf::SRSliceProxy* slc) -> double {
    int protonShowerIndex = varTruthProtonMatchedShowerIndex(slc);
    if(protonShowerIndex>=0){
      return slc->reco.shw.at(protonShowerIndex).pfp.trackScore;
    }
    else{
      return -999;
    }

  });

  //====   For a given true proton (truth_index), find a reco track whose best-matched is this muon

  const Var varTruthProtonMatchedTrackIndex([](const caf::SRSliceProxy* slc) -> double {

    int truth_idx = varProtonTruthIndex(slc);

    return GetMatchedRecoTrackIndex(slc, truth_idx);

  });

  const Var varTruthProtonMatchedTrackContainedness([](const caf::SRSliceProxy* slc) -> int {
    int ProtonTrackIndex = varTruthProtonMatchedTrackIndex(slc);
    //==== -1 : No track found
    //==== 0 : Exiting
    //==== 1 : Contained
    if(ProtonTrackIndex>=0){

      auto const& trk_Proton = slc->reco.trk.at(ProtonTrackIndex);
      bool isContained = fv.isContained(trk_Proton.end.x, trk_Proton.end.y, trk_Proton.end.z);
      if(isContained) return 1;
      else return 0;
    }
    else{
      return -1;
    }

  });

  const Var varTruthProtonMatchedTrackChi2Proton([](const caf::SRSliceProxy* slc) -> double {
    int protonTrackIndex = varTruthProtonMatchedTrackIndex(slc);
    if(protonTrackIndex>=0){

      auto const& trk = slc->reco.trk.at(protonTrackIndex);
      float Chi2Proton(-999.);
      if(trk.bestplane == 0){
        Chi2Proton = trk.chi2pid[0].chi2_proton;
      }
      else if(trk.bestplane == 1){
        Chi2Proton = trk.chi2pid[1].chi2_proton;
      }
      else{
        Chi2Proton = trk.chi2pid[2].chi2_proton;
      }
      if(isnan(Chi2Proton)) return -999.;
      else return Chi2Proton;
/*
      auto const& trk = slc->reco.trk.at(protonTrackIndex);
      int bp = trk.bestplane;
      float Chi2Proton = trk.chi2pid[bp].chi2_proton;
      if(isnan(Chi2Proton)) return -999.;
      else return Chi2Proton;
*/
    }
    else{
      return 999999;
    }
  });

  const Var varTruthProtonMatchedTrackChi2ProtonCollection([](const caf::SRSliceProxy* slc) -> double {
    int protonTrackIndex = varTruthProtonMatchedTrackIndex(slc);
    if(protonTrackIndex>=0){

      auto const& trk = slc->reco.trk.at(protonTrackIndex);
      float Chi2Proton = trk.chi2pid[2].chi2_proton;
      if(isnan(Chi2Proton)) return -999.;
      else return Chi2Proton;
/*
      auto const& trk = slc->reco.trk.at(protonTrackIndex);
      int bp = trk.bestplane;
      float Chi2Proton = trk.chi2pid[bp].chi2_proton;
      if(isnan(Chi2Proton)) return -999.;
      else return Chi2Proton;
*/
    }
    else{
      return 999999;
    }
  });

  const Var varTruthProtonMatchedTrackReducedChi2Proton([](const caf::SRSliceProxy* slc) -> double {
    int protonTrackIndex = varTruthProtonMatchedTrackIndex(slc);
    if(protonTrackIndex>=0){

      auto const& trk = slc->reco.trk.at(protonTrackIndex);
      int bp = trk.bestplane;
      float Chi2Proton = trk.chi2pid[bp].chi2_proton;
      int Chi2Ndof = trk.chi2pid[bp].pid_ndof;
      if( isnan(Chi2Proton) || isnan(Chi2Ndof) ) return -999.;
      else if(Chi2Ndof==0) return -999.;
      else return Chi2Proton/double(Chi2Ndof);

    }
    else{
      return 999999;
    }
  });

  const Var varTruthProtonMatchedTrackChi2Muon([](const caf::SRSliceProxy* slc) -> double {
    int protonTrackIndex = varTruthProtonMatchedTrackIndex(slc);
    if(protonTrackIndex>=0){

      auto const& trk = slc->reco.trk.at(protonTrackIndex);
      int bp = trk.bestplane;
      float Chi2Muon = trk.chi2pid[bp].chi2_muon;
      if(isnan(Chi2Muon)) return -999.;
      else return Chi2Muon;

    }
    else{
      return 999999;
    }

  });

  const Var varTruthProtonMatchedTrackReducedChi2Muon([](const caf::SRSliceProxy* slc) -> double {
    int protonTrackIndex = varTruthProtonMatchedTrackIndex(slc);
    if(protonTrackIndex>=0){

      auto const& trk = slc->reco.trk.at(protonTrackIndex);
      int bp = trk.bestplane;
      float Chi2Muon = trk.chi2pid[bp].chi2_muon;
      int Chi2Ndof = trk.chi2pid[bp].pid_ndof;
      if(isnan(Chi2Muon)) return -999.;
      else if(Chi2Ndof==0) return -999.;
      else return Chi2Muon/double(Chi2Ndof);

    }
    else{
      return 999999;
    }

  });

  const Var varTruthProtonMatchedTrackRangeP([](const caf::SRSliceProxy* slc) -> double {
    int ProtonTrackIndex = varTruthProtonMatchedTrackIndex(slc);
    if(ProtonTrackIndex>=0){
      auto const& trk_Proton = slc->reco.trk.at(ProtonTrackIndex);
      return trk_Proton.rangeP.p_proton;
    }
    else{
      return -999;
    }

  });
  const Var varTruthProtonMatchedTrackMCSP([](const caf::SRSliceProxy* slc) -> double {
    int ProtonTrackIndex = varTruthProtonMatchedTrackIndex(slc);
    if(ProtonTrackIndex>=0){
      auto const& trk_Proton = slc->reco.trk.at(ProtonTrackIndex);
      if(isnan(trk_Proton.mcsP.fwdP_proton)) return -999.;
      else return trk_Proton.mcsP.fwdP_proton;
    }
    else{
      return -999;
    }

  });
  const Var varTruthProtonMatchedTrackCombinedP([](const caf::SRSliceProxy* slc) -> double {
    int ProtonTrackIndex = varTruthProtonMatchedTrackIndex(slc);
    if(ProtonTrackIndex>=0){
      int mtc = varTruthProtonMatchedTrackContainedness(slc); // 0 : exiting, 1 : contained
      if(mtc==0) return varTruthProtonMatchedTrackMCSP(slc);
      else if(mtc==1) return varTruthProtonMatchedTrackRangeP(slc);
      else{
        cout << "[varTruthProtonMatchedTrackCombinedP] wtf?" << endl;
        exit(EXIT_FAILURE);
        return varTruthProtonMatchedTrackMCSP(slc); //
      }
    }
    else{
      return -999;
    }

  });

  const Var varTruthProtonMatchedTrackLength([](const caf::SRSliceProxy* slc) -> double {
    int ProtonTrackIndex = varTruthProtonMatchedTrackIndex(slc);
    if(ProtonTrackIndex>=0){
      auto const& trk_Proton = slc->reco.trk.at(ProtonTrackIndex);
      return trk_Proton.len;
    }
    else{
      return -999;
    }

  });



  //====   For a given true proton (truth_index), find a reco stub whose best-matched is this muon

  const Var varTruthProtonMatchedStubIndex([](const caf::SRSliceProxy* slc) -> double {

    int truth_idx = varProtonTruthIndex(slc);

    return GetMatchedRecoStubIndex(slc, truth_idx);

  });

  const Var varTruthProtonMatchedStubE([](const caf::SRSliceProxy* slc) -> double {
    int ProtonStubIndex = varTruthProtonMatchedStubIndex(slc);
    if(ProtonStubIndex>=0){
      auto const& stub_Proton = slc->reco.stub.at(ProtonStubIndex);

      for(unsigned ip=0; ip<stub_Proton.planes.size(); ip++){
        double sumQ(0.);
        for(unsigned ih=0; ih<stub_Proton.planes[ip].hits.size(); ih++){
          sumQ += stub_Proton.planes[ip].hits.at(ih).charge;
        }
      }
      return 1.;
    }
    else{
      return -999;
    }

  });


  //==== Reco variables

  //====   Muon


  //====   All muon-like tracks
  const MultiVar varAllMuonTrackIndices([](const caf::SRSliceProxy* slc) -> std::vector<double> {

    std::vector<double> out;

    //==== The (dis)qualification of a slice is based upon the track level information.
    float Atslc;
    float Chi2Proton, Chi2Muon;
    bool AtSlice, Contained, MaybeMuonExiting, MaybeMuonContained;

    for(std::size_t i(0); i < slc->reco.trk.size(); ++i){

      auto const& trk = slc->reco.trk.at(i);
      //==== First we calculate the distance of each track to the slice vertex.
      Atslc = sqrt( pow( slc->vertex.x - trk.start.x, 2.0 ) +
              pow( slc->vertex.y - trk.start.y, 2.0 ) +
              pow( slc->vertex.z - trk.start.z, 2.0 ) );
      //==== We require that the distance of the track from the slice is less than
      //==== 10 cm and that the parent of the track has been marked as the primary.
      AtSlice = ( Atslc < 10.0 );//&& trk.parent_is_primary);

      int bp = trk.bestplane;
      Chi2Proton = trk.chi2pid[bp].chi2_proton;
      Chi2Muon = trk.chi2pid[bp].chi2_muon;
      bool PassChi2 = Chi2Proton > 60 && Chi2Muon < 30; // MuonSel__MuonChi2LT30_and_ProtonChi2GT60

      Contained = fv.isContained(trk.end.x, trk.end.y, trk.end.z);

      MaybeMuonExiting = ( !Contained && trk.len > 100);
      MaybeMuonContained = ( Contained && PassChi2 && trk.len > 50. );

      if( AtSlice && ( MaybeMuonExiting || MaybeMuonContained ) ){
        out.push_back(i);
      }

    }

    return out;
  });

  //====   Longest one only
  const Var varMuonTrackInd([](const caf::SRSliceProxy* slc) -> int {

    //==== The (dis)qualification of a slice is based upon the track level information.
    float Atslc, Longest(0);
    float Chi2Proton, Chi2Muon;
    bool AtSlice, Contained, MaybeMuonExiting, MaybeMuonContained;
    int PTrackInd(-1);

    for(std::size_t i(0); i < slc->reco.trk.size(); ++i){

      auto const& trk = slc->reco.trk.at(i);
      //==== First we calculate the distance of each track to the slice vertex.
      Atslc = sqrt( pow( slc->vertex.x - trk.start.x, 2.0 ) +
              pow( slc->vertex.y - trk.start.y, 2.0 ) +
              pow( slc->vertex.z - trk.start.z, 2.0 ) );
      //===== We require that the distance of the track from the slice is less than
      //==== 10 cm and that the parent of the track has been marked as the primary.
      AtSlice = ( Atslc < 10.0 );//&& trk.parent_is_primary);

      int bp = trk.bestplane;
      Chi2Proton = trk.chi2pid[bp].chi2_proton;
      Chi2Muon = trk.chi2pid[bp].chi2_muon;
      bool PassChi2 = Chi2Proton > 60 && Chi2Muon < 30; // MuonSel__MuonChi2LT30_and_ProtonChi2GT60

      Contained = fv.isContained(trk.end.x, trk.end.y, trk.end.z);

      MaybeMuonExiting = ( !Contained && trk.len > 100);
      MaybeMuonContained = ( Contained && PassChi2 && trk.len > 50. );

      if( AtSlice && ( MaybeMuonExiting || MaybeMuonContained ) && trk.len > Longest ){
        Longest = trk.len;
        PTrackInd = i;
      }


    }

    return PTrackInd;

  });

  const Var varMuonRecoP([](const caf::SRSliceProxy* slc) -> float {
    float p(-5.f);
    bool Contained(false);
    
    if( varMuonTrackInd(slc) >= 0 ){
      auto const& trk = slc->reco.trk.at(varMuonTrackInd(slc));
      Contained = fv.isContained(trk.end.x, trk.end.y, trk.end.z);
      if(Contained){
        p = trk.rangeP.p_muon;
      }
      else{
        p = trk.mcsP.fwdP_muon;
      }
    }
    return p;
  });

  const Var varMuonCaloPlane0P([](const caf::SRSliceProxy* slc) -> float {
    float p(-5.f);

    if( varMuonTrackInd(slc) >= 0 ){
      auto const& trk = slc->reco.trk.at(varMuonTrackInd(slc));
      double KE = trk.calo[0].ke/1000.; // Note ke for now is in MeV
      double E_Muon = KE+M_MUON;
      p = sqrt( E_Muon*E_Muon - M_MUON*M_MUON );
    }
    return p;
  });

  const Var varMuonCaloPlane1P([](const caf::SRSliceProxy* slc) -> float {
    float p(-5.f);

    if( varMuonTrackInd(slc) >= 0 ){
      auto const& trk = slc->reco.trk.at(varMuonTrackInd(slc));
      double KE = trk.calo[1].ke/1000.; // Note ke for now is in MeV
      double E_Muon = KE+M_MUON;
      p = sqrt( E_Muon*E_Muon - M_MUON*M_MUON );
    }
    return p;
  });

  const Var varMuonCaloPlane2P([](const caf::SRSliceProxy* slc) -> float {
    float p(-5.f);

    if( varMuonTrackInd(slc) >= 0 ){
      auto const& trk = slc->reco.trk.at(varMuonTrackInd(slc));
      double KE = trk.calo[2].ke/1000.; // Note ke for now is in MeV
      double E_Muon = KE+M_MUON;
      p = sqrt( E_Muon*E_Muon - M_MUON*M_MUON );
    }
    return p;
  });

  const Var varMuonLength([](const caf::SRSliceProxy* slc) -> float {
    if( varMuonTrackInd(slc) >= 0 ){
      auto const& trk = slc->reco.trk.at(varMuonTrackInd(slc));
      return trk.len;
    }
    else{
      return -999.;
    }
  });

  const Var varMuonChi2Muon([](const caf::SRSliceProxy* slc) -> float {
    if( varMuonTrackInd(slc) >= 0 ){
      auto const& trk = slc->reco.trk.at(varMuonTrackInd(slc));
      int bp = trk.bestplane;
      return trk.chi2pid[bp].chi2_muon;
    }
    else{
      return -999.;
    }
  });

  const Var varMuonChi2Proton([](const caf::SRSliceProxy* slc) -> float {
    if( varMuonTrackInd(slc) >= 0 ){
      auto const& trk = slc->reco.trk.at(varMuonTrackInd(slc));
      int bp = trk.bestplane;
      return trk.chi2pid[bp].chi2_proton;
    }
    else{
      return -999.;
    }
  });

  const Var varMuonReducedChi2Muon([](const caf::SRSliceProxy* slc) -> float {
    if( varMuonTrackInd(slc) >= 0 ){
      auto const& trk = slc->reco.trk.at(varMuonTrackInd(slc));
      int bp = trk.bestplane;
      double Chi2Muon = trk.chi2pid[bp].chi2_muon;
      int ndof = trk.chi2pid[bp].pid_ndof;
      if(ndof==0) return -999.;
      else return Chi2Muon/double(ndof);
    }
    else{
      return -999.;
    }
  });

  const Var varMuonRecoStartX([](const caf::SRSliceProxy* slc) -> float {
    if( varMuonTrackInd(slc) >= 0 ){
      auto const& trk = slc->reco.trk.at(varMuonTrackInd(slc));
      return trk.start.x;
    }
    return -9999999999.;
  });
  const Var varMuonRecoStartY([](const caf::SRSliceProxy* slc) -> float {
    if( varMuonTrackInd(slc) >= 0 ){
      auto const& trk = slc->reco.trk.at(varMuonTrackInd(slc));
      return trk.start.y;
    }
    return -9999999999.;
  });
  const Var varMuonRecoStartZ([](const caf::SRSliceProxy* slc) -> float {
    if( varMuonTrackInd(slc) >= 0 ){
      auto const& trk = slc->reco.trk.at(varMuonTrackInd(slc));
      return trk.start.z;
    }
    return -9999999999.;
  });

  const Var varMuonTrackFromVertex([](const caf::SRSliceProxy* slc) -> float {
    if( varMuonTrackInd(slc) >= 0 ){

      auto const& trk = slc->reco.trk.at(varMuonTrackInd(slc));
      return sqrt( pow( slc->vertex.x - trk.start.x, 2.0 ) +
              pow( slc->vertex.y - trk.start.y, 2.0 ) +
              pow( slc->vertex.z - trk.start.z, 2.0 ) );

    }
    return -9999999999.;
  });

  //==== Start from a reco Track, and look at its best match gen-particle.
  //==== This means that the get-particle may not be a true muon
  const Var varMuonBestmatchP([](const caf::SRSliceProxy* slc) -> float {
    float p(-5.f);
    
    if( varMuonTrackInd(slc) >= 0 ){
      auto const& trk = slc->reco.trk.at(varMuonTrackInd(slc));
      p = sqrt( pow(trk.truth.p.genp.x, 2) + 
                pow(trk.truth.p.genp.y, 2) + 
                pow(trk.truth.p.genp.z, 2) );
    }
    return p;
  });

  const Var varMuonBestmatchPDG([](const caf::SRSliceProxy* slc) -> int {
    int pdg(-9999999);

    if( varMuonTrackInd(slc) >= 0 ){
      auto const& trk = slc->reco.trk.at(varMuonTrackInd(slc));
      pdg = trk.truth.p.pdg;
    }
    return pdg;
  });

  const Var varMuonPResidual([](const caf::SRSliceProxy* slc) -> float {
    float pr(-999);

    if( varMuonTrackInd(slc) >= 0 ){
      double RecoP = varMuonRecoP(slc);
      double TrueP = varMuonBestmatchP(slc);
      pr = RecoP-TrueP;
    }
    return pr;
  });

  const Var varMuonPResidualFraction([](const caf::SRSliceProxy* slc) -> float {
    float prf(-999);

    if( varMuonTrackInd(slc) >= 0 ){
      double RecoP = varMuonRecoP(slc);
      double TrueP = varMuonBestmatchP(slc);
      if(TrueP>0.){
        prf = (RecoP-TrueP)/TrueP;
      }
    }
    return prf;
  });

  const Var varMuonCaloPlane0PResidualFraction([](const caf::SRSliceProxy* slc) -> float {
    float prf(-999);

    if( varMuonTrackInd(slc) >= 0 ){
      double RecoP = varMuonCaloPlane0P(slc);
      double TrueP = varMuonBestmatchP(slc);
      if(TrueP>0.){
        prf = (RecoP-TrueP)/TrueP;
      }
    }
    return prf;
  });

  const Var varMuonCaloPlane1PResidualFraction([](const caf::SRSliceProxy* slc) -> float {
    float prf(-999);

    if( varMuonTrackInd(slc) >= 0 ){
      double RecoP = varMuonCaloPlane1P(slc);
      double TrueP = varMuonBestmatchP(slc);
      if(TrueP>0.){
        prf = (RecoP-TrueP)/TrueP;
      }
    }
    return prf;
  });

  const Var varMuonCaloPlane2PResidualFraction([](const caf::SRSliceProxy* slc) -> float {
    float prf(-999);

    if( varMuonTrackInd(slc) >= 0 ){
      double RecoP = varMuonCaloPlane2P(slc);
      double TrueP = varMuonBestmatchP(slc);
      if(TrueP>0.){
        prf = (RecoP-TrueP)/TrueP;
      }
    }
    return prf;
  });

  const Var varMuonRecoCosineTheta([](const caf::SRSliceProxy* slc) -> float {
    float costh(-5.f);

    if( varMuonTrackInd(slc) >= 0 ){
      auto const& trk = slc->reco.trk.at(varMuonTrackInd(slc));
      costh = trk.costh;
    }
    return costh;
  });

  const Var varMuonRecoNuMICosineTheta([](const caf::SRSliceProxy* slc) -> float {
    float costh(-5.f);

    if( varMuonTrackInd(slc) >= 0 ){
      auto const& trk = slc->reco.trk.at(varMuonTrackInd(slc));
      TVector3 v3(trk.dir.x, trk.dir.y, trk.dir.z);
      double angleNuMI = v3.Angle(NuDirection_NuMI);
      costh = TMath::Cos(angleNuMI);
    }
    return costh;
  });

  const Var varMuonBestmatchCosineTheta([](const caf::SRSliceProxy* slc) -> float {
    float costh(-5.f);

    if( varMuonTrackInd(slc) >= 0 ){
      auto const& trk = slc->reco.trk.at(varMuonTrackInd(slc));
      TVector3 v3(trk.truth.p.genp.x, trk.truth.p.genp.y, trk.truth.p.genp.z);
      costh = v3.CosTheta();
    }
    return costh;
  });

  const Var varMuonCosineThetaResidual([](const caf::SRSliceProxy* slc) -> float {
    float costhr(-999);

    if( varMuonTrackInd(slc) >= 0 ){
      double RecoCosineTheta = varMuonRecoCosineTheta(slc);
      double TrueCosineTheta = varMuonBestmatchCosineTheta(slc);
      costhr = RecoCosineTheta-TrueCosineTheta;
    }
    return costhr;
  });

  const Var varMuonCosineThetaResidualFraction([](const caf::SRSliceProxy* slc) -> float {
    float costhrf(-999);

    if( varMuonTrackInd(slc) >= 0 ){
      double RecoCosineTheta = varMuonRecoCosineTheta(slc);
      double TrueCosineTheta = varMuonBestmatchCosineTheta(slc);
      costhrf = (RecoCosineTheta-TrueCosineTheta)/TrueCosineTheta;
    }
    return costhrf;
  });

  //==== Proton

  //====   Track-based

  const Var varNProtonCandTrack([](const caf::SRSliceProxy* slc) -> int {
    const int muTrackInd = varMuonTrackInd(slc);

    //==== The (dis)qualification of a slice is based upon the track level information.
    int nProtonTrackCand(0);

    for(std::size_t i(0); i < slc->reco.trk.size(); ++i){

      if((int)i==muTrackInd) continue;

      auto const& trk = slc->reco.trk.at(i);
      bool Contained = fv.isContained(trk.end.x, trk.end.y, trk.end.z);

      if(Contained){
        nProtonTrackCand++;
      }

    }

    return nProtonTrackCand;

  });

  const Var varNProtonCandMatched([](const caf::SRSliceProxy* slc) -> int {
    const int muTrackInd = varMuonTrackInd(slc);

    //==== The (dis)qualification of a slice is based upon the track level information.
    int nProtonTrackCand(0);

    for(std::size_t i(0); i < slc->reco.trk.size(); ++i){

      if((int)i==muTrackInd) continue;

      auto const& trk = slc->reco.trk.at(i);
      bool Contained = fv.isContained(trk.end.x, trk.end.y, trk.end.z);

      if(Contained){
        if(abs(trk.truth.p.pdg)==2212) nProtonTrackCand++;
      }

    }

    return nProtonTrackCand;

  });

  const Var varProtonTrackInd([](const caf::SRSliceProxy* slc) -> int {

    //==== The (dis)qualification of a slice is based upon the track level information.
    float Atslc, Longest(0);
    float Chi2Proton;//, Chi2Muon;
    bool AtSlice, Contained;
    int PTrackInd(-1);

    const int muTrackInd = varMuonTrackInd(slc);

    for(std::size_t i(0); i < slc->reco.trk.size(); ++i){

      if((int)i==muTrackInd) continue;

      auto const& trk = slc->reco.trk.at(i);
      //==== First we calculate the distance of each track to the slice vertex.
      Atslc = sqrt( pow( slc->vertex.x - trk.start.x, 2.0 ) +
              pow( slc->vertex.y - trk.start.y, 2.0 ) +
              pow( slc->vertex.z - trk.start.z, 2.0 ) );
      //===== We require that the distance of the track from the slice is less than
      //==== 10 cm and that the parent of the track has been marked as the primary.
      AtSlice = ( Atslc < 10.0 );//&& trk.parent_is_primary);

      Contained = fv.isContained(trk.end.x, trk.end.y, trk.end.z);

      int bp = trk.bestplane;
      Chi2Proton = trk.chi2pid[bp].chi2_proton;
      //Chi2Muon = trk.chi2pid[bp].chi2_muon;
      bool PassChi2 = Chi2Proton < 100;

      //==== TEST
      PassChi2 = true; // not appling chi2

      if( AtSlice && Contained && PassChi2 && trk.len > Longest ){
        Longest = trk.len;
        PTrackInd = i;
      }

    }

    return PTrackInd;

  });

  const Var varProtonCaloP([](const caf::SRSliceProxy* slc) -> float {
    float p(-5.f);

    if( varProtonTrackInd(slc) >= 0 ){
      auto const& trk = slc->reco.trk.at(varProtonTrackInd(slc));
      int bp = trk.bestplane;
      double best_calo = trk.calo[bp].ke;
      best_calo = best_calo/1000.; // Note ke for now is in MeV
      double E_Proton = best_calo+M_PROTON;
      p = sqrt( E_Proton*E_Proton - M_PROTON*M_PROTON );
    }
    return p;

  });

  const Var varProtonRecoP([](const caf::SRSliceProxy* slc) -> float {
    float p(-5.f);
    bool Contained(false);
    
    if( varProtonTrackInd(slc) >= 0 ){
      auto const& trk = slc->reco.trk.at(varProtonTrackInd(slc));

      Contained = fv.isContained(trk.end.x, trk.end.y, trk.end.z);
      if(Contained){
        p = trk.rangeP.p_proton;
      }
      else{
        p = trk.mcsP.fwdP_proton;
      }

    }
    return p;
  });

  const Var varProtonLength([](const caf::SRSliceProxy* slc) -> float {
    if( varProtonTrackInd(slc) >= 0 ){
      auto const& trk = slc->reco.trk.at(varProtonTrackInd(slc));
      return trk.len;
    }
    else{
      return -999.;
    }
  });

  const Var varProtonChi2Proton([](const caf::SRSliceProxy* slc) -> float {
    if( varProtonTrackInd(slc) >= 0 ){
      auto const& trk = slc->reco.trk.at(varProtonTrackInd(slc));
      int bp = trk.bestplane;
      return trk.chi2pid[bp].chi2_proton;
    }
    else{
      return -999.;
    }
  });

  const Var varProtonReducedChi2Proton([](const caf::SRSliceProxy* slc) -> float {
    if( varProtonTrackInd(slc) >= 0 ){
      auto const& trk = slc->reco.trk.at(varProtonTrackInd(slc));
      int bp = trk.bestplane;
      double Chi2Proton = trk.chi2pid[bp].chi2_proton;
      int ndof = trk.chi2pid[bp].pid_ndof;
      if(ndof==0) return -999.;
      else return Chi2Proton/double(ndof);

    }
    else{
      return -999.;
    }
  });

  const Var varProtonBestmatchP([](const caf::SRSliceProxy* slc) -> float {
    float p(-5.f);
    
    if( varProtonTrackInd(slc) >= 0 ){
      auto const& trk = slc->reco.trk.at(varProtonTrackInd(slc));
      if( !isnan(trk.truth.p.genp.x) ){
        p = sqrt( pow(trk.truth.p.genp.x, 2) + 
                  pow(trk.truth.p.genp.y, 2) + 
                  pow(trk.truth.p.genp.z, 2) );
      }
    }
    return p;
  });

  const Var varProtonBestmatchPDG([](const caf::SRSliceProxy* slc) -> int {
    int pdg(-9999999);

    if( varProtonTrackInd(slc) >= 0 ){
      auto const& trk = slc->reco.trk.at(varProtonTrackInd(slc));
      pdg = trk.truth.p.pdg;
    }
    return pdg;
  });

  const Var varProtonPResidual([](const caf::SRSliceProxy* slc) -> float {
    float pr(-999);

    if( varProtonTrackInd(slc) >= 0 ){
      double RecoP = varProtonRecoP(slc);
      double TrueP = varProtonBestmatchP(slc);
      pr = RecoP-TrueP;
    }
    return pr;
  });

  const Var varProtonPResidualFraction([](const caf::SRSliceProxy* slc) -> float {
    float prf(-999);

    if( varProtonTrackInd(slc) >= 0 ){
      double RecoP = varProtonRecoP(slc);
      double TrueP = varProtonBestmatchP(slc);
      if(TrueP>0.){
        prf = (RecoP-TrueP)/TrueP;
      }
    }
    return prf;
  });

  const Var varProtonRecoCosineTheta([](const caf::SRSliceProxy* slc) -> float {
    float costh(-5.f);

    if( varProtonTrackInd(slc) >= 0 ){
      auto const& trk = slc->reco.trk.at(varProtonTrackInd(slc));
      costh = trk.costh;
    }
    return costh;
  });

  const Var varProtonRecoNuMICosineTheta([](const caf::SRSliceProxy* slc) -> float {
    float costh(-5.f);

    if( varProtonTrackInd(slc) >= 0 ){
      auto const& trk = slc->reco.trk.at(varProtonTrackInd(slc));
      TVector3 v3(trk.dir.x, trk.dir.y, trk.dir.z);
      double angleNuMI = v3.Angle(NuDirection_NuMI);
      costh = TMath::Cos(angleNuMI);
    }
    return costh;
  });

  const Var varProtonBestmatchCosineTheta([](const caf::SRSliceProxy* slc) -> float {
    float costh(-5.f);

    if( varProtonTrackInd(slc) >= 0 ){
      auto const& trk = slc->reco.trk.at(varProtonTrackInd(slc));
      if( !isnan(trk.truth.p.genp.x) ){
        TVector3 v3(trk.truth.p.genp.x, trk.truth.p.genp.y, trk.truth.p.genp.z);
        costh = v3.CosTheta();
      }
    }
    return costh;
  });

  const Var varProtonCosineThetaResidual([](const caf::SRSliceProxy* slc) -> float {

    float costhr(-999);

    if( varProtonTrackInd(slc) >= 0 ){
      double RecoCosineTheta = varProtonRecoCosineTheta(slc);
      double TrueCosineTheta = varProtonBestmatchCosineTheta(slc);
      costhr = RecoCosineTheta-TrueCosineTheta;
    }
    return costhr;
  });

  const Var varProtonCosineThetaResidualFraction([](const caf::SRSliceProxy* slc) -> float {

    float costhrf(-999);

    if( varProtonTrackInd(slc) >= 0 ){
      double RecoCosineTheta = varProtonRecoCosineTheta(slc);
      double TrueCosineTheta = varProtonBestmatchCosineTheta(slc);
      costhrf = (RecoCosineTheta-TrueCosineTheta)/TrueCosineTheta;
    }
    return costhrf;
  });

  //====   Stub-based

  const Var varNStub([](const caf::SRSliceProxy* slc) -> int {

    return slc->reco.stub.size();

  });

  //==== Neutrino

  const Var varNeutrinoCombinedEnergy([](const caf::SRSliceProxy* slc) -> double {

    double E_Nu(-5.f);
    if( varMuonTrackInd(slc) >= 0 && varProtonTrackInd(slc) >= 0 ){

      double P_Muon = varMuonRecoP(slc); // comdined
      double E_Muon = sqrt( P_Muon*P_Muon + M_MUON*M_MUON );

      double P_Proton = varProtonRecoP(slc);
      double E_Proton = sqrt( P_Proton*P_Proton + M_PROTON*M_PROTON );
      double KE_Proton = E_Proton - M_PROTON;

      E_Nu = E_Muon + KE_Proton + E_EffNuclB;

    }
    return E_Nu;

  });

  const Var varNeutrinoCombinedEnergyResidual([](const caf::SRSliceProxy* slc) -> double {

    double E_Nur(-999);
    if( varMuonTrackInd(slc) >= 0 && varProtonTrackInd(slc) >= 0 ){

      double RecoE = varNeutrinoCombinedEnergy(slc);
      double TrueE = varNeutrinoTruthE(slc);

      E_Nur = RecoE-TrueE;

    }
    return E_Nur;

  });

  const Var varNeutrinoCombinedEnergyResidualFraction([](const caf::SRSliceProxy* slc) -> double {

    double E_Nurf(-999);
    if( varMuonTrackInd(slc) >= 0 && varProtonTrackInd(slc) >= 0 ){

      double RecoE = varNeutrinoCombinedEnergy(slc);
      double TrueE = varNeutrinoTruthE(slc);
      if(TrueE>0.){
        E_Nurf = (RecoE-TrueE)/TrueE;
      }

    }
    return E_Nurf;

  });

  //==== https://s3.cern.ch/inspire-prod-files-9/93642a13c46438d97680971700e2013c
  const Var varNeutrinoQE([](const caf::SRSliceProxy* slc) -> double {

    int muonTrackIndex = varMuonTrackInd(slc);
    if(muonTrackIndex>=0){
      auto const& trk_Muon = slc->reco.trk.at(muonTrackIndex);

      double mu_costh = trk_Muon.costh;
      double mu_p = varMuonRecoP(slc);
      double mu_E = sqrt(mu_p*mu_p+M_MUON*M_MUON);

      double EQE_num = M_PROTON*M_PROTON - (M_NEUTRON-E_EffNuclB)*(M_NEUTRON-E_EffNuclB) - M_MUON*M_MUON + 2.*(M_NEUTRON-E_EffNuclB)*mu_E;
      double EQE_den = 2.*(M_NEUTRON - E_EffNuclB - mu_E + mu_p * mu_costh);

      return EQE_num/EQE_den;

    }
    else{
      return -999.;
    }

  });

  const Var varNeutrinoQEResidual([](const caf::SRSliceProxy* slc) -> double {

    int muonTrackIndex = varMuonTrackInd(slc);
    if(muonTrackIndex>=0){

      double RecoE = varNeutrinoQE(slc);
      double TrueE = varNeutrinoTruthE(slc);

      double E_Nur = RecoE-TrueE;

      return E_Nur;

    }
    else{
      return -999.;
    }

  });

  const Var varNeutrinoQEResidualFraction([](const caf::SRSliceProxy* slc) -> double {

    int muonTrackIndex = varMuonTrackInd(slc);
    double E_Nurf(-999.);
    if(muonTrackIndex>=0){

      double RecoE = varNeutrinoQE(slc);
      double TrueE = varNeutrinoTruthE(slc);
      if(TrueE>0.){
        E_Nurf = (RecoE-TrueE)/TrueE;
      }

    }

    return E_Nurf;

  });

  const Var varNeutrinoTestEnergy([](const caf::SRSliceProxy* slc) -> double {

    double E_Nu(-5.f);
    if( varMuonTrackInd(slc) >= 0 && varProtonTrackInd(slc) >= 0 ){

      double P_Muon = varMuonRecoP(slc); // comdined
      double E_Muon = sqrt( P_Muon*P_Muon + M_MUON*M_MUON );

      double P_Proton = varProtonCaloP(slc);
      double E_Proton = sqrt( P_Proton*P_Proton + M_PROTON*M_PROTON );
      double KE_Proton = E_Proton - M_PROTON;

      E_Nu = E_Muon + KE_Proton + E_EffNuclB;

    }
    return E_Nu;

  });


}

