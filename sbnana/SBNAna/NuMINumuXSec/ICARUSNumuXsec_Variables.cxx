#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_Variables.h"

using namespace std;
using namespace ana;

namespace ICARUSNumuXsec{

  //==== Spill variable

  const SpillVar spillvarCountSpill([](const caf::SRSpillProxy *sr) ->int {
    return 0.;
  });

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

  const SpillMultiVar spillvarMCNeutrinoE([](const caf::SRSpillProxy *sr)
  {
    std::vector<double> rets;
    for(const auto& _nu : sr->mc.nu){
      rets.push_back(_nu.E);
    }
    return rets;
  });

  const SpillVar spillvarNSlice([](const caf::SRSpillProxy *sr)
  {
    return sr->nslc;
  });
  const SpillVar spillvarNTrack([](const caf::SRSpillProxy *sr)
  { 
    return sr->reco.ntrk;
  });
  const SpillMultiVar spillNuDirectionX([](const caf::SRSpillProxy *sr)
  {
    std::vector<double> rets;
    for(const auto& _nu : sr->mc.nu){
      double this_x = _nu.prod_vtx.x/100.;
      double this_y = _nu.prod_vtx.y/100.;
      double this_z = _nu.prod_vtx.z/100.;
      TVector3 this_coord = nct.GetICARUSCoord(this_x, this_y, this_z);
      rets.push_back(this_coord.X());
    }
    return rets;
  });
  const SpillMultiVar spillNuDirectionY([](const caf::SRSpillProxy *sr)
  {
    std::vector<double> rets;
    for(const auto& _nu : sr->mc.nu){
      double this_x = _nu.prod_vtx.x/100.;
      double this_y = _nu.prod_vtx.y/100.;
      double this_z = _nu.prod_vtx.z/100.;
      TVector3 this_coord = nct.GetICARUSCoord(this_x, this_y, this_z);
      rets.push_back(this_coord.Y());
    }
    return rets;
  });
  const SpillMultiVar spillNuDirectionZ([](const caf::SRSpillProxy *sr)
  {
    std::vector<double> rets;
    for(const auto& _nu : sr->mc.nu){
      double this_x = _nu.prod_vtx.x/100.;
      double this_y = _nu.prod_vtx.y/100.;
      double this_z = _nu.prod_vtx.z/100.;
      TVector3 this_coord = nct.GetICARUSCoord(this_x, this_y, this_z);
      rets.push_back(this_coord.Z());
    }
    return rets;
  });

  //====   Spill-based Truth neutrino information
  const SpillMultiVar spillvarNeutrinoQEEnergyResidual([](const caf::SRSpillProxy *sr)
  {
    std::vector<double> rets;
    for(const auto& _nu : sr->mc.nu){

      double max_E(-999);
      int truth_idx(-1);
      for(std::size_t i(0); i < _nu.prim.size(); ++i){
        if( abs(_nu.prim.at(i).pdg)== 13 ){
          if(isnan(_nu.prim.at(i).genE)) continue;
          double this_E = _nu.prim.at(i).genE;
          if(this_E>max_E){
            max_E = this_E;
            truth_idx = i;
          }
        }
      }

      if(truth_idx>=0){

        auto const& pMuon = _nu.prim.at(truth_idx);

        TVector3 pMuon_Momentum(pMuon.genp.x, pMuon.genp.y, pMuon.genp.z);

        double mu_costh = pMuon_Momentum.CosTheta();
        double mu_p = pMuon_Momentum.Mag();
        double mu_E = sqrt(mu_p*mu_p+M_MUON*M_MUON);

        double EQE_num = M_PROTON*M_PROTON - (M_NEUTRON-E_EffNuclB)*(M_NEUTRON-E_EffNuclB) - M_MUON*M_MUON + 2.*(M_NEUTRON-E_EffNuclB)*mu_E;
        double EQE_den = 2.*(M_NEUTRON - E_EffNuclB - mu_E + mu_p * mu_costh);

        double EQE = EQE_num/EQE_den;
        rets.push_back( EQE-_nu.E );

      }

    }
    return rets;
  });

  const SpillVar spillTEST([](const caf::SRSpillProxy *sr)
  {
/*
    double twg = sr->hdr.triggerinfo.trigger_within_gate;
    std::cout << "[spillTEST] twg = " << twg << std::endl;
    return twg;
*/
    double a(1.);
    for(const auto& it: sr->hdr.triggerinfo){
      std::cout << "[spillTEST] time = " << it.global_trigger_det_time << std::endl;
    }
    return a;
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

  const MultiVar varAllTrackEndPositionX([](const caf::SRSliceProxy* slc) -> std::vector<double> {
    std::vector<double> out;
    for (auto const& trk : slc->reco.trk) {
      out.push_back(trk.end.x);
    }
    return out;
  });

  const MultiVar varAllTrackEndPositionY([](const caf::SRSliceProxy* slc) -> std::vector<double> {
    std::vector<double> out;
    for (auto const& trk : slc->reco.trk) {
      out.push_back(trk.end.y);
    }
    return out;
  });

  const MultiVar varAllTrackEndPositionZ([](const caf::SRSliceProxy* slc) -> std::vector<double> {
    std::vector<double> out;
    for (auto const& trk : slc->reco.trk) {
      out.push_back(trk.end.z);
    }
    return out;
  });
  const MultiVar varAllTrackDirectionX([](const caf::SRSliceProxy* slc) -> std::vector<double> {
    std::vector<double> out;
    for (auto const& trk : slc->reco.trk) {
      out.push_back(trk.dir.x);
    }
    return out;
  });
  const MultiVar varAllTrackDirectionY([](const caf::SRSliceProxy* slc) -> std::vector<double> {
    std::vector<double> out;
    for (auto const& trk : slc->reco.trk) {
      out.push_back(trk.dir.y);
    }
    return out;
  });
  const MultiVar varAllTrackDirectionZ([](const caf::SRSliceProxy* slc) -> std::vector<double> {
    std::vector<double> out;
    for (auto const& trk : slc->reco.trk) {
      out.push_back(trk.dir.z);
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
  //====   Scattering
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
  //====   Vertex
  const Var varTruthVtxX([](const caf::SRSliceProxy* slc) -> double {
    if(isnan(slc->truth.position.x)) return -999999.;
    else{
      return slc->truth.position.x;
    }
  });
  const Var varTruthVtxY([](const caf::SRSliceProxy* slc) -> double {
    if(isnan(slc->truth.position.y)) return -999999.;
    else{
      return slc->truth.position.y;
    }
  });
  const Var varTruthVtxZ([](const caf::SRSliceProxy* slc) -> double {
    if(isnan(slc->truth.position.z)) return -999999.;
    else{
      return slc->truth.position.z;
    }
  });
  //====     Vertex residual
  const Var varVtxResidualX([](const caf::SRSliceProxy* slc) -> double {
    if(isnan(slc->vertex.x) || isnan(slc->truth.position.x)){
      return -999999999.;
    }
    else{
      return slc->vertex.x-slc->truth.position.x;
    }
  });
  const Var varVtxResidualY([](const caf::SRSliceProxy* slc) -> double {
    if(isnan(slc->vertex.y) || isnan(slc->truth.position.y)){
      return -999999999.;
    }
    else{
      return slc->vertex.y-slc->truth.position.y;
    }
  });
  const Var varVtxResidualZ([](const caf::SRSliceProxy* slc) -> double {
    if(isnan(slc->vertex.z) || isnan(slc->truth.position.z)){
      return -999999999.;
    }
    else{
      return slc->vertex.z-slc->truth.position.z;
    }
  });


  //====   Number of particles
  const Var varTruthNNeutron = SIMPLEVAR(truth.nneutron);
  const Var varTruthNPiMinus = SIMPLEVAR(truth.npiminus);
  const Var varTruthNPiPlus = SIMPLEVAR(truth.npiplus);
  const Var varTruthNChargedPion = varTruthNPiMinus+varTruthNPiPlus;
  const Var varTruthNPiZero = SIMPLEVAR(truth.npizero);
  const Var varTruthNProton = SIMPLEVAR(truth.nproton);

  //====   Truth Muon

  const Var varMuonTruthIndex([](const caf::SRSliceProxy* slc) -> int {
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

  const Var varMuonTruthLength([](const caf::SRSliceProxy* slc) -> double {

    int truth_idx = varMuonTruthIndex(slc);
    if(truth_idx>=0){
      return slc->truth.prim.at(truth_idx).length;
    }
    else{
      return -999.;
    }

  });

  //====   Truth Proton

  const Var varProtonTruthIndex([](const caf::SRSliceProxy* slc) -> int {

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

  //====   Truth Charged Pion

  const Var varChargedPionTruthIndex([](const caf::SRSliceProxy* slc) -> int {

    double max_E(-999);
    int truth_idx(-1);
    for(std::size_t i(0); i < slc->truth.prim.size(); ++i){
      if( abs(slc->truth.prim.at(i).pdg)== 211 ){
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

  const Var varChargedPionTruthP([](const caf::SRSliceProxy* slc) -> double {

    int truth_idx = varChargedPionTruthIndex(slc);
    if(truth_idx>=0){
      double genE = slc->truth.prim.at(truth_idx).genE;
      double genP2 = genE*genE-M_CHARGEDPION*M_CHARGEDPION;
      if(genP2>=0) return sqrt(genP2);
      else return 0;
    }
    else{
      return -999.;
    }

  });

  const Var varChargedPionTruthT([](const caf::SRSliceProxy* slc) -> double {

    int truth_idx = varChargedPionTruthIndex(slc);
    if(truth_idx>=0){
      double genE = slc->truth.prim.at(truth_idx).genE;
      return (genE-M_CHARGEDPION);
    }
    else{
      return -999.;
    }

  });

  const Var varChargedPionTruthLength([](const caf::SRSliceProxy* slc) -> double {
    
    int truth_idx = varChargedPionTruthIndex(slc);
    if(truth_idx>=0){
      return slc->truth.prim.at(truth_idx).length;
    }
    else{
      return -999.;
    }

  });

  //====   Truth Muon+Proton

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

  //==== For a given true muon (truth_index), find a reco whose best-matched is this muon

  const Var varTruthMuonMatchedTrackIndex([](const caf::SRSliceProxy* slc) -> int {

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
      bool isContained = fv_track.isContained(trk.end.x, trk.end.y, trk.end.z);
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

  const Var varTruthMuonMatchedTrackChi2Pion([](const caf::SRSliceProxy* slc) -> double {
    int muonTrackIndex = varTruthMuonMatchedTrackIndex(slc);
    if(muonTrackIndex>=0){

      auto const& trk = slc->reco.trk.at(muonTrackIndex);
      int bp = trk.bestplane;
      float Chi2Pion = trk.chi2pid[bp].chi2_pion;
      if(isnan(Chi2Pion)) return -999.;
      else return Chi2Pion;

    }
    else{
      return 999999;
    }

  });

  const Var varTruthMuonMatchedTrackReducedChi2Pion([](const caf::SRSliceProxy* slc) -> double {
    int muonTrackIndex = varTruthMuonMatchedTrackIndex(slc);
    if(muonTrackIndex>=0){

      auto const& trk = slc->reco.trk.at(muonTrackIndex);
      int bp = trk.bestplane;
      float Chi2Pion = trk.chi2pid[bp].chi2_pion;
      int Chi2Ndof = trk.chi2pid[bp].pid_ndof;
      if(isnan(Chi2Pion)) return -999.;
      else if(Chi2Ndof==0) return -999.;
      else return Chi2Pion/double(Chi2Ndof);

    }
    else{
      return 999999;
    }

  });

  //====   Matched reco track positions
  const Var varTruthMuonMatchedTrackEndPositionX([](const caf::SRSliceProxy* slc) -> double {
    int muonTrackIndex = varTruthMuonMatchedTrackIndex(slc);
    if(muonTrackIndex>=0){
      auto const& trk = slc->reco.trk.at(muonTrackIndex);
      return trk.end.x;
    }
    else{
      return -999;
    }

  });
  const Var varTruthMuonMatchedTrackEndPositionY([](const caf::SRSliceProxy* slc) -> double {
    int muonTrackIndex = varTruthMuonMatchedTrackIndex(slc);
    if(muonTrackIndex>=0){
      auto const& trk = slc->reco.trk.at(muonTrackIndex);
      return trk.end.y;
    }
    else{
      return -999;
    }

  });
  const Var varTruthMuonMatchedTrackEndPositionZ([](const caf::SRSliceProxy* slc) -> double {
    int muonTrackIndex = varTruthMuonMatchedTrackIndex(slc);
    if(muonTrackIndex>=0){
      auto const& trk = slc->reco.trk.at(muonTrackIndex);
      return trk.end.z;
    }
    else{
      return -999;
    }

  });

  //====   Matched reco angle

  const Var varTruthMuonMatchedTrackNuMICosineTheta([](const caf::SRSliceProxy* slc) -> double {
    float costh(-5.f);
    int muonTrackIndex = varTruthMuonMatchedTrackIndex(slc);
    if(muonTrackIndex>=0){
      auto const& trk = slc->reco.trk.at(muonTrackIndex);
      TVector3 v3(trk.dir.x, trk.dir.y, trk.dir.z);
      double angleNuMI = v3.Angle(NuDirection_NuMI);
      costh = TMath::Cos(angleNuMI);
    }
    return costh;
  });

  //====   Matched reco length

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

  //====   Matched reco momentum

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

  //====   Matched reco momentum residual

  const Var varTruthMuonMatchedTrackRangePResidual([](const caf::SRSliceProxy* slc) -> double {
    int muonTrackIndex = varTruthMuonMatchedTrackIndex(slc);
    if(muonTrackIndex>=0){
      double TruthP = varMuonTruthP(slc);
      double RecoP = varTruthMuonMatchedTrackRangeP(slc);
      return (RecoP-TruthP);
    }
    else{
      return -999;
    }

  });

  const Var varTruthMuonMatchedTrackMCSPResidual([](const caf::SRSliceProxy* slc) -> double {
    int muonTrackIndex = varTruthMuonMatchedTrackIndex(slc);
    if(muonTrackIndex>=0){
      double TruthP = varMuonTruthP(slc);
      double RecoP = varTruthMuonMatchedTrackMCSP(slc);
      return (RecoP-TruthP);
    }
    else{
      return -999;
    }

  });
  const Var varTruthMuonMatchedTrackCombinedPResidual([](const caf::SRSliceProxy* slc) -> double {
    int muonTrackIndex = varTruthMuonMatchedTrackIndex(slc);
    if(muonTrackIndex>=0){
      double TruthP = varMuonTruthP(slc);
      double RecoP = varTruthMuonMatchedTrackCombinedP(slc);
      return (RecoP-TruthP);
    }
    else{
      return -999;
    }

  });

  //====   Matched reco momentum residual fraction

  const Var varTruthMuonMatchedTrackRangePResidualFraction([](const caf::SRSliceProxy* slc) -> double {
    int muonTrackIndex = varTruthMuonMatchedTrackIndex(slc);
    if(muonTrackIndex>=0){
      double TruthP = varMuonTruthP(slc);
      double RecoP = varTruthMuonMatchedTrackRangeP(slc);
      return (RecoP-TruthP)/TruthP;
    }
    else{
      return -999;
    }

  });

  const Var varTruthMuonMatchedTrackMCSPResidualFraction([](const caf::SRSliceProxy* slc) -> double {
    int muonTrackIndex = varTruthMuonMatchedTrackIndex(slc);
    if(muonTrackIndex>=0){
      double TruthP = varMuonTruthP(slc);
      double RecoP = varTruthMuonMatchedTrackMCSP(slc);
      return (RecoP-TruthP)/TruthP;
    }
    else{
      return -999;
    }

  });

  const Var varTruthMuonMatchedTrackCombinedPResidualFraction([](const caf::SRSliceProxy* slc) -> double {
    int muonTrackIndex = varTruthMuonMatchedTrackIndex(slc);
    if(muonTrackIndex>=0){
      double TruthP = varMuonTruthP(slc);
      double RecoP = varTruthMuonMatchedTrackCombinedP(slc);
      return (RecoP-TruthP)/TruthP;
    }
    else{
      return -999;
    }

  });

  //====   Matched reco track other properties

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

  //====   dEdX vs rr

  const MultiVar varTruthMuonMatchedTrackEnddedx([](const caf::SRSliceProxy* slc) -> std::vector<double> {
    int muonTrackIndex = varTruthMuonMatchedTrackIndex(slc);
    vector<double> ret;
    if(muonTrackIndex>=0){
      auto const& trk = slc->reco.trk.at(muonTrackIndex);
      int bp = trk.bestplane;
      for( auto const& pt : trk.calo[bp].points ){
        ret.push_back(pt.dedx);
      }
    }
    return ret;
  });

  const MultiVar varTruthMuonMatchedTrackEndrr([](const caf::SRSliceProxy* slc) -> std::vector<double> {
    int muonTrackIndex = varTruthMuonMatchedTrackIndex(slc);
    vector<double> ret;
    if(muonTrackIndex>=0){
      auto const& trk = slc->reco.trk.at(muonTrackIndex);
      int bp = trk.bestplane;
      for( auto const& pt : trk.calo[bp].points ){
        ret.push_back(pt.rr);
      }
    }
    return ret;

  });

  const MultiVar varTruthMuonMatchedTrackFrontdedx([](const caf::SRSliceProxy* slc) -> std::vector<double> {
    int muonTrackIndex = varTruthMuonMatchedTrackIndex(slc);
    vector<double> ret;
    if(muonTrackIndex>=0){
      auto const& trk = slc->reco.trk.at(muonTrackIndex);
      int bp = trk.bestplane;
      for( auto const& pt : trk.calo[bp].points ){
        ret.push_back(pt.dedx);
      }
    }
    return ret;
  });

  const MultiVar varTruthMuonMatchedTrackFrontdedxTemplate([](const caf::SRSliceProxy* slc) -> std::vector<double> {
    int muonTrackIndex = varTruthMuonMatchedTrackIndex(slc);
    vector<double> ret;
    if(muonTrackIndex>=0){
      auto const& trk = slc->reco.trk.at(muonTrackIndex);
      int bp = trk.bestplane;
      vector<float> rrs;
      for( auto const& pt : trk.calo[bp].points ){
        ret.push_back( dedxtempt.GetdEdX(pt.rr, 3) );
      }
    }
    return ret;
  });

  const MultiVar varTruthMuonMatchedTrackFrontdedxDiff([](const caf::SRSliceProxy* slc) -> std::vector<double> {
    int muonTrackIndex = varTruthMuonMatchedTrackIndex(slc);
    vector<double> ret;
    if(muonTrackIndex>=0){
      auto const& trk = slc->reco.trk.at(muonTrackIndex);
      int bp = trk.bestplane;
      vector<float> rrs;
      for( auto const& pt : trk.calo[bp].points ){
        ret.push_back( pt.dedx - dedxtempt.GetdEdX(pt.rr, 3) );
      }
    }
    return ret;
  });

  const MultiVar varTruthMuonMatchedTrackFrontrr([](const caf::SRSliceProxy* slc) -> std::vector<double> {
    int muonTrackIndex = varTruthMuonMatchedTrackIndex(slc);
    vector<double> ret;
    if(muonTrackIndex>=0){
      auto const& trk = slc->reco.trk.at(muonTrackIndex);
      int bp = trk.bestplane;
      vector<float> rrs;
      for( auto const& pt : trk.calo[bp].points ){
        rrs.push_back(pt.rr);
      }
      float rrmax = !rrs.empty() ? *std::max_element(rrs.begin(), rrs.end()) : 0.;
      for( auto const& rr : rrs ){
        ret.push_back(rrmax-rr);
      }
    }
    return ret;

  });

  const Var varTruthMuonMatchedTrackFrontLargedEdX([](const caf::SRSliceProxy* slc) -> int {
    int muonTrackIndex = varTruthMuonMatchedTrackIndex(slc);
    if(muonTrackIndex>=0){

      int nLargedEdX(0);
      auto const& trk = slc->reco.trk.at(muonTrackIndex);
      int bp = trk.bestplane;
      vector<float> rrs, dedxs;
      for( auto const& pt : trk.calo[bp].points ){
        rrs.push_back(pt.rr);
        dedxs.push_back(pt.dedx);
      }
      float rrmax = !rrs.empty() ? *std::max_element(rrs.begin(), rrs.end()) : 0.;
      for(unsigned int i=0; i<rrs.size(); i++){
        float dedx = dedxs.at(i);
        double this_frontrr = rrmax-rrs.at(i);
        if(this_frontrr>1. && this_frontrr<25.){
          if(dedx>100.){
            nLargedEdX++;
          }
        }
      }
      return nLargedEdX;
    }
    else{
      return 0.;
    }

  });


  //====   Michel study

  const Var varTruthMuonMatchedTrackStitchedTrackIndex([](const caf::SRSliceProxy* slc) -> int {
    int muonTrackIndex = varTruthMuonMatchedTrackIndex(slc);
    if(muonTrackIndex>=0){

      auto const& trk = slc->reco.trk.at(muonTrackIndex);
      TVector3 trk_end(trk.end.x, trk.end.y, trk.end.z);

      int sttrk_ind(-1);
      double dist(99999999.);
      for(std::size_t i(0); i < slc->reco.trk.size(); ++i){
        if(i==(std::size_t)muonTrackIndex) continue;
        auto const& sttrk = slc->reco.trk.at(i);
        TVector3 sttrk_start(sttrk.start.x, sttrk.start.y, sttrk.start.z);
        double this_dist = (trk_end-sttrk_start).Mag();
        if(this_dist<dist){
          dist = this_dist;
          sttrk_ind = i;
        }
      }

      return sttrk_ind;

    }
    else{
      return -999;
    }

  });
  const Var varTruthMuonMatchedTrackStitchedTrackDistance([](const caf::SRSliceProxy* slc) -> float {
    int muonStitchedTrackIndex = varTruthMuonMatchedTrackStitchedTrackIndex(slc);
    if(muonStitchedTrackIndex>=0){
      auto const& trk = slc->reco.trk.at(varTruthMuonMatchedTrackIndex(slc));
      auto const& sttrk = slc->reco.trk.at(muonStitchedTrackIndex);
      TVector3 trk_end(trk.end.x, trk.end.y, trk.end.z);
      TVector3 sttrk_start(sttrk.start.x, sttrk.start.y, sttrk.start.z);
      return (trk_end-sttrk_start).Mag();
    }
    else{
      return -9999999999.;
    }
  });
  const Var varTruthMuonMatchedTrackStitchedTrackLength([](const caf::SRSliceProxy* slc) -> int {
    int muonStitchedTrackIndex = varTruthMuonMatchedTrackStitchedTrackIndex(slc);
    if(muonStitchedTrackIndex>=0){
      return slc->reco.trk.at(muonStitchedTrackIndex).len;
    }
    else{
      return -99999;
    }
  });


  const Var varTruthMuonMatchedTrackStitchedTrackPDG([](const caf::SRSliceProxy* slc) -> int {
    int muonStitchedTrackIndex = varTruthMuonMatchedTrackStitchedTrackIndex(slc);
    if(muonStitchedTrackIndex>=0){
      return slc->reco.trk.at(muonStitchedTrackIndex).truth.p.pdg;
    }
    else{
      return -99999;
    }
  });

  const Var varTruthMuonMatchedTrackStitchedTrackCosine([](const caf::SRSliceProxy* slc) -> float {
    int muonStitchedTrackIndex = varTruthMuonMatchedTrackStitchedTrackIndex(slc);
    if(muonStitchedTrackIndex>=0){
      auto const& trk = slc->reco.trk.at(varTruthMuonMatchedTrackIndex(slc));
      auto const& sttrk = slc->reco.trk.at(muonStitchedTrackIndex);
      TVector3 trk_end_dir(trk.dir_end.x, trk.dir_end.y, trk.dir_end.z);
      TVector3 sttrk_dir(sttrk.dir.x, sttrk.dir.y, sttrk.dir.z);
      double angle = trk_end_dir.Angle(sttrk_dir);
      return TMath::Cos(angle);
    }
    else{
      return -9999999999.;
    }
  });

  const Var varTruthMuonMatchedTrackStitchedTrackChi2Muon([](const caf::SRSliceProxy* slc) -> float {
    int muonStitchedTrackIndex = varTruthMuonMatchedTrackStitchedTrackIndex(slc);
    if(muonStitchedTrackIndex>=0){
      auto const& sttrk = slc->reco.trk.at(muonStitchedTrackIndex);
      int bp = sttrk.bestplane;
      return sttrk.chi2pid[bp].chi2_muon;
    }
    else{
      return -9999999999.;
    }
  });

  const Var varTruthMuonMatchedTrackStitchedTrackChi2Proton([](const caf::SRSliceProxy* slc) -> float {
    int muonStitchedTrackIndex = varTruthMuonMatchedTrackStitchedTrackIndex(slc);
    if(muonStitchedTrackIndex>=0){
      auto const& sttrk = slc->reco.trk.at(muonStitchedTrackIndex);
      int bp = sttrk.bestplane;
      return sttrk.chi2pid[bp].chi2_proton;
    }
    else{
      return -9999999999.;
    }
  });

  const Var varTruthMuonMatchedTrackStitchedTrackChi2Pion([](const caf::SRSliceProxy* slc) -> float {
    int muonStitchedTrackIndex = varTruthMuonMatchedTrackStitchedTrackIndex(slc);
    if(muonStitchedTrackIndex>=0){
      auto const& sttrk = slc->reco.trk.at(muonStitchedTrackIndex);
      int bp = sttrk.bestplane;
      return sttrk.chi2pid[bp].chi2_pion;
    }
    else{
      return -9999999999.;
    }
  });

  const Var varTruthMuonMatchedTrackStitchedShowerIndex([](const caf::SRSliceProxy* slc) -> int {
    int muonTrackIndex = varTruthMuonMatchedTrackIndex(slc);
    if(muonTrackIndex>=0){

      auto const& trk = slc->reco.trk.at(muonTrackIndex);
      TVector3 trk_end(trk.end.x, trk.end.y, trk.end.z);

      int stshw_ind(-1);
      double dist(99999999.);
      for(std::size_t i(0); i < slc->reco.shw.size(); ++i){
        auto const& stshw = slc->reco.shw.at(i);
        TVector3 stshw_start(stshw.start.x, stshw.start.y, stshw.start.z);
        double this_dist = (trk_end-stshw_start).Mag();
        if(this_dist<dist){
          dist = this_dist;
          stshw_ind = i;
        }
      }

      return stshw_ind;

    }
    else{
      return -999;
    }

  });
  const Var varTruthMuonMatchedTrackStitchedShowerDistance([](const caf::SRSliceProxy* slc) -> float {
    int muonStitchedShowerIndex = varTruthMuonMatchedTrackStitchedShowerIndex(slc);
    if(muonStitchedShowerIndex>=0){
      auto const& trk = slc->reco.trk.at(varTruthMuonMatchedTrackIndex(slc));
      auto const& stshw = slc->reco.shw.at(muonStitchedShowerIndex);
      TVector3 trk_end(trk.end.x, trk.end.y, trk.end.z);
      TVector3 stshw_start(stshw.start.x, stshw.start.y, stshw.start.z);
      return (trk_end-stshw_start).Mag();
    }
    else{
      return -9999999999.;
    }
  });
  const Var varTruthMuonMatchedTrackStitchedShowerPDG([](const caf::SRSliceProxy* slc) -> int {
    int muonStitchedShowerIndex = varTruthMuonMatchedTrackStitchedShowerIndex(slc);
    if(muonStitchedShowerIndex>=0){
      return slc->reco.shw.at(muonStitchedShowerIndex).truth.p.pdg;
    }
    else{
      return -99999;
    }
  });

  const Var varTruthMuonMatchedTrackStitchedShowerEnergy([](const caf::SRSliceProxy* slc) -> float {
    int muonStitchedShowerIndex = varTruthMuonMatchedTrackStitchedShowerIndex(slc);
    if(muonStitchedShowerIndex>=0){
      auto const& stshw = slc->reco.shw.at(muonStitchedShowerIndex);
      return stshw.bestplane_energy;
    }
    else{
      return -9999999999.;
    }
  });
  const Var varTruthMuonMatchedTrackStitchedShowerCosine([](const caf::SRSliceProxy* slc) -> float {
    int muonStitchedShowerIndex = varTruthMuonMatchedTrackStitchedShowerIndex(slc);
    if(muonStitchedShowerIndex>=0){
      auto const& trk = slc->reco.trk.at(varTruthMuonMatchedTrackIndex(slc));
      auto const& stshw = slc->reco.shw.at(muonStitchedShowerIndex);
      TVector3 trk_end_dir(trk.dir_end.x, trk.dir_end.y, trk.dir_end.z);
      TVector3 stshw_dir(stshw.dir.x, stshw.dir.y, stshw.dir.z);
      double angle = trk_end_dir.Angle(stshw_dir);
      return TMath::Cos(angle);
    }
    else{
      return -9999999999.;
    }
  });
  //====   Count all close-enought track/shower daughters
  const Var varTruthMuonMatchedTrackNDaughterTracks([](const caf::SRSliceProxy* slc) -> int {
    int muonTrackIndex = varTruthMuonMatchedTrackIndex(slc);
    if(muonTrackIndex>=0){
      return GetNDaughterTracks(slc, muonTrackIndex);
    }
    else{
      return 0;
    }
  });
  const Var varTruthMuonMatchedTrackNDaughterShowers([](const caf::SRSliceProxy* slc) -> int { 
    int muonTrackIndex = varTruthMuonMatchedTrackIndex(slc);
    if(muonTrackIndex>=0){
      return GetNDaughterShowers(slc, muonTrackIndex);
    }
    else{
      return 0;
    }
  });

  //====   Matched reco shower

  const Var varTruthMuonMatchedShowerIndex([](const caf::SRSliceProxy* slc) -> int {

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

  //==== For a given true proton (truth_index), find a reco whose best-matched is this proton

  const Var varTruthProtonMatchedTrackIndex([](const caf::SRSliceProxy* slc) -> int {

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
      bool isContained = fv_track.isContained(trk_Proton.end.x, trk_Proton.end.y, trk_Proton.end.z);
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

  const Var varTruthProtonMatchedTrackChi2Pion([](const caf::SRSliceProxy* slc) -> double {
    int protonTrackIndex = varTruthProtonMatchedTrackIndex(slc);
    if(protonTrackIndex>=0){

      auto const& trk = slc->reco.trk.at(protonTrackIndex);
      int bp = trk.bestplane;
      float Chi2Pion = trk.chi2pid[bp].chi2_pion;
      if(isnan(Chi2Pion)) return -999.;
      else return Chi2Pion;

    }
    else{
      return 999999;
    }

  });

  const Var varTruthProtonMatchedTrackReducedChi2Pion([](const caf::SRSliceProxy* slc) -> double {
    int protonTrackIndex = varTruthProtonMatchedTrackIndex(slc);
    if(protonTrackIndex>=0){

      auto const& trk = slc->reco.trk.at(protonTrackIndex);
      int bp = trk.bestplane;
      float Chi2Pion = trk.chi2pid[bp].chi2_pion;
      int Chi2Ndof = trk.chi2pid[bp].pid_ndof;
      if(isnan(Chi2Pion)) return -999.;
      else if(Chi2Ndof==0) return -999.;
      else return Chi2Pion/double(Chi2Ndof);

    }
    else{
      return 999999;
    }

  });

  //====   Matched reco length

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

  //====   Matched reco momentum

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

  //====   Matched reco momentum residual

  const Var varTruthProtonMatchedTrackRangePResidual([](const caf::SRSliceProxy* slc) -> double {
    int ProtonTrackIndex = varTruthProtonMatchedTrackIndex(slc);
    if(ProtonTrackIndex>=0){
      double TruthP = varProtonTruthP(slc);
      double RecoP = varTruthProtonMatchedTrackRangeP(slc);
      return (RecoP-TruthP);
    }
    else{
      return -999;
    }

  });

  const Var varTruthProtonMatchedTrackMCSPResidual([](const caf::SRSliceProxy* slc) -> double {
    int ProtonTrackIndex = varTruthProtonMatchedTrackIndex(slc);
    if(ProtonTrackIndex>=0){
      double TruthP = varProtonTruthP(slc);
      double RecoP = varTruthProtonMatchedTrackMCSP(slc);
      return (RecoP-TruthP);
    }
    else{
      return -999;
    }

  });
  const Var varTruthProtonMatchedTrackCombinedPResidual([](const caf::SRSliceProxy* slc) -> double {
    int ProtonTrackIndex = varTruthProtonMatchedTrackIndex(slc);
    if(ProtonTrackIndex>=0){
      double TruthP = varProtonTruthP(slc);
      double RecoP = varTruthProtonMatchedTrackCombinedP(slc);
      return (RecoP-TruthP);
    }
    else{
      return -999;
    }

  });

  //====   Matched reco momentum residual fraction

  const Var varTruthProtonMatchedTrackRangePResidualFraction([](const caf::SRSliceProxy* slc) -> double {
    int ProtonTrackIndex = varTruthProtonMatchedTrackIndex(slc);
    if(ProtonTrackIndex>=0){
      double TruthP = varProtonTruthP(slc);
      double RecoP = varTruthProtonMatchedTrackRangeP(slc);
      return (RecoP-TruthP)/TruthP;
    }
    else{
      return -999;
    }

  });

  const Var varTruthProtonMatchedTrackMCSPResidualFraction([](const caf::SRSliceProxy* slc) -> double {
    int ProtonTrackIndex = varTruthProtonMatchedTrackIndex(slc);
    if(ProtonTrackIndex>=0){
      double TruthP = varProtonTruthP(slc);
      double RecoP = varTruthProtonMatchedTrackMCSP(slc);
      return (RecoP-TruthP)/TruthP;
    }
    else{
      return -999;
    }

  });

  const Var varTruthProtonMatchedTrackCombinedPResidualFraction([](const caf::SRSliceProxy* slc) -> double {
    int ProtonTrackIndex = varTruthProtonMatchedTrackIndex(slc);
    if(ProtonTrackIndex>=0){
      double TruthP = varProtonTruthP(slc);
      double RecoP = varTruthProtonMatchedTrackCombinedP(slc);
      return (RecoP-TruthP)/TruthP;
    }
    else{
      return -999;
    }

  });

  //====   Matched reco track other properties

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

  //====   dEdX vs rr

  const MultiVar varTruthProtonMatchedTrackEnddedx([](const caf::SRSliceProxy* slc) -> std::vector<double> {
    int protonTrackIndex = varTruthProtonMatchedTrackIndex(slc);
    vector<double> ret;
    if(protonTrackIndex>=0){
      auto const& trk = slc->reco.trk.at(protonTrackIndex);
      int bp = trk.bestplane;
      for( auto const& pt : trk.calo[bp].points ){
        ret.push_back(pt.dedx);
      }
    }
    return ret;
  });

  const MultiVar varTruthProtonMatchedTrackEndrr([](const caf::SRSliceProxy* slc) -> std::vector<double> {
    int protonTrackIndex = varTruthProtonMatchedTrackIndex(slc);
    vector<double> ret;
    if(protonTrackIndex>=0){
      auto const& trk = slc->reco.trk.at(protonTrackIndex);
      int bp = trk.bestplane;
      for( auto const& pt : trk.calo[bp].points ){
        ret.push_back(pt.rr);
      }
    }
    return ret;

  });




  //====   Matched reco shower

  const Var varTruthProtonMatchedShowerIndex([](const caf::SRSliceProxy* slc) -> int {

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

  //====   Matched stub

  const Var varTruthProtonMatchedStubIndex([](const caf::SRSliceProxy* slc) -> int {

    int truth_idx = varProtonTruthIndex(slc);

    return GetMatchedRecoStubIndex(slc, truth_idx);

  });

  const Var varTruthProtonMatchedStubE([](const caf::SRSliceProxy* slc) -> double {
    int ProtonStubIndex = varTruthProtonMatchedStubIndex(slc);
    if(ProtonStubIndex>=0){
      auto const& stub_Proton = slc->reco.stub.at(ProtonStubIndex);

      int nHitMax(-999);
      double ret(-999.);
      for(unsigned ip=0; ip<stub_Proton.planes.size(); ip++){

        double sumQ(0.);

        if(stub_Proton.planes[ip].p!=2) continue;

        for(unsigned ih=0; ih<stub_Proton.planes[ip].hits.size(); ih++){
          sumQ += stub_Proton.planes[ip].hits.at(ih).charge;
        }

        if((int)stub_Proton.planes[ip].hits.size() >= nHitMax){
          nHitMax = stub_Proton.planes[ip].hits.size();
          ret = sumQ;
        }

      }
      //std::cout << "[varTruthProtonMatchedStubE] nHitMax = " << nHitMax << ", sumQ = " << ret << std::endl;
      static double p0 = -0.00236892;
      static double p1 = 0.383538;
      if(ret<=0.) return -999.;
      else return (ret*23.6e-9-p0)/p1; // to GeV
    }
    else{
      return -999;
    }

  });

  const Var varTruthProtonMatchedStubLength([](const caf::SRSliceProxy* slc) -> double {
    int ProtonStubIndex = varTruthProtonMatchedStubIndex(slc);
    if(ProtonStubIndex>=0){
      auto const& stub_Proton = slc->reco.stub.at(ProtonStubIndex);
      TVector3 v3_Pos_start(stub_Proton.vtx.x, stub_Proton.vtx.y, stub_Proton.vtx.z);
      TVector3 v3_Pos_end(stub_Proton.end.x, stub_Proton.end.y, stub_Proton.end.z);

      return (v3_Pos_start-v3_Pos_end).Mag();

    }
    else{
      return -999;
    }
  });

  //==== For a given true charged pion (truth_index), find a reco track whose best-matched is this charged pion

  const Var varTruthChargedPionMatchedTrackIndex([](const caf::SRSliceProxy* slc) -> int {

    int truth_idx = varChargedPionTruthIndex(slc);

    return GetMatchedRecoTrackIndex(slc, truth_idx);

  });

  const Var varTruthChargedPionMatchedTrackContainedness([](const caf::SRSliceProxy* slc) -> int {
    int cpionTrackIndex = varTruthChargedPionMatchedTrackIndex(slc);
    //==== -1 : No track found
    //==== 0 : Exiting
    //==== 1 : Contained
    if(cpionTrackIndex>=0){

      auto const& trk = slc->reco.trk.at(cpionTrackIndex);
      bool isContained = fv_track.isContained(trk.end.x, trk.end.y, trk.end.z);
      if(isContained) return 1;
      else return 0;
    }
    else{
      return -1;
    }

  });

  const Var varTruthChargedPionMatchedTrackChi2Proton([](const caf::SRSliceProxy* slc) -> double {
    int cpionTrackIndex = varTruthChargedPionMatchedTrackIndex(slc);
    if(cpionTrackIndex>=0){

      auto const& trk = slc->reco.trk.at(cpionTrackIndex);
      int bp = trk.bestplane;
      float Chi2Proton = trk.chi2pid[bp].chi2_proton;
      if(isnan(Chi2Proton)) return -999.;
      else return Chi2Proton;

    }
    else{
      return 999999;
    }
  });

  const Var varTruthChargedPionMatchedTrackReducedChi2Proton([](const caf::SRSliceProxy* slc) -> double {
    int cpionTrackIndex = varTruthChargedPionMatchedTrackIndex(slc);
    if(cpionTrackIndex>=0){

      auto const& trk = slc->reco.trk.at(cpionTrackIndex);
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

  const Var varTruthChargedPionMatchedTrackChi2Muon([](const caf::SRSliceProxy* slc) -> double {
    int cpionTrackIndex = varTruthChargedPionMatchedTrackIndex(slc);
    if(cpionTrackIndex>=0){

      auto const& trk = slc->reco.trk.at(cpionTrackIndex);
      int bp = trk.bestplane;
      float Chi2Muon = trk.chi2pid[bp].chi2_muon;
      if(isnan(Chi2Muon)) return -999.;
      else return Chi2Muon;

    }
    else{
      return 999999;
    }
  });

  const Var varTruthChargedPionMatchedTrackReducedChi2Muon([](const caf::SRSliceProxy* slc) -> double {
    int cpionTrackIndex = varTruthChargedPionMatchedTrackIndex(slc);
    if(cpionTrackIndex>=0){

      auto const& trk = slc->reco.trk.at(cpionTrackIndex);
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

  const Var varTruthChargedPionMatchedTrackChi2Pion([](const caf::SRSliceProxy* slc) -> double {
    int cpionTrackIndex = varTruthChargedPionMatchedTrackIndex(slc);
    if(cpionTrackIndex>=0){

      auto const& trk = slc->reco.trk.at(cpionTrackIndex);
      int bp = trk.bestplane;
      float Chi2Pion = trk.chi2pid[bp].chi2_pion;
      if(isnan(Chi2Pion)) return -999.;
      else return Chi2Pion;

    }
    else{
      return 999999;
    }
  });

  const Var varTruthChargedPionMatchedTrackReducedChi2Pion([](const caf::SRSliceProxy* slc) -> double {
    int cpionTrackIndex = varTruthChargedPionMatchedTrackIndex(slc);
    if(cpionTrackIndex>=0){

      auto const& trk = slc->reco.trk.at(cpionTrackIndex);
      int bp = trk.bestplane;
      float Chi2Pion = trk.chi2pid[bp].chi2_pion;
      int Chi2Ndof = trk.chi2pid[bp].pid_ndof;
      if(isnan(Chi2Pion)) return -999.;
      else if(Chi2Ndof==0) return -999.;
      else return Chi2Pion/double(Chi2Ndof);

    }
    else{
      return 999999;
    }
  });

  const Var varTruthChargedPionMatchedTrackLength([](const caf::SRSliceProxy* slc) -> double {
    int cpionTrackIndex = varTruthChargedPionMatchedTrackIndex(slc);
    if(cpionTrackIndex>=0){
      auto const& trk = slc->reco.trk.at(cpionTrackIndex);
      return trk.len;
    }
    else{
      return -999;
    }

  });

  //====   Matched reco momentum

  const Var varTruthChargedPionMatchedTrackRangeP([](const caf::SRSliceProxy* slc) -> double {
    int cpionTrackIndex = varTruthChargedPionMatchedTrackIndex(slc);
    if(cpionTrackIndex>=0){
      auto const& trk = slc->reco.trk.at(cpionTrackIndex);
      return trk.rangeP.p_pion;
    }
    else{
      return -999;
    }

  });

  //====   Matched reco momentum residual fraction

  const Var varTruthChargedPionMatchedTrackRangePResidualFraction([](const caf::SRSliceProxy* slc) -> double {
    int cpionTrackIndex = varTruthChargedPionMatchedTrackIndex(slc);
    if(cpionTrackIndex>=0){
      double TruthP = varChargedPionTruthP(slc);
      double RecoP = varTruthChargedPionMatchedTrackRangeP(slc);
      return (RecoP-TruthP)/TruthP;
    }
    else{
      return -999;
    }

  });

  //====   dEdX vs rr

  const MultiVar varTruthChargedPionMatchedTrackEnddedx([](const caf::SRSliceProxy* slc) -> std::vector<double> {
    int cpionTrackIndex = varTruthChargedPionMatchedTrackIndex(slc);
    vector<double> ret;
    if(cpionTrackIndex>=0){
      auto const& trk = slc->reco.trk.at(cpionTrackIndex);
      int bp = trk.bestplane;
      for( auto const& pt : trk.calo[bp].points ){
        ret.push_back(pt.dedx);
      }
    }
    return ret;
  });

  const MultiVar varTruthChargedPionMatchedTrackEndrr([](const caf::SRSliceProxy* slc) -> std::vector<double> {
    int cpionTrackIndex = varTruthChargedPionMatchedTrackIndex(slc);
    vector<double> ret;
    if(cpionTrackIndex>=0){
      auto const& trk = slc->reco.trk.at(cpionTrackIndex);
      int bp = trk.bestplane;
      for( auto const& pt : trk.calo[bp].points ){
        ret.push_back(pt.rr);
      }
    }
    return ret;

  });

  const MultiVar varTruthChargedPionMatchedTrackFrontdedx([](const caf::SRSliceProxy* slc) -> std::vector<double> {
    int cpionTrackIndex = varTruthChargedPionMatchedTrackIndex(slc);
    vector<double> ret;
    if(cpionTrackIndex>=0){
      auto const& trk = slc->reco.trk.at(cpionTrackIndex);
      int bp = trk.bestplane;
      for( auto const& pt : trk.calo[bp].points ){
        ret.push_back(pt.dedx);
      }
    }
    return ret;
  });

  const MultiVar varTruthChargedPionMatchedTrackFrontdedxTemplate([](const caf::SRSliceProxy* slc) -> std::vector<double> {
    int cpionTrackIndex = varTruthChargedPionMatchedTrackIndex(slc);
    vector<double> ret;
    if(cpionTrackIndex>=0){
      auto const& trk = slc->reco.trk.at(cpionTrackIndex);
      int bp = trk.bestplane;
      for( auto const& pt : trk.calo[bp].points ){
        ret.push_back( dedxtempt.GetdEdX(pt.rr, 2) );
      }
    }
    return ret;
  });

  const MultiVar varTruthChargedPionMatchedTrackFrontdedxDiff([](const caf::SRSliceProxy* slc) -> std::vector<double> {
    int cpionTrackIndex = varTruthChargedPionMatchedTrackIndex(slc);
    vector<double> ret;
    if(cpionTrackIndex>=0){
      auto const& trk = slc->reco.trk.at(cpionTrackIndex);
      int bp = trk.bestplane;
      for( auto const& pt : trk.calo[bp].points ){
        ret.push_back( pt.dedx - dedxtempt.GetdEdX(pt.rr, 2) );
      }
    }
    return ret;
  });

  const MultiVar varTruthChargedPionMatchedTrackFrontrr([](const caf::SRSliceProxy* slc) -> std::vector<double> {
    int cpionTrackIndex = varTruthChargedPionMatchedTrackIndex(slc);
    vector<double> ret;
    if(cpionTrackIndex>=0){
      auto const& trk = slc->reco.trk.at(cpionTrackIndex);
      int bp = trk.bestplane;
      vector<float> rrs;
      for( auto const& pt : trk.calo[bp].points ){
        //std::cout << pt.rr << "\t" << pt.dedx << std::endl;
        rrs.push_back(pt.rr);
      }
      float rrmax = !rrs.empty() ? *std::max_element(rrs.begin(), rrs.end()) : 0.;
      for( auto const& rr : rrs ){
        ret.push_back(rrmax-rr);
      }
    }
    return ret;

  });

  const Var varTruthChargedPionMatchedTrackFrontLargedEdX([](const caf::SRSliceProxy* slc) -> int {
    int cpionTrackIndex = varTruthChargedPionMatchedTrackIndex(slc);
    if(cpionTrackIndex>=0){
      
      int nLargedEdX(0);
      auto const& trk = slc->reco.trk.at(cpionTrackIndex);
      int bp = trk.bestplane;
      vector<float> rrs, dedxs;
      for( auto const& pt : trk.calo[bp].points ){
        rrs.push_back(pt.rr);
        dedxs.push_back(pt.dedx);
      }
      float rrmax = !rrs.empty() ? *std::max_element(rrs.begin(), rrs.end()) : 0.;
      for(unsigned int i=0; i<rrs.size(); i++){
        float dedx = dedxs.at(i);
        double this_frontrr = rrmax-rrs.at(i);
        if(this_frontrr>1. && this_frontrr<25.){
          if(dedx>100.){
            nLargedEdX++;
          }
        }
      }
      return nLargedEdX;

    }
    else{
      return 0.;
    }

  });

  //====   Decay study

  const Var varTruthChargedPionMatchedTrackStitchedTrackIndex([](const caf::SRSliceProxy* slc) -> int {
    int cpionTrackIndex = varTruthChargedPionMatchedTrackIndex(slc);
    if(cpionTrackIndex>=0){

      auto const& trk = slc->reco.trk.at(cpionTrackIndex);
      TVector3 trk_end(trk.end.x, trk.end.y, trk.end.z);

      int sttrk_ind(-1);
      double dist(99999999.);
      for(std::size_t i(0); i < slc->reco.trk.size(); ++i){
        if(i==(std::size_t)cpionTrackIndex) continue;
        auto const& sttrk = slc->reco.trk.at(i);
        TVector3 sttrk_start(sttrk.start.x, sttrk.start.y, sttrk.start.z);
        double this_dist = (trk_end-sttrk_start).Mag();
        if(this_dist<dist){
          dist = this_dist;
          sttrk_ind = i;
        }
      }

      return sttrk_ind;

    }
    else{
      return -999;
    }

  });
  const Var varTruthChargedPionMatchedTrackStitchedTrackDistance([](const caf::SRSliceProxy* slc) -> float {
    int cpionStitchedTrackIndex = varTruthChargedPionMatchedTrackStitchedTrackIndex(slc);
    if(cpionStitchedTrackIndex>=0){
      auto const& trk = slc->reco.trk.at(varTruthChargedPionMatchedTrackIndex(slc));
      auto const& sttrk = slc->reco.trk.at(cpionStitchedTrackIndex);
      TVector3 trk_end(trk.end.x, trk.end.y, trk.end.z);
      TVector3 sttrk_start(sttrk.start.x, sttrk.start.y, sttrk.start.z);
      return (trk_end-sttrk_start).Mag();
    }
    else{
      return -9999999999.;
    }
  });
  const Var varTruthChargedPionMatchedTrackStitchedTrackPDG([](const caf::SRSliceProxy* slc) -> int {
    int cpionStitchedTrackIndex = varTruthChargedPionMatchedTrackStitchedTrackIndex(slc);
    if(cpionStitchedTrackIndex>=0){
      return slc->reco.trk.at(cpionStitchedTrackIndex).truth.p.pdg;
    }
    else{
      return -99999;
    }
  });

  const Var varTruthChargedPionMatchedTrackStitchedTrackLength([](const caf::SRSliceProxy* slc) -> float {

    int cpionStitchedTrackIndex = varTruthChargedPionMatchedTrackStitchedTrackIndex(slc);
    if(cpionStitchedTrackIndex>=0){
      return slc->reco.trk.at(cpionStitchedTrackIndex).len;
    }
    else{
      return -99999.;
    }
  });

  const Var varTruthChargedPionMatchedTrackStitchedTrackCosine([](const caf::SRSliceProxy* slc) -> float {
    int cpionStitchedTrackIndex = varTruthChargedPionMatchedTrackStitchedTrackIndex(slc);
    if(cpionStitchedTrackIndex>=0){
      auto const& trk = slc->reco.trk.at(varTruthChargedPionMatchedTrackIndex(slc));
      auto const& sttrk = slc->reco.trk.at(cpionStitchedTrackIndex);
      TVector3 trk_end_dir(trk.dir_end.x, trk.dir_end.y, trk.dir_end.z);
      TVector3 sttrk_dir(sttrk.dir.x, sttrk.dir.y, sttrk.dir.z);
      double angle = trk_end_dir.Angle(sttrk_dir);
      return TMath::Cos(angle);
    }
    else{
      return -9999999999.;
    }
  });

  const Var varTruthChargedPionMatchedTrackStitchedTrackChi2Muon([](const caf::SRSliceProxy* slc) -> float {
    int cpionStitchedTrackIndex = varTruthChargedPionMatchedTrackStitchedTrackIndex(slc);
    if(cpionStitchedTrackIndex>=0){
      auto const& sttrk = slc->reco.trk.at(cpionStitchedTrackIndex);
      int bp = sttrk.bestplane;
      return sttrk.chi2pid[bp].chi2_muon;
    }
    else{
      return -9999999999.;
    }
  });

  const Var varTruthChargedPionMatchedTrackStitchedTrackChi2Proton([](const caf::SRSliceProxy* slc) -> float {
    int cpionStitchedTrackIndex = varTruthChargedPionMatchedTrackStitchedTrackIndex(slc);
    if(cpionStitchedTrackIndex>=0){
      auto const& sttrk = slc->reco.trk.at(cpionStitchedTrackIndex);
      int bp = sttrk.bestplane;
      return sttrk.chi2pid[bp].chi2_proton;
    }
    else{
      return -9999999999.;
    }
  });

  const Var varTruthChargedPionMatchedTrackStitchedTrackChi2Pion([](const caf::SRSliceProxy* slc) -> float {
    int cpionStitchedTrackIndex = varTruthChargedPionMatchedTrackStitchedTrackIndex(slc);
    if(cpionStitchedTrackIndex>=0){
      auto const& sttrk = slc->reco.trk.at(cpionStitchedTrackIndex);
      int bp = sttrk.bestplane;
      return sttrk.chi2pid[bp].chi2_pion;
    }
    else{
      return -9999999999.;
    }
  });




  const Var varTruthChargedPionMatchedTrackStitchedShowerIndex([](const caf::SRSliceProxy* slc) -> int {
    int cpionTrackIndex = varTruthChargedPionMatchedTrackIndex(slc);
    if(cpionTrackIndex>=0){

      auto const& trk = slc->reco.trk.at(cpionTrackIndex);
      TVector3 trk_end(trk.end.x, trk.end.y, trk.end.z);

      int stshw_ind(-1);
      double dist(99999999.);
      for(std::size_t i(0); i < slc->reco.shw.size(); ++i){
        auto const& stshw = slc->reco.shw.at(i);
        TVector3 stshw_start(stshw.start.x, stshw.start.y, stshw.start.z);
        double this_dist = (trk_end-stshw_start).Mag();
        if(this_dist<dist){
          dist = this_dist;
          stshw_ind = i;
        }
      }

      return stshw_ind;

    }
    else{
      return -999;
    }

  });
  const Var varTruthChargedPionMatchedTrackStitchedShowerDistance([](const caf::SRSliceProxy* slc) -> float {
    int cpionStitchedShowerIndex = varTruthChargedPionMatchedTrackStitchedShowerIndex(slc);
    if(cpionStitchedShowerIndex>=0){
      auto const& trk = slc->reco.trk.at(varTruthChargedPionMatchedTrackIndex(slc));
      auto const& stshw = slc->reco.shw.at(cpionStitchedShowerIndex);
      TVector3 trk_end(trk.end.x, trk.end.y, trk.end.z);
      TVector3 stshw_start(stshw.start.x, stshw.start.y, stshw.start.z);
      return (trk_end-stshw_start).Mag();
    }
    else{
      return -9999999999.;
    }
  });
  const Var varTruthChargedPionMatchedTrackStitchedShowerPDG([](const caf::SRSliceProxy* slc) -> int {
    int cpionStitchedShowerIndex = varTruthChargedPionMatchedTrackStitchedShowerIndex(slc);
    if(cpionStitchedShowerIndex>=0){
      return slc->reco.shw.at(cpionStitchedShowerIndex).truth.p.pdg;
    }
    else{
      return -99999;
    }
  });

  const Var varTruthChargedPionMatchedTrackStitchedShowerEnergy([](const caf::SRSliceProxy* slc) -> float {
    int cpionStitchedShowerIndex = varTruthChargedPionMatchedTrackStitchedShowerIndex(slc);
    if(cpionStitchedShowerIndex>=0){
      auto const& stshw = slc->reco.shw.at(cpionStitchedShowerIndex);
      return stshw.bestplane_energy;
    }
    else{
      return -9999999999.;
    }
  });
  const Var varTruthChargedPionMatchedTrackStitchedShowerCosine([](const caf::SRSliceProxy* slc) -> float {
    int cpionStitchedShowerIndex = varTruthChargedPionMatchedTrackStitchedShowerIndex(slc);
    if(cpionStitchedShowerIndex>=0){
      auto const& trk = slc->reco.trk.at(varTruthChargedPionMatchedTrackIndex(slc));
      auto const& stshw = slc->reco.shw.at(cpionStitchedShowerIndex);
      TVector3 trk_end_dir(trk.dir_end.x, trk.dir_end.y, trk.dir_end.z);
      TVector3 stshw_dir(stshw.dir.x, stshw.dir.y, stshw.dir.z);
      double angle = trk_end_dir.Angle(stshw_dir);
      return TMath::Cos(angle);
    }
    else{
      return -9999999999.;
    }
  });
  //====   Count all close-enought track/shower daughters
  const Var varTruthChargedPionMatchedTrackNDaughterTracks([](const caf::SRSliceProxy* slc) -> int {
    int cpionTrackIndex = varTruthChargedPionMatchedTrackIndex(slc);
    if(cpionTrackIndex>=0){
      return GetNDaughterTracks(slc, cpionTrackIndex);
    }
    else{
      return 0;
    }
  });
  const Var varTruthChargedPionMatchedTrackNDaughterShowers([](const caf::SRSliceProxy* slc) -> int {
    int cpionTrackIndex = varTruthChargedPionMatchedTrackIndex(slc);
    if(cpionTrackIndex>=0){
      return GetNDaughterShowers(slc, cpionTrackIndex);
    }
    else{
      return 0;
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

      Contained = fv_track.isContained(trk.end.x, trk.end.y, trk.end.z);

      MaybeMuonExiting = ( !Contained && trk.len > 100);
      MaybeMuonContained = ( Contained && PassChi2 && trk.len > 50. );

      if( AtSlice && ( MaybeMuonExiting || MaybeMuonContained ) ){
        out.push_back(i);
      }

    }

    return out;
  });

  const Var varNMuonCandTrack([](const caf::SRSliceProxy* slc) -> int {
    const int nMuonCandTrack = varAllMuonTrackIndices(slc).size();
    return nMuonCandTrack;
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

      Contained = fv_track.isContained(trk.end.x, trk.end.y, trk.end.z);

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
      Contained = fv_track.isContained(trk.end.x, trk.end.y, trk.end.z);
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

  const Var varMuonChi2Pion([](const caf::SRSliceProxy* slc) -> float {
    if( varMuonTrackInd(slc) >= 0 ){
      auto const& trk = slc->reco.trk.at(varMuonTrackInd(slc));
      int bp = trk.bestplane;
      return trk.chi2pid[bp].chi2_pion;
    }
    else{
      return -999.;
    }
  });

  const Var varMuonChi2MuonReEval([](const caf::SRSliceProxy* slc) -> float {
    if( varMuonTrackInd(slc) >= 0 ){
      auto const& trk = slc->reco.trk.at(varMuonTrackInd(slc));
      int bp = trk.bestplane;

      int npt(0), counter(0);
      double chi2mu(0.);
      for( auto const& pt : trk.calo[bp].points ){
        double rr = pt.rr;
        double dedx = pt.dedx;
        double temp_dedx = dedxtempt.GetdEdX(rr, 3);
        double temp_dedx_err = dedxtempt.GetdEdXErr(rr, 3);
        double errdedx = 0.04231+0.0001783*dedx*dedx; //resolution on dE/dx
        errdedx *= dedx;
        if( counter!=0 && rr<26 && dedx<=1000. ){
          double this_chi2 = pow((dedx-temp_dedx)/std::sqrt(pow(temp_dedx_err,2)+pow(errdedx,2)),2);
          std::cout << "[JSKIMDEBUG] rr = " << rr << ", dedx = " << dedx << ", temp_dedx = " << temp_dedx << ", temp_dedx_err = " << temp_dedx_err << ", errdedx = " << errdedx << " -> this_chi2 = " << this_chi2 << std::endl;
          chi2mu += pow((dedx-temp_dedx)/std::sqrt(pow(temp_dedx_err,2)+pow(errdedx,2)),2);
          npt++;
        }
        counter++;
      }
      std::cout << "[JSKIMDEBUG] -> chi2mu = " << chi2mu << ", chi2mu/ndf = " << chi2mu/(double)npt << std::endl;
      return chi2mu/(double)npt;
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
  const Var varMuonRecoEndX([](const caf::SRSliceProxy* slc) -> float {
    if( varMuonTrackInd(slc) >= 0 ){
      auto const& trk = slc->reco.trk.at(varMuonTrackInd(slc));
      return trk.end.x;
    }
    return -9999999999.;
  });
  const Var varMuonRecoEndY([](const caf::SRSliceProxy* slc) -> float {
    if( varMuonTrackInd(slc) >= 0 ){
      auto const& trk = slc->reco.trk.at(varMuonTrackInd(slc));
      return trk.end.y;
    }
    return -9999999999.;
  });
  const Var varMuonRecoEndZ([](const caf::SRSliceProxy* slc) -> float {
    if( varMuonTrackInd(slc) >= 0 ){
      auto const& trk = slc->reco.trk.at(varMuonTrackInd(slc));
      return trk.end.z;
    }
    return -9999999999.;
  });

  const Var varMuonRecoTrackFromVertex([](const caf::SRSliceProxy* slc) -> float {
    if( varMuonTrackInd(slc) >= 0 ){

      auto const& trk = slc->reco.trk.at(varMuonTrackInd(slc));
      return sqrt( pow( slc->vertex.x - trk.start.x, 2.0 ) +
              pow( slc->vertex.y - trk.start.y, 2.0 ) +
              pow( slc->vertex.z - trk.start.z, 2.0 ) );

    }
    return -9999999999.;
  });

  const Var varMuonRecoDirX([](const caf::SRSliceProxy* slc) -> float {
    if( varMuonTrackInd(slc) >= 0 ){
      auto const& trk = slc->reco.trk.at(varMuonTrackInd(slc));
      return trk.dir.x;
    }
    return -9999999999.;
  });
  const Var varMuonRecoDirY([](const caf::SRSliceProxy* slc) -> float {
    if( varMuonTrackInd(slc) >= 0 ){
      auto const& trk = slc->reco.trk.at(varMuonTrackInd(slc));
      return trk.dir.y;
    }
    return -9999999999.;
  });
  const Var varMuonRecoDirZ([](const caf::SRSliceProxy* slc) -> float {
    if( varMuonTrackInd(slc) >= 0 ){
      auto const& trk = slc->reco.trk.at(varMuonTrackInd(slc));
      return trk.dir.z;
    }
    return -9999999999.;
  });

  const Var varMuonRecoForceDownDirX([](const caf::SRSliceProxy* slc) -> float {
    if( varMuonTrackInd(slc) >= 0 ){
      auto const& trk = slc->reco.trk.at(varMuonTrackInd(slc));
      double flip = trk.dir.y>0 ? -1. : +1.; // if original track is upward(>0), flip it
      return trk.dir.x*flip;
    }
    return -9999999999.;
  });
  const Var varMuonRecoForceDownDirY([](const caf::SRSliceProxy* slc) -> float {
    if( varMuonTrackInd(slc) >= 0 ){
      auto const& trk = slc->reco.trk.at(varMuonTrackInd(slc));
      double flip = trk.dir.y>0 ? -1. : +1.; // if original track is upward(>0), flip it
      return trk.dir.y*flip;
    }
    return -9999999999.;
  });
  const Var varMuonRecoForceDownDirZ([](const caf::SRSliceProxy* slc) -> float {
    if( varMuonTrackInd(slc) >= 0 ){
      auto const& trk = slc->reco.trk.at(varMuonTrackInd(slc));
      double flip = trk.dir.y>0 ? -1. : +1.; // if original track is upward(>0), flip it
      return trk.dir.z*flip;
    }
    return -9999999999.;
  });
  //==== dedx 
  const MultiVar varMuonTrackCalodedx([](const caf::SRSliceProxy* slc) -> std::vector<double> {

    std::vector<double> out;
    if( varMuonTrackInd(slc) >= 0 ){
      auto const& trk = slc->reco.trk.at(varMuonTrackInd(slc));
      int bp = trk.bestplane;
      for( auto const& pt : trk.calo[bp].points ){
        out.push_back(pt.dedx);
      }
    }

    return out;

  });
  const MultiVar varMuonTrackCalorr([](const caf::SRSliceProxy* slc) -> std::vector<double> {

    std::vector<double> out;
    if( varMuonTrackInd(slc) >= 0 ){
      auto const& trk = slc->reco.trk.at(varMuonTrackInd(slc));
      int bp = trk.bestplane;
      std::cout << "[varMuonTrackCalorr] Printing rr" << std::endl;
      for( auto const& pt : trk.calo[bp].points ){
        out.push_back(pt.rr);
        std::cout << pt.rr << std::endl;
      }
    }

    return out;

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

  //==== Proton

  //====   Track-based

  const Var varNProtonCandTrack([](const caf::SRSliceProxy* slc) -> int {
    const int muTrackInd = varMuonTrackInd(slc);

    //==== The (dis)qualification of a slice is based upon the track level information.
    int nProtonTrackCand(0);

    for(std::size_t i(0); i < slc->reco.trk.size(); ++i){

      if((int)i==muTrackInd) continue;

      auto const& trk = slc->reco.trk.at(i);
      bool Contained = fv_track.isContained(trk.end.x, trk.end.y, trk.end.z);

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
      bool Contained = fv_track.isContained(trk.end.x, trk.end.y, trk.end.z);

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

      Contained = fv_track.isContained(trk.end.x, trk.end.y, trk.end.z);

      int bp = trk.bestplane;
      Chi2Proton = trk.chi2pid[bp].chi2_proton;
      //Chi2Muon = trk.chi2pid[bp].chi2_muon;
      bool PassChi2 = Chi2Proton < 100;

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

      Contained = fv_track.isContained(trk.end.x, trk.end.y, trk.end.z);
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

  //====   Stub-based

  const Var varNStub([](const caf::SRSliceProxy* slc) -> int {

    return slc->reco.stub.size();

  });

  //==== Muon+Proton

  const Var varMuonProtonCosineTheta([](const caf::SRSliceProxy* slc) -> double {

    if( varMuonTrackInd(slc) >= 0 && varProtonTrackInd(slc) >= 0 ){

      auto const& trk_Muon = slc->reco.trk.at(varMuonTrackInd(slc));
      auto const& trk_Proton = slc->reco.trk.at(varProtonTrackInd(slc));

      TVector3 v3_Muon(trk_Muon.dir.x, trk_Muon.dir.y, trk_Muon.dir.z);
      TVector3 v3_Proton(trk_Proton.dir.x, trk_Proton.dir.y, trk_Proton.dir.z);

      return v3_Muon.Unit().Dot( v3_Proton.Unit() );

    }
    else{
      return -9999.;
    }

  });

  //==== Neutrino

  //====   Energy by Emu+Kp+Eb
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

    double E_Nu(-5.f);
    if( varMuonTrackInd(slc) >= 0 && varProtonTrackInd(slc) >= 0 ){
      double TruthE = varNeutrinoTruthE(slc);
      double RecoE = varNeutrinoCombinedEnergy(slc);
      return (RecoE-TruthE);
    }
    return E_Nu;

  });
  const Var varNeutrinoCombinedEnergyResidualFraction([](const caf::SRSliceProxy* slc) -> double {

    double E_Nu(-5.f);
    if( varMuonTrackInd(slc) >= 0 && varProtonTrackInd(slc) >= 0 ){
      double TruthE = varNeutrinoTruthE(slc);
      double RecoE = varNeutrinoCombinedEnergy(slc);
      return (RecoE-TruthE)/TruthE;
    }
    return E_Nu;

  });

  //==== https://s3.cern.ch/inspire-prod-files-9/93642a13c46438d97680971700e2013c
  const Var varNeutrinoQEEnergy([](const caf::SRSliceProxy* slc) -> double {

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

  const Var varNeutrinoQEEnergyResidual([](const caf::SRSliceProxy* slc) -> double {

    double E_Nu(-5.f);
    if( varMuonTrackInd(slc) >= 0 ){
      double TruthE = varNeutrinoTruthE(slc);
      double RecoE = varNeutrinoQEEnergy(slc);
      return (RecoE-TruthE);
    }
    return E_Nu;

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

