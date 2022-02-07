#pragma once

#include "sbnana/SBNAna/Cuts/Cuts.h"
#include "sbnana/SBNAna/Vars/Vars.h"
#include "sbnana/SBNAna/Cuts/VolumeDefinitions.h"

#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_Contants.h"
#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_Utilities.h"
#include "sbnana/SBNAna/Cuts/NumuCutsIcarus202106.h"
#include "TVector3.h"
#include <iostream>

using namespace std;
using namespace ana;

namespace ICARUSNumuXsec{

  //==== Slice variables

  const Var varCountSlice([](const caf::SRSliceProxy* slc) ->int {
    return 0.;
  });

  const Var varFMScore([](const caf::SRSliceProxy* slc) -> double {
    if(isnan(slc->fmatch.score)) return -2.;
    else if(slc->fmatch.score<0) return -1.;
    else return slc->fmatch.score;
  });

  const Var varFMTime([](const caf::SRSliceProxy* slc) -> double {
    if(isnan(slc->fmatch.time)) return -62.;
    else if(slc->fmatch.time<-50) return -61.;
    else return slc->fmatch.time;
  });

  const Var varNuScore([](const caf::SRSliceProxy* slc) -> double {
    if(isnan(slc->nu_score)) return -0.1;
    else if(slc->nu_score<0) return -0.05;
    else return slc->nu_score;
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


  //==== GENIE interaction code
  //==== https://internal.dunescience.org/doxygen/namespacesimb.html#a2cce734d1b71408bbc7d98d148ac4360
  const Var varGENIEIntCode([](const caf::SRSliceProxy* slc) -> int {
    //if(!isnan(slc->truth.genie_mode)){
      if(slc->truth.genie_mode<0) return -1;
      else if(slc->truth.genie_mode>13) return 14;
      else return slc->truth.genie_mode;
    //}
    //return -2;
  });

  //==== Truth variables

  //const Var varNeutrinoTruthE = SIMPLEVAR(truth.E);
  //const Var varTruthQ2 = SIMPLEVAR(truth.Q2);

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

  const Var varMuonTruthP([](const caf::SRSliceProxy* slc) -> double {

    double max_E(-999);
    for(std::size_t i(0); i < slc->truth.prim.size(); ++i){
      if( abs(slc->truth.prim.at(i).pdg)== 13 ){
        if(isnan(slc->truth.prim.at(i).genE)) continue;
        double this_E = slc->truth.prim.at(i).genE;
        if(this_E>max_E){
          max_E = this_E;
        }
      }
    }
    double max_P = sqrt(max_E*max_E-M_MUON*M_MUON);
    return max_P;

  });

  const Var varMuonTruthCosineTheta([](const caf::SRSliceProxy* slc) -> double {

    double max_E(-999);
    float costh(-5.f);

    for(std::size_t i(0); i < slc->truth.prim.size(); ++i){
      if( abs(slc->truth.prim.at(i).pdg)== 13 ){
        if(isnan(slc->truth.prim.at(i).genE)) continue;
        double this_E = slc->truth.prim.at(i).genE;
        if(this_E>max_E){
          max_E = this_E;
          TVector3 v3(slc->truth.prim.at(i).genp.x, slc->truth.prim.at(i).genp.y, slc->truth.prim.at(i).genp.z);
          costh = v3.CosTheta();
        }
      }
    }
    return costh;

  });

  const Var varMuonTruthNuMICosineTheta([](const caf::SRSliceProxy* slc) -> double {

    double max_E(-999);
    float costh(-5.f);

    for(std::size_t i(0); i < slc->truth.prim.size(); ++i){
      if( abs(slc->truth.prim.at(i).pdg)== 13 ){
        if(isnan(slc->truth.prim.at(i).genE)) continue;
        double this_E = slc->truth.prim.at(i).genE;
        if(this_E>max_E){
          max_E = this_E;
          TVector3 v3(slc->truth.prim.at(i).genp.x, slc->truth.prim.at(i).genp.y, slc->truth.prim.at(i).genp.z);
          double angleNuMI = v3.Angle(NuDirection_NuMI);
          costh = TMath::Cos(angleNuMI);
        }
      }
    }
    return costh;

  });

  const Var varProtonTruthP([](const caf::SRSliceProxy* slc) -> double {

    double max_E(-999);
    //if(isnan(slc->truth.prim.size())) return max_E;

    for(std::size_t i(0); i < slc->truth.prim.size(); ++i){
      if( abs(slc->truth.prim.at(i).pdg)== 2212 ){
        if(isnan(slc->truth.prim.at(i).genE)) continue;
        double this_E = slc->truth.prim.at(i).genE;
        if(this_E>max_E){
          max_E = this_E;
        }
      }
    }
    double max_P = max_E*max_E-M_PROTON*M_PROTON; // before sqrt, check if positive
    if(max_P<0) return 0.;
    else return sqrt(max_P);

  });

  const Var varProtonTruthCosineTheta([](const caf::SRSliceProxy* slc) -> double {

    double max_E(-999);
    float costh(-5.f);

    for(std::size_t i(0); i < slc->truth.prim.size(); ++i){
      if( abs(slc->truth.prim.at(i).pdg)== 2212 ){
        if(isnan(slc->truth.prim.at(i).genE)) continue;
        double this_E = slc->truth.prim.at(i).genE;
        if(this_E>max_E){
          max_E = this_E;
          TVector3 v3(slc->truth.prim.at(i).genp.x, slc->truth.prim.at(i).genp.y, slc->truth.prim.at(i).genp.z);
          costh = v3.CosTheta();
        }
      }
    }
    return costh;

  });

  const Var varProtonTruthNuMICosineTheta([](const caf::SRSliceProxy* slc) -> double {

    double max_E(-999);
    float costh(-5.f);

    for(std::size_t i(0); i < slc->truth.prim.size(); ++i){
      if( abs(slc->truth.prim.at(i).pdg)== 2212 ){
        if(isnan(slc->truth.prim.at(i).genE)) continue;
        double this_E = slc->truth.prim.at(i).genE;
        if(this_E>max_E){
          max_E = this_E;
          TVector3 v3(slc->truth.prim.at(i).genp.x, slc->truth.prim.at(i).genp.y, slc->truth.prim.at(i).genp.z);
          double angleNuMI = v3.Angle(NuDirection_NuMI);
          costh = TMath::Cos(angleNuMI);
        }
      }
    }
    return costh;

  });

  const Var varTruthMuonProtonCosineTheta([](const caf::SRSliceProxy* slc) -> double {

    double max_E(-999);
    TVector3 v3_Muon;
    //==== Muon
    bool MuonFound = false;
    for(std::size_t i(0); i < slc->truth.prim.size(); ++i){
      if( abs(slc->truth.prim.at(i).pdg)== 13 ){
        if(isnan(slc->truth.prim.at(i).genE)) continue;
        double this_E = slc->truth.prim.at(i).genE;
        if(this_E>max_E){
          max_E = this_E;
          TVector3 v3(slc->truth.prim.at(i).genp.x, slc->truth.prim.at(i).genp.y, slc->truth.prim.at(i).genp.z);
          v3_Muon = v3;
          MuonFound = true;
        }
      }
    }

    max_E = 999.;
    TVector3 v3_Proton;
    //==== Proton
    bool ProtonFound = false;
    for(std::size_t i(0); i < slc->truth.prim.size(); ++i){
      if( abs(slc->truth.prim.at(i).pdg)== 2212 ){
        if(isnan(slc->truth.prim.at(i).genE)) continue;
        double this_E = slc->truth.prim.at(i).genE;
        if(this_E>max_E){
          max_E = this_E;
          TVector3 v3(slc->truth.prim.at(i).genp.x, slc->truth.prim.at(i).genp.y, slc->truth.prim.at(i).genp.z);
          v3_Proton = v3;
          ProtonFound = true;
        }
      }
    }

    if(MuonFound&&ProtonFound){
      return TMath::Cos( v3_Muon.Angle(v3_Proton) );
    }
    else{
      return -999.;
    }


  });

  //==== For a given truth, find the matching Reco track

  //====   For a given true muon (truth_index), find a reco track whose best-matched is this muon

  const Var varTruthMuonMatchedTrackIndex([](const caf::SRSliceProxy* slc) -> double {

    //cout << "[varGetMatchedTrackForProton] called" << endl;

    int truth_index(-1);
    double max_E(-999);
    for(std::size_t i(0); i < slc->truth.prim.size(); ++i){
      if( abs(slc->truth.prim.at(i).pdg)== 13 ){
        if(isnan(slc->truth.prim.at(i).genE)) continue;
        double this_E = slc->truth.prim.at(i).genE;
        if(this_E>max_E){
          max_E = this_E;
          truth_index = i;
        }
      }
    }

    if(truth_index<0) return -999.;

    //cout << "[varGetMatchedTrackForProton] Truth proton genE = " << slc->truth.prim.at(truth_index).genE << endl;
    double truth_E = slc->truth.prim.at(truth_index).genE;

    int PTrackInd(-1);
    for(std::size_t i(0); i < slc->reco.trk.size(); ++i){
      auto const& trk = slc->reco.trk.at(i);
      if(isnan(trk.truth.p.genE)) continue;
      double this_matched_E = trk.truth.p.genE;
      if( abs(trk.truth.p.pdg)==13 && fabs(truth_E-this_matched_E)/truth_E < 0.01 ){
        PTrackInd = i;
      }
    }
    if(PTrackInd<0) return -999.;
    //cout << "[varGetMatchedTrackForProton] Truth proton genE = " << slc->truth.prim.at(truth_index).genE << endl;
    //cout << "[varGetMatchedTrackForProton] Found track's matched genE = " << slc->reco.trk.at(PTrackInd).truth.p.genE << endl;

    return PTrackInd;

  });

  const Var varTruthMuonMatchedTrackContainedness([](const caf::SRSliceProxy* slc) -> int {
    int muonTrackIndex = varTruthMuonMatchedTrackIndex(slc);
    //==== -1 : No track found
    //==== 0 : Exiting
    //==== 1 : Contained
    if(muonTrackIndex>=0){

      auto const& trk_Muon = slc->reco.trk.at(muonTrackIndex);
/*
      double XMargin = 25.;
      double YMargin = 25.;
      double ZMarginUp = 30.;
      double ZMarginDown = 50.;
*/

      bool isContained = fv.isContained(trk_Muon.end.x, trk_Muon.end.y, trk_Muon.end.z);
/*
      //==== Cryo0
      if( trk_Muon.end.x < 0 ){
        isContained = ( !isnan(trk_Muon.end.x) &&
                      ( trk_Muon.end.x < -71.1 - XMargin && trk_Muon.end.x > -369.33 + XMargin ) &&
                      !isnan(trk_Muon.end.y) &&
                      ( trk_Muon.end.y > -181.7 + YMargin && trk_Muon.end.y < 134.8 - YMargin ) &&
                      !isnan(trk_Muon.end.z) &&
                      ( trk_Muon.end.z > -895.95 + ZMarginUp && trk_Muon.end.z < 895.95 - ZMarginDown ) );
      }
      //==== Cryo1
      else{
        isContained = ( !isnan(trk_Muon.end.x) &&
                      ( trk_Muon.end.x > 71.1 + XMargin && trk_Muon.end.x < 369.33 - XMargin ) &&
                      !isnan(trk_Muon.end.y) &&
                      ( trk_Muon.end.y > -181.7 + YMargin && trk_Muon.end.y < 134.8 - YMargin ) &&
                      !isnan(trk_Muon.end.z) &&
                      ( trk_Muon.end.z > -895.95 + ZMarginUp && trk_Muon.end.z < 895.95 - ZMarginDown ) );
      }
*/
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

  const Var varTruthMuonMatchedTrackRangeP([](const caf::SRSliceProxy* slc) -> double {
    int muonTrackIndex = varTruthMuonMatchedTrackIndex(slc);
    if(muonTrackIndex>=0){
      auto const& trk_Muon = slc->reco.trk.at(muonTrackIndex);
      return trk_Muon.rangeP.p_muon;
    }
    else{
      return -999;
    }

  });
  const Var varTruthMuonMatchedTrackMCSP([](const caf::SRSliceProxy* slc) -> double {
    int muonTrackIndex = varTruthMuonMatchedTrackIndex(slc);
    if(muonTrackIndex>=0){
      auto const& trk_Muon = slc->reco.trk.at(muonTrackIndex);
      double P(-999.);
      if(!isnan(trk_Muon.mcsP.fwdP_muon)){
        P = trk_Muon.mcsP.fwdP_muon;
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

  //====   For a given true proton (truth_index), find a reco track whose best-matched is this muon

  const Var varTruthProtonMatchedTrackIndex([](const caf::SRSliceProxy* slc) -> double {

    //cout << "[varGetMatchedTrackForProton] called" << endl;

    int truth_index(-1);
    double max_E(-999);
    for(std::size_t i(0); i < slc->truth.prim.size(); ++i){
      if( abs(slc->truth.prim.at(i).pdg)== 2212 ){
        if(isnan(slc->truth.prim.at(i).genE)) continue;
        double this_E = slc->truth.prim.at(i).genE;
        if(this_E>max_E){
          max_E = this_E;
          truth_index = i;
        }
      }
    }

    if(truth_index<0) return -999.;

    //cout << "[varGetMatchedTrackForProton] Truth proton genE = " << slc->truth.prim.at(truth_index).genE << endl;
    double truth_E = slc->truth.prim.at(truth_index).genE;

    int PTrackInd(-1);
    for(std::size_t i(0); i < slc->reco.trk.size(); ++i){
      auto const& trk = slc->reco.trk.at(i);
      if(isnan(trk.truth.p.genE)) continue;
      double this_matched_E = trk.truth.p.genE;
      if( abs(trk.truth.p.pdg)==2212 && fabs(truth_E-this_matched_E)/truth_E < 0.01 ){
        PTrackInd = i;
      }
    }
    if(PTrackInd<0) return -999.;
    //cout << "[varGetMatchedTrackForProton] Truth proton genE = " << slc->truth.prim.at(truth_index).genE << endl;
    //cout << "[varGetMatchedTrackForProton] Found track's matched genE = " << slc->reco.trk.at(PTrackInd).truth.p.genE << endl;

    return PTrackInd;

  });

  const Var varTruthProtonMatchedTrackContainedness([](const caf::SRSliceProxy* slc) -> int {
    int ProtonTrackIndex = varTruthProtonMatchedTrackIndex(slc);
    //==== -1 : No track found
    //==== 0 : Exiting
    //==== 1 : Contained
    if(ProtonTrackIndex>=0){

      auto const& trk_Proton = slc->reco.trk.at(ProtonTrackIndex);
/*
      double XMargin = 25.;
      double YMargin = 25.;
      double ZMarginUp = 30.;
      double ZMarginDown = 50.;
*/
      bool isContained = fv.isContained(trk_Proton.end.x, trk_Proton.end.y, trk_Proton.end.z);
/*
      //==== Cryo0
      if( trk_Proton.end.x < 0 ){
        isContained = ( !isnan(trk_Proton.end.x) &&
                      ( trk_Proton.end.x < -71.1 - XMargin && trk_Proton.end.x > -369.33 + XMargin ) &&
                      !isnan(trk_Proton.end.y) &&
                      ( trk_Proton.end.y > -181.7 + YMargin && trk_Proton.end.y < 134.8 - YMargin ) &&
                      !isnan(trk_Proton.end.z) &&
                      ( trk_Proton.end.z > -895.95 + ZMarginUp && trk_Proton.end.z < 895.95 - ZMarginDown ) );
      }
      //==== Cryo1
      else{
        isContained = ( !isnan(trk_Proton.end.x) &&
                      ( trk_Proton.end.x > 71.1 + XMargin && trk_Proton.end.x < 369.33 - XMargin ) &&
                      !isnan(trk_Proton.end.y) &&
                      ( trk_Proton.end.y > -181.7 + YMargin && trk_Proton.end.y < 134.8 - YMargin ) &&
                      !isnan(trk_Proton.end.z) &&
                      ( trk_Proton.end.z > -895.95 + ZMarginUp && trk_Proton.end.z < 895.95 - ZMarginDown ) );
      }
*/
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


  //==== Reco variables

  //====   Muon

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
/*
      Contained = ( !isnan(trk.end.x) &&
                  ( trk.end.x < -71.1 - 25 && trk.end.x > -369.33 + 25 ) &&
                  !isnan(trk.end.y) &&
                  ( trk.end.y > -181.7 + 25 && trk.end.y < 134.8 - 25 ) &&
                  !isnan(trk.end.z) &&
                  ( trk.end.z > -895.95 + 30 && trk.end.z < 895.95 - 50 ) );
*/

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
/*
      Contained = ( !isnan(trk.end.x) &&
                  ( trk.end.x < -71.1 - 25 && trk.end.x > -369.33 + 25 ) &&
                  !isnan(trk.end.y) &&
                  ( trk.end.y > -181.7 + 25 && trk.end.y < 134.8 - 25 ) &&
                  !isnan(trk.end.z) &&
                  ( trk.end.z > -895.95 + 30 && trk.end.z < 895.95 - 50 ) );
*/
      if(Contained){
        p = trk.rangeP.p_muon;
      }
      else{
        p = trk.mcsP.fwdP_muon;
      }
/*
      if(p<0.1){
        cout << "[varMuonRecoP] Contained = " << Contained << endl;
        cout << "[varMuonRecoP] trk.rangeP.p_muon = " << trk.rangeP.p_muon << endl;
        cout << "[varMuonRecoP] trk.mcsP.fwdP_muon = " << trk.mcsP.fwdP_muon << endl;
        cout << "[varMuonRecoP] -> p = " << p << endl;
      }
*/

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

  const Var varProtonTrackInd([](const caf::SRSliceProxy* slc) -> int {

    //==== The (dis)qualification of a slice is based upon the track level information.
    float Atslc, Longest(0);
    float Chi2Proton;//, Chi2Muon;
    bool AtSlice, Contained;
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

      Contained = fv.isContained(trk.end.x, trk.end.y, trk.end.z);
/*
      Contained = ( !isnan(trk.end.x) &&
                  ( trk.end.x < -71.1 - 25 && trk.end.x > -369.33 + 25 ) &&
                  !isnan(trk.end.y) &&
                  ( trk.end.y > -181.7 + 25 && trk.end.y < 134.8 - 25 ) &&
                  !isnan(trk.end.z) &&
                  ( trk.end.z > -895.95 + 30 && trk.end.z < 895.95 - 50 ) );
*/
      int bp = trk.bestplane;
      Chi2Proton = trk.chi2pid[bp].chi2_proton;
      //Chi2Muon = trk.chi2pid[bp].chi2_muon;
      bool PassChi2 = Chi2Proton < 60;

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
/*
      Contained = ( !isnan(trk.end.x) &&
                  ( trk.end.x < -71.1 - 25 && trk.end.x > -369.33 + 25 ) &&
                  !isnan(trk.end.y) &&
                  ( trk.end.y > -181.7 + 25 && trk.end.y < 134.8 - 25 ) &&
                  !isnan(trk.end.z) &&
                  ( trk.end.z > -895.95 + 30 && trk.end.z < 895.95 - 50 ) );
*/
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

