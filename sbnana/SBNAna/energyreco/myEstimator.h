#ifndef myEstimator_h
#define myEstimator_h

#include "sbnana/SBNAna/energyreco/myMuonSelection.h"
#include "sbnana/SBNAna/energyreco/myProtonSelection.h"
#include "TVector3.h"

#define bindingE 0.040

using namespace std;

//==== Muon estimator

//====   Reco

const Var varMuonTrackRangeP([](const caf::SRSliceProxy* slc) -> double {
  int muonTrackIndex = varMuonTrackIndex(slc);
  if(muonTrackIndex>=0){
    auto const& trk_Muon = slc->reco.trk.at(muonTrackIndex);
    return trk_Muon.rangeP.p_muon;
  }
  else{
    return -999;
  }
});

const Var varMuonTrackMCSP([](const caf::SRSliceProxy* slc) -> double {
  int muonTrackIndex = varMuonTrackIndex(slc);
  if(muonTrackIndex>=0){
    auto const& trk_Muon = slc->reco.trk.at(muonTrackIndex);
    return trk_Muon.mcsP.fwdP_muon;
  }   
  else{
    return -999;
  }   
}); 

const Var varMuonTrackCombinedP([](const caf::SRSliceProxy* slc) -> double {
  int muonTrackIndex = varMuonTrackIndex(slc);
  if(muonTrackIndex>=0){
    int mtc = varMuonTrackContainedness(slc); // 0 : exiting, 1 : contained
    if(mtc==0) return varMuonTrackMCSP(slc);
    else if(mtc==1) return varMuonTrackRangeP(slc);
    else{
      cout << "[varMuonTrackCombinedP] wtf?" << endl;
      exit(EXIT_FAILURE);
      return varMuonTrackMCSP(slc); //
    }
  }   
  else{
    return -999;
  }   
}); 

const Var varMuonTrackCaloP([](const caf::SRSliceProxy* slc) -> double {
  int muonTrackIndex = varMuonTrackIndex(slc);
  if(muonTrackIndex>=0){
    auto const& trk_Muon = slc->reco.trk.at(muonTrackIndex);

    double best_calo = 0.;
    if(trk_Muon.bestplane == 0){
      best_calo = trk_Muon.calo0.ke;
    }
    else if(trk_Muon.bestplane == 1){
      best_calo = trk_Muon.calo1.ke;
    }
    else{
      best_calo = trk_Muon.calo2.ke;
    }

    double KE_Muon = best_calo/1000.;
    double E_Muon = KE_Muon+M_MUON;
    double P_Muon = sqrt( E_Muon*E_Muon - M_MUON*M_MUON );
    return P_Muon;

  }
  else{
    return -999;
  }
});

const Var varMuonTrackPlane0CaloP([](const caf::SRSliceProxy* slc) -> double {
  int muonTrackIndex = varMuonTrackIndex(slc);
  if(muonTrackIndex>=0){
    auto const& trk_Muon = slc->reco.trk.at(muonTrackIndex);

    double best_calo = trk_Muon.calo0.ke;
    double KE_Muon = best_calo/1000.;
    double E_Muon = KE_Muon+M_MUON;
    double P_Muon = sqrt( E_Muon*E_Muon - M_MUON*M_MUON );
    return P_Muon;

  }
  else{
    return -999;
  }
});

const Var varMuonTrackPlane1CaloP([](const caf::SRSliceProxy* slc) -> double {
  int muonTrackIndex = varMuonTrackIndex(slc);
  if(muonTrackIndex>=0){ 
    auto const& trk_Muon = slc->reco.trk.at(muonTrackIndex);
      
    double best_calo = trk_Muon.calo1.ke;
    double KE_Muon = best_calo/1000.;
    double E_Muon = KE_Muon+M_MUON;
    double P_Muon = sqrt( E_Muon*E_Muon - M_MUON*M_MUON );
    return P_Muon;

  }
  else{
    return -999;
  }
});

const Var varMuonTrackPlane2CaloP([](const caf::SRSliceProxy* slc) -> double {
  int muonTrackIndex = varMuonTrackIndex(slc);
  if(muonTrackIndex>=0){ 
    auto const& trk_Muon = slc->reco.trk.at(muonTrackIndex);
      
    double best_calo = trk_Muon.calo2.ke;
    double KE_Muon = best_calo/1000.;
    double E_Muon = KE_Muon+M_MUON;
    double P_Muon = sqrt( E_Muon*E_Muon - M_MUON*M_MUON );
    return P_Muon;

  }
  else{
    return -999;
  }
});

const Var varMuonTrackCaloKE([](const caf::SRSliceProxy* slc) -> double {
  int muonTrackIndex = varMuonTrackIndex(slc);
  if(muonTrackIndex>=0){
    auto const& trk_Muon = slc->reco.trk.at(muonTrackIndex);

    double best_calo = 0.;
    if(trk_Muon.bestplane == 0){
      best_calo = trk_Muon.calo0.ke;
    }
    else if(trk_Muon.bestplane == 1){
      best_calo = trk_Muon.calo1.ke;
    }
    else{
      best_calo = trk_Muon.calo2.ke;
    }

    double KE_Muon = best_calo/1000.;
    return KE_Muon;

  }
  else{
    return -999;
  }
});

//====   Truth
//====     track.p has the information of the best match SRTrueParticle
//====     https://github.com/SBNSoftware/sbnanaobj/blob/733b185944bf4565ddb451c9aee3ade8715d4114/sbnanaobj/StandardRecord/SRTrackTruth.h#L33

const Var varMuonTrackMatchedTruthE([](const caf::SRSliceProxy* slc) -> double {
  int muonTrackIndex = varMuonTrackIndex(slc);
  if(muonTrackIndex>=0){
    auto const& trk_Muon = slc->reco.trk.at(muonTrackIndex);
    double E_Muon = trk_Muon.truth.p.genE;
    return E_Muon;
  }
  else{
    return -999.;
  }
});

const Var varMuonTrackMatchedTruthP([](const caf::SRSliceProxy* slc) -> double {
  int muonTrackIndex = varMuonTrackIndex(slc);
  if(muonTrackIndex>=0){ 
    auto const& trk_Muon = slc->reco.trk.at(muonTrackIndex);
    double E_Muon = trk_Muon.truth.p.genE;
    double P_Muon = sqrt( E_Muon*E_Muon - M_MUON*M_MUON );
    return P_Muon;
  }
  else{
    return -999.;
  }
});

const Var varMuonTrackMatchedTruthPDG([](const caf::SRSliceProxy* slc) -> double {
  int muonTrackIndex = varMuonTrackIndex(slc);
  if(muonTrackIndex>=0){
    return slc->reco.trk.at(muonTrackIndex).truth.p.pdg;
  }
  else{
    return -999.;
  }
});

//====   Residual study

const Var varMuonTrackRangePResidual([](const caf::SRSliceProxy* slc) -> double {
  double P_Muon_Reco = varMuonTrackRangeP(slc);
  double P_Muon_Truth = varMuonTrackMatchedTruthP(slc);
  double P_Muon_Diff = P_Muon_Reco-P_Muon_Truth;
  return P_Muon_Diff;
});

const Var varMuonTrackMCSPResidual([](const caf::SRSliceProxy* slc) -> double {
  double P_Muon_Reco = varMuonTrackMCSP(slc);
  double P_Muon_Truth = varMuonTrackMatchedTruthP(slc);
  double P_Muon_Diff = P_Muon_Reco-P_Muon_Truth;
  return P_Muon_Diff;
});

const Var varMuonTrackCombinedPResidual([](const caf::SRSliceProxy* slc) -> double {
  double P_Muon_Reco = varMuonTrackCombinedP(slc);
  double P_Muon_Truth = varMuonTrackMatchedTruthP(slc);
  double P_Muon_Diff = P_Muon_Reco-P_Muon_Truth;
  return P_Muon_Diff;
});

const Var varMuonTrackCaloPResidual([](const caf::SRSliceProxy* slc) -> double {
  double P_Muon_Reco = varMuonTrackCaloP(slc);
  double P_Muon_Truth = varMuonTrackMatchedTruthP(slc);
  double P_Muon_Diff = P_Muon_Reco-P_Muon_Truth;
  return P_Muon_Diff;
});

const Var varMuonTrackPlane0CaloPResidual([](const caf::SRSliceProxy* slc) -> double {
  double P_Muon_Reco = varMuonTrackPlane0CaloP(slc);
  double P_Muon_Truth = varMuonTrackMatchedTruthP(slc);
  double P_Muon_Diff = P_Muon_Reco-P_Muon_Truth;
  return P_Muon_Diff;
});

const Var varMuonTrackPlane1CaloPResidual([](const caf::SRSliceProxy* slc) -> double {
  double P_Muon_Reco = varMuonTrackPlane1CaloP(slc);
  double P_Muon_Truth = varMuonTrackMatchedTruthP(slc);
  double P_Muon_Diff = P_Muon_Reco-P_Muon_Truth;
  return P_Muon_Diff;
});

const Var varMuonTrackPlane2CaloPResidual([](const caf::SRSliceProxy* slc) -> double {
  double P_Muon_Reco = varMuonTrackPlane2CaloP(slc);
  double P_Muon_Truth = varMuonTrackMatchedTruthP(slc);
  double P_Muon_Diff = P_Muon_Reco-P_Muon_Truth;
  return P_Muon_Diff;
});

const Var varMuonTrackRangePResidualFraction([](const caf::SRSliceProxy* slc) -> double {
  double P_Muon_Reco = varMuonTrackRangeP(slc);
  double P_Muon_Truth = varMuonTrackMatchedTruthP(slc);
  double P_Muon_Diff = P_Muon_Reco-P_Muon_Truth;
  return P_Muon_Diff/P_Muon_Truth;
});

const Var varMuonTrackMCSPResidualFraction([](const caf::SRSliceProxy* slc) -> double {
  double P_Muon_Reco = varMuonTrackMCSP(slc);
  double P_Muon_Truth = varMuonTrackMatchedTruthP(slc);
  double P_Muon_Diff = P_Muon_Reco-P_Muon_Truth;
  return P_Muon_Diff/P_Muon_Truth;
});

const Var varMuonTrackCombinedPResidualFraction([](const caf::SRSliceProxy* slc) -> double {
  double P_Muon_Reco = varMuonTrackCombinedP(slc);
  double P_Muon_Truth = varMuonTrackMatchedTruthP(slc);
  double P_Muon_Diff = P_Muon_Reco-P_Muon_Truth;
  return P_Muon_Diff/P_Muon_Truth;
});

const Var varMuonTrackCaloPResidualFraction([](const caf::SRSliceProxy* slc) -> double {
  double P_Muon_Reco = varMuonTrackCaloP(slc);
  double P_Muon_Truth = varMuonTrackMatchedTruthP(slc);
  double P_Muon_Diff = P_Muon_Reco-P_Muon_Truth;
  return P_Muon_Diff/P_Muon_Truth;
});

const Var varMuonTrackPlane0CaloPResidualFraction([](const caf::SRSliceProxy* slc) -> double {
  double P_Muon_Reco = varMuonTrackPlane0CaloP(slc);
  double P_Muon_Truth = varMuonTrackMatchedTruthP(slc);
  double P_Muon_Diff = P_Muon_Reco-P_Muon_Truth;
  return P_Muon_Diff/P_Muon_Truth;
});

const Var varMuonTrackPlane1CaloPResidualFraction([](const caf::SRSliceProxy* slc) -> double {
  double P_Muon_Reco = varMuonTrackPlane1CaloP(slc);
  double P_Muon_Truth = varMuonTrackMatchedTruthP(slc);
  double P_Muon_Diff = P_Muon_Reco-P_Muon_Truth;
  return P_Muon_Diff/P_Muon_Truth;
});

const Var varMuonTrackPlane2CaloPResidualFraction([](const caf::SRSliceProxy* slc) -> double {
  double P_Muon_Reco = varMuonTrackPlane2CaloP(slc);
  double P_Muon_Truth = varMuonTrackMatchedTruthP(slc);
  double P_Muon_Diff = P_Muon_Reco-P_Muon_Truth;
  return P_Muon_Diff/P_Muon_Truth;
});

//==== Proton

//====   Reco

const Var varProtonTrackRangeP([](const caf::SRSliceProxy* slc) -> double {
  int protonTrackIndex = varProtonTrackIndex(slc);
  if(protonTrackIndex>=0){
    auto const& trk_Proton = slc->reco.trk.at(protonTrackIndex);
    return trk_Proton.rangeP.p_proton;
  }
  else{
    return -999;
  }
});

const Var varProtonTrackMCSP([](const caf::SRSliceProxy* slc) -> double {
  int protonTrackIndex = varProtonTrackIndex(slc);
  if(protonTrackIndex>=0){
    auto const& trk_Proton = slc->reco.trk.at(protonTrackIndex);
    return trk_Proton.mcsP.fwdP_proton;
  }
  else{
    return -999;
  }
});

const Var varProtonTrackCombinedP([](const caf::SRSliceProxy* slc) -> double {
  int protonTrackIndex = varProtonTrackIndex(slc);
  if(protonTrackIndex>=0){
    int mtc = varProtonTrackContainedness(slc); // 0 : exiting, 1 : contained
    if(mtc==0) return varProtonTrackMCSP(slc);
    else if(mtc==1) return varProtonTrackRangeP(slc);
    else{
      cout << "[varProtonTrackCombinedP] wtf?" << endl;
      exit(EXIT_FAILURE);
      return varProtonTrackMCSP(slc); //
    }
  }
  else{
    return -999;
  }
});

const Var varProtonTrackCaloP([](const caf::SRSliceProxy* slc) -> double {
  int protonTrackIndex = varProtonTrackIndex(slc);
  if(protonTrackIndex>=0){
    auto const& trk_Proton = slc->reco.trk.at(protonTrackIndex);

    double best_calo = 0.;
    if(trk_Proton.bestplane == 0){
      best_calo = trk_Proton.calo0.ke;
    }
    else if(trk_Proton.bestplane == 1){
      best_calo = trk_Proton.calo1.ke;
    }
    else{
      best_calo = trk_Proton.calo2.ke;
    }

    double KE_Proton = best_calo/1000.;
    double E_Proton = KE_Proton+M_PROTON;
    double p_proton = sqrt( E_Proton*E_Proton - M_PROTON*M_PROTON );
    return p_proton;

  }
  else{
    return -999;
  }
});

const Var varProtonTrackPlane0CaloP([](const caf::SRSliceProxy* slc) -> double {
  int protonTrackIndex = varProtonTrackIndex(slc);
  if(protonTrackIndex>=0){
    auto const& trk_Proton = slc->reco.trk.at(protonTrackIndex);

    double best_calo = trk_Proton.calo0.ke;
    double KE_Proton = best_calo/1000.;
    double E_Proton = KE_Proton+M_PROTON;
    double p_proton = sqrt( E_Proton*E_Proton - M_PROTON*M_PROTON );
    return p_proton;

  }
  else{
    return -999;
  }
});

const Var varProtonTrackPlane1CaloP([](const caf::SRSliceProxy* slc) -> double {
  int protonTrackIndex = varProtonTrackIndex(slc);
  if(protonTrackIndex>=0){
    auto const& trk_Proton = slc->reco.trk.at(protonTrackIndex);

    double best_calo = trk_Proton.calo1.ke;
    double KE_Proton = best_calo/1000.;
    double E_Proton = KE_Proton+M_PROTON;
    double p_proton = sqrt( E_Proton*E_Proton - M_PROTON*M_PROTON );
    return p_proton;

  }
  else{
    return -999;
  }
});

const Var varProtonTrackPlane2CaloP([](const caf::SRSliceProxy* slc) -> double {
  int protonTrackIndex = varProtonTrackIndex(slc);
  if(protonTrackIndex>=0){
    auto const& trk_Proton = slc->reco.trk.at(protonTrackIndex);

    double best_calo = trk_Proton.calo2.ke;
    double KE_Proton = best_calo/1000.;
    double E_Proton = KE_Proton+M_PROTON;
    double p_proton = sqrt( E_Proton*E_Proton - M_PROTON*M_PROTON );
    return p_proton;

  }
  else{
    return -999;
  }
});

const Var varProtonTrackCaloKE([](const caf::SRSliceProxy* slc) -> double {
  int protonTrackIndex = varProtonTrackIndex(slc);
  if(protonTrackIndex>=0){
    auto const& trk_Proton = slc->reco.trk.at(protonTrackIndex);

    double best_calo = 0.;
    if(trk_Proton.bestplane == 0){
      best_calo = trk_Proton.calo0.ke;
    }
    else if(trk_Proton.bestplane == 1){
      best_calo = trk_Proton.calo1.ke;
    }
    else{
      best_calo = trk_Proton.calo2.ke;
    }

    double KE_Proton = best_calo/1000.;
    return KE_Proton;

  }
  else{
    return -999;
  }
});

//====   Truth
//====     track.p has the information of the best match SRTrueParticle
//====     https://github.com/SBNSoftware/sbnanaobj/blob/733b185944bf4565ddb451c9aee3ade8715d4114/sbnanaobj/StandardRecord/SRTrackTruth.h#L33

const Var varProtonTrackMatchedTruthE([](const caf::SRSliceProxy* slc) -> double {
  int protonTrackIndex = varProtonTrackIndex(slc);
  if(protonTrackIndex>=0){
    auto const& trk_Proton = slc->reco.trk.at(protonTrackIndex);
    double E_Proton = trk_Proton.truth.p.genE;
    return E_Proton;
  }
  else{
    return -999.;
  }
});

const Var varProtonTrackMatchedTruthP([](const caf::SRSliceProxy* slc) -> double {
  int protonTrackIndex = varProtonTrackIndex(slc);
  if(protonTrackIndex>=0){
    auto const& trk_Proton = slc->reco.trk.at(protonTrackIndex);
    double E_Proton = trk_Proton.truth.p.genE;
    double p_proton = sqrt( E_Proton*E_Proton - M_PROTON*M_PROTON );
    return p_proton;
  }
  else{
    return -999.;
  }
});

const Var varProtonTrackMatchedTruthPDG([](const caf::SRSliceProxy* slc) -> int {
  int protonTrackIndex = varProtonTrackIndex(slc);
  if(protonTrackIndex>=0){
    return slc->reco.trk.at(protonTrackIndex).truth.p.pdg;
  }
  else{
    return -999999;
  }
});

//====   Residual study

const Var varProtonTrackRangePResidual([](const caf::SRSliceProxy* slc) -> double {
  double p_proton_Reco = varProtonTrackRangeP(slc);
  double p_proton_Truth = varProtonTrackMatchedTruthP(slc);
  double p_proton_Diff = p_proton_Reco-p_proton_Truth;
  return p_proton_Diff;
});

const Var varProtonTrackMCSPResidual([](const caf::SRSliceProxy* slc) -> double {
  double p_proton_Reco = varProtonTrackMCSP(slc);
  double p_proton_Truth = varProtonTrackMatchedTruthP(slc);
  double p_proton_Diff = p_proton_Reco-p_proton_Truth;
  return p_proton_Diff;
});

const Var varProtonTrackCombinedPResidual([](const caf::SRSliceProxy* slc) -> double {
  double p_proton_Reco = varProtonTrackCombinedP(slc);
  double p_proton_Truth = varProtonTrackMatchedTruthP(slc);
  double p_proton_Diff = p_proton_Reco-p_proton_Truth;
  return p_proton_Diff;
});

const Var varProtonTrackCaloPResidual([](const caf::SRSliceProxy* slc) -> double {
  double p_proton_Reco = varProtonTrackCaloP(slc);
  double p_proton_Truth = varProtonTrackMatchedTruthP(slc);
  double p_proton_Diff = p_proton_Reco-p_proton_Truth;
  return p_proton_Diff;
});

const Var varProtonTrackPlane0CaloPResidual([](const caf::SRSliceProxy* slc) -> double {
  double p_proton_Reco = varProtonTrackPlane0CaloP(slc);
  double p_proton_Truth = varProtonTrackMatchedTruthP(slc);
  double p_proton_Diff = p_proton_Reco-p_proton_Truth;
  return p_proton_Diff;
});

const Var varProtonTrackPlane1CaloPResidual([](const caf::SRSliceProxy* slc) -> double {
  double p_proton_Reco = varProtonTrackPlane1CaloP(slc);
  double p_proton_Truth = varProtonTrackMatchedTruthP(slc);
  double p_proton_Diff = p_proton_Reco-p_proton_Truth;
  return p_proton_Diff;
});

const Var varProtonTrackPlane2CaloPResidual([](const caf::SRSliceProxy* slc) -> double {
  double p_proton_Reco = varProtonTrackPlane2CaloP(slc);
  double p_proton_Truth = varProtonTrackMatchedTruthP(slc);
  double p_proton_Diff = p_proton_Reco-p_proton_Truth;
  return p_proton_Diff;
});


const Var varProtonTrackRangePResidualFraction([](const caf::SRSliceProxy* slc) -> double {
  double p_proton_Reco = varProtonTrackRangeP(slc);
  double p_proton_Truth = varProtonTrackMatchedTruthP(slc);
  double p_proton_Diff = p_proton_Reco-p_proton_Truth;
  return p_proton_Diff/p_proton_Truth;
});

const Var varProtonTrackMCSPResidualFraction([](const caf::SRSliceProxy* slc) -> double {
  double p_proton_Reco = varProtonTrackMCSP(slc);
  double p_proton_Truth = varProtonTrackMatchedTruthP(slc);
  double p_proton_Diff = p_proton_Reco-p_proton_Truth;
  return p_proton_Diff/p_proton_Truth;
});

const Var varProtonTrackCombinedPResidualFraction([](const caf::SRSliceProxy* slc) -> double {
  double p_proton_Reco = varProtonTrackCombinedP(slc);
  double p_proton_Truth = varProtonTrackMatchedTruthP(slc);
  double p_proton_Diff = p_proton_Reco-p_proton_Truth;
  return p_proton_Diff/p_proton_Truth;
});

const Var varProtonTrackCaloPResidualFraction([](const caf::SRSliceProxy* slc) -> double {
  double p_proton_Reco = varProtonTrackCaloP(slc);
  double p_proton_Truth = varProtonTrackMatchedTruthP(slc);
  double p_proton_Diff = p_proton_Reco-p_proton_Truth;
  return p_proton_Diff/p_proton_Truth;
});

const Var varProtonTrackPlane0CaloPResidualFraction([](const caf::SRSliceProxy* slc) -> double {
  double p_proton_Reco = varProtonTrackPlane0CaloP(slc);
  double p_proton_Truth = varProtonTrackMatchedTruthP(slc);
  double p_proton_Diff = p_proton_Reco-p_proton_Truth;
  return p_proton_Diff/p_proton_Truth;
});

const Var varProtonTrackPlane1CaloPResidualFraction([](const caf::SRSliceProxy* slc) -> double {
  double p_proton_Reco = varProtonTrackPlane1CaloP(slc);
  double p_proton_Truth = varProtonTrackMatchedTruthP(slc);
  double p_proton_Diff = p_proton_Reco-p_proton_Truth;
  return p_proton_Diff/p_proton_Truth;
});

const Var varProtonTrackPlane2CaloPResidualFraction([](const caf::SRSliceProxy* slc) -> double {
  double p_proton_Reco = varProtonTrackPlane2CaloP(slc);
  double p_proton_Truth = varProtonTrackMatchedTruthP(slc);
  double p_proton_Diff = p_proton_Reco-p_proton_Truth;
  return p_proton_Diff/p_proton_Truth;
});

//==== Neutrino

//====   Reco

const Var varNeutrinoCaloEnergy([](const caf::SRSliceProxy* slc) -> double {

  double KE_Muon = varMuonTrackCaloKE(slc); // calo-based
  //==== calo-based
  double E_Muon = KE_Muon+M_MUON;

  double KE_Proton = varProtonTrackCaloKE(slc);

  double Eb = bindingE; // TODO binding energy

  double E_Nu = E_Muon + KE_Proton + Eb;

  return E_Nu;

});

const Var varNeutrinoCombinedEnergy([](const caf::SRSliceProxy* slc) -> double {

  double P_Muon = varMuonTrackCombinedP(slc); // comdinbed
  //==== combined
  double E_Muon = sqrt( P_Muon*P_Muon + M_MUON*M_MUON );

  double KE_Proton = varProtonTrackCaloKE(slc);

  double Eb = bindingE; // TODO binding energy will be added later

  double E_Nu = E_Muon + KE_Proton + Eb;

  return E_Nu;

});

//==== https://s3.cern.ch/inspire-prod-files-9/93642a13c46438d97680971700e2013c
const Var varNeutrinoQE([](const caf::SRSliceProxy* slc) -> double {

  int muonTrackIndex = varMuonTrackIndex(slc);
  if(muonTrackIndex>=0){
    auto const& trk_Muon = slc->reco.trk.at(muonTrackIndex);

    double Eb = bindingE; // TODO binding energy will be added later
    double mu_costh = trk_Muon.costh;
    double mu_p = varMuonTrackCombinedP(slc);
    double mu_E = sqrt(mu_p*mu_p+M_MUON*M_MUON);

    double EQE_num = M_PROTON*M_PROTON - (M_NEUTRON-Eb)*(M_NEUTRON-Eb) - M_MUON*M_MUON + 2.*(M_NEUTRON-Eb)*mu_E;
    double EQE_den = 2.*(M_NEUTRON - Eb - mu_E + mu_p * mu_costh);

    return EQE_num/EQE_den;

  }
  else{
    return -999.;
  }

});

//====   Fake-Reco, using truth

const Var varNeutrinoFakeRecoEnergy([](const caf::SRSliceProxy* slc) -> double {

  double E_Muon = varMuonTrackMatchedTruthE(slc);

  double E_Proton = varProtonTrackMatchedTruthE(slc);
  double KE_Proton = E_Proton - M_PROTON;

  double Eb = bindingE; // TODO binding energy will be added later

  double E_Nu = E_Muon + KE_Proton + Eb;

  return E_Nu;

});

//====   Truth

const Var varNeutrinoTruthE = SIMPLEVAR(truth.E);

//====   Residual

const Var varNeutrinoCaloEnergyResidual([](const caf::SRSliceProxy* slc) -> double {

  double E_Neutrino_Reco = varNeutrinoCaloEnergy(slc);
  double E_Neutrino_Truth = varNeutrinoTruthE(slc);
  double E_Neutrino_Diff = E_Neutrino_Reco-E_Neutrino_Truth;
  return E_Neutrino_Diff;

});

const Var varNeutrinoCombinedEnergyResidual([](const caf::SRSliceProxy* slc) -> double {
  
  double E_Neutrino_Reco = varNeutrinoCombinedEnergy(slc);
  double E_Neutrino_Truth = varNeutrinoTruthE(slc);
  double E_Neutrino_Diff = E_Neutrino_Reco-E_Neutrino_Truth;
  return E_Neutrino_Diff;

});

const Var varNeutrinoQEResidual([](const caf::SRSliceProxy* slc) -> double {
  
  double E_Neutrino_Reco = varNeutrinoQE(slc);
  double E_Neutrino_Truth = varNeutrinoTruthE(slc);
  double E_Neutrino_Diff = E_Neutrino_Reco-E_Neutrino_Truth;
  return E_Neutrino_Diff;

});

const Var varNeutrinoCaloEnergyResidualFraction([](const caf::SRSliceProxy* slc) -> double {

  double E_Neutrino_Reco = varNeutrinoCaloEnergy(slc);
  double E_Neutrino_Truth = varNeutrinoTruthE(slc);
  double E_Neutrino_Diff = E_Neutrino_Reco-E_Neutrino_Truth;
  return E_Neutrino_Diff/E_Neutrino_Truth;

});

const Var varNeutrinoCombinedEnergyResidualFraction([](const caf::SRSliceProxy* slc) -> double {

  double E_Neutrino_Reco = varNeutrinoCombinedEnergy(slc);
  double E_Neutrino_Truth = varNeutrinoTruthE(slc);
  double E_Neutrino_Diff = E_Neutrino_Reco-E_Neutrino_Truth;
  return E_Neutrino_Diff/E_Neutrino_Truth;

});

const Var varNeutrinoQEResidualFraction([](const caf::SRSliceProxy* slc) -> double {

  double E_Neutrino_Reco = varNeutrinoQE(slc);
  double E_Neutrino_Truth = varNeutrinoTruthE(slc);
  double E_Neutrino_Diff = E_Neutrino_Reco-E_Neutrino_Truth;
  return E_Neutrino_Diff/E_Neutrino_Truth;

});

const Var varMuonTruthP([](const caf::SRSliceProxy* slc) -> double {

  double max_E(-999);
  for(std::size_t i(0); i < slc->truth.prim.size(); ++i){
    if( abs(slc->truth.prim.at(i).pdg)== 13 ){
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

const Var varProtonTruthP([](const caf::SRSliceProxy* slc) -> double {

  double max_E(-999);
  for(std::size_t i(0); i < slc->truth.prim.size(); ++i){
    if( abs(slc->truth.prim.at(i).pdg)== 2212 ){
      double this_E = slc->truth.prim.at(i).genE;
      if(this_E>max_E){
        max_E = this_E;
      }
    }
  }
  double max_P = sqrt(max_E*max_E-M_PROTON*M_PROTON);
  return max_P;

});

const Var varProtonTruthCosineTheta([](const caf::SRSliceProxy* slc) -> double {

  double max_E(-999);
  float costh(-5.f);

  for(std::size_t i(0); i < slc->truth.prim.size(); ++i){
    if( abs(slc->truth.prim.at(i).pdg)== 2212 ){
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


//==== TEST

//==== For a given true muon (truth_index), find a reco track whose best-matched is this muon
const Var varTruthMuonMatchedTrackIndex([](const caf::SRSliceProxy* slc) -> double {

  //cout << "[varGetMatchedTrackForProton] called" << endl;

  int truth_index(-1);
  double max_E(-999);
  for(std::size_t i(0); i < slc->truth.prim.size(); ++i){
    if( abs(slc->truth.prim.at(i).pdg)== 13 ){
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

    double XMargin = 25.;
    double YMargin = 25.;
    double ZMarginUp = 30.;
    double ZMarginDown = 50.;

    bool isContained;
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
    float Chi2Proton;
    if(trk.bestplane == 0){
      Chi2Proton = trk.chi2pid0.chi2_proton;
    }
    else if(trk.bestplane == 1){
      Chi2Proton = trk.chi2pid1.chi2_proton;
    }
    else{
      Chi2Proton = trk.chi2pid2.chi2_proton;
    }

    return Chi2Proton;

  }
  else{
    return 999999;
  }
});

const Var varTruthMuonMatchedTrackNormalizedChi2Proton([](const caf::SRSliceProxy* slc) -> double {
  int muonTrackIndex = varTruthMuonMatchedTrackIndex(slc);
  if(muonTrackIndex>=0){

    auto const& trk = slc->reco.trk.at(muonTrackIndex);
    float Chi2Proton;
    int Chi2Ndof;
    if(trk.bestplane == 0){
      Chi2Proton = trk.chi2pid0.chi2_proton;
      Chi2Ndof = trk.chi2pid0.pid_ndof;
    }
    else if(trk.bestplane == 1){
      Chi2Proton = trk.chi2pid1.chi2_proton;
      Chi2Ndof = trk.chi2pid1.pid_ndof;
    }
    else{
      Chi2Proton = trk.chi2pid2.chi2_proton;
      Chi2Ndof = trk.chi2pid2.pid_ndof;
    }

    return Chi2Proton/Chi2Ndof;

  }
  else{
    return 999999;
  }
});

const Var varTruthMuonMatchedTrackChi2Muon([](const caf::SRSliceProxy* slc) -> double {
  int muonTrackIndex = varTruthMuonMatchedTrackIndex(slc);
  if(muonTrackIndex>=0){
    
    auto const& trk = slc->reco.trk.at(muonTrackIndex);
    float Chi2Muon;
    if(trk.bestplane == 0){
      Chi2Muon = trk.chi2pid0.chi2_muon;
    }
    else if(trk.bestplane == 1){
      Chi2Muon = trk.chi2pid1.chi2_muon;
    }
    else{
      Chi2Muon = trk.chi2pid2.chi2_muon;
    }

    return Chi2Muon;

  }
  else{
    return 999999;
  }

});

const Var varTruthMuonMatchedTrackNormalizedChi2Muon([](const caf::SRSliceProxy* slc) -> double {
  int muonTrackIndex = varTruthMuonMatchedTrackIndex(slc);
  if(muonTrackIndex>=0){

    auto const& trk = slc->reco.trk.at(muonTrackIndex);
    float Chi2Muon;
    int Chi2Ndof;
    if(trk.bestplane == 0){
      Chi2Muon = trk.chi2pid0.chi2_muon;
      Chi2Ndof = trk.chi2pid0.pid_ndof;
    }
    else if(trk.bestplane == 1){
      Chi2Muon = trk.chi2pid1.chi2_muon;
      Chi2Ndof = trk.chi2pid1.pid_ndof;
    }
    else{
      Chi2Muon = trk.chi2pid2.chi2_muon;
      Chi2Ndof = trk.chi2pid2.pid_ndof;
    }

    return Chi2Muon/Chi2Ndof;

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
    return trk_Muon.mcsP.fwdP_muon;
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

const Var varTruthProtonMatchedTrackIndex([](const caf::SRSliceProxy* slc) -> double {

  //cout << "[varGetMatchedTrackForProton] called" << endl;

  int truth_index(-1);
  double max_E(-999);
  for(std::size_t i(0); i < slc->truth.prim.size(); ++i){
    if( abs(slc->truth.prim.at(i).pdg)== 2212 ){
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

    double XMargin = 25.;
    double YMargin = 25.;
    double ZMarginUp = 30.;
    double ZMarginDown = 50.;

    bool isContained;
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
    float Chi2Proton;
    if(trk.bestplane == 0){
      Chi2Proton = trk.chi2pid0.chi2_proton;
    }
    else if(trk.bestplane == 1){
      Chi2Proton = trk.chi2pid1.chi2_proton;
    }
    else{
      Chi2Proton = trk.chi2pid2.chi2_proton;
    }

    return Chi2Proton;

  }
  else{
    return 999999;
  }
});

const Var varTruthProtonMatchedTrackNormalizedChi2Proton([](const caf::SRSliceProxy* slc) -> double {
  int protonTrackIndex = varTruthProtonMatchedTrackIndex(slc);
  if(protonTrackIndex>=0){

    auto const& trk = slc->reco.trk.at(protonTrackIndex);
    float Chi2Proton;
    int Chi2Ndof;
    if(trk.bestplane == 0){
      Chi2Proton = trk.chi2pid0.chi2_proton;
      Chi2Ndof = trk.chi2pid0.pid_ndof;
    }
    else if(trk.bestplane == 1){
      Chi2Proton = trk.chi2pid1.chi2_proton;
      Chi2Ndof = trk.chi2pid1.pid_ndof;
    }
    else{
      Chi2Proton = trk.chi2pid2.chi2_proton;
      Chi2Ndof = trk.chi2pid2.pid_ndof;
    }

    return Chi2Proton/Chi2Ndof;

  }
  else{
    return 999999;
  }
});

const Var varTruthProtonMatchedTrackChi2Muon([](const caf::SRSliceProxy* slc) -> double {
  int protonTrackIndex = varTruthProtonMatchedTrackIndex(slc);
  if(protonTrackIndex>=0){
    
    auto const& trk = slc->reco.trk.at(protonTrackIndex);
    float Chi2Muon;
    if(trk.bestplane == 0){
      Chi2Muon = trk.chi2pid0.chi2_muon;
    }
    else if(trk.bestplane == 1){
      Chi2Muon = trk.chi2pid1.chi2_muon;
    }
    else{
      Chi2Muon = trk.chi2pid2.chi2_muon;
    }

    return Chi2Muon;

  }
  else{
    return 999999;
  }

});

const Var varTruthProtonMatchedTrackNormalizedChi2Muon([](const caf::SRSliceProxy* slc) -> double {
  int protonTrackIndex = varTruthProtonMatchedTrackIndex(slc);
  if(protonTrackIndex>=0){

    auto const& trk = slc->reco.trk.at(protonTrackIndex);
    float Chi2Muon;
    int Chi2Ndof;
    if(trk.bestplane == 0){
      Chi2Muon = trk.chi2pid0.chi2_muon;
      Chi2Ndof = trk.chi2pid0.pid_ndof;
    }
    else if(trk.bestplane == 1){
      Chi2Muon = trk.chi2pid1.chi2_muon;
      Chi2Ndof = trk.chi2pid1.pid_ndof;
    }
    else{
      Chi2Muon = trk.chi2pid2.chi2_muon;
      Chi2Ndof = trk.chi2pid2.pid_ndof;
    }

    return Chi2Muon/Chi2Ndof;

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
    return trk_Proton.mcsP.fwdP_proton;
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

#endif
