#pragma once

#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_Variables.h"

namespace ICARUSNumuXsec{

  //==== NC
  const Cut cutIsNuMuNC([](const caf::SRSliceProxy* slc) {
      return ( kIsNuSlice(slc) && slc->truth.isnc && ( slc->truth.pdg == 14 || slc->truth.pdg == -14 ) );
    });

  const Cut cutTruthNoPiZero([](const caf::SRSliceProxy *slc){
    return (varTruthNPiZero(slc)==0);
  });

  const Cut cutTruthNoChargedPion([](const caf::SRSliceProxy *slc){
    return (varTruthNChargedPion(slc)==0);
  });

  const Cut cutTruthNoNeutron([](const caf::SRSliceProxy *slc){
    return (varTruthNNeutron(slc)==0);
  });

  //==== FV

  const Cut cutRFiducial([](const caf::SRSliceProxy* slc) {

    double XMargin = 25.;
    double YMargin = 25.;
    double ZMarginUp = 30.;
    double ZMarginDown = 50.;

    return ( !isnan(slc->vertex.x) &&
           ( ( slc->vertex.x < -71.1 - XMargin && slc->vertex.x > -369.33 + XMargin ) ||
           ( slc->vertex.x > 71.1 + XMargin && slc->vertex.x < 369.33 - XMargin ) ) &&
           !isnan(slc->vertex.y) &&
           ( slc->vertex.y > -181.7 + YMargin && slc->vertex.y < 134.8 - YMargin ) &&
           !isnan(slc->vertex.z) &&
           ( slc->vertex.z > -895.95 + ZMarginUp && slc->vertex.z < 895.95 - ZMarginDown ) );
  });

  //==== FMScore
  const Cut cutFMScore([](const caf::SRSliceProxy* slc) {
      return ( !isnan(slc->fmatch.score) && slc->fmatch.score < 6.0 );
    });

  const Cut cutHasMuon([](const caf::SRSliceProxy* slc) {
    return (varMuonTrackInd(slc) >= 0);
  });

  const Cut cutMuonContained([](const caf::SRSliceProxy* slc) {

    bool Contained(false);
    if( varMuonTrackInd(slc) >= 0 ){
      auto const& trk = slc->reco.trk.at(varMuonTrackInd(slc));
      Contained = ( !isnan(trk.end.x) &&
                  ( trk.end.x < -71.1 - 25 && trk.end.x > -369.33 + 25 ) &&
                  !isnan(trk.end.y) &&
                  ( trk.end.y > -181.7 + 25 && trk.end.y < 134.8 - 25 ) &&
                  !isnan(trk.end.z) &&
                  ( trk.end.z > -895.95 + 30 && trk.end.z < 895.95 - 50 ) );
      return Contained;
    }
    return Contained;

  });


  const Cut cutMuonMatchedToMuon([](const caf::SRSliceProxy* slc) {
    return ( abs(varMuonTruePDG(slc)) == 13 );
  });

  const Cut cutMuonMatchedToProton([](const caf::SRSliceProxy* slc) {
    return ( abs(varMuonTruePDG(slc)) == 2212 );
  });


  const Cut cutHasProton([](const caf::SRSliceProxy* slc) {
    return (varProtonTrackInd(slc) >= 0);
  });

  const Cut cutProtonContained([](const caf::SRSliceProxy* slc) {

    bool Contained(false);
    if( varProtonTrackInd(slc) >= 0 ){
      auto const& trk = slc->reco.trk.at(varProtonTrackInd(slc));
      Contained = ( !isnan(trk.end.x) &&
                  ( trk.end.x < -71.1 - 25 && trk.end.x > -369.33 + 25 ) &&
                  !isnan(trk.end.y) &&
                  ( trk.end.y > -181.7 + 25 && trk.end.y < 134.8 - 25 ) &&
                  !isnan(trk.end.z) &&
                  ( trk.end.z > -895.95 + 30 && trk.end.z < 895.95 - 50 ) );
      return Contained;
    }
    return Contained;

  });


  const Cut cutProtonMatchedToMuon([](const caf::SRSliceProxy* slc) {
    return ( abs(varProtonTruePDG(slc)) == 13 );
  });

  const Cut cutProtonMatchedToProton([](const caf::SRSliceProxy* slc) {
    return ( abs(varProtonTruePDG(slc)) == 2212 );
  });


}

