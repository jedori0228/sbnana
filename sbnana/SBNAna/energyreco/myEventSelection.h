#ifndef myEventSelection_h
#define myEventSelection_h

#include "sbnana/SBNAna/Cuts/NumuCutsIcarus202106.h"
#include "sbnana/SBNAna/Cuts/Cuts.h"
#include "sbnana/SBNAna/energyreco/myMuonSelection.h"
#include "sbnana/SBNAna/energyreco/myProtonSelection.h"
#include "sbnana/SBNAna/energyreco/myTruth.h"

//==== Truth fiducial
const Cut cutmyTFiducial([](const caf::SRSliceProxy* slc) {

  double XMargin = 25.;
  double YMargin = 25.;
  double ZMarginUp = 30.;
  double ZMarginDown = 50.;

  return ( !isnan(slc->truth.position.x) &&
         ( ( slc->truth.position.x < -71.1 - XMargin && slc->truth.position.x > -369.33 + XMargin ) ||
         ( slc->truth.position.x > 71.1 + XMargin && slc->truth.position.x < 369.33 - XMargin ) ) &&
         !isnan(slc->truth.position.y) &&
         ( slc->truth.position.y > -181.7 + YMargin && slc->truth.position.y < 134.8 - YMargin ) &&
         !isnan(slc->truth.position.z) &&
         ( slc->truth.position.z > -895.95 + ZMarginUp && slc->truth.position.z < 895.95 - ZMarginDown ) );
});
//==== Reco fiducial
const Cut cutmyRFiducial([](const caf::SRSliceProxy* slc) {

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


//==== fmatch
const Cut cutmyFMScore([](const caf::SRSliceProxy* slc) {
    return ( !isnan(slc->fmatch.score) && slc->fmatch.score < 10.0 );
  });

const Cut cutHasMuonTrack([](const caf::SRSliceProxy *slc){
  return ( varMuonTrackIndex(slc) >= 0);
});

const Cut cutHasProtonTrack([](const caf::SRSliceProxy *slc){
  return ( varProtonTrackIndex(slc) >= 0);
});

const Cut cutIsMuonTrackLong([](const caf::SRSliceProxy *slc){
  int muonTrackIndex = varMuonTrackIndex(slc);
  if(muonTrackIndex>=0){
    auto const& trk_Muon = slc->reco.trk.at(muonTrackIndex);

    if( trk_Muon.len>50. ){
      return true;
    }
    else{
      return false;
    }

  }
  else{
    return false;
  }

});

//==== Truth information

const Cut cutTruthIsCC([](const caf::SRSliceProxy *slc){
  return (!isnan(slc->truth.E) && slc->truth.iscc);
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

#endif
