// BH: note that for simplicity with our current p, mu candidate selection and shower rejection, we use 0.45 as the separator in the selection.
// To make it like it was before, use 0.5 for the track shower separation instead!!!

#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

using namespace ana;

bool IsValidTrkIdx( const caf::SRSliceProxy* slice, const unsigned int idxTrk ) {
  return slice->reco.npfp > idxTrk;
}

bool IsTracklikeTrack( const caf::SRSliceProxy* slice, const unsigned int idxTrk ) {
  return (!std::isnan(slice->reco.pfp.at(idxTrk).trackScore) && slice->reco.pfp.at(idxTrk).trackScore > 0.45);
}

bool IsShowerlike( const caf::SRSliceProxy* slice, const unsigned int idxShw ) {
  return (!std::isnan(slice->reco.pfp.at(idxShw).trackScore) && slice->reco.pfp.at(idxShw).trackScore > 0. && slice->reco.pfp.at(idxShw).trackScore <= 0.45 );
}

bool IsPrimaryPFP( const caf::SRSliceProxy* slice, const unsigned int idxTrk ) {
  return slice->reco.pfp.at(idxTrk).parent_is_primary;
}
