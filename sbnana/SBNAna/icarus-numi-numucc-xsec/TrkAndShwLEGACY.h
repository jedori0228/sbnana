// BH: This keeps shower-like definition at 0.45 (right now is only used in the photon cut...), BUT uses 0.5 as trackScore for the tracklike cut used for proton and muon candidate selections

#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

using namespace ana;

bool IsValidTrkIdx( const caf::SRSliceProxy* slice, const unsigned int idxTrk ) {
  return slice->reco.npfp > idxTrk;
}

bool IsTracklikeTrack( const caf::SRSliceProxy* slice, const unsigned int idxTrk ) {
  return (!std::isnan(slice->reco.pfp.at(idxTrk).trackScore) && slice->reco.pfp.at(idxTrk).trackScore > 0.5);
}

bool IsShowerlike( const caf::SRSliceProxy* slice, const unsigned int idxShw ) {
  return (!std::isnan(slice->reco.pfp.at(idxShw).trackScore) && slice->reco.pfp.at(idxShw).trackScore > 0. && slice->reco.pfp.at(idxShw).trackScore <= 0.45 );
}

bool IsPrimaryPFP( const caf::SRSliceProxy* slice, const unsigned int idxTrk ) {
  return slice->reco.pfp.at(idxTrk).parent_is_primary;
}
