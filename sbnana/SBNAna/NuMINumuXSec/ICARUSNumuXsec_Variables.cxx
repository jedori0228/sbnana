#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_Variables.h"

using namespace std;
using namespace ana;

namespace ICARUSNumuXsec{

  // SpillVar
  const SpillVar spillvarCountSpill([](const caf::SRSpillProxy *sr) ->int {
    return 0.;
  });
  // - Test
  const SpillVar spillvarTest([](const caf::SRSpillProxy *sr) ->int {
    double ret = 0.;
    return ret;
  });
  // - PMT-CRT matching
  const SpillMultiVar spillvarOpFlashTime([](const caf::SRSpillProxy *sr)
  {
    std::vector<double> rets;
    for(const auto& opflash : sr->opflashes){
      rets.push_back(opflash.firsttime);
    }
    return rets;
  });
  const SpillMultiVar spillvarValidOpFlashTime([](const caf::SRSpillProxy *sr)
  {
    std::vector<double> rets;
    for(const auto& opflash : sr->opflashes){
      if(opflash.onbeamtime) rets.push_back(opflash.firsttime); // TODO I'm using a hacked version of onbeamtime..
    }
    return rets;
  });
  const SpillMultiVar spillvarInTimeOpFlashTime([](const caf::SRSpillProxy *sr)
  {
    vector<double> validTimes = spillvarValidOpFlashTime(sr);
    std::vector<double> rets;
    for(const auto& opt : validTimes){
      if( cpmt.IsInTime(opt) ) rets.push_back(opt);
    }
    return rets;
  });
  const SpillMultiVar spillvarCRTPMTTime([](const caf::SRSpillProxy *sr)
  {
    vector<double> intimeTimes = spillvarInTimeOpFlashTime(sr);
    std::vector<double> rets;
    for(const auto& opt : intimeTimes){
      int crtHitIdx = cpmt.GetMatchedCRTHitIndex(opt, sr->crt_hits, 0);
      if(crtHitIdx>=0){
        rets.push_back( sr->crt_hits.at(crtHitIdx).t1 - opt );
      }
    }
    return rets;
  });
  const SpillMultiVar spillvarCRTPMTMatchingID([](const caf::SRSpillProxy *sr)
  {
    vector<double> intimeTimes = spillvarInTimeOpFlashTime(sr);
    std::vector<double> rets;
    if(intimeTimes.size()>0){
      for(const auto& opt : intimeTimes){
        rets.push_back( cpmt.GetMatchID(opt, sr->crt_hits) );
      }
    }
    return rets;
  });

  const SpillVar spillvarCRTPMTMatchingEventID([](const caf::SRSpillProxy *sr)
  {
/*
    vector<double> intimeTimes = spillvarInTimeOpFlashTime(sr);
    if(intimeTimes.size()>0){
      int CountNoMatching = 0;
      int CountEnteringTrack = 0;
      int CountExitingTrack = 0;
      int CountUnknown = 0
      for(const auto& opt : intimeTimes){
        int this_ID = cpmt.GetMatchID(opt, sr->crt_hits);
        if(this_ID==0) CountNoMatching++;
        if(this_ID==1 || this_ID==2 || this_ID==3 || this_ID==6 || this_ID==7) CountEnteringTrack++;
        if(this_ID==4 || this_ID==5) CountExitingTrack++;
        if(this_ID==8) CountUnknown++;
      }
    }
    else{
      return 9;
    }
*/
    return 0;
  });

  // Var
  // - Flash matching
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
  // - Longest track
  //   - index
  const Var varLongestTrackIndex([](const caf::SRSliceProxy* slc) -> int {
    int ret(-1);
    double lmax(-999.);

    for(std::size_t i(0); i < slc->reco.trk.size(); ++i){

      const auto& trk = slc->reco.trk.at(i);

      bool pass = false;
      if( fv_track.isContained(trk.end.x, trk.end.y, trk.end.z) ){
        pass = trk.len>50.;
      }
      else{
        pass = trk.len>100.;
      }
      if(!pass) continue;

      if(trk.len>lmax){
        lmax = trk.len;
        ret = i;
      }
    }

    return ret;
  });
  //   - length
  const Var varLongestTrackLength([](const caf::SRSliceProxy* slc) -> double {
    int ltidx = varLongestTrackIndex(slc);
    if(ltidx>=0){
      const auto& trk = slc->reco.trk.at(ltidx);
      return trk.len;
    }
    else{
      return -999.;
    }
  });
  //   - direction
  const Var varLongestTrackDirectionX([](const caf::SRSliceProxy* slc) -> double {
    int ltidx = varLongestTrackIndex(slc);
    if(ltidx>=0){
      const auto& trk = slc->reco.trk.at(ltidx);
      return trk.dir.x;
    }
    else{
      return -999.;
    }
  });
  const Var varLongestTrackDirectionY([](const caf::SRSliceProxy* slc) -> double {
    int ltidx = varLongestTrackIndex(slc);
    if(ltidx>=0){
      const auto& trk = slc->reco.trk.at(ltidx);
      return trk.dir.y;
    }
    else{
      return -999.;
    }
  });
  const Var varLongestTrackDirectionZ([](const caf::SRSliceProxy* slc) -> double {
    int ltidx = varLongestTrackIndex(slc);
    if(ltidx>=0){
      const auto& trk = slc->reco.trk.at(ltidx);
      return trk.dir.z;
    }
    else{
      return -999.;
    }
  });
  const Var varLongestTrackDirectionXZ([](const caf::SRSliceProxy* slc) -> double {
    int ltidx = varLongestTrackIndex(slc);
    if(ltidx>=0){
      const auto& trk = slc->reco.trk.at(ltidx);
      return sqrt(trk.dir.x*trk.dir.x+trk.dir.z*trk.dir.z);
    }
    else{
      return -999.;
    }
  });
  const Var varLongestTrackForceDownDirectionX([](const caf::SRSliceProxy* slc) -> double {
    int ltidx = varLongestTrackIndex(slc);
    if(ltidx>=0){
      const auto& trk = slc->reco.trk.at(ltidx);
      double flip = trk.dir.y>0 ? -1. : +1.; // if original track is upward(>0), flip it
      return trk.dir.x*flip;
    }
    else{
      return -999.;
    }
  });
  const Var varLongestTrackForceDownDirectionY([](const caf::SRSliceProxy* slc) -> double {
    int ltidx = varLongestTrackIndex(slc);
    if(ltidx>=0){
      const auto& trk = slc->reco.trk.at(ltidx);
      double flip = trk.dir.y>0 ? -1. : +1.; // if original track is upward(>0), flip it
      return trk.dir.y*flip;
    }
    else{
      return -999.;
    }
  });
  const Var varLongestTrackForceDownDirectionZ([](const caf::SRSliceProxy* slc) -> double {
    int ltidx = varLongestTrackIndex(slc);
    if(ltidx>=0){
      const auto& trk = slc->reco.trk.at(ltidx);
      double flip = trk.dir.y>0 ? -1. : +1.; // if original track is upward(>0), flip it
      return trk.dir.z*flip;
    }
    else{
      return -999.;
    }
  });
  // - Reco muon
  //   - index
  const Var varMuonTrackInd([](const caf::SRSliceProxy* slc) -> int {

    // The (dis)qualification of a slice is based upon the track level information.
    float Atslc, Longest(0);
    float Chi2Proton, Chi2Muon;
    bool AtSlice, Contained, MaybeMuonExiting, MaybeMuonContained;
    int PTrackInd(-1);

    for(std::size_t i(0); i < slc->reco.trk.size(); ++i){

      auto const& trk = slc->reco.trk.at(i);
      // First we calculate the distance of each track to the slice vertex.
      Atslc = sqrt( pow( slc->vertex.x - trk.start.x, 2.0 ) +
              pow( slc->vertex.y - trk.start.y, 2.0 ) +
              pow( slc->vertex.z - trk.start.z, 2.0 ) );
      // We require that the distance of the track from the slice is less than
      // 10 cm and that the parent of the track has been marked as the primary.
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
  //   - length
  const Var varMuonLength([](const caf::SRSliceProxy* slc) -> float {
    if( varMuonTrackInd(slc) >= 0 ){
      auto const& trk = slc->reco.trk.at(varMuonTrackInd(slc));
      return trk.len;
    }
    else{
      return -999.;
    }
  });
  //   - direction
  const Var varMuonRecoDirectionX([](const caf::SRSliceProxy* slc) -> float {
    if( varMuonTrackInd(slc) >= 0 ){
      auto const& trk = slc->reco.trk.at(varMuonTrackInd(slc));
      return trk.dir.x;
    }
    return -9999999999.;
  });
  const Var varMuonRecoDirectionY([](const caf::SRSliceProxy* slc) -> float {
    if( varMuonTrackInd(slc) >= 0 ){
      auto const& trk = slc->reco.trk.at(varMuonTrackInd(slc));
      return trk.dir.y;
    }
    return -9999999999.;
  });
  const Var varMuonRecoDirectionZ([](const caf::SRSliceProxy* slc) -> float {
    if( varMuonTrackInd(slc) >= 0 ){
      auto const& trk = slc->reco.trk.at(varMuonTrackInd(slc));
      return trk.dir.z;
    }
    return -9999999999.;
  });
  const Var varMuonRecoForceDownDirectionX([](const caf::SRSliceProxy* slc) -> float {
    if( varMuonTrackInd(slc) >= 0 ){
      auto const& trk = slc->reco.trk.at(varMuonTrackInd(slc));
      double flip = trk.dir.y>0 ? -1. : +1.; // if original track is upward(>0), flip it
      return trk.dir.x*flip;
    }
    return -9999999999.;
  });
  const Var varMuonRecoForceDownDirectionY([](const caf::SRSliceProxy* slc) -> float {
    if( varMuonTrackInd(slc) >= 0 ){
      auto const& trk = slc->reco.trk.at(varMuonTrackInd(slc));
      double flip = trk.dir.y>0 ? -1. : +1.; // if original track is upward(>0), flip it
      return trk.dir.y*flip;
    }
    return -9999999999.;
  });
  const Var varMuonRecoForceDownDirectionZ([](const caf::SRSliceProxy* slc) -> float {
    if( varMuonTrackInd(slc) >= 0 ){
      auto const& trk = slc->reco.trk.at(varMuonTrackInd(slc));
      double flip = trk.dir.y>0 ? -1. : +1.; // if original track is upward(>0), flip it
      return trk.dir.z*flip;
    }
    return -9999999999.;
  });
  //   - angle
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
  //   - pid
  extern const Var varMuonChi2Muon;
  extern const Var varMuonChi2Proton;
  extern const Var varMuonChi2Pion;
  //   - momentum
  extern const Var varMuonRecoP;


}
