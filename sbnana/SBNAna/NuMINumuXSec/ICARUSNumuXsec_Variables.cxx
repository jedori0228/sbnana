#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_Variables.h"

using namespace std;
using namespace ana;

namespace ICARUSNumuXsec{


  // - PMT
  const SpillMultiVar OpFlashFirstTime([](const caf::SRSpillProxy *sr) -> vector<double> {

    std::vector<double> rets;

    double manual_shift = 0.;
    if(!sr->hdr.ismc) manual_shift -= 4.;

    for(const auto& opflash : sr->opflashes){
      rets.push_back( opflash.firsttime + manual_shift );
    }

    return rets;

  });
  const SpillMultiVar OpFlashTime([](const caf::SRSpillProxy *sr) -> vector<double> {

    std::vector<double> rets;

    double manual_shift = 0.;
    if(!sr->hdr.ismc) manual_shift -= 4.;

    for(const auto& opflash : sr->opflashes){
      rets.push_back( opflash.time + manual_shift );
    }
    return rets;

  });

  const SpillMultiVar OpFlashTimeAfterSignalSelection([](const caf::SRSpillProxy *sr) -> vector<double> {

    std::vector<double> rets;

    // First check if there is a slice that pass signal selection
    bool HasSlicePassSelection = false;
    for(std::size_t i(0); i < sr->slc.size(); ++i){
      const auto& slc = sr->slc.at(i);
      if( kNuMISelection_1muNp0pi(&slc) ){
        HasSlicePassSelection = true;
        break;
      }
    }

    if(HasSlicePassSelection){

      double manual_shift = 0.;
      if(!sr->hdr.ismc) manual_shift -= 4.;

      for(const auto& opflash : sr->opflashes){
        rets.push_back( opflash.time + manual_shift );
      }
    }
    return rets;

  });

  // - Recalc chi2

  const Var kNuMIRecoMuonChi2MuonPlusMichel5cmShift([](const caf::SRSliceProxy* slc) -> double {
    if( kNuMIMuonCandidateIdx(slc) >= 0 ){
      auto const& trk = slc->reco.pfp.at(kNuMIMuonCandidateIdx(slc)).trk;
      const bool Contained = isContainedVol(trk.end.x,trk.end.y,trk.end.z);
      if(Contained){
        return dedxtempt.CalculateChi2MuonPlusMichel(trk.calo[2], 5., true);
      }
      else return -5;
    }
    else{
      return -1;
    }
  });
  const Var kNuMIRecoMuonChi2MuonPlusMichel5cmNoShift([](const caf::SRSliceProxy* slc) -> double {
    if( kNuMIMuonCandidateIdx(slc) >= 0 ){
      auto const& trk = slc->reco.pfp.at(kNuMIMuonCandidateIdx(slc)).trk;
      const bool Contained = isContainedVol(trk.end.x,trk.end.y,trk.end.z);
      if(Contained){
        return dedxtempt.CalculateChi2MuonPlusMichel(trk.calo[2], 5., false);
      }
      else return -5;
    }
    else{
      return -1;
    }
  });

  const Var kNuMIRecoMuonChi2MuonPlusMichel10cmShift([](const caf::SRSliceProxy* slc) -> double {
    if( kNuMIMuonCandidateIdx(slc) >= 0 ){
      auto const& trk = slc->reco.pfp.at(kNuMIMuonCandidateIdx(slc)).trk;
      const bool Contained = isContainedVol(trk.end.x,trk.end.y,trk.end.z);
      if(Contained){
        return dedxtempt.CalculateChi2MuonPlusMichel(trk.calo[2], 10., true);
      }
      else return -5;
    }
    else{
      return -1;
    }
  });
  const Var kNuMIRecoMuonChi2MuonPlusMichel10cmNoShift([](const caf::SRSliceProxy* slc) -> double {
    if( kNuMIMuonCandidateIdx(slc) >= 0 ){
      auto const& trk = slc->reco.pfp.at(kNuMIMuonCandidateIdx(slc)).trk;
      const bool Contained = isContainedVol(trk.end.x,trk.end.y,trk.end.z);
      if(Contained){
        return dedxtempt.CalculateChi2MuonPlusMichel(trk.calo[2], 10., false);
      }
      else return -5;
    }
    else{
      return -1;
    }
  });

  const Var kNuMIRecoMuonFloatChi2MuonPlusMichel1cm([](const caf::SRSliceProxy* slc) -> double {
    if( kNuMIMuonCandidateIdx(slc) >= 0 ){
      auto const& trk = slc->reco.pfp.at(kNuMIMuonCandidateIdx(slc)).trk;
      const bool Contained = isContainedVol(trk.end.x,trk.end.y,trk.end.z);
      if(Contained){
        return dedxtempt.CalculateFloatingChi2MuonPlusMichel(trk.calo[2], 1.);
      }
      else return -5;
    }
    else{
      return -1;
    }
  });
  const Var kNuMIRecoMuonFloatChi2MuonPlusMichel2cm([](const caf::SRSliceProxy* slc) -> double {
    if( kNuMIMuonCandidateIdx(slc) >= 0 ){
      auto const& trk = slc->reco.pfp.at(kNuMIMuonCandidateIdx(slc)).trk;
      const bool Contained = isContainedVol(trk.end.x,trk.end.y,trk.end.z);
      if(Contained){
        return dedxtempt.CalculateFloatingChi2MuonPlusMichel(trk.calo[2], 2.);
      }
      else return -5;
    }
    else{
      return -1;
    }
  });
  const Var kNuMIRecoMuonFloatChi2MuonPlusMichel3cm([](const caf::SRSliceProxy* slc) -> double {
    if( kNuMIMuonCandidateIdx(slc) >= 0 ){
      auto const& trk = slc->reco.pfp.at(kNuMIMuonCandidateIdx(slc)).trk;
      const bool Contained = isContainedVol(trk.end.x,trk.end.y,trk.end.z);
      if(Contained){
        return dedxtempt.CalculateFloatingChi2MuonPlusMichel(trk.calo[2], 3.);
      }
      else return -5;
    }
    else{
      return -1;
    }
  });
  const Var kNuMIRecoMuonFloatChi2MuonPlusMichel4cm([](const caf::SRSliceProxy* slc) -> double {
    if( kNuMIMuonCandidateIdx(slc) >= 0 ){
      auto const& trk = slc->reco.pfp.at(kNuMIMuonCandidateIdx(slc)).trk;
      const bool Contained = isContainedVol(trk.end.x,trk.end.y,trk.end.z);
      if(Contained){
        return dedxtempt.CalculateFloatingChi2MuonPlusMichel(trk.calo[2], 4.);
      }
      else return -5;
    }
    else{
      return -1;
    }
  });
  const Var kNuMIRecoMuonFloatChi2MuonPlusMichel5cm([](const caf::SRSliceProxy* slc) -> double {
    if( kNuMIMuonCandidateIdx(slc) >= 0 ){
      auto const& trk = slc->reco.pfp.at(kNuMIMuonCandidateIdx(slc)).trk;
      const bool Contained = isContainedVol(trk.end.x,trk.end.y,trk.end.z);
      if(Contained){
        return dedxtempt.CalculateFloatingChi2MuonPlusMichel(trk.calo[2], 5.);
      }
      else return -5;
    }
    else{
      return -1;
    }
  });
  const Var kNuMIRecoMuonFloatChi2MuonPlusMichel10cm([](const caf::SRSliceProxy* slc) -> double {
    if( kNuMIMuonCandidateIdx(slc) >= 0 ){
      auto const& trk = slc->reco.pfp.at(kNuMIMuonCandidateIdx(slc)).trk;
      const bool Contained = isContainedVol(trk.end.x,trk.end.y,trk.end.z);
      if(Contained){
        //return dedxtempt.CalculateFloatingChi2MuonPlusMichel(trk.calo[2], 10.);


        double new_chi2 = dedxtempt.CalculateFloatingChi2MuonPlusMichel(trk.calo[2], 10.);
        double old_chi2 = trk.chi2pid[2].chi2_muon;
        printf("[kNuMIRecoMuonFloatChi2MuonPlusMichel10cm] chi2 (original, fitted) = (%1.3f, %1.3f)\n", old_chi2, new_chi2);
        return new_chi2;

      }
      else return -5;
    }
    else{
      return -1;
    }
  });
  const Var kNuMIRecoMuonFloatChi2MuonPlusMichel15cm([](const caf::SRSliceProxy* slc) -> double {
    if( kNuMIMuonCandidateIdx(slc) >= 0 ){
      auto const& trk = slc->reco.pfp.at(kNuMIMuonCandidateIdx(slc)).trk;
      const bool Contained = isContainedVol(trk.end.x,trk.end.y,trk.end.z);
      if(Contained){
        return dedxtempt.CalculateFloatingChi2MuonPlusMichel(trk.calo[2], 15.);
      }
      else return -5;
    }
    else{
      return -1;
    }
  });
  const Var kNuMIRecoMuonFloatChi2MuonPlusMichelDebug([](const caf::SRSliceProxy* slc) -> double {
    if( kNuMIMuonCandidateIdx(slc) >= 0 ){
      auto const& trk = slc->reco.pfp.at(kNuMIMuonCandidateIdx(slc)).trk;
      const bool Contained = isContainedVol(trk.end.x,trk.end.y,trk.end.z);
      if(Contained){
        return dedxtempt.CalculateFloatingChi2MuonPlusMichel(trk.calo[2], -1.);
      }
      else return -5;
    }
    else{
      return -1;
    }
  });
  const Var kNuMIRecoMuonFloatChi2MuonPlusMichel15cmDelta([](const caf::SRSliceProxy* slc) -> double {
    if( kNuMIMuonCandidateIdx(slc) >= 0 ){
      auto const& trk = slc->reco.pfp.at(kNuMIMuonCandidateIdx(slc)).trk;
      const bool Contained = isContainedVol(trk.end.x,trk.end.y,trk.end.z);
      if(Contained){
        double chi2muon_original = trk.chi2pid[2].chi2_muon;
        double chi2muon_fited = dedxtempt.CalculateFloatingChi2MuonPlusMichel(trk.calo[2], 15.);
        return chi2muon_original-chi2muon_fited;
      }
      else return -5;
    }
    else{
      return -1;
    }
  });
  const Var kNuMIRecoMuonFloatChi2MuonPlusMichelDebugDelta([](const caf::SRSliceProxy* slc) -> double {
    if( kNuMIMuonCandidateIdx(slc) >= 0 ){
      auto const& trk = slc->reco.pfp.at(kNuMIMuonCandidateIdx(slc)).trk;
      const bool Contained = isContainedVol(trk.end.x,trk.end.y,trk.end.z);
      if(Contained){
        double chi2muon_original = trk.chi2pid[2].chi2_muon;
        double chi2muon_fited = dedxtempt.CalculateFloatingChi2MuonPlusMichel(trk.calo[2], -1.); 
        return chi2muon_original-chi2muon_fited;
      }
      else return -5;
    }
    else{
      return -1;
    }
  });

  const Var kNuMILeadingChargedPionCandidateChi2MuonRecalc0p5([](const caf::SRSliceProxy* slc) -> double {
    int chargedpion_index = kNuMILeadingChargedPionCandidateInd(slc);
    double ret = -5.f;
    if(chargedpion_index>=0){
      auto const& trk = slc->reco.pfp.at(chargedpion_index).trk;
      ret = dedxtempt.CalculateChi2(trk.calo[2], 3, 0.5);
    }
    return ret;
  });
  const Var kNuMILeadingChargedPionCandidateChi2MuonRecalc1p0([](const caf::SRSliceProxy* slc) -> double {
    int chargedpion_index = kNuMILeadingChargedPionCandidateInd(slc);
    double ret = -5.f;
    if(chargedpion_index>=0){
      auto const& trk = slc->reco.pfp.at(chargedpion_index).trk;
      ret = dedxtempt.CalculateChi2(trk.calo[2], 3, 1.0);
    }
    return ret;
  });



  // - Test
  const SpillMultiVar spillvarTest([](const caf::SRSpillProxy *sr) -> vector<double> {

    std::vector<double> rets;

    for(std::size_t i(0); i < sr->slc.size(); ++i){
      const auto& slc = sr->slc.at(i);

      double MuonRecoP = kNuMIMuonCandidateRecoP(&slc);
      double MuonTrueP = kNuMIMuonTrueP(&slc);
      if(MuonRecoP<0 || MuonTrueP<0) continue;

      if( !kNuMIMuonCandidateContained(&slc) ) continue;

      int MuonIdx = kNuMIMuonCandidateIdx(&slc);
      const auto& trk = slc.reco.pfp[MuonIdx].trk;

/*
      printf("\n(run, subrun, event) = (%d, %d, %d)", sr->hdr.run.GetValue(), sr->hdr.subrun.GetValue(), sr->hdr.evt.GetValue());
      printf("- Slice index = %ld\n", i);
      printf("- Muon P (reco, true) = (%1.1f, %1.1f)\n", MuonRecoP, MuonTrueP);
      printf("- Muon track start (x,y,z) = (%1.2f, %1.2f, %1.2f)\n", trk.start.x.GetValue(), trk.start.y.GetValue(), trk.start.z.GetValue());
      printf("- Muon track end (x,y,z) = (%1.2f, %1.2f, %1.2f)\n", trk.end.x.GetValue(), trk.end.y.GetValue(), trk.end.z.GetValue());
*/

      double FracDiff = (MuonRecoP-MuonTrueP)/MuonTrueP;
      if( FracDiff<-0.4 ){

        double trk_end_x = fabs(trk.end.x);
        double dCathode = fabs( trk_end_x - 210.21500 );

        if( dCathode < 5 ){

          printf("\n(run, subrun, event) = (%d, %d, %d)", sr->hdr.run.GetValue(), sr->hdr.subrun.GetValue(), sr->hdr.evt.GetValue());
          printf("- Slice index = %ld\n", i);
          printf("- Muon P (reco, true) = (%1.1f, %1.1f)\n", MuonRecoP, MuonTrueP);
          printf("- Muon track start (x,y,z) = (%1.2f, %1.2f, %1.2f)\n", trk.start.x.GetValue(), trk.start.y.GetValue(), trk.start.z.GetValue());
          printf("- Muon track end (x,y,z) = (%1.2f, %1.2f, %1.2f)\n", trk.end.x.GetValue(), trk.end.y.GetValue(), trk.end.z.GetValue());

          printf("- Neutrino time = %1.2f\n", slc.truth.time.GetValue());

        }
      }


    }

    return rets;


  });

  // SpillVar
  const SpillVar spillvarCountSpill([](const caf::SRSpillProxy *sr) -> double {
    return 0.;
  });
  // Trigger
  const SpillVar TriggerWithinGate([](const caf::SRSpillProxy *sr) -> double {
    if(sr->hdr.ismc){
      return sr->hdr.triggerinfo.trigger_within_gate;
    }
    else{
      return sr->hdr.triggerinfo.trigger_within_gate-4.;
    }
  });
  const SpillVar TriggerInfoTriggerType([](const caf::SRSpillProxy *sr) -> double {
    if(sr->hdr.ismc){
      return -1.;
    }
    else{
      return sr->hdr.triggerinfo.trigger_type;
    }
  });
  const SpillVar TriggerInfoSourceType([](const caf::SRSpillProxy *sr) -> double {
    if(sr->hdr.ismc){
      return -1.;
    }
    else{
      return sr->hdr.triggerinfo.source_type;
    }
  });
  const SpillVar spillvarNTrack([](const caf::SRSpillProxy *sr) -> double {
    int nTrk=0;
    for(std::size_t i(0); i < sr->slc.size(); ++i){
      const auto& slc = sr->slc.at(i);
      for(const auto& pfp: slc.reco.pfp){
        bool IsTrack = IsPFPTrack(pfp);
        if( IsTrack ) nTrk++;
      }
    }
    return nTrk;
  });
  const SpillVar spillvarNShower([](const caf::SRSpillProxy *sr) ->int {
    int nShw=0;
    for(std::size_t i(0); i < sr->slc.size(); ++i){
      const auto& slc = sr->slc.at(i);
      for(const auto& pfp: slc.reco.pfp){
        bool IsShower = IsPFPShower(pfp);
        if( IsShower ) nShw++;
      }
    }
    return nShw;
  });
  // - CRT Hit
  const SpillMultiVar spillvarSideCRTHitPe([](const caf::SRSpillProxy *sr){
    std::vector<double> rets;
    for(const auto& hit : sr->crt_hits){
      if(hit.plane>=40&&hit.plane<=49){
        rets.push_back(hit.pe);
      }
    }
    return rets;
  });
  const SpillMultiVar spillvarTopCRTHitPe([](const caf::SRSpillProxy *sr){
    std::vector<double> rets;
    for(const auto& hit : sr->crt_hits){
      if(hit.plane>=30&&hit.plane<=39){
        rets.push_back(hit.pe);
      }
    }
    return rets;
  });
  // - Pos
  const SpillMultiVar spillvarEastWestCRTHitPosX([](const caf::SRSpillProxy *sr){
    std::vector<double> rets;
    for(const auto& hit : sr->crt_hits){
      if(hit.plane>=40&&hit.plane<=45){
        rets.push_back(hit.position.x);
      }
    }
    return rets;
  });
  const SpillMultiVar spillvarFlashPosX([](const caf::SRSpillProxy *sr){
    std::vector<double> rets;
    for(const auto& opflash : sr->opflashes){
      rets.push_back(opflash.center.x);
    }
    return rets;
  });
  // - OpFlash
  const SpillMultiVar spillvarOpFlashPeakToFirstTime([](const caf::SRSpillProxy *sr)
  {
    std::vector<double> rets;
    for(const auto& opflash : sr->opflashes){
      rets.push_back( 1000.*(opflash.time - opflash.firsttime) );
    }
    return rets;
  });
  // - CRTHit
  const SpillMultiVar spillvarTopCRTHitTime([](const caf::SRSpillProxy *sr)
  {
    std::vector<double> rets;
    for(const auto& hit : sr->crt_hits){
      if(hit.plane>=30&&hit.plane<=39){
        double this_crttime = sr->hdr.ismc ? hit.t0 : hit.t1;
        rets.push_back(this_crttime);
      }
    }
    return rets;
  });
  const SpillMultiVar spillvarSideCRTHitTime([](const caf::SRSpillProxy *sr)
  {
    std::vector<double> rets;
    for(const auto& hit : sr->crt_hits){
      if(hit.plane>=40&&hit.plane<=49){
        double this_crttime = sr->hdr.ismc ? hit.t0 : hit.t1;
        rets.push_back(this_crttime);
      }
    }
    return rets;
  });
  // - PMT-CRT matching
  const SpillMultiVar spillvarTopCRTPMTTime([](const caf::SRSpillProxy *sr)
  {
    vector<double> intimeTimes = ICARUSCRTPMTMatching::spillvarInTimeOpFlashTime(sr);
    std::vector<double> rets;
    for(const auto& opt : intimeTimes){
      std::vector<int> crtHitIdices = ICARUSCRTPMTMatching::cpmt.GetMatchedCRTHitIndex(opt, sr->crt_hits, 0);
      for(const auto& crtHitIdx:crtHitIdices){
        double this_crttime = sr->hdr.ismc ? sr->crt_hits.at(crtHitIdx).t0 : sr->crt_hits.at(crtHitIdx).t1;
        rets.push_back( this_crttime - opt );
      }
    }
    return rets;
  });
  const SpillMultiVar spillvarSideCRTPMTTime([](const caf::SRSpillProxy *sr)
  {
    vector<double> intimeTimes = ICARUSCRTPMTMatching::spillvarInTimeOpFlashTime(sr);
    std::vector<double> rets;
    for(const auto& opt : intimeTimes){
      std::vector<int> crtHitIdices = ICARUSCRTPMTMatching::cpmt.GetMatchedCRTHitIndex(opt, sr->crt_hits, 1);
      for(const auto& crtHitIdx:crtHitIdices){
        double this_crttime = sr->hdr.ismc ? sr->crt_hits.at(crtHitIdx).t0 : sr->crt_hits.at(crtHitIdx).t1;
        rets.push_back( this_crttime - opt );
      }
    }
    return rets;
  });
  // Nu from spillvar
  const SpillVar TruthFirstNuEnergy([](const caf::SRSpillProxy *sr) -> double {
    if(sr->mc.nu.size()==0){
      return -999.;
    }
    else{
      return sr->mc.nu[0].E;
    }
  });

  // Var
  const Var varNuDirectionX([](const caf::SRSliceProxy* slc) ->int {
      double this_x = slc->truth.prod_vtx.x/100.;
      double this_y = slc->truth.prod_vtx.y/100.;
      double this_z = slc->truth.prod_vtx.z/100.;
      TVector3 this_coord = nct.GetICARUSCoord(this_x, this_y, this_z).Unit();
      this_coord *= -1.;
      return this_coord.X();
  });
  const Var varNuDirectionY([](const caf::SRSliceProxy* slc) ->int {
      double this_x = slc->truth.prod_vtx.x/100.;
      double this_y = slc->truth.prod_vtx.y/100.;
      double this_z = slc->truth.prod_vtx.z/100.;
      TVector3 this_coord = nct.GetICARUSCoord(this_x, this_y, this_z).Unit();
      this_coord *= -1.;
      return this_coord.Y();
  });
  const Var varNuDirectionZ([](const caf::SRSliceProxy* slc) ->int {
      double this_x = slc->truth.prod_vtx.x/100.;
      double this_y = slc->truth.prod_vtx.y/100.;
      double this_z = slc->truth.prod_vtx.z/100.;
      TVector3 this_coord = nct.GetICARUSCoord(this_x, this_y, this_z).Unit();
      this_coord *= -1.;
      return this_coord.Z();
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
  const Var varTruthW([](const caf::SRSliceProxy* slc) -> double {
    if(isnan(slc->truth.w)) return -999.;
    else return slc->truth.w;
  });
  // - Truth n particle
  const Var varTruthNProton([](const caf::SRSliceProxy* slc) -> double {
    ICARUSNumuXsec::InteractionTool::NParticles nptls = intt.GetNParticles(slc);
    return nptls.NProton;
  });
  const Var varTruthNNeutron([](const caf::SRSliceProxy* slc) -> double {
    ICARUSNumuXsec::InteractionTool::NParticles nptls = intt.GetNParticles(slc);
    return nptls.NNeutron;
  });
  const Var varTruthNPip([](const caf::SRSliceProxy* slc) -> double {
    ICARUSNumuXsec::InteractionTool::NParticles nptls = intt.GetNParticles(slc);
    return nptls.NPip;
  });
  const Var varTruthNPim([](const caf::SRSliceProxy* slc) -> double {
    ICARUSNumuXsec::InteractionTool::NParticles nptls = intt.GetNParticles(slc);
    return nptls.NPim;
  });
  const Var varTruthNPi0([](const caf::SRSliceProxy* slc) -> double {
    ICARUSNumuXsec::InteractionTool::NParticles nptls = intt.GetNParticles(slc);
    return nptls.NPi0;
  });

  // - Slice var
  const Var varCountSlice([](const caf::SRSliceProxy* slc) ->int {
    return 0.;
  });
  const Var varIsClearCosmic([](const caf::SRSliceProxy* slc) ->int {
    if(slc->is_clear_cosmic) return 1;
    else return 0;
  });
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
  // - NuID
  const Var varCRLongestTrackDirY([](const caf::SRSliceProxy* slc) -> double {
    return slc->nuid.crlongtrkdiry;
  });
  // - Long-enough tracks
  const MultiVar varLongTrackDirectionY([](const caf::SRSliceProxy* slc) -> std::vector<double> {
    std::vector<double> rets;
    for(std::size_t i(0); i < slc->reco.pfp.size(); ++i){

      bool IsTrack = IsPFPTrack( slc->reco.pfp.at(i) );
      if(!IsTrack) continue;
      const auto& trk = slc->reco.pfp.at(i).trk;
      if(isnan(trk.len)) continue;
      if(trk.len>50.) rets.push_back(trk.dir.y);

    }
    return rets;
  });
  // - Primary tracks
  const MultiVar PrimaryTrackIndices([](const caf::SRSliceProxy* slc) -> std::vector<double> {
    std::vector<double> rets;
    for(unsigned int i_pfp=0; i_pfp<slc->reco.pfp.size(); i_pfp++){
      const auto& pfp = slc->reco.pfp.at(i_pfp);
      if( !IsPFPTrack(pfp) ) continue; // TODO
      const auto& trk = pfp.trk;
      if(isnan(trk.start.x)) continue;
      const float Atslc = std::hypot(slc->vertex.x - trk.start.x,
                                     slc->vertex.y - trk.start.y,
                                     slc->vertex.z - trk.start.z);
      if (Atslc < 10. && pfp.parent_is_primary){
        rets.push_back(i_pfp);
      }
    }
    return rets;
  });
  const Var NPrimaryTracks([](const caf::SRSliceProxy* slc) -> double {
    return PrimaryTrackIndices(slc).size();
  });

  // - Longest track
  //   - index
  const Var LongestTrackIndex([](const caf::SRSliceProxy* slc) -> int {
    int ret(-1);
    double lmax(-999.);

    for(std::size_t i(0); i < slc->reco.pfp.size(); ++i){

      bool IsTrack = IsPFPTrack( slc->reco.pfp.at(i) );
      if(!IsTrack) continue;

      const auto& trk = slc->reco.pfp.at(i).trk;

      if(isnan(trk.len)) continue;
      if(isnan(trk.end.x)) continue;

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
  const Var LongestTrackLength([](const caf::SRSliceProxy* slc) -> double {
    int ltidx = LongestTrackIndex(slc);
    if(ltidx>=0){
      const auto& trk = slc->reco.pfp.at(ltidx).trk;
      return trk.len;
    }
    else{
      return -999.;
    }
  });
  //   - direction
  const Var LongestTrackDirectionX([](const caf::SRSliceProxy* slc) -> double {
    int ltidx = LongestTrackIndex(slc);
    if(ltidx>=0){
      const auto& trk = slc->reco.pfp.at(ltidx).trk;
      return trk.dir.x;
    }
    else{
      return -999.;
    }
  });
  const Var LongestTrackDirectionY([](const caf::SRSliceProxy* slc) -> double {
    int ltidx = LongestTrackIndex(slc);
    if(ltidx>=0){
      const auto& trk = slc->reco.pfp.at(ltidx).trk;
      return trk.dir.y;
    }
    else{
      return -999.;
    }
  });
  const Var LongestTrackDirectionZ([](const caf::SRSliceProxy* slc) -> double {
    int ltidx = LongestTrackIndex(slc);
    if(ltidx>=0){
      const auto& trk = slc->reco.pfp.at(ltidx).trk;
      return trk.dir.z;
    }
    else{
      return -999.;
    }
  });
  const Var LongestTrackDirectionXZ([](const caf::SRSliceProxy* slc) -> double {
    int ltidx = LongestTrackIndex(slc);
    if(ltidx>=0){
      const auto& trk = slc->reco.pfp.at(ltidx).trk;
      return sqrt(trk.dir.x*trk.dir.x+trk.dir.z*trk.dir.z);
    }
    else{
      return -999.;
    }
  });
  const Var LongestTrackForceDownDirectionX([](const caf::SRSliceProxy* slc) -> double {
    int ltidx = LongestTrackIndex(slc);
    if(ltidx>=0){
      const auto& trk = slc->reco.pfp.at(ltidx).trk;
      double flip = trk.dir.y>0 ? -1. : +1.; // if original track is upward(>0), flip it
      return trk.dir.x*flip;
    }
    else{
      return -999.;
    }
  });
  const Var LongestTrackForceDownDirectionY([](const caf::SRSliceProxy* slc) -> double {
    int ltidx = LongestTrackIndex(slc);
    if(ltidx>=0){
      const auto& trk = slc->reco.pfp.at(ltidx).trk;
      double flip = trk.dir.y>0 ? -1. : +1.; // if original track is upward(>0), flip it
      return trk.dir.y*flip;
    }
    else{
      return -999.;
    }
  });
  const Var LongestTrackForceDownDirectionZ([](const caf::SRSliceProxy* slc) -> double {
    int ltidx = LongestTrackIndex(slc);
    if(ltidx>=0){
      const auto& trk = slc->reco.pfp.at(ltidx).trk;
      double flip = trk.dir.y>0 ? -1. : +1.; // if original track is upward(>0), flip it
      return trk.dir.z*flip;
    }
    else{
      return -999.;
    }
  });
  // - Reco muon
  // - Stub
  const Var NStubs([](const caf::SRSliceProxy* slc) -> double {
    return slc->reco.stub.size();
  });
  const MultiVar StubCollectionCharges([](const caf::SRSliceProxy* slc) -> std::vector<double> {
    vector<double> rets;
    for(const auto& stub: slc->reco.stub){
      double ret = 0.;
      for(const auto& stub_plane: stub.planes){
        if(stub_plane.p==caf::kCollection){
          for(const auto& stub_hit: stub_plane.hits){
            ret += stub_hit.charge/1e5;
          }
        }
      }
      if(ret>0.) rets.push_back(ret);
      else rets.push_back(-999.);

    }
    return rets;
  });

  // - In case of neutrino overlapped with cosmic
  const Var NNuPFP([](const caf::SRSliceProxy* slc) -> double {
    int ret = 0;
    for(const auto& pfp: slc->reco.pfp){
      if(pfp.trk.truth.p.interaction_id!=-1 && pfp.trk.truth.p.interaction_id!=INT_MIN){
        ret++;
/*
        if(kIsCosmic(slc)){
          std::cout << "[JSKIIMDEBUG][NNuPFP] This is a cosmic slice, but I found a pfp with in_id = " << pfp.trk.truth.p.interaction_id << std::endl;
          std::cout << "[JSKIIMDEBUG][NNuPFP] - pdg = " << (pfp.trk.truth.p.pdg) << std::endl;
        }
*/
      }
    }
    return ret;
  });
  const Var NCosmicPFP([](const caf::SRSliceProxy* slc) -> double {
    int ret = 0;
    for(const auto& pfp: slc->reco.pfp){
      if(pfp.trk.truth.p.interaction_id==-1) ret++;
    }
    return ret;
  });

  // - Test
  const Var SliceTestVar([](const caf::SRSliceProxy* slc) -> double {
/*
    int chargedpion_index = kNuMILeadingChargedPionCandidateInd(slc);
    int ret = -5.; // PDG can be negative using this instead of -5
    if(chargedpion_index>=0){
      const auto& trk = slc->reco.pfp[chargedpion_index].trk;
      double x = trk.truth.p.genp.x;
      double y = trk.truth.p.genp.y;
      double z = trk.truth.p.genp.z;
      if(isnan(x) || isnan(y) || isnan(z)) return -4.;
      printf("gen p (x,y,z) = (%1.3f, %1.3f, %1.3f)\n",x,y,z);
      TVector3 vec_genp(x,y,z);
      ret = vec_genp.Unit().X();
      printf("-> ret = %1.3f\n"%(ret));
    }
    return ret;
*/
return 0.;
  });

  // - For trigger eff study
  const SpillVar NuMuSliceLongestTrackLenForTriggerEff([](const caf::SRSpillProxy *sr) -> double {

    int nNuMuCCSlice = 0;
    // Require a single NuMuCC matched slice
    for(std::size_t i(0); i < sr->slc.size(); ++i){
      const auto& slc = sr->slc.at(i);

      // NuMuCC matched slice
      bool isNuMu = kIsNuSlice(&slc) && ( slc.truth.pdg == 14 || slc.truth.pdg == -14 );
      bool isCC = kIsNuSlice(&slc) && slc.truth.iscc;
      bool isNuMuCC = isNuMu && isCC;
      if(!isNuMuCC) continue;
      nNuMuCCSlice += 1;

      if(nNuMuCCSlice==2) return -999.;
      else return LongestTrackLength(&slc);

    } // END loop slice

    return -999.;

  });
  const SpillVar InTimeCosmicSliceLongestTrackLenForTriggerEff([](const caf::SRSpillProxy *sr) -> double {

    int nIntimeCosmicMuSlice = 0;
    // Require a single NuMuCC matched slice
    for(std::size_t i(0); i < sr->slc.size(); ++i){
      const auto& slc = sr->slc.at(i);

      int ltidx = LongestTrackIndex(&slc);
      if(ltidx<0) continue;

      const auto& trk = slc.reco.pfp.at(ltidx).trk;
      bool isCosmicMatch = (trk.truth.p.interaction_id == -1);
      if(!isCosmicMatch) continue;

      bool isInTime = (0.1<trk.truth.p.genT) && (trk.truth.p.genT<9.5);
      if(!isInTime) continue;

      nIntimeCosmicMuSlice += 1;

      if(nIntimeCosmicMuSlice==2) return -999.;
      else return trk.len;

    } // END loop slice

    return -999.;

  });

  const Var Pass_VtxInFV([](const caf::SRSliceProxy* slc) -> int {
    bool Pass = kNuMIVertexInFV(slc);
    if(Pass) return 1;
    else return 0;
  });
  const Var Pass_NotClearCosmic([](const caf::SRSliceProxy* slc) -> int {
    bool Pass = kNuMINotClearCosmic(slc);
    if(Pass) return 1;
    else return 0;
  });
  const Var Pass_HasMuon([](const caf::SRSliceProxy* slc) -> int {
    bool Pass = kNuMIHasMuonCandidate(slc);
    if(Pass) return 1;
    else return 0;
  });
  const Var Pass_HasProton([](const caf::SRSliceProxy* slc) -> int {
    bool Pass = kNuMIHasProtonCandidate(slc);
    if(Pass) return 1;
    else return 0;
  });
  const Var Pass_ProtonPCut([](const caf::SRSliceProxy* slc) -> int {
    bool Pass = kNuMIProtonCandidateRecoPTreshold(slc);
    if(Pass) return 1;
    else return 0;
  });
  const Var Pass_PrimaryHadronContained([](const caf::SRSliceProxy* slc) -> int {
    bool Pass = kNuMIAllPrimaryHadronsContained(slc);
    if(Pass) return 1;
    else return 0;
  });
  const Var Pass_NoChargedPionTrack([](const caf::SRSliceProxy* slc) -> int {
    bool Pass = kNuMINoSecondPrimaryMuonlikeTracks(slc);
    if(Pass) return 1;
    else return 0;
  });
  const Var Pass_NoNeutralPionShower([](const caf::SRSliceProxy* slc) -> int {
    bool Pass = kNuMICutPhotons(slc);
    if(Pass) return 1;
    else return 0;
  });
  const Var Pass_MuonContained([](const caf::SRSliceProxy* slc) -> int {
    bool Pass = kNuMIMuonCandidateContained(slc);
    if(Pass) return 1;
    else return 0;
  });

  const TruthVar kTruth_BNBDefaultWeight([](const caf::SRTrueInteractionProxy *nu) -> double {
    // For NuMI, we have a CV correction, PPFX
    // If BNB has simialr thing, use that here
    // For now we are returning 1.
    return 1.;
  });

}
