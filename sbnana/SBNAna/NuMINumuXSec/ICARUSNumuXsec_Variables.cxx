#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_Variables.h"
#include "sbnana/SBNAna/Vars/BeamExposureVars.h"

using namespace std;
using namespace ana;

namespace ICARUSNumuXsec{


  // - PMT
  const SpillMultiVar OpFlashFirstTime([](const caf::SRSpillProxy *sr) -> vector<double> {

    std::vector<double> rets;

    double manual_shift = 0.;
    if(!sr->hdr.ismc){
      if(sr->hdr.run<9300){
        manual_shift -= 4.;
      }
    }

    for(const auto& opflash : sr->opflashes){
      rets.push_back( opflash.firsttime + manual_shift );
    }

    return rets;

  });
  const SpillMultiVar OpFlashTime([](const caf::SRSpillProxy *sr) -> vector<double> {

    std::vector<double> rets;

    double manual_shift = 0.;
    if(!sr->hdr.ismc){
      if(sr->hdr.run<9300){
        manual_shift -= 4.;
      }
    }

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
      if(!sr->hdr.ismc){
        if(sr->hdr.run<9300){
          manual_shift -= 4.;
        }
      }

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
  const SpillMultiVar spillvarTest([](const caf::SRSpillProxy *sr) -> std::vector<double> {

    std::vector<double> ret;

    int SigSelSliceIdx = kNuMI_SignalSelectionSlice_Idx(sr);
    if(SigSelSliceIdx<0) return ret;
    const auto& slc = sr->slc[SigSelSliceIdx];
    int this_cuttype = kNuMISliceSignalType(&slc);

    if(this_cuttype!=5){

      double spillTriggerTime = kNuMISpillTriggerTime(sr);

      printf("\n(run, subrun, event) = (%d, %d, %d), nSlice = %ld\n", sr->hdr.run.GetValue(), sr->hdr.subrun.GetValue(), sr->hdr.evt.GetValue(),sr->slc.size());
      printf("Trigger time = %f\n", spillTriggerTime);
      printf("- Number of true_particles = %ld\n", sr->true_particles.size());

      for(unsigned int i_tp=0; i_tp<sr->true_particles.size(); i_tp++){

        const auto& prim = sr->true_particles[i_tp];

        double this_genT = prim.genT;
        double this_time_diff = fabs( spillTriggerTime - this_genT );

        double start_x = prim.start.x;
        if( (start_x>-9998) && (this_time_diff < 0.3) ){
        //if( (this_time_diff < 1.0)  ){
          printf("  - i_tp = %d\n", i_tp);
          printf("    - pdg = %d\n", prim.pdg.GetValue());
          printf("    - genT = %f (genT - Trigger = %f)\n", this_genT, this_genT - spillTriggerTime);
          printf("    - interaction_id = %d\n", prim.interaction_id.GetValue());
          printf("    - start = (%1.2f, %1.2f, %1.2f)\n", prim.start.x.GetValue(), prim.start.y.GetValue(), prim.start.z.GetValue());
          printf("    - end = (%1.2f, %1.2f, %1.2f)\n", prim.end.x.GetValue(), prim.end.y.GetValue(), prim.end.z.GetValue());
          const float dist = std::hypot(prim.end.x - prim.start.x, prim.end.y - prim.start.y, prim.end.z - prim.start.z);
          printf("    - |end-start| = %1.2f\n", dist);
          printf("    - genp = (%1.2f, %1.2f, %1.2f)\n", prim.genp.x.GetValue(), prim.genp.y.GetValue(), prim.genp.z.GetValue());
          printf("    - startE = %1.3f\n", prim.startE.GetValue());
          printf("    - end_process = %d\n", prim.end_process.GetValue());
/*
          int parent_id = prim.parent;

          const auto& prim_parent = sr->true_particles[parent_id];

          printf("    - parent info:\n");

          printf("      - pdg = %d\n", prim_parent.pdg.GetValue());
          printf("      - genT = %f (genT - Trigger = %f)\n", prim_parent.genT.GetValue(), prim_parent.genT.GetValue() - spillTriggerTime);
          printf("      - start = (%1.2f, %1.2f, %1.2f)\n", prim_parent.start.x.GetValue(), prim_parent.start.y.GetValue(), prim_parent.start.z.GetValue());
          printf("      - end = (%1.2f, %1.2f, %1.2f)\n", prim_parent.end.x.GetValue(), prim_parent.end.y.GetValue(), prim_parent.end.z.GetValue());
          const float dist_parent = std::hypot(prim_parent.end.x - prim_parent.start.x, prim_parent.end.y - prim_parent.start.y, prim_parent.end.z - prim_parent.start.z);
          printf("      - |end-start| = %1.2f\n", dist_parent);
          printf("      - genp = (%1.2f, %1.2f, %1.2f)\n", prim_parent.genp.x.GetValue(), prim_parent.genp.y.GetValue(), prim_parent.genp.z.GetValue());
          printf("      - startE = %1.3f\n", prim_parent.startE.GetValue());
          printf("      - end_process = %d\n", prim_parent.end_process.GetValue());
*/

        }
      }


      printf("- Number of truth nu = %ld\n",sr->mc.nu.size());
      for(unsigned int i_nu=0; i_nu<sr->mc.nu.size(); i_nu++){
        const auto& nu = sr->mc.nu[i_nu];
        printf("  - i_nu = %d\n", i_nu);
        printf("    - GENIE mode = %d\n", nu.genie_mode.GetValue());
        printf("    - time = %1.3f\n",nu.time.GetValue());
        printf("    - E = %f\n", nu.E.GetValue());
        printf("    - nu pos = (%1.2f, %1.2f, %1.2f)\n", nu.position.x.GetValue(), nu.position.y.GetValue(), nu.position.z.GetValue());
        printf("    - nu.prim.size() = %ld\n", nu.prim.size());
/*
        for(std::size_t j(0); j < nu.prim.size(); ++j){
          const auto& prim = nu.prim[j];
          printf("  - %ld-th prim\n",j);
          printf("    - pdg = %d\n", prim.pdg.GetValue());
          const double this_mass = ptlt.GetMass(prim.pdg.GetValue());
          printf("    - Mass = %1.3f\n", this_mass);
          printf("    - GStatus = %d\n", prim.gstatus.GetValue());
          printf("    - Parent ID = %d\n", prim.parent.GetValue());
          printf("    - ndaughters = %ld\n", prim.daughters.size());
          printf("    - start = (%1.2f, %1.2f, %1.2f)\n", prim.start.x.GetValue(), prim.start.y.GetValue(), prim.start.z.GetValue());
          printf("    - end = (%1.2f, %1.2f, %1.2f)\n", prim.end.x.GetValue(), prim.end.y.GetValue(), prim.end.z.GetValue());
          const float dist = std::hypot(prim.end.x - prim.start.x, prim.end.y - prim.start.y, prim.end.z - prim.start.z);
          printf("    - |end-start| = %1.2f\n", dist);
          printf("    - genp = (%1.2f, %1.2f, %1.2f)\n", prim.genp.x.GetValue(), prim.genp.y.GetValue(), prim.genp.z.GetValue());
          printf("    - startE = %1.3f\n", prim.startE.GetValue());
          printf("    - startE-Mass = %1.3f\n", prim.startE.GetValue()-this_mass);
          printf("    - startE-endE = %1.3f\n", prim.startE.GetValue()-prim.endE.GetValue());
          printf("    - end_process = %d\n", prim.end_process.GetValue());

        }
*/

      }

    }
    
    return ret;


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


  const Var IsCosmicSlice([](const caf::SRSliceProxy* slc) -> int {
    if( slc->truth.index < 0 ) return 1;
    else return 0;
  });

  const TruthVar DummyTruthVar([](const caf::SRTrueInteractionProxy *nu) -> int {
    return 1;
  });

  const Var kNuMINUANCECode([](const caf::SRSliceProxy* slc) -> int {
    if( slc->truth.index < 0 ) return -1;
    else return slc->truth.genie_inttype;
  });

  const SpillVar kNuMIG3ChaseSpillWeightByClosesetNu( [](const caf::SRSpillProxy *sr) -> double {

    if(!sr->hdr.ismc) return 1.;

    // Find the triggering neutrino. If none, skip this slice.
    double spillTime = kNuMISpillTriggerTime(sr);
    if ( !(spillTime > -5.) ) return 1.;
    double theMinDeltaT = 999999.;
    unsigned int idxMinDeltaT = 0;
    for ( unsigned int idx=0; idx < sr->mc.nu.size(); ++idx ) {
      if ( fabs(sr->mc.nu[idx].time - spillTime) < theMinDeltaT ) {
        theMinDeltaT = fabs(sr->mc.nu[idx].time - spillTime);
        idxMinDeltaT = idx;
      }
    }

    double wtForTriggeringNu = kGetTruthNuMIFluxWeightG3Chase(&sr->mc.nu[idxMinDeltaT]);

    return wtForTriggeringNu;

  });

  // LifetimeVariation
  const Var kNuMI_MuonMatchedTrack_TrackScore([](const caf::SRSliceProxy* slc) -> double {

    int TrueMuonIndex = kTruth_MuonIndex(&slc->truth);
    if(TrueMuonIndex<0) return -1.;

    const auto& prim_Muon = slc->truth.prim[TrueMuonIndex];

    for(const auto& pfp: slc->reco.pfp){

      if( pfp.trk.truth.p.G4ID == prim_Muon.G4ID ){
        return pfp.trackScore;
      }

    }


    return -2.;

  });

  const Var kNuMI_MuonMatchedTrack_Chi2Muon([](const caf::SRSliceProxy* slc) -> double {

    int TrueMuonIndex = kTruth_MuonIndex(&slc->truth);
    if(TrueMuonIndex<0) return -1.;

    const auto& prim_Muon = slc->truth.prim[TrueMuonIndex];

    for(const auto& pfp: slc->reco.pfp){

      if( pfp.trk.truth.p.G4ID == prim_Muon.G4ID ){
        return pfp.trk.chi2pid[2].chi2_muon;
      }

    }


    return -2.;

  });

  const Var kNuMI_MuonMatchedTrack_Chi2Proton([](const caf::SRSliceProxy* slc) -> double {
    
    int TrueMuonIndex = kTruth_MuonIndex(&slc->truth);
    if(TrueMuonIndex<0) return -1.;

    const auto& prim_Muon = slc->truth.prim[TrueMuonIndex];

    for(const auto& pfp: slc->reco.pfp){

      if( pfp.trk.truth.p.G4ID == prim_Muon.G4ID ){
        return pfp.trk.chi2pid[2].chi2_proton;
      }

    }


    return -2.;

  });

  const Var kNuMI_ProtonMatchedTrack_TrackScore([](const caf::SRSliceProxy* slc) -> double {

    int TrueProtonIndex = kTruth_ProtonIndex(&slc->truth);
    if(TrueProtonIndex<0) return -1.;

    const auto& prim_Proton = slc->truth.prim[TrueProtonIndex];

    for(const auto& pfp: slc->reco.pfp){

      if( pfp.trk.truth.p.G4ID == prim_Proton.G4ID ){
        return pfp.trackScore;
      }

    }


    return -2.;

  });

  const Var kNuMI_ProtonMatchedTrack_Chi2Muon([](const caf::SRSliceProxy* slc) -> double {

    int TrueProtonIndex = kTruth_ProtonIndex(&slc->truth);
    if(TrueProtonIndex<0) return -1.;

    const auto& prim_Proton = slc->truth.prim[TrueProtonIndex];

    for(const auto& pfp: slc->reco.pfp){

      if( pfp.trk.truth.p.G4ID == prim_Proton.G4ID ){
        return pfp.trk.chi2pid[2].chi2_muon;
      }

    }


    return -2.;

  });

  const Var kNuMI_ProtonMatchedTrack_Chi2Proton([](const caf::SRSliceProxy* slc) -> double {
    
    int TrueProtonIndex = kTruth_ProtonIndex(&slc->truth);
    if(TrueProtonIndex<0) return -1.;

    const auto& prim_Proton = slc->truth.prim[TrueProtonIndex];

    for(const auto& pfp: slc->reco.pfp){

      if( pfp.trk.truth.p.G4ID == prim_Proton.G4ID ){
        return pfp.trk.chi2pid[2].chi2_proton;
      }

    }


    return -2.;

  });

  // Rock study
  const SpillCut kNuMI_HasSignalSelectionSlice ( [](const caf::SRSpillProxy *sr) {

    bool HasSlicePassSelection = false;
    for(std::size_t i(0); i < sr->slc.size(); ++i){
      const auto& slc = sr->slc.at(i);
      if( kNuMISelection_1muNp0pi(&slc) ){
        HasSlicePassSelection = true;
        break;
      }
    }

    double spillTriggerTime = kNuMISpillTriggerTime(sr);

    return HasSlicePassSelection && (spillTriggerTime > -0.1 && spillTriggerTime < 10.1);

  });
  const SpillVar kNuMI_SignalSelectionSlice_CutType ( [](const caf::SRSpillProxy *sr) -> int {

    int ret = 0;
    for(std::size_t i(0); i < sr->slc.size(); ++i){
      const auto& slc = sr->slc.at(i);
      if( kNuMISelection_1muNp0pi(&slc) ){
        ret = kNuMISliceSignalType(&slc);
        break;
      }
    }
    return ret;

  });

  const SpillMultiVar kNuMI_SignalSelectionSlices_Others_MuonTrackMatched_genT ( [](const caf::SRSpillProxy *sr) {

    std::vector<double> rets;
    for(std::size_t i(0); i < sr->slc.size(); ++i){
      const auto& slc = sr->slc.at(i);
      if( kNuMISelection_1muNp0pi(&slc) ){
        if( kNuMISliceSignalType(&slc)==5 ){

          const auto& MuonIdx = kNuMIMuonCandidateIdx(&slc);
          const auto& MuonTrk = slc.reco.pfp[MuonIdx].trk;

          int MuonTrk_Truth_G4ID = MuonTrk.truth.p.G4ID;

          for(unsigned int i_p=0; i_p<sr->true_particles.size(); i_p++){
            if( sr->true_particles[i_p].G4ID==MuonTrk_Truth_G4ID ){
              rets.push_back( sr->true_particles[i_p].genT );
            }
          }

        }

      }
    }

    return rets;

  });

  const SpillMultiVar kNuMI_SignalSelectionSlices_Others_MuonTrackMatched_pdg ( [](const caf::SRSpillProxy *sr) {

    std::vector<double> rets;
    for(std::size_t i(0); i < sr->slc.size(); ++i){
      const auto& slc = sr->slc.at(i);
      if( kNuMISelection_1muNp0pi(&slc) ){
        if( kNuMISliceSignalType(&slc)==5 ){

          const auto& MuonIdx = kNuMIMuonCandidateIdx(&slc);
          const auto& MuonTrk = slc.reco.pfp[MuonIdx].trk;

          int MuonTrk_Truth_G4ID = MuonTrk.truth.p.G4ID;

          for(unsigned int i_p=0; i_p<sr->true_particles.size(); i_p++){
            if( sr->true_particles[i_p].G4ID==MuonTrk_Truth_G4ID ){
              rets.push_back( sr->true_particles[i_p].pdg );
            }
          }

        }

      }
    }

    return rets;

  });

  const SpillMultiVar kNuMI_SignalSelectionSlices_Others_MuonTrackMatched_interaction_id ( [](const caf::SRSpillProxy *sr) {

    std::vector<double> rets;
    for(std::size_t i(0); i < sr->slc.size(); ++i){
      const auto& slc = sr->slc.at(i);
      if( kNuMISelection_1muNp0pi(&slc) ){
        if( kNuMISliceSignalType(&slc)==5 ){

          const auto& MuonIdx = kNuMIMuonCandidateIdx(&slc);
          const auto& MuonTrk = slc.reco.pfp[MuonIdx].trk;

          int MuonTrk_Truth_G4ID = MuonTrk.truth.p.G4ID;

          for(unsigned int i_p=0; i_p<sr->true_particles.size(); i_p++){
            if( sr->true_particles[i_p].G4ID==MuonTrk_Truth_G4ID ){
              rets.push_back( sr->true_particles[i_p].interaction_id );
            }
          }

        }

      }
    }

    return rets;

  });

  const SpillMultiVar kNuMI_TrueNeutrino_PosX ( [](const caf::SRSpillProxy *sr) -> std::vector<double> {

    std::vector<double> rets;
    for(std::size_t i(0); i < sr->mc.nu.size(); ++i){
      rets.push_back(sr->mc.nu[i].position.x);
    }

    return rets;

  });

  const SpillMultiVar kNuMI_TrueNeutrino_PosY ( [](const caf::SRSpillProxy *sr) -> std::vector<double> {

    std::vector<double> rets;
    for(std::size_t i(0); i < sr->mc.nu.size(); ++i){
      rets.push_back(sr->mc.nu[i].position.y);
    }

    return rets;

  });

  const SpillMultiVar kNuMI_TrueNeutrino_PosZ ( [](const caf::SRSpillProxy *sr) -> std::vector<double> {

    std::vector<double> rets;
    for(std::size_t i(0); i < sr->mc.nu.size(); ++i){
      rets.push_back(sr->mc.nu[i].position.z);
    }

    return rets;

  });

  const SpillVar kNuMI_trigger_within_gate ( [](const caf::SRSpillProxy *sr) -> double {
    return sr->hdr.triggerinfo.trigger_within_gate;
  });
  const SpillVar kNuMI_NumberOfNeutrinos ( [](const caf::SRSpillProxy *sr) -> int {
    return sr->mc.nu.size();
  });
  const SpillVar kNuMI_TriggerNeutrino_Idx ( [](const caf::SRSpillProxy *sr) -> int {

    double spillTriggerTime = kNuMISpillTriggerTime(sr);
    if( ! (spillTriggerTime > -0.1 && spillTriggerTime < 10.1) ) return -1;

    double ClosestTime = 9999.;
    int TrigNuIdx = -2;
    for(std::size_t i(0); i < sr->mc.nu.size(); ++i){
      double this_time_diff = fabs( spillTriggerTime - sr->mc.nu[i].time );
      if( this_time_diff < ClosestTime ){
        ClosestTime = this_time_diff;
        TrigNuIdx = i;
      }
    }

    return TrigNuIdx;


  });
  const SpillVar kNuMI_TriggerNeutrino_time ( [](const caf::SRSpillProxy *sr) -> double {

    double spillTriggerTime = kNuMISpillTriggerTime(sr);
    if( ! (spillTriggerTime > -0.1 && spillTriggerTime < 10.1) ) return 99999999.;

    int TrigNuIdx = kNuMI_TriggerNeutrino_Idx(sr);
    if(TrigNuIdx<0) return 99999999.;

    const auto& nu = sr->mc.nu[TrigNuIdx];
    return nu.time;

  });
  const SpillVar kNuMI_TriggerNeutrino_PosX ( [](const caf::SRSpillProxy *sr) -> double {

    double spillTriggerTime = kNuMISpillTriggerTime(sr);
    if( ! (spillTriggerTime > -0.1 && spillTriggerTime < 10.1) ) return 99999999.;

    int TrigNuIdx = kNuMI_TriggerNeutrino_Idx(sr);
    if(TrigNuIdx<0) return 99999999.;

    const auto& nu = sr->mc.nu[TrigNuIdx];
    return nu.position.x;

  });
  const SpillVar kNuMI_TriggerNeutrino_PosY ( [](const caf::SRSpillProxy *sr) -> double {

    double spillTriggerTime = kNuMISpillTriggerTime(sr);
    if( ! (spillTriggerTime > -0.1 && spillTriggerTime < 10.1) ) return 99999999.;

    int TrigNuIdx = kNuMI_TriggerNeutrino_Idx(sr);
    if(TrigNuIdx<0) return 99999999.;

    const auto& nu = sr->mc.nu[TrigNuIdx];
    return nu.position.y;

  });
  const SpillVar kNuMI_TriggerNeutrino_PosZ ( [](const caf::SRSpillProxy *sr) -> double {

    double spillTriggerTime = kNuMISpillTriggerTime(sr);
    if( ! (spillTriggerTime > -0.1 && spillTriggerTime < 10.1) ) return 99999999.;

    int TrigNuIdx = kNuMI_TriggerNeutrino_Idx(sr);
    if(TrigNuIdx<0) return 99999999.;

    const auto& nu = sr->mc.nu[TrigNuIdx];
    return nu.position.z;

  });
  const SpillVar kNuMI_HasIntimeCosmic ( [](const caf::SRSpillProxy *sr) -> int {

    bool HasIntimeCosmic = false;
    for(unsigned int i_p=0; i_p<sr->true_particles.size(); i_p++){
      //if( sr->true_particles[i_p].interaction_id==-1 && abs(sr->true_particles[i_p].pdg)==13 ){
      if( sr->true_particles[i_p].interaction_id==-1 ){
        if( sr->true_particles[i_p].genT>-0.1 && sr->true_particles[i_p].genT<10.1 ){
          HasIntimeCosmic = true;
          break;
        }
      }
    }

    if(HasIntimeCosmic) return 1;
    else return 0;

  });
  const SpillVar kNuMI_NIntimeCosmic ( [](const caf::SRSpillProxy *sr) -> int {

    int NIntimeCosmic = 0;
    for(unsigned int i_p=0; i_p<sr->true_particles.size(); i_p++){
      //if( sr->true_particles[i_p].interaction_id==-1 && abs(sr->true_particles[i_p].pdg)==13 ){
      if( sr->true_particles[i_p].interaction_id==-1 ){
        if( sr->true_particles[i_p].genT>-0.1 && sr->true_particles[i_p].genT<10.1 ){
          NIntimeCosmic++;
        }
      }
    }

    return NIntimeCosmic;

  });
  const SpillVar kNuMI_IntimeCosmicClosestTime ( [](const caf::SRSpillProxy *sr) -> double {

    double spillTriggerTime = kNuMISpillTriggerTime(sr);
    double ClosestTime = 9999.;
    for(unsigned int i_p=0; i_p<sr->true_particles.size(); i_p++){
      //if( sr->true_particles[i_p].interaction_id==-1 && abs(sr->true_particles[i_p].pdg)==13 ){
      if( sr->true_particles[i_p].interaction_id==-1 ){
        if( sr->true_particles[i_p].genT>-0.1 && sr->true_particles[i_p].genT<10.1 ){
          double this_timediff = fabs(sr->true_particles[i_p].genT  - spillTriggerTime);
          if(this_timediff<ClosestTime){
            ClosestTime = this_timediff;
          }
        }
      }
    }

    return ClosestTime;

  });
  const SpillVar kNuMI_TriggerNeutrino_PrimaryMuonIdx ( [](const caf::SRSpillProxy *sr) -> int {

    double spillTriggerTime = kNuMISpillTriggerTime(sr);
    if( ! (spillTriggerTime > -0.1 && spillTriggerTime < 10.1) ) return -1;

    int TrigNuIdx = kNuMI_TriggerNeutrino_Idx(sr);
    if(TrigNuIdx<0) return -2;

    const auto& nu = sr->mc.nu[TrigNuIdx];
    int MuonInd = -1;
    for(unsigned int i_prim=0; i_prim<nu.prim.size(); i_prim++){

      const auto& prim = nu.prim[i_prim];
      if(abs(prim.pdg)==13){
        MuonInd = i_prim;
        break;
      }

    }
    if(MuonInd==-1){

      for(unsigned int i_prim=0; i_prim<nu.prim.size(); i_prim++){

        const auto& prim = nu.prim[i_prim];
        if(abs(prim.pdg)==211){
          MuonInd = i_prim;
          break;
        }

      }

    }

    return MuonInd;

  });
  const SpillVar kNuMI_TriggerNeutrino_PrimaryMuon_pdg ( [](const caf::SRSpillProxy *sr) -> double {

    double spillTriggerTime = kNuMISpillTriggerTime(sr);
    if( ! (spillTriggerTime > -0.1 && spillTriggerTime < 10.1) ) return -1;

    int TrigNuIdx = kNuMI_TriggerNeutrino_Idx(sr);
    if(TrigNuIdx<0) return -2;

    int PrimMuIdx = kNuMI_TriggerNeutrino_PrimaryMuonIdx(sr);
    if(PrimMuIdx<0) return 9999999.;

    const auto& mu = sr->mc.nu[TrigNuIdx].prim[PrimMuIdx];

    return mu.pdg;

  });
  const SpillVar kNuMI_TriggerNeutrino_PrimaryMuon_genE ( [](const caf::SRSpillProxy *sr) -> double {

    double spillTriggerTime = kNuMISpillTriggerTime(sr);
    if( ! (spillTriggerTime > -0.1 && spillTriggerTime < 10.1) ) return -1;

    int TrigNuIdx = kNuMI_TriggerNeutrino_Idx(sr);
    if(TrigNuIdx<0) return -2;

    int PrimMuIdx = kNuMI_TriggerNeutrino_PrimaryMuonIdx(sr);
    if(PrimMuIdx<0) return 9999999.;

    const auto& mu = sr->mc.nu[TrigNuIdx].prim[PrimMuIdx];

    return mu.genE;

  });
  const SpillVar kNuMI_TriggerNeutrino_PrimaryMuon_GenX ( [](const caf::SRSpillProxy *sr) -> double {

    double spillTriggerTime = kNuMISpillTriggerTime(sr);
    if( ! (spillTriggerTime > -0.1 && spillTriggerTime < 10.1) ) return -1;

    int TrigNuIdx = kNuMI_TriggerNeutrino_Idx(sr);
    if(TrigNuIdx<0) return -2;

    int PrimMuIdx = kNuMI_TriggerNeutrino_PrimaryMuonIdx(sr);
    if(PrimMuIdx<0) return 9999999.;

    const auto& mu = sr->mc.nu[TrigNuIdx].prim[PrimMuIdx];

    return mu.gen.x;

  });
  const SpillVar kNuMI_TriggerNeutrino_PrimaryMuon_GenY ( [](const caf::SRSpillProxy *sr) -> double {

    double spillTriggerTime = kNuMISpillTriggerTime(sr);
    if( ! (spillTriggerTime > -0.1 && spillTriggerTime < 10.1) ) return -1;

    int TrigNuIdx = kNuMI_TriggerNeutrino_Idx(sr);
    if(TrigNuIdx<0) return -2;

    int PrimMuIdx = kNuMI_TriggerNeutrino_PrimaryMuonIdx(sr);
    if(PrimMuIdx<0) return 9999999.;

    const auto& mu = sr->mc.nu[TrigNuIdx].prim[PrimMuIdx];

    return mu.gen.y;

  });
  const SpillVar kNuMI_TriggerNeutrino_PrimaryMuon_GenZ ( [](const caf::SRSpillProxy *sr) -> double {

    double spillTriggerTime = kNuMISpillTriggerTime(sr);
    if( ! (spillTriggerTime > -0.1 && spillTriggerTime < 10.1) ) return -1;

    int TrigNuIdx = kNuMI_TriggerNeutrino_Idx(sr);
    if(TrigNuIdx<0) return -2;

    int PrimMuIdx = kNuMI_TriggerNeutrino_PrimaryMuonIdx(sr);
    if(PrimMuIdx<0) return 9999999.;

    const auto& mu = sr->mc.nu[TrigNuIdx].prim[PrimMuIdx];

    return mu.gen.z;

  });
  const SpillVar kNuMI_TriggerNeutrino_PrimaryMuon_StartX ( [](const caf::SRSpillProxy *sr) -> double {
    
    double spillTriggerTime = kNuMISpillTriggerTime(sr); 
    if( ! (spillTriggerTime > -0.1 && spillTriggerTime < 10.1) ) return -1;
    
    int TrigNuIdx = kNuMI_TriggerNeutrino_Idx(sr);
    if(TrigNuIdx<0) return -2;
    
    int PrimMuIdx = kNuMI_TriggerNeutrino_PrimaryMuonIdx(sr);
    if(PrimMuIdx<0) return 9999999.;

    const auto& mu = sr->mc.nu[TrigNuIdx].prim[PrimMuIdx];

    return mu.start.x;

  });
  const SpillVar kNuMI_TriggerNeutrino_PrimaryMuon_StartY ( [](const caf::SRSpillProxy *sr) -> double {

    double spillTriggerTime = kNuMISpillTriggerTime(sr);
    if( ! (spillTriggerTime > -0.1 && spillTriggerTime < 10.1) ) return -1;

    int TrigNuIdx = kNuMI_TriggerNeutrino_Idx(sr);
    if(TrigNuIdx<0) return -2;

    int PrimMuIdx = kNuMI_TriggerNeutrino_PrimaryMuonIdx(sr);
    if(PrimMuIdx<0) return 9999999.;

    const auto& mu = sr->mc.nu[TrigNuIdx].prim[PrimMuIdx];

    return mu.start.y;

  });
  const SpillVar kNuMI_TriggerNeutrino_PrimaryMuon_StartZ ( [](const caf::SRSpillProxy *sr) -> double {

    double spillTriggerTime = kNuMISpillTriggerTime(sr);
    if( ! (spillTriggerTime > -0.1 && spillTriggerTime < 10.1) ) return -1;

    int TrigNuIdx = kNuMI_TriggerNeutrino_Idx(sr);
    if(TrigNuIdx<0) return -2;

    int PrimMuIdx = kNuMI_TriggerNeutrino_PrimaryMuonIdx(sr);
    if(PrimMuIdx<0) return 9999999.;

    const auto& mu = sr->mc.nu[TrigNuIdx].prim[PrimMuIdx];

    return mu.start.z;

  });
  const SpillVar kNuMI_TriggerNeutrino_PrimaryMuon_EndX ( [](const caf::SRSpillProxy *sr) -> double {

    double spillTriggerTime = kNuMISpillTriggerTime(sr);
    if( ! (spillTriggerTime > -0.1 && spillTriggerTime < 10.1) ) return -1;

    int TrigNuIdx = kNuMI_TriggerNeutrino_Idx(sr);
    if(TrigNuIdx<0) return -2;

    int PrimMuIdx = kNuMI_TriggerNeutrino_PrimaryMuonIdx(sr);
    if(PrimMuIdx<0) return 9999999.;

    const auto& mu = sr->mc.nu[TrigNuIdx].prim[PrimMuIdx];

    return mu.end.x;

  });
  const SpillVar kNuMI_TriggerNeutrino_PrimaryMuon_EndY ( [](const caf::SRSpillProxy *sr) -> double {

    double spillTriggerTime = kNuMISpillTriggerTime(sr);
    if( ! (spillTriggerTime > -0.1 && spillTriggerTime < 10.1) ) return -1;

    int TrigNuIdx = kNuMI_TriggerNeutrino_Idx(sr);
    if(TrigNuIdx<0) return -2;

    int PrimMuIdx = kNuMI_TriggerNeutrino_PrimaryMuonIdx(sr);
    if(PrimMuIdx<0) return 9999999.;

    const auto& mu = sr->mc.nu[TrigNuIdx].prim[PrimMuIdx];

    return mu.end.y;

  });
  const SpillVar kNuMI_TriggerNeutrino_PrimaryMuon_EndZ ( [](const caf::SRSpillProxy *sr) -> double {

    double spillTriggerTime = kNuMISpillTriggerTime(sr);
    if( ! (spillTriggerTime > -0.1 && spillTriggerTime < 10.1) ) return -1;

    int TrigNuIdx = kNuMI_TriggerNeutrino_Idx(sr);
    if(TrigNuIdx<0) return -2;

    int PrimMuIdx = kNuMI_TriggerNeutrino_PrimaryMuonIdx(sr);
    if(PrimMuIdx<0) return 9999999.;

    const auto& mu = sr->mc.nu[TrigNuIdx].prim[PrimMuIdx];

    return mu.end.z;

  });

  const SpillVar kNuMI_TriggerTrueParticle_Idx ( [](const caf::SRSpillProxy *sr) -> int {

    double spillTriggerTime = kNuMISpillTriggerTime(sr);
    if( ! (spillTriggerTime > -0.1 && spillTriggerTime < 10.1) ) return -1;

    double ClosestTime = 9999.;
    int TrigTPIdx = -2;
    for(unsigned int i_p=0; i_p<sr->true_particles.size(); i_p++){
      double this_time_diff = fabs( spillTriggerTime - sr->true_particles[i_p].genT );
      if( this_time_diff < ClosestTime ){
        ClosestTime = this_time_diff;
        TrigTPIdx = i_p;
      }
    }

    return TrigTPIdx;


  });
  const SpillVar kNuMI_TriggerTrueParticle_pdg ( [](const caf::SRSpillProxy *sr) -> int {
    
    double spillTriggerTime = kNuMISpillTriggerTime(sr); 
    if( ! (spillTriggerTime > -0.1 && spillTriggerTime < 10.1) ) return -1;
    
    int TrigTPIdx = kNuMI_TriggerTrueParticle_Idx(sr);
    if(TrigTPIdx<0) return 99999999;

    else return sr->true_particles[TrigTPIdx].pdg;


  });
  const SpillVar kNuMI_SignalSelectionSlice_Idx( [](const caf::SRSpillProxy *sr) -> int {
    int ret = -1;
    for(std::size_t i(0); i < sr->slc.size(); ++i){
      const auto& slc = sr->slc.at(i);
      if( kNuMISelection_1muNp0pi(&slc) ){
        ret = i;
        break;
      }
    }
    return ret;
  });
  const SpillVar kNuMI_SignalSelectionSlice_Truth_time( [](const caf::SRSpillProxy *sr) -> double {

    int SigSelSliceIdx = kNuMI_SignalSelectionSlice_Idx(sr);
    if(SigSelSliceIdx<0) return -99999.;
    const auto& slc = sr->slc[SigSelSliceIdx];

    return slc.truth.time;

  });
  const SpillVar kNuMI_SignalSelectionSlice_Truth_vtx_x( [](const caf::SRSpillProxy *sr) -> double {
    
    int SigSelSliceIdx = kNuMI_SignalSelectionSlice_Idx(sr);
    if(SigSelSliceIdx<0) return -99999.;
    const auto& slc = sr->slc[SigSelSliceIdx];

    return slc.truth.prod_vtx.x;

  });
  const SpillVar kNuMI_SignalSelectionSlice_Truth_vtx_y( [](const caf::SRSpillProxy *sr) -> double {

    int SigSelSliceIdx = kNuMI_SignalSelectionSlice_Idx(sr);
    if(SigSelSliceIdx<0) return -99999.;
    const auto& slc = sr->slc[SigSelSliceIdx];

    return slc.truth.prod_vtx.y;

  });
  const SpillVar kNuMI_SignalSelectionSlice_Truth_vtx_z( [](const caf::SRSpillProxy *sr) -> double {

    int SigSelSliceIdx = kNuMI_SignalSelectionSlice_Idx(sr);
    if(SigSelSliceIdx<0) return -99999.;
    const auto& slc = sr->slc[SigSelSliceIdx];

    return slc.truth.prod_vtx.z;

  });
  const SpillVar kNuMI_MuonTrackMatchedTPIdx( [](const caf::SRSpillProxy *sr) -> int {

    int SigSelSliceIdx = kNuMI_SignalSelectionSlice_Idx(sr);
    if(SigSelSliceIdx<0) return -2;
    const auto& slc = sr->slc[SigSelSliceIdx];

    const auto& MuonIdx = kNuMIMuonCandidateIdx(&slc);
    const auto& MuonTrk = slc.reco.pfp[MuonIdx].trk;
    int MuonTrk_Truth_G4ID = MuonTrk.truth.p.G4ID;

    int TPIdx = -1;
    for(unsigned int i_p=0; i_p<sr->true_particles.size(); i_p++){
      if( sr->true_particles[i_p].G4ID==MuonTrk_Truth_G4ID ){
        TPIdx = i_p;
        break;
      }
    }

    return TPIdx;

  });
  const SpillVar kNuMI_MuonTrackMatchedTP_genT( [](const caf::SRSpillProxy *sr) -> double {

    int TPIdx = kNuMI_MuonTrackMatchedTPIdx(sr);
    if(TPIdx<0) return -999999.;

    const auto& tp = sr->true_particles[TPIdx];
    return tp.genT;

  });
  const SpillVar kNuMI_MuonTrackMatchedTP_pdg( [](const caf::SRSpillProxy *sr) -> int {

    int TPIdx = kNuMI_MuonTrackMatchedTPIdx(sr);
    if(TPIdx<0) return -999999.;

    const auto& tp = sr->true_particles[TPIdx];
    return tp.pdg;

  });
  const SpillVar kNuMI_MuonTrackMatchedTP_genE( [](const caf::SRSpillProxy *sr) -> double {

    int TPIdx = kNuMI_MuonTrackMatchedTPIdx(sr);
    if(TPIdx<0) return -999999.;

    const auto& tp = sr->true_particles[TPIdx];
    return tp.genE;

  });
  const SpillVar kNuMI_ProtonTrackMatchedTPIdx( [](const caf::SRSpillProxy *sr) -> int {
    
    int SigSelSliceIdx = kNuMI_SignalSelectionSlice_Idx(sr);
    if(SigSelSliceIdx<0) return -2;
    const auto& slc = sr->slc[SigSelSliceIdx];
    
    const auto& ProtonIdx = kNuMIProtonCandidateIdx(&slc);
    const auto& ProtonTrk = slc.reco.pfp[ProtonIdx].trk;
    int ProtonTrk_Truth_G4ID = ProtonTrk.truth.p.G4ID;
    
    int TPIdx = -1;
    for(unsigned int i_p=0; i_p<sr->true_particles.size(); i_p++){
      if( sr->true_particles[i_p].G4ID==ProtonTrk_Truth_G4ID ){
        TPIdx = i_p;
        break;
      }
    }

    return TPIdx;

  });
  const SpillVar kNuMI_ProtonTrackMatchedTP_genT( [](const caf::SRSpillProxy *sr) -> double {

    int TPIdx = kNuMI_ProtonTrackMatchedTPIdx(sr);
    if(TPIdx<0) return -999999.;

    const auto& tp = sr->true_particles[TPIdx];
    return tp.genT;

  });
  const SpillVar kNuMI_ProtonTrackMatchedTP_pdg( [](const caf::SRSpillProxy *sr) -> int {

    int TPIdx = kNuMI_ProtonTrackMatchedTPIdx(sr);
    if(TPIdx<0) return -999999.;

    const auto& tp = sr->true_particles[TPIdx];
    return tp.pdg;

  }); 
  const SpillVar kNuMI_ProtonTrackMatchedTP_genE( [](const caf::SRSpillProxy *sr) -> double {
  
    int TPIdx = kNuMI_ProtonTrackMatchedTPIdx(sr);
    if(TPIdx<0) return -999999.;

    const auto& tp = sr->true_particles[TPIdx];
    return tp.genE;
    
  });




  const SpillMultiVar kNuMI_IntimeCosmics_genT ( [](const caf::SRSpillProxy *sr) -> std::vector<double> {

    std::vector<double> rets;
    for(unsigned int i_p=0; i_p<sr->true_particles.size(); i_p++){
      if( sr->true_particles[i_p].interaction_id==-1 && abs(sr->true_particles[i_p].pdg)==13 ){
        if( sr->true_particles[i_p].genT>-0.1 && sr->true_particles[i_p].genT<10.1 ){
          rets.push_back( sr->true_particles[i_p].genT );
        }
      }
    }

    return rets;

  });

  const SpillMultiVar kNuMI_IntimeCosmics_TimeFromTrig ( [](const caf::SRSpillProxy *sr) -> std::vector<double> {

    std::vector<double> rets;
    double spillTriggerTime = kNuMISpillTriggerTime(sr);
    for(unsigned int i_p=0; i_p<sr->true_particles.size(); i_p++){
      if( sr->true_particles[i_p].interaction_id==-1 && abs(sr->true_particles[i_p].pdg)==13 ){
        if( sr->true_particles[i_p].genT>-0.1 && sr->true_particles[i_p].genT<10.1 ){
          rets.push_back( sr->true_particles[i_p].genT - spillTriggerTime );
        }
      }
    }

    return rets;

  });


  const SpillCut kNuMI_Spill_NoIntimeCosmic ( [](const caf::SRSpillProxy *sr) {
    return ( kNuMI_HasIntimeCosmic(sr)== 0);
  });
  const SpillCut kNuMI_Spill_HasIntimeNu ( [](const caf::SRSpillProxy *sr) {
    double TrigTime = kNuMI_trigger_within_gate(sr);
    double TrigNuTime = kNuMI_TriggerNeutrino_time(sr);

    if(TrigNuTime>9999) return false;
    else{
      return ( fabs(TrigNuTime - TrigTime) < 1.0 );
    }

  });

}
