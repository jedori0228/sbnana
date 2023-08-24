#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_Variables.h"

using namespace std;
using namespace ana;

namespace ICARUSNumuXsec{

  // - Test
  const SpillMultiVar spillvarTest([](const caf::SRSpillProxy *sr) -> vector<double> {

    std::vector<double> rets;

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
  // - Some truth
  const Var varTruePDG([](const caf::SRSliceProxy* slc) -> int {
    return kNuMITruePDG(slc);
  });

  const Var varTrueTarget([](const caf::SRSliceProxy* slc) -> int {
    return kNuMITrueTarget(slc);
  });
  // - GENIE interaction code
  // - https://internal.dunescience.org/doxygen/namespacesimb.html#a2cce734d1b71408bbc7d98d148ac4360
  const Var varGENIEIntCode([](const caf::SRSliceProxy* slc) -> int {
    return kNuMITrueMode(slc);
  });
  // - Truth interaction
  const Var varNeutrinoTruthE([](const caf::SRSliceProxy* slc) -> double {
    return kNuMITrueNuE(slc);
  });
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
    std::cout << "[JSKIMDEBUG][cutIsQELike] called" << std::endl;
    printf("- genie_mode: %d\n",slc->truth.genie_mode.GetValue());
    int n_proton=0;
    int n_neutron=0;
    int n_pip=0;
    int n_pim=0;
    int n_pi0=0;
    for(unsigned i=0; i<slc->truth.ghepptl.size(); i++){
      const auto& ghepptl = slc->truth.ghepptl.at(i);
      if(ghepptl.gstatus!=1) continue;
      const int& pdg = ghepptl.pdg;
      if( pdg==2212 ) n_proton++;
      else if( pdg==2112 ) n_neutron++;
      else if( pdg==211 ) n_pip++;
      else if( pdg==-211 ) n_pim++;
      else if( pdg==111 ) n_pi0++;
      else{

      }

    } // END ghep loop

    printf("(p, n, pi+, pi-, pi0) = (%d, %d, %d, %d, %d)\n", n_proton, n_neutron, n_pip, n_pim, n_pi0);
*/

    std::cout << "[JSKIMDEBUG][cutIsQELike] called" << std::endl;
    printf("- genie_mode: %d\n",slc->truth.genie_mode.GetValue());
    int n_muon=0;
    int n_proton=0;
    int n_neutron=0;
    int n_pip=0;
    int n_pim=0;
    int n_pi0=0;

    for(const auto& prm: slc->truth.prim){
      const int pdg = prm.pdg;
      const int apdg = abs(pdg);
      if(apdg==13) n_muon++;

      if( pdg==2212 ) n_proton++;
      else if( pdg==2112 ) n_neutron++;
      else if( pdg==211 ) n_pip++;
      else if( pdg==-211 ) n_pim++;
      else if( pdg==111 ) n_pi0++;

    } // END prim loop

    printf("(mu, p, n, pi+, pi-, pi0) = (%d, %d, %d, %d, %d, %d)\n", n_muon, n_proton, n_neutron, n_pip, n_pim, n_pi0);


    return 0;
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

}
