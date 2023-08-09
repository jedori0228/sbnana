// Selections and variables to plot -- the main stuff

#include "sbnana/CAFAna/Core/Spectrum.h"
#include "sbnana/CAFAna/Core/Binning.h"
#include "sbnana/CAFAna/Core/Var.h"

#include "sbnanaobj/StandardRecord/SRTrueInteraction.h"
#include "sbnana/SBNAna/Vars/Vars.h"

#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

// ROOT
#include "TMath.h"
#include "TVector3.h"
#include "TGraphAsymmErrors.h"
#include "TRandom3.h"
#include "TPaveText.h"
#include "TColor.h"
#include "TLegend.h"
#include "TH1D.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TGraph.h"

// C++
#include <vector>
#include <string>
#include <map>
#include <utility>
#include <cassert>

// The defs file
#include "sbnana/SBNAna/icarus-numi-numucc-xsec/helper_defs.h"

// Trk + Shw dual fits
#include "sbnana/SBNAna/icarus-numi-numucc-xsec/TrkAndShw.h"

using namespace ana;

////////////////////////////////////////////////////////////////
// The 1muNp0pi selection:
////////////////////////////////////////////////////////////////

const SpillVar kSpillTriggerTime ( [](const caf::SRSpillProxy *sr) -> double {
    double triggerTime = 0.;

    double foundTriggerTime = sr->hdr.triggerinfo.trigger_within_gate;
    if ( !std::isnan(foundTriggerTime) && !std::isinf(foundTriggerTime) && foundTriggerTime < -15. ) triggerTime = -15.;
    else if ( std::isnan(foundTriggerTime) || std::isinf(foundTriggerTime) ) triggerTime = -16.;
    else if ( foundTriggerTime > 30. ) triggerTime = 30.;
    else triggerTime = foundTriggerTime;

    return triggerTime;
  });

const SpillCut kIsValidTriggerNuMI ( [](const caf::SRSpillProxy *sr) {
    double spillTriggerTime = kSpillTriggerTime(sr);
    return spillTriggerTime > -0.1 && spillTriggerTime < 9.7;
  });

const Cut kRFiducialNew([](const caf::SRSliceProxy* slc) {
    if ( std::isnan(slc->vertex.x) || std::isnan(slc->vertex.y) || std::isnan(slc->vertex.z) ) return false;
    return isInFV(slc->vertex.x,slc->vertex.y,slc->vertex.z);
  });

const Cut kNotClearCosmic([](const caf::SRSliceProxy* slc) {
    return !slc->is_clear_cosmic;
  });

// Muon candidate - both with track-like tracks and all tracks:
// Redo the slice vars from NumuVars in SBNAna to use 10cm for containment
//   (and to consider the other cryostat), and now even more so to deal with shower-like tracks...
const Var kPTrackIndNew([](const caf::SRSliceProxy* slc) -> int {
    float Longest(0);
    int PTrackInd(-1);

    unsigned int idxTrk = 0;
    while ( IsValidTrkIdx(slc, idxTrk) ) {
      if ( !IsTracklikeTrack(slc, idxTrk) ) { 
        idxTrk+=1;
        continue;
      }
      auto const& trk = slc->reco.pfp.at(idxTrk).trk;
      unsigned int thisIdx = idxTrk;
      idxTrk+=1;

      if ( std::isnan(trk.start.x) || std::isnan(trk.len) || trk.len <= 0. ) continue;
      if ( std::isnan(slc->vertex.x) || std::isnan(slc->vertex.y) || std::isnan(slc->vertex.z) ) return -1;
      const float Atslc = std::hypot(slc->vertex.x - trk.start.x,
                                     slc->vertex.y - trk.start.y,
                                     slc->vertex.z - trk.start.z);
      const bool isPrimCandidate = (Atslc < 10. && IsPrimaryPFP(slc,thisIdx));

      if ( !isPrimCandidate || trk.calo[2].nhit < 5 ) continue;
      const bool Contained = isContainedVol(trk.end.x,trk.end.y,trk.end.z);
      const float Chi2Proton = trk.chi2pid[2].chi2_proton;
      const float Chi2Muon = trk.chi2pid[2].chi2_muon;
      if ( (!Contained && trk.len > 100.) || (Contained && trk.len > 50. && Chi2Proton > 60. && Chi2Muon < 30.) ) {
        if ( trk.len <= Longest ) continue;
        Longest = trk.len;
        PTrackInd = thisIdx;
      }
    }

    return PTrackInd;
  });

const Cut kPTrackNew([](const caf::SRSliceProxy* slc) {
    return ( kPTrackIndNew(slc) >= 0 );
  });

const Var kProtonTrackIndNew([](const caf::SRSliceProxy* slc) -> int {
    int primaryInd = kPTrackIndNew(slc);

    float Longest(0);
    int PTrackInd(-1);

    unsigned int idxTrk = 0;
    while ( IsValidTrkIdx(slc, idxTrk) ) {
      int thisIdxInt = idxTrk;
      if ( thisIdxInt == primaryInd ) {
        idxTrk+=1;
        continue; // skip the particle which is the muon candidate!
      }
      if ( !IsTracklikeTrack(slc, idxTrk) ) { 
        idxTrk+=1;
        continue;
      }
      auto const& trk = slc->reco.pfp.at(idxTrk).trk; //GetTrack( slc, idxTrk );
      unsigned int thisIdx = idxTrk;
      idxTrk+=1;

      if ( std::isnan(trk.start.x) || std::isnan(trk.len) || trk.len <= 0. ) continue;
      if ( std::isnan(slc->vertex.x) || std::isnan(slc->vertex.y) || std::isnan(slc->vertex.z) ) return -1;
      const float Atslc = std::hypot(slc->vertex.x - trk.start.x,
                                     slc->vertex.y - trk.start.y,
                                     slc->vertex.z - trk.start.z);
      const bool isPrimCandidate = (Atslc < 10. && IsPrimaryPFP(slc,thisIdx));

      if ( !isPrimCandidate || trk.calo[2].nhit < 5 ) continue;
      const bool Contained = isContainedVol(trk.end.x,trk.end.y,trk.end.z);
      const float Chi2Proton = trk.chi2pid[2].chi2_proton;
      const float Chi2Muon = trk.chi2pid[2].chi2_muon;
      if ( Contained && Chi2Proton < 90. && Chi2Muon > 30. ) {
        if ( trk.len <= Longest ) continue;
        Longest = trk.len;
        PTrackInd = thisIdx;
      }
    }

    return PTrackInd;
  });

const Cut kProtonTrack([](const caf::SRSliceProxy* slc) {
    return ( kProtonTrackIndNew(slc) >= 0 );
  });

const Cut kAllPrimaryHadronsContained([](const caf::SRSliceProxy* slc) {
    // Considers all track fits, not just track-like PFPs
    int primaryInd = kPTrackIndNew(slc);
    if ( primaryInd < 0 ) return false;
    unsigned int idxPrim = (unsigned int)primaryInd;

    unsigned int idxTrk = 0;
    while ( IsValidTrkIdx(slc, idxTrk) ) {
      int thisIdxInt = idxTrk;
      if ( thisIdxInt == primaryInd ) {
        idxTrk+=1;
        continue; // skip the particle which is the muon candidate!
      }
      auto const& trk = slc->reco.pfp.at(idxTrk).trk;
      unsigned int thisIdx = idxTrk;
      idxTrk+=1;

      if ( std::isnan(trk.start.x) || std::isnan(trk.len) || trk.len <= 0. ) continue;
      if ( std::isnan(slc->vertex.x) || std::isnan(slc->vertex.y) || std::isnan(slc->vertex.z) ) return false;
      const float Atslc = std::hypot(slc->vertex.x - trk.start.x,
                                     slc->vertex.y - trk.start.y,
                                     slc->vertex.z - trk.start.z);
      const bool isPrimCandidate = (Atslc < 10. && IsPrimaryPFP(slc,thisIdx));

      if ( !isPrimCandidate ) continue;
      const bool Contained = isContainedVol(trk.end.x,trk.end.y,trk.end.z);
      if ( !Contained ) return false;
    }

    return true;
  });

const Cut kNoSecondPrimaryMuonlikeTracks([](const caf::SRSliceProxy* slc) {
    int primaryInd = kPTrackIndNew(slc);
    if ( primaryInd < 0 ) return false;
    unsigned int idxPrim = (unsigned int)primaryInd;

    unsigned int idxTrk = 0;
    while ( IsValidTrkIdx(slc, idxTrk) ) {
      int thisIdxInt = idxTrk;
      if ( thisIdxInt == primaryInd ) {
        idxTrk+=1;
        continue; // skip the particle which is the muon candidate!
      }
      auto const& trk = slc->reco.pfp.at(idxTrk).trk;
      unsigned int thisIdx = idxTrk;
      idxTrk+=1;

      if ( std::isnan(trk.start.x) || std::isnan(trk.len) || trk.len <= 0. ) continue;
      if ( std::isnan(slc->vertex.x) || std::isnan(slc->vertex.y) || std::isnan(slc->vertex.z) ) return false;
      const float Atslc = std::hypot(slc->vertex.x - trk.start.x,
                                     slc->vertex.y - trk.start.y,
                                     slc->vertex.z - trk.start.z);
      const bool isPrimCandidate = (Atslc < 10. && IsPrimaryPFP(slc,thisIdx));

      if ( !isPrimCandidate || trk.calo[2].nhit < 5 ) continue;
      const float Chi2Proton = trk.chi2pid[2].chi2_proton;
      const float Chi2Muon = trk.chi2pid[2].chi2_muon;
      if ( Chi2Proton > 60. && Chi2Muon < 30. ) return false;
    }

    return true;
  });

// NOTE: this selection does NOT include the MIP-without-Bragg variable.

const Cut kRecoProtonPInKinematicTreshold([](const caf::SRSliceProxy* slc) {
    float p(-5.f);

    if ( kProtonTrackIndNew(slc) >= 0 )
    {
      auto const& trk = slc->reco.pfp.at(kProtonTrackIndNew(slc)).trk;
      p = trk.rangeP.p_proton;
    }

    return (p > 0.4 && p < 1.0);
  });

// And the attempt to cut photons (without rejecting extra stuff considered protons and muons... needs more looking...)
const Cut kCutPhotons([](const caf::SRSliceProxy* slc) {
    int primaryInd = kPTrackIndNew(slc);
    if ( primaryInd < 0 ) return false;
    unsigned int idxPrim = (unsigned int)primaryInd;

    int primaryProtonInd = kProtonTrackIndNew(slc);
    if ( primaryProtonInd < 0 ) return false;
    unsigned int idxPrimProton = (unsigned int)primaryProtonInd;

    unsigned int idxTrk = 0;
    while ( IsValidTrkIdx(slc, idxTrk) ) {
      int thisIdxInt = idxTrk;
      if ( thisIdxInt == primaryInd || thisIdxInt == primaryProtonInd ) {
        idxTrk+=1;
        continue; // skip the particle which is the muon or leading proton candidate!
      }
      if ( !IsShowerlike(slc, idxTrk) ) { 
        idxTrk+=1;
        continue; // skip things with track score > 0.45
      }
      auto const& shw = slc->reco.pfp.at(idxTrk).shw;
      unsigned int thisIdx = idxTrk;
      idxTrk+=1;

      // Check if shower fit even seems kind-of valid:
      if ( std::isnan(shw.start.x) || (shw.start.x > -5.5 && shw.start.x < -4.5) ||
           std::isnan(shw.len) || shw.len <= 0. ) continue;

      // if it meets this then we're not going to cut on it...
      if ( std::isnan(shw.plane[2].energy) || std::isinf(shw.plane[2].energy) || shw.plane[2].energy <= 0.04 ) continue;

      // and... if it meets then then we're not going to cut on it...
      if ( std::isnan(shw.conversion_gap) || std::isinf(shw.conversion_gap) || shw.conversion_gap <= 5. ) continue;

      // if we got here, then it should be the case that the fit seems valid and:
      // shwE > 0.040 GeV
      // trackScore < 0.45 (technically <= 0.45)
      // conversionGap > 5. cm
      return false;
    }

    // guess we're not cutting anything
    return true;
  });


// The combined cuts
const Cut kNuMISelection_1muNp0pi = kRFiducialNew && kNotClearCosmic && /*Preselection*/
                                    kPTrackNew && kProtonTrack && kRecoProtonPInKinematicTreshold && /*Mu, P candidates*/
                                    kAllPrimaryHadronsContained && kNoSecondPrimaryMuonlikeTracks && /*Esp. chg pi rejection*/
                                    kCutPhotons;                                                     /*Esp. pi0 rejection*/

// Version from plot approvals (USE WITH "Legacy_TrkAndShw.h" as the one determining IsTracklikeTrack !!)
const Cut kNuMISelection_1muNp0pi_Legacy = kRFiducialNew && kNotClearCosmic && /*Preselection*/
                                           kPTrackNew && kProtonTrack && /*Mu, P candidates*/
                                           kAllPrimaryHadronsContained && kNoSecondPrimaryMuonlikeTracks; /*Esp. chg pi rejection*/

////////////////////////////////////////////////////////////////
// Additional potentially useful selection cuts if one is interested in contained muon candidates:
////////////////////////////////////////////////////////////////

const Cut kIsMuonCandidateContained([](const caf::SRSliceProxy* slc) {
    int muonCandidate = kPTrackIndNew(slc);

    if ( muonCandidate < 0 ) return false;
    else {
      auto const& trk = slc->reco.pfp.at(kPTrackIndNew(slc)).trk;
      const bool Contained = isContainedVol(trk.end.x,trk.end.y,trk.end.z);
      return Contained;
    }

    return false;
  });

const Cut kRejectSplitMuons([](const caf::SRSliceProxy* slc) {
    int primaryInd = kPTrackIndNew(slc);
    if ( primaryInd < 0 ) return false;
    unsigned int idxPrim = (unsigned int)primaryInd;

    unsigned int idxTrk = 0;
    while ( IsValidTrkIdx(slc, idxTrk) ) {
      int thisIdxInt = idxTrk;
      if ( thisIdxInt == primaryInd ) {
        idxTrk+=1;
        continue; // skip the particle which is the muon candidate!
      }
      auto const& trk = slc->reco.pfp.at(idxTrk).trk;
      unsigned int thisIdx = idxTrk;
      idxTrk+=1;

      if ( std::isnan(trk.start.x) || std::isnan(trk.len) || trk.len <= 0. ) continue;
      if ( std::isnan(slc->vertex.x) || std::isnan(slc->vertex.y) || std::isnan(slc->vertex.z) ) return false;
      const float Atslc = std::hypot(slc->vertex.x - trk.start.x,
                                     slc->vertex.y - trk.start.y,
                                     slc->vertex.z - trk.start.z);
      const bool isSplitCandidate = (Atslc > 10.);

      if ( !isSplitCandidate || trk.calo[2].nhit < 5 ) continue;
      const float Chi2Proton = trk.chi2pid[2].chi2_proton;
      const float Chi2Muon = trk.chi2pid[2].chi2_muon;
      if ( !(Chi2Proton > 60. && Chi2Muon < 30.) ) continue;

      if ( trk.len > 10. ) return false;
    }

    return true;
  });

////////////////////////////////////////////////////////////////
// X-SEC related Vars:
////////////////////////////////////////////////////////////////

const Var kRecoMuonPNew([](const caf::SRSliceProxy* slc) -> float {
    float p(-5.f);

    if ( kPTrackIndNew(slc) >= 0 )
    {
      auto const& trk = slc->reco.pfp.at(kPTrackIndNew(slc)).trk;
      const bool Contained = isContainedVol(trk.end.x,trk.end.y,trk.end.z);
      if(Contained) p = trk.rangeP.p_muon;
      else p = trk.mcsP.fwdP_muon;
    }
    return p;
  });

// True version
const Var kTrueMuonPNew([](const caf::SRSliceProxy* slc) -> float {
    float p(-5.f);
    if ( slc->truth.index < 0 ) return p;

    for ( auto const& prim : slc->truth.prim ) {
      if ( prim.start_process != 0 ) continue; // to handle fact that primaries now include more particles.
      double momentum = sqrt( (prim.startp.x*prim.startp.x) + (prim.startp.y*prim.startp.y) + (prim.startp.z*prim.startp.z) );
      if ( abs(prim.pdg) == 13 && momentum > p ) {
        p = momentum;
      }
    }

    return p;
  });

const Var kRecoProtonPNew([](const caf::SRSliceProxy* slc) -> float {
    float p(-5.f);

    if ( kProtonTrackIndNew(slc) >= 0 )
    {
      auto const& trk = slc->reco.pfp.at(kProtonTrackIndNew(slc)).trk;
      p = trk.rangeP.p_proton;
    }
    return p;
  });

// True version
const Var kTrueProtonPNew([](const caf::SRSliceProxy* slc) -> float {
    float p(-5.f);
    if ( slc->truth.index < 0 ) return p;

    for ( auto const& prim : slc->truth.prim ) {
      if ( prim.start_process != 0 ) continue; // to handle fact that primaries now include more particles.
      double momentum = sqrt( (prim.startp.x*prim.startp.x) + (prim.startp.y*prim.startp.y) + (prim.startp.z*prim.startp.z) );
      if ( abs(prim.pdg) == 2212 && prim.contained && momentum > p ) {
        p = momentum;
      }
    }

    return p;
  });

const Var kRecoMuonCosThNuMI([](const caf::SRSliceProxy* slc) -> float {
    float costh(-5.f);

    if ( kPTrackIndNew(slc) >= 0 ) {
      auto const& mutrk = slc->reco.pfp.at(kPTrackIndNew(slc)).trk;

      TVector3 muDir(mutrk.dir.x, mutrk.dir.y, mutrk.dir.z);
      muDir = muDir.Unit();

      costh = TMath::Cos( muDir.Angle(rFromNuMI) );
    }

    return costh;
  });

// True version
const Var kTrueMuonCosThNuMI([](const caf::SRSliceProxy* slc) -> float {
    float p(-5.f);
    float cosThNuMI(-5.f);
    if ( slc->truth.index < 0 ) return cosThNuMI;

    for ( auto const& prim : slc->truth.prim ) {
      if ( prim.start_process != 0 ) continue; // to handle fact that primaries now include more particles.
      double momentum = sqrt( (prim.startp.x*prim.startp.x) + (prim.startp.y*prim.startp.y) + (prim.startp.z*prim.startp.z) );
      if ( abs(prim.pdg) == 13 && momentum > p ) {
        p = momentum;

        TVector3 muDir(prim.end.x - prim.start.x, prim.end.y - prim.start.y, prim.end.z - prim.start.z);
        muDir = muDir.Unit();
        cosThNuMI = TMath::Cos( muDir.Angle(rFromNuMI) );
      }
    }

    return cosThNuMI;
  });

const Var kRecoMuonProtonCosTh([](const caf::SRSliceProxy* slc) -> float {
    float costh(-5.f);

    if ( kPTrackIndNew(slc) >= 0 && kProtonTrackIndNew(slc) >= 0 ) {
      auto const& mutrk = slc->reco.pfp.at(kPTrackIndNew(slc)).trk;
      auto const& ptrk = slc->reco.pfp.at(kProtonTrackIndNew(slc)).trk;

      TVector3 muDir(mutrk.dir.x, mutrk.dir.y, mutrk.dir.z);
      muDir = muDir.Unit();
      TVector3 pDir(ptrk.dir.x, ptrk.dir.y, ptrk.dir.z);
      pDir = pDir.Unit();

      costh = TMath::Cos( muDir.Angle(pDir) );
    }

    return costh;
  });

// True version
const Var kTrueMuonProtonCosTh([](const caf::SRSliceProxy* slc) -> float {
    float pM(-5.f);
    float pP(-5.f);
    float pX(0.f);
    float pY(0.f);
    float pZ(0.f);
    float mX(0.f);
    float mY(0.f);
    float mZ(0.f);
    if ( slc->truth.index < 0 ) return -5.f;

    for ( auto const& prim : slc->truth.prim ) {
      if ( prim.start_process != 0 ) continue; // to handle fact that primaries now include more particles.
      double momentum = sqrt( (prim.startp.x*prim.startp.x) + (prim.startp.y*prim.startp.y) + (prim.startp.z*prim.startp.z) );
      if ( abs(prim.pdg) == 13 && momentum > pM ) {
        pM = momentum;

        TVector3 muDir(prim.end.x - prim.start.x, prim.end.y - prim.start.y, prim.end.z - prim.start.z);
        muDir = muDir.Unit();
        mX = muDir.X();
        mY = muDir.Y();
        mZ = muDir.Z();
      }
      if ( abs(prim.pdg) == 2212 && prim.contained && momentum > pP ) {
        pP = momentum;

        TVector3 pDir(prim.end.x - prim.start.x, prim.end.y - prim.start.y, prim.end.z - prim.start.z);
        pDir = pDir.Unit();
        pX = pDir.X();
        pY = pDir.Y();
        pZ = pDir.Z();
      }
    }

    TVector3 muonDir(mX, mY, mZ);
    muonDir = muonDir.Unit();
    TVector3 protonDir(pX, pY, pZ);
    protonDir = protonDir.Unit();
    return TMath::Cos( muonDir.Angle(protonDir) );
  });
