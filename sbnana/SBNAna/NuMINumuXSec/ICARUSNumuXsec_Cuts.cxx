#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_Cuts.h"

using namespace ana;
using namespace std;

namespace ICARUSNumuXsec{

  //==== Neutrino flavor
  const Cut cutIsNuMu([](const caf::SRSliceProxy* slc) {
    return (kIsNuSlice(slc) && ( slc->truth.pdg == 14 || slc->truth.pdg == -14 ));
    });
  const Cut cutIsNuE([](const caf::SRSliceProxy* slc) {
    return (kIsNuSlice(slc) && ( slc->truth.pdg == 12 || slc->truth.pdg == -12 ));
    });

  //==== GENIE Interaction code
  const Cut cutIsQE([](const caf::SRSliceProxy* slc) {
    return varGENIEIntCode(slc)==0;
    });
  const Cut cutIsRes([](const caf::SRSliceProxy* slc) {
    return varGENIEIntCode(slc)==1;
    });
  const Cut cutIsDIS([](const caf::SRSliceProxy* slc) {
    return varGENIEIntCode(slc)==2;
    });
  const Cut cutIsCoh([](const caf::SRSliceProxy* slc) {
    return varGENIEIntCode(slc)==3;
    });
  const Cut cutIsCohElastic([](const caf::SRSliceProxy* slc) {
    return varGENIEIntCode(slc)==4;
    });
  const Cut cutIsElectronScattering([](const caf::SRSliceProxy* slc) {
    return varGENIEIntCode(slc)==5;
    });
  const Cut cutIsIMDAnnihilation([](const caf::SRSliceProxy* slc) {
    return varGENIEIntCode(slc)==6;
    });
  const Cut cutIsInverseBetaDecay([](const caf::SRSliceProxy* slc) {
    return varGENIEIntCode(slc)==7;
    });
  const Cut cutIsGlashowResonance([](const caf::SRSliceProxy* slc) {
    return varGENIEIntCode(slc)==8;
    });
  const Cut cutIsAMNuGamma([](const caf::SRSliceProxy* slc) {
    return varGENIEIntCode(slc)==9;
    });
  const Cut cutIsMEC([](const caf::SRSliceProxy* slc) {
    return varGENIEIntCode(slc)==10;
    });
  const Cut cutIsDiffractive([](const caf::SRSliceProxy* slc) {
    return varGENIEIntCode(slc)==11;
    });
  const Cut cutIsEM([](const caf::SRSliceProxy* slc) {
    return varGENIEIntCode(slc)==12;
    });
  const Cut cutIsWeakMix([](const caf::SRSliceProxy* slc) {
    return varGENIEIntCode(slc)==13;
    });
  const Cut cutIsUnknownInteractionType1([](const caf::SRSliceProxy* slc) {
    return varGENIEIntCode(slc)==-1;
    });
  const Cut cutIsUnknownInteractionType2([](const caf::SRSliceProxy* slc) {
    return varGENIEIntCode(slc)>13;
    });
  const Cut cutIsUnknownInteractionType3([](const caf::SRSliceProxy* slc) {
    return varGENIEIntCode(slc)<-1;
    });

  //==== CC vs NC

  const Cut cutIsCC([](const caf::SRSliceProxy* slc) {
      return ( kIsNuSlice(slc) && slc->truth.iscc );
    });
  const Cut cutIsNC([](const caf::SRSliceProxy* slc) {
      return ( kIsNuSlice(slc) && slc->truth.isnc );
    });

  //==== Number of truth particles

  const Cut cutHasTruthMuon([](const caf::SRSliceProxy* slc) {
      return ( varMuonTruthIndex(slc)>=0 );
    });
  const Cut cutHasTruthMuonHasRecoTrack([](const caf::SRSliceProxy* slc) {
    return ( varTruthMuonMatchedTrackIndex(slc)>=0 );
    });

  const Cut cutHasTruthProton([](const caf::SRSliceProxy* slc) {
      return ( varProtonTruthIndex(slc)>=0 );
    });
  const Cut cutHasTruthProtonHasRecoTrack([](const caf::SRSliceProxy* slc) {
    return ( varTruthProtonMatchedTrackIndex(slc)>=0 );
    });
  const Cut cutHasTruthChargedPion([](const caf::SRSliceProxy* slc) {
      return ( varChargedPionTruthIndex(slc)>=0 );
    });
  const Cut cutHasTruthChargedPionHasRecoTrack([](const caf::SRSliceProxy* slc) {
    return ( varTruthChargedPionMatchedTrackIndex(slc)>=0 );
    });

  const Cut cutTruthNoPiZero = (varTruthNPiZero==0);
  const Cut cutTruthNoChargedPion = (varTruthNChargedPion==0);
  const Cut cutTruthNoNeutron = (varTruthNNeutron==0);

  //==== NuMu-CC categories

  const Cut cutIsNuMuCC([](const caf::SRSliceProxy* slc) {
      return ( cutIsNuMu(slc) && cutIsCC(slc) );
    });

  const Cut cutIsNuMuCCQE([](const caf::SRSliceProxy* slc) {
      return ( cutIsNuMuCC(slc) && cutIsQE(slc) );
    });
  const Cut cutIsNuMuCCRes([](const caf::SRSliceProxy* slc) {
      return ( cutIsNuMuCC(slc) && cutIsRes(slc) );
    });
  const Cut cutIsNuMuCCMEC([](const caf::SRSliceProxy* slc) {
      return ( cutIsNuMuCC(slc) && cutIsMEC(slc) );
    });
  const Cut cutIsNuMuCCDIS([](const caf::SRSliceProxy* slc) {
      return ( cutIsNuMuCC(slc) && cutIsDIS(slc) );
    });

  //==== NuMu-NC categories

  const Cut cutIsNuMuNC([](const caf::SRSliceProxy* slc) {
      return ( kIsNuSlice(slc) && slc->truth.isnc && ( slc->truth.pdg == 14 || slc->truth.pdg == -14 ) );
    });

  //==== NuE

  const Cut cutIsNuECC([](const caf::SRSliceProxy* slc) {
      return ( cutIsNuE(slc) && cutIsCC(slc) );
    });

  //==== Truth kinematic cuts

  const Cut cutTFiducial([](const caf::SRSliceProxy* slc) {
    if( !isnan(slc->truth.position.x) ) return fv.isContained(slc->truth.position.x, slc->truth.position.y, slc->truth.position.z);
    else return false;
  });

  const Cut cutTruthMuonTCut([](const caf::SRSliceProxy* slc) {

    double truthMuonT = varMuonTruthT(slc);
    static const double truthMuonTThreshold_Contained = 0.150;
    static const double truthMuonTThreshold_Exiting = 0.250;
    bool passMuonTCut = cutTruthMuonContained(slc) ? (truthMuonT > truthMuonTThreshold_Contained) : (truthMuonT > truthMuonTThreshold_Exiting);

    return passMuonTCut;

    });

  const Cut cutTruthProtonTCut([](const caf::SRSliceProxy* slc) {

    double truthProtonT = varProtonTruthT(slc);
    static const double truthProtonTThreshold = 0.050;
    bool passProtonTCut = (truthProtonT > truthProtonTThreshold);

    return passProtonTCut;

    });

  //==== NuMu+P signal signal definition using truth

  const Cut cutIsNuMuCCSignalDef([](const caf::SRSliceProxy* slc) {

    //==== 1. NuMu-CC
    //==== 2. Truth fiducial
    //==== 3. Muon T cut
    //==== 4. Proton T cut

    return ( cutIsNuMuCC(slc) && cutTFiducial(slc) && cutTruthMuonTCut(slc) && cutTruthProtonTCut(slc) );

    });

  //==== FV

  const Cut cutRFiducial([](const caf::SRSliceProxy* slc) {

    if( !isnan(slc->vertex.x) ) return fv.isContained(slc->vertex.x, slc->vertex.y, slc->vertex.z);
    else return false;

  });

  //==== FMScore

  const Cut cutFMScore([](const caf::SRSliceProxy* slc) {
/*
    bool passCut = !isnan(slc->fmatch.score) && slc->fmatch.score < 12.0 && slc->fmatch.score >= 0;
    std::cout << "slc->fmatch.score = " << slc->fmatch.score;
    if(passCut) std::cout << " -> Pass" << std::endl;
    else std::cout << " -> Fail" << std::endl;
*/
    return ( !isnan(slc->fmatch.score) && slc->fmatch.score < 12.0 && slc->fmatch.score >= 0 );
  });

  const Cut kSlcFlashMatchDataCut([](const caf::SRSliceProxy *slc)
				{
				  return (kSlcHasFlashMatch(slc) && slc->fmatch.score>0. && slc->fmatch.score<12.);
				});

  const Cut cutFMTime([](const caf::SRSliceProxy* slc) {
    //return ( !isnan(slc->fmatch.time) && slc->fmatch.time>=0 && slc->fmatch.time<=16 );
    return ( !isnan(slc->fmatch.time) && slc->fmatch.time>=-0.2 && slc->fmatch.time<=9.9 );
  });

  //==== NuScore

  const Cut cutNuScore([](const caf::SRSliceProxy* slc) {
    return ( !isnan(slc->nu_score) && slc->nu_score > 0.4 );
  });

  const Cut cutSliceCRLongestTrackDirY([](const caf::SRSliceProxy* slc) {
    if(isnan(slc->nuid.crlongtrkdiry)){
      return false;
    }
    else{
      return slc->nuid.crlongtrkdiry>-0.7;
    }
  });

  const Cut cutSliceNuVertexYTop([](const caf::SRSliceProxy* slc) {
    if(isnan(slc->nuid.nuvtxy)){
      return false;
    }
    else{
      return slc->nuid.nuvtxy>80.;
    }
  });


  //==== 220811 Vertex test
  const Cut cutVertexYPos([](const caf::SRSliceProxy* slc) {
    double vtxy = varVertexRecoY(slc);
    return ( vtxy > 0.1 && vtxy < 60. );
  });
  const Cut cutVertexYNeg([](const caf::SRSliceProxy* slc) {
    double vtxy = varVertexRecoY(slc);
    return ( vtxy < -0.1 && vtxy > -60. );
  });

  //==== Muon related

  const Cut cutHasMuon([](const caf::SRSliceProxy* slc) {
    return (varMuonTrackInd(slc) >= 0);
  });

  const Cut cutMuonContained([](const caf::SRSliceProxy* slc) {

    if( varMuonTrackInd(slc) >= 0 ){
      auto const& trk = slc->reco.trk.at(varMuonTrackInd(slc));
      return fv_track.isContained(trk.end.x, trk.end.y, trk.end.z);
    }
    else{
      return false;
    }

  });

  const Cut cutTruthMuonContained([](const caf::SRSliceProxy* slc) {

    double max_E(-999);
    int muonidx(-1);
    for(std::size_t i(0); i < slc->truth.prim.size(); ++i){
      if( abs(slc->truth.prim.at(i).pdg)== 13 ){
        if(isnan(slc->truth.prim.at(i).genE)) continue;
        double this_E = slc->truth.prim.at(i).genE;
        if(this_E>max_E){
          max_E = this_E;
          muonidx = i;
        }
      }
    }
    if(muonidx>=0){
      return bool(slc->truth.prim.at(muonidx).contained);
    }
    else{
      return false;
    }

  });

  const Cut cutRecoMuonTruthContained([](const caf::SRSliceProxy* slc) {

    if( varMuonTrackInd(slc) >= 0 ){
      auto const& trk = slc->reco.trk.at(varMuonTrackInd(slc));
      return bool(trk.truth.p.contained);
    }
    else{
      return false;
    }

  });

  const Cut cutTruthMuonMatchedTrackContained([](const caf::SRSliceProxy* slc) {
    if( varTruthMuonMatchedTrackIndex(slc) >= 0 ){
      return (varTruthMuonMatchedTrackContainedness(slc)==1);
    }
    else{
      return false;
    }

  });

  const Cut cutMuonMatchedToMuon([](const caf::SRSliceProxy* slc) {
    return ( abs(varMuonBestmatchPDG(slc)) == 13 );
  });

  const Cut cutMuonMatchedToProton([](const caf::SRSliceProxy* slc) {
    return ( abs(varMuonBestmatchPDG(slc)) == 2212 );
  });

  //==== Michel study

  const Cut cutTruthMuonMatchedTrackHasStitchedTrack([](const caf::SRSliceProxy* slc) {
    return ( varTruthMuonMatchedTrackStitchedTrackIndex(slc) >= 0 );
  });

  const Cut cutTruthMuonMatchedTrackStitchedTrackMatched([](const caf::SRSliceProxy* slc) {
    int muonTrackIndex = varTruthMuonMatchedTrackIndex(slc);
    int muonStitchedTrackIndex = varTruthMuonMatchedTrackStitchedTrackIndex(slc);
    if(muonStitchedTrackIndex>=0){

      auto const& motherTrack = slc->reco.trk.at(muonTrackIndex);
      int motherTrack_G4ID = motherTrack.truth.p.G4ID;

      auto const& stTrack = slc->reco.trk.at(muonStitchedTrackIndex);
      int stTrack_parent_G4ID = stTrack.truth.p.parent;

      if(motherTrack_G4ID==stTrack_parent_G4ID){
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

  const Cut cutTruthMuonMatchedTrackHasStitchedShower([](const caf::SRSliceProxy* slc) {
    return ( varTruthMuonMatchedTrackStitchedShowerIndex(slc) >= 0 );
  });

  const Cut cutTruthMuonMatchedTrackStitchedShowerMatched([](const caf::SRSliceProxy* slc) {
    int muonTrackIndex = varTruthMuonMatchedTrackIndex(slc);
    int muonStitchedShowerIndex = varTruthMuonMatchedTrackStitchedShowerIndex(slc);
    if(muonStitchedShowerIndex>=0){

      auto const& motherTrack = slc->reco.trk.at(muonTrackIndex);
      int motherTrack_G4ID = motherTrack.truth.p.G4ID;
  
      auto const& stShower = slc->reco.shw.at(muonStitchedShowerIndex);
      int stShower_parent_G4ID = stShower.truth.p.parent;
  
      if(motherTrack_G4ID==stShower_parent_G4ID){
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

  //==== Michel tagging for muon

  const Cut cutTruthMuonMatchedTrackIsMichelTagged([](const caf::SRSliceProxy* slc) {

    int muonTrackIndex = varTruthMuonMatchedTrackIndex(slc);
    if( muonTrackIndex >= 0 ){
      int shwind = GetMichelShowerIndex(slc, muonTrackIndex);
      return (shwind>=0);
    }
    else{
      return false;
    }

  });

 //==== Michel tagging for muon

  const Cut cutMuonIsMichelTagged([](const caf::SRSliceProxy* slc) {

    int muonTrackIndex = varMuonTrackInd(slc);
    if( muonTrackIndex >= 0 ){
      int shwind = GetMichelShowerIndex(slc, muonTrackIndex);
      return (shwind>=0);
    }
    else{
      return false;
    }

  });

  //==== Pion tagging for muon

  const Cut cutMuonIsPionTagged([](const caf::SRSliceProxy* slc) {

    int muonTrackIndex = varMuonTrackInd(slc);
    if( muonTrackIndex >= 0 ){
      auto const& trk = slc->reco.trk.at(muonTrackIndex);
      if( fv_track.isContained(trk.end.x, trk.end.y, trk.end.z) ){
        return IsPionTagged(slc, muonTrackIndex);
      }
      else{
        return false;
      }
    }
    else{
      return false;
    }

  });


  //==== Proton related

  const Cut cutHasProton([](const caf::SRSliceProxy* slc) {
    return (varProtonTrackInd(slc) >= 0);
  });

  const Cut cutProtonContained([](const caf::SRSliceProxy* slc) {

    if( varProtonTrackInd(slc) >= 0 ){
      auto const& trk = slc->reco.trk.at(varProtonTrackInd(slc));
      return fv_track.isContained(trk.end.x, trk.end.y, trk.end.z);
    }
    else{
      return false;
    }

  });

  const Cut cutTruthProtonContained([](const caf::SRSliceProxy* slc) {

    double max_E(-999);
    int muonidx(-1);
    for(std::size_t i(0); i < slc->truth.prim.size(); ++i){
      if( abs(slc->truth.prim.at(i).pdg)==2212 ){
        if(isnan(slc->truth.prim.at(i).genE)) continue;
        double this_E = slc->truth.prim.at(i).genE;
        if(this_E>max_E){
          max_E = this_E;
          muonidx = i;
        }
      }
    }
    if(muonidx>=0){
      return bool(slc->truth.prim.at(muonidx).contained);
    }
    else{
      return false;
    }

  });

  const Cut cutProtonHighMomentum([](const caf::SRSliceProxy* slc) {

    if( varProtonTrackInd(slc) >= 0 ){
      return (varProtonRecoP(slc)>1.);
    }
    else{
      return false;
    }

  });

  const Cut cutProtonMatchedToMuon([](const caf::SRSliceProxy* slc) {
    return ( abs(varProtonBestmatchPDG(slc)) == 13 );
  });

  const Cut cutProtonMatchedToProton([](const caf::SRSliceProxy* slc) {
    return ( abs(varProtonBestmatchPDG(slc)) == 2212 );
  });

  const Cut cutTruthProtonLargePResidualFraction([](const caf::SRSliceProxy* slc) {
    return ( varTruthProtonMatchedTrackRangePResidualFraction(slc) >= 0.2 );
  });

  //==== Charged pion related

  const Cut cutTruthChargedPionContained([](const caf::SRSliceProxy* slc) {

    double max_E(-999);
    int muonidx(-1);
    for(std::size_t i(0); i < slc->truth.prim.size(); ++i){
      if( abs(slc->truth.prim.at(i).pdg)== 211 ){
        if(isnan(slc->truth.prim.at(i).genE)) continue;
        double this_E = slc->truth.prim.at(i).genE;
        if(this_E>max_E){
          max_E = this_E;
          muonidx = i;
        }
      }
    }
    if(muonidx>=0){
      return bool(slc->truth.prim.at(muonidx).contained);
    }
    else{
      return false;
    }

  });

  const Cut cutTruthChargedPionMatchedTrackContained([](const caf::SRSliceProxy* slc) {
    if( varTruthChargedPionMatchedTrackIndex(slc) >= 0 ){
      return (varTruthChargedPionMatchedTrackContainedness(slc)==1);
    }
    else{
      return false;
    }

  });

  //==== Michel study

  const Cut cutTruthChargedPionMatchedTrackHasStitchedTrack([](const caf::SRSliceProxy* slc) {
    return ( varTruthChargedPionMatchedTrackStitchedTrackIndex(slc) >= 0 );
  });

  const Cut cutTruthChargedPionMatchedTrackStitchedTrackMatched([](const caf::SRSliceProxy* slc) {
    int cpionTrackIndex = varTruthChargedPionMatchedTrackIndex(slc);
    int cpionStitchedTrackIndex = varTruthChargedPionMatchedTrackStitchedTrackIndex(slc);
    if(cpionStitchedTrackIndex>=0){

      auto const& motherTrack = slc->reco.trk.at(cpionTrackIndex);
      int motherTrack_G4ID = motherTrack.truth.p.G4ID;

      auto const& stTrack = slc->reco.trk.at(cpionStitchedTrackIndex);
      int stTrack_parent_G4ID = stTrack.truth.p.parent;

      if(motherTrack_G4ID==stTrack_parent_G4ID){
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

  const Cut cutTruthChargedPionMatchedTrackHasStitchedShower([](const caf::SRSliceProxy* slc) {
    return ( varTruthChargedPionMatchedTrackStitchedShowerIndex(slc) >= 0 );
  });

  const Cut cutTruthChargedPionMatchedTrackStitchedShowerMatched([](const caf::SRSliceProxy* slc) {
    int cpionTrackIndex = varTruthChargedPionMatchedTrackIndex(slc);
    int cpionStitchedShowerIndex = varTruthChargedPionMatchedTrackStitchedShowerIndex(slc);
    if(cpionStitchedShowerIndex>=0){

      auto const& motherTrack = slc->reco.trk.at(cpionTrackIndex);
      int motherTrack_G4ID = motherTrack.truth.p.G4ID;
  
      auto const& stShower = slc->reco.shw.at(cpionStitchedShowerIndex);
      int stShower_parent_G4ID = stShower.truth.p.parent;
  
      if(motherTrack_G4ID==stShower_parent_G4ID){
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

  //==== CRT 

  const SpillCut spillcutFDTopCRTHitVeto(
      [](const caf::SRSpillProxy* sr){
        for (auto const& crtHit: sr->crt_hits){
          if(crtHit.plane>=30 && crtHit.plane<=34){
            if (crtHit.t0 > 0.1 && crtHit.t0 < 9.0 && crtHit.pe > 100) // CRTNuMIWindow
              return false;
          }
        }
        return true;
      }
      );

  const SpillCut spillcutFDSideCRTHitVeto(
      [](const caf::SRSpillProxy* sr){
        for (auto const& crtHit: sr->crt_hits){
          if(crtHit.plane>=40 && crtHit.plane<=47){
            if (crtHit.t0 > 0.1 && crtHit.t0 < 9.0 && crtHit.pe > 100) // CRTNuMIWindow
              return false;
          }
        }
        return true;
      }
      );

  const SpillCut spillcutFDTopCRTHitVetoTestMatching(
      [](const caf::SRSpillProxy* sr){

        vector<TVector3> muonTrackStarts;
        for(auto const& slc : sr->slc){
          if(!cutRFiducial(&slc)) continue;
          if(!kNotClearCosmic(&slc)) continue;
          if(!cutNuScore(&slc)) continue;
          if(!cutFMScore(&slc)) continue;
          if(!cutHasMuon(&slc)) continue;

          vector<double> muIndices = varAllMuonTrackIndices(&slc);

          for(auto const& muInd: muIndices){
            auto const& trk = slc.reco.trk.at(muInd);
            muonTrackStarts.emplace_back(trk.start.x, trk.start.y, trk.start.z);
          }

        }

        if(muonTrackStarts.size()==0) return true;

        //std::cout << "[spillcutFDTopCRTHitVetoTestMatching] Muon track found, and looping over CRTHits" << std::endl;
        //std::cout << "[spillcutFDTopCRTHitVetoTestMatching] muonTrackStarts.size() = " << muonTrackStarts.size() << std::endl;

        for (auto const& crtHit: sr->crt_hits){

          if(crtHit.plane>=30 && crtHit.plane<=34){

            if(crtHit.t0 > 0.1 && crtHit.t0 < 9.0 && crtHit.pe > 100){
/*
              std::cout << "[spillcutFDTopCRTHitVetoTestMatching] In-time CRTHit found:" << std::endl;
              std::cout << "[spillcutFDTopCRTHitVetoTestMatching] x = " << crtHit.position.x << std::endl;
              std::cout << "[spillcutFDTopCRTHitVetoTestMatching] y = " << crtHit.position.y << std::endl;
              std::cout << "[spillcutFDTopCRTHitVetoTestMatching] z = " << crtHit.position.z << std::endl;
*/
              for(auto& muonTrackStart: muonTrackStarts){
/*
                std::cout << "[spillcutFDTopCRTHitVetoTestMatching] Muon track information" << std::endl;
                std::cout << "[spillcutFDTopCRTHitVetoTestMatching] x = " << muonTrackStart.X() << std::endl;
                std::cout << "[spillcutFDTopCRTHitVetoTestMatching] y = " << muonTrackStart.Y() << std::endl;
                std::cout << "[spillcutFDTopCRTHitVetoTestMatching] z = " << muonTrackStart.Z() << std::endl;
*/
                //==== Cryo using sign of X
                if( crtHit.position.x * muonTrackStart.X() > 0 ) return false;
                //==== TPC using sign of Z
                if( crtHit.position.z * muonTrackStart.Z() > 0 ) return false;

              } // END Loop vertes 

            } // END If intime-ed CRT hit exist

          } // END If top crt

        } // END Loop CRTHit


        return true;

      }
      );


  const SpillCut spillcutFDSideCRTHitVetoTestMatching(
      [](const caf::SRSpillProxy* sr){

        vector<TVector3> muonTrackStarts;
        for(auto const& slc : sr->slc){
          if(!cutRFiducial(&slc)) continue;
          if(!kNotClearCosmic(&slc)) continue;
          if(!cutNuScore(&slc)) continue;
          if(!cutFMScore(&slc)) continue;
          if(!cutHasMuon(&slc)) continue;

          vector<double> muIndices = varAllMuonTrackIndices(&slc);

          for(auto const& muInd: muIndices){
            auto const& trk = slc.reco.trk.at(muInd);
            muonTrackStarts.emplace_back(trk.start.x, trk.start.y, trk.start.z);
          }
        }

        if(muonTrackStarts.size()==0) return true;

        for (auto const& crtHit: sr->crt_hits){

          if(crtHit.plane>=40 && crtHit.plane<=47){

            if(crtHit.t0 > 0.1 && crtHit.t0 < 9.0 && crtHit.pe > 100){

              for(auto& muonTrackStart: muonTrackStarts){
                //==== Cryo using sign of X
                if( crtHit.position.x * muonTrackStart.X() > 0 ) return false;
                //==== TPC using sign of Z
                if( crtHit.position.z * muonTrackStart.Z() > 0 ) return false;
              } // END Loop vertes 

            } // END If intime-ed CRT hit exist

          } // END If top crt

        } // END Loop CRTHit


        return true;

      }
      );

  const Cut cutNominal_ContainedMuon([](const caf::SRSliceProxy* slc) {
    return ( cutRFiducial(slc) && 
             cutFMScore(slc) && 
             cutFMTime(slc) && 
             cutSliceCRLongestTrackDirY(slc) && 
             cutHasMuon(slc) && 
             cutMuonContained(slc)
           );
    });

  const Cut cutNominal_ExitingMuon([](const caf::SRSliceProxy* slc) {
    return ( cutRFiducial(slc) && 
             cutFMScore(slc) && 
             cutFMTime(slc) && 
             cutSliceCRLongestTrackDirY(slc) && 
             cutHasMuon(slc) && 
             !cutMuonContained(slc) 
           );
    });

  const Cut cutMuonProtonCosineTheta([](const caf::SRSliceProxy* slc) {
    return ( varMuonProtonCosineTheta(slc)>-0.9 );
    });



}

