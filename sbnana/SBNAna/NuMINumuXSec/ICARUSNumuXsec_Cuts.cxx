#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_Cuts.h"

using namespace ana;
using namespace std;

namespace ICARUSNumuXsec{

  const Cut cutIsNuMu([](const caf::SRSliceProxy* slc) {
    return (kIsNuSlice(slc) && ( slc->truth.pdg == 14 || slc->truth.pdg == -14 ));
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

  const Cut cutIsCC([](const caf::SRSliceProxy* slc) {
      return ( kIsNuSlice(slc) && slc->truth.iscc );
    });
  const Cut cutIsNC([](const caf::SRSliceProxy* slc) {
      return ( kIsNuSlice(slc) && slc->truth.isnc );
    });

  //==== Nue

  const Cut cutIsNuECC([](const caf::SRSliceProxy* slc) {
      return ( kIsNuSlice(slc) && slc->truth.iscc && ( slc->truth.pdg == 12 || slc->truth.pdg == -12 ) );
    });

  //==== NC

  const Cut cutIsNuMuNC([](const caf::SRSliceProxy* slc) {
      return ( kIsNuSlice(slc) && slc->truth.isnc && ( slc->truth.pdg == 14 || slc->truth.pdg == -14 ) );
    });

  const Cut cutTruthNoPiZero = (varTruthNPiZero==0);
  const Cut cutTruthNoChargedPion = (varTruthNChargedPion==0);
  const Cut cutTruthNoNeutron = (varTruthNNeutron==0);

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
    passProtonTCut = true; //TODO

    return passProtonTCut;

    });

  const Cut cutIsNuMuCCSignalDef([](const caf::SRSliceProxy* slc) {

    //==== 1. NuMu-CC
    //==== 2. Truth fiducial
    //==== 3. Muon T cut
    //==== 4. Proton T cut

    return ( kIsNuMuCC(slc) && cutTFiducial(slc) && cutTruthMuonTCut(slc) && cutTruthProtonTCut(slc) );

    });

  //==== FV

  const Cut cutRFiducial([](const caf::SRSliceProxy* slc) {

    if( !isnan(slc->vertex.x) ) return fv.isContained(slc->vertex.x, slc->vertex.y, slc->vertex.z);
    else return false;

  });

  //==== FMScore

  const Cut cutFMScore([](const caf::SRSliceProxy* slc) {
      //return ( !isnan(slc->fmatch.score) && slc->fmatch.score < 6.0 );
      //return ( !isnan(slc->fmatch.score) && slc->fmatch.score < 7.0 && (slc->fmatch.time>-0.2 && slc->fmatch.time<9.9) );
      return ( !isnan(slc->fmatch.score) && slc->fmatch.score < 12.0 && slc->fmatch.score >= 0 );
      //return true;
    });

  //==== NuScore

  const Cut cutNuScore([](const caf::SRSliceProxy* slc) {
    //return ( !isnan(slc->nu_score) && slc->nu_score > 0.4 );
    //return ( !isnan(slc->nu_score) && slc->nu_score > 0.2 );
    return true;
  });

  //==== Muon related

  const Cut cutHasMuon([](const caf::SRSliceProxy* slc) {
    return (varMuonTrackInd(slc) >= 0);
  });

  const Cut cutMuonContained([](const caf::SRSliceProxy* slc) {

    if( varMuonTrackInd(slc) >= 0 ){
      auto const& trk = slc->reco.trk.at(varMuonTrackInd(slc));
      return fv.isContained(trk.end.x, trk.end.y, trk.end.z);
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

  const Cut cutZeroMomentum([](const caf::SRSliceProxy* slc) {

    if( varMuonTrackInd(slc) >= 0 ){
      return (varMuonRecoP(slc)==0.);
/*
      //TODO TEST
      if(varMuonRecoP(slc)==0.){
        auto const& trk = slc->reco.trk.at(varMuonTrackInd(slc));
        int trackStartCryo = fv.containedCryo(trk.start.x, trk.start.y, trk.start.z);
        std::cout << "[cutZeroMomentum] ====================================" << std::endl;
        std::cout << "[cutZeroMomentum] trackStartCryo = " << trackStartCryo << std::endl;
        std::cout << "[cutZeroMomentum] Start = " << trk.start.x << ", " << trk.start.y << ", " << trk.start.z << std::endl;
        std::cout << "[cutZeroMomentum] End = " <<  trk.end.x << ", " << trk.end.y << ", " << trk.end.z << std::endl;
        return true;
      }
      else{
        return false;
      }
*/


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


  const Cut cutHasProton([](const caf::SRSliceProxy* slc) {
    return (varProtonTrackInd(slc) >= 0);
  });

  const Cut cutProtonContained([](const caf::SRSliceProxy* slc) {

    if( varProtonTrackInd(slc) >= 0 ){
      auto const& trk = slc->reco.trk.at(varProtonTrackInd(slc));
      return fv.isContained(trk.end.x, trk.end.y, trk.end.z);
/*
      Contained = ( !isnan(trk.end.x) &&
                  ( trk.end.x < -71.1 - 25 && trk.end.x > -369.33 + 25 ) &&
                  !isnan(trk.end.y) &&
                  ( trk.end.y > -181.7 + 25 && trk.end.y < 134.8 - 25 ) &&
                  !isnan(trk.end.z) &&
                  ( trk.end.z > -895.95 + 30 && trk.end.z < 895.95 - 50 ) );
*/
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
          auto const& trk = slc.reco.trk.at(varMuonTrackInd(&slc));
          muonTrackStarts.emplace_back(trk.start.x, trk.start.y, trk.start.z);
        }

        for (auto const& crtHit: sr->crt_hits){

          if(crtHit.plane>=30 && crtHit.plane<=34){

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


  const SpillCut spillcutFDSideCRTHitVetoTestMatching(
      [](const caf::SRSpillProxy* sr){

        vector<TVector3> muonTrackStarts;
        for(auto const& slc : sr->slc){
          if(!cutRFiducial(&slc)) continue;
          if(!kNotClearCosmic(&slc)) continue;
          if(!cutNuScore(&slc)) continue;
          if(!cutFMScore(&slc)) continue;
          if(!cutHasMuon(&slc)) continue;
          auto const& trk = slc.reco.trk.at(varMuonTrackInd(&slc));
          muonTrackStarts.emplace_back(trk.start.x, trk.start.y, trk.start.z);
        }

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


}

