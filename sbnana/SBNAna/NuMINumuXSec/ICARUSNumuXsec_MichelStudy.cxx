#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_MichelStudy.h"

using namespace std;
using namespace ana;
using namespace ICARUSNumuXsec::TruthMatch;

namespace ICARUSNumuXsec{

namespace MichelStudy{

  const SpillMultiVar SliceIndiciesWithMichelElectron([](const caf::SRSpillProxy* sr) -> vector<double> {

    vector<double> rets;

    for(unsigned int i_slc=0; i_slc<sr->slc.size(); i_slc++){

      const auto& slc = sr->slc.at(i_slc);

      int cpionTrackIndex = TruthChargedPionMatchedTrackIndex(&slc);
      if(cpionTrackIndex>=0){

        // matched reco trk
        const auto& trk = slc.reco.pfp.at(cpionTrackIndex).trk;

        int idx_daughter_muon = -1;
        for(const auto& d_id: trk.truth.p.daughters){
          for(unsigned int i_tp=0; i_tp<sr->true_particles.size(); i_tp++){
            const auto& tp = sr->true_particles.at(i_tp);
            if( (int)tp.G4ID==(int)d_id ){
              if(abs(tp.pdg)==13){
                idx_daughter_muon = i_tp;
                break;
              }
            }
          }
        }

        if(idx_daughter_muon>=0){

          int idx_michel_electron = -1;
          const auto& tp_daughter_muon = sr->true_particles.at(idx_daughter_muon);
          for(const auto& d_id: tp_daughter_muon.daughters){
            for(unsigned int i_tp=0; i_tp<sr->true_particles.size(); i_tp++){
              const auto& tp = sr->true_particles.at(i_tp); 
              if( (int)tp.G4ID==(int)d_id ){
                if(abs(tp.pdg)==11){
                  idx_michel_electron = i_tp;
                  break;
                }
              }
            }
          }

          if(idx_michel_electron>=0){
            rets.push_back(i_slc);
          }

        }


      }

    }

    return rets;

  });


  const SpillMultiVar TruthChargedPionWithMichel_TruthChargedPionKEs([](const caf::SRSpillProxy* sr) -> vector<double> {

    vector<double> rets;

    vector<double> sliceIndices = SliceIndiciesWithMichelElectron(sr);
    for(const auto slice_idx: sliceIndices){
      const auto& slc = sr->slc.at(slice_idx);
      double pion_ke = TruthChargedPionKE(&slc);
      if(pion_ke>=0.) rets.push_back(pion_ke); 
    }

    return rets;

  });
  const SpillMultiVar TruthChargedPionWithMichel_TruthMichelElectronKEs([](const caf::SRSpillProxy* sr) -> vector<double> {

    vector<double> rets;

    for(unsigned int i_slc=0; i_slc<sr->slc.size(); i_slc++){

      const auto& slc = sr->slc.at(i_slc);

      int cpionTrackIndex = TruthChargedPionMatchedTrackIndex(&slc);
      if(cpionTrackIndex>=0){

        // matched reco trk
        const auto& trk = slc.reco.pfp.at(cpionTrackIndex).trk;

        int idx_daughter_muon = -1;
        for(const auto& d_id: trk.truth.p.daughters){
          for(unsigned int i_tp=0; i_tp<sr->true_particles.size(); i_tp++){
            const auto& tp = sr->true_particles.at(i_tp);
            if( (int)tp.G4ID==(int)d_id ){
              if(abs(tp.pdg)==13){
                idx_daughter_muon = i_tp;
                break;
              }
            }
          }
        }
        if(idx_daughter_muon>=0){

          int idx_michel_electron = -1;
          const auto& tp_daughter_muon = sr->true_particles.at(idx_daughter_muon);
          for(const auto& d_id: tp_daughter_muon.daughters){
            for(unsigned int i_tp=0; i_tp<sr->true_particles.size(); i_tp++){
              const auto& tp = sr->true_particles.at(i_tp); 
              if( (int)tp.G4ID==(int)d_id ){
                if(abs(tp.pdg)==11){
                  idx_michel_electron = i_tp;
                  break;
                }
              }
            }
          }

          if(idx_michel_electron>=0){
            const auto& tp_michel_electron = sr->true_particles.at(idx_michel_electron);
            rets.push_back( tp_michel_electron.genE-M_ELECTRON );
          }

        }

      }

    }

    return rets;

  });

  const SpillMultiVar TruthChargedPionWithMichel_TruthMichelElectronMatchedRecoShowerEnergies([](const caf::SRSpillProxy* sr) -> vector<double> {

    vector<double> rets;

    for(unsigned int i_slc=0; i_slc<sr->slc.size(); i_slc++){

      const auto& slc = sr->slc.at(i_slc);

      int cpionTrackIndex = TruthChargedPionMatchedTrackIndex(&slc);
      if(cpionTrackIndex>=0){

        // matched reco trk
        const auto& trk = slc.reco.pfp.at(cpionTrackIndex).trk;

        int idx_daughter_muon = -1;
        for(const auto& d_id: trk.truth.p.daughters){
          for(unsigned int i_tp=0; i_tp<sr->true_particles.size(); i_tp++){
            const auto& tp = sr->true_particles.at(i_tp);
            if( (int)tp.G4ID==(int)d_id ){
              if(abs(tp.pdg)==13){
                idx_daughter_muon = i_tp;
                break;
              }
            }
          }
        }
        if(idx_daughter_muon>=0){

          int idx_michel_electron = -1;
          const auto& tp_daughter_muon = sr->true_particles.at(idx_daughter_muon);
          for(const auto& d_id: tp_daughter_muon.daughters){
            for(unsigned int i_tp=0; i_tp<sr->true_particles.size(); i_tp++){
              const auto& tp = sr->true_particles.at(i_tp); 
              if( (int)tp.G4ID==(int)d_id ){
                if(abs(tp.pdg)==11){
                  idx_michel_electron = i_tp;
                  break;
                }
              }
            }
          }

          if(idx_michel_electron>=0){

            // now loop over reco pfps in this slice, and find anything matched to michel electron

            for(const auto& pfp: slc.reco.pfp){

              const auto& shw = pfp.shw;
              double matcehd_shower_energy = -0.5;
              if((int)shw.truth.p.G4ID==(int)idx_michel_electron){
                matcehd_shower_energy = shw.plane[2].energy;
              }
              rets.push_back( matcehd_shower_energy );

            }

          }

        }

      }

    }

    return rets;

  });

  const Var MichelElectronMatchedRecoShowerIndex([](const caf::SRSliceProxy* slc) -> double {
    int cpionindex = TruthChargedPionIndex(slc);
    double ret = -999.;
    if(cpionindex>=0){
      const auto& prim_cpion = slc->truth.prim.at(cpionindex);
      for(unsigned int i_pfp=0; i_pfp<slc->reco.pfp.size(); i_pfp++){
        const auto& pfp = slc->reco.pfp.at(i_pfp);
        const auto& shw = pfp.shw;

        bool isElectron = abs(shw.truth.p.pdg)==11;
        const auto& parentG4ID = shw.truth.p.parent;

        if(isElectron && ((int)parentG4ID==(int)prim_cpion.G4ID)){
          ret = i_pfp;
          break;
        }
      }
    }
    return ret;
  });

  const Var MichelTest([](const caf::SRSliceProxy* slc) -> double {

    if(slc->truth.prim.size()==0) return -1;
    std::cout << "[JSKIMDEBUG][MichelTest] Printing primaries of this slice" << std::endl;
    for(unsigned int i=0; i<slc->truth.prim.size(); i++){
      const auto& prim = slc->truth.prim.at(i);
      std::cout << "[JSKIMDEBUG][MichelTest] - index: " << i << std::endl;
      std::cout << "[JSKIMDEBUG][MichelTest]   - pdg = " << prim.pdg << std::endl;
      std::cout << "[JSKIMDEBUG][MichelTest]   - status = " << prim.gstatus << std::endl;
      std::cout << "[JSKIMDEBUG][MichelTest]   - start_process = " << prim.start_process << std::endl;

    }
    return 1.;

  });


} // end namespace TwoTrack

} // end namespace ICARUSNumuXsec
