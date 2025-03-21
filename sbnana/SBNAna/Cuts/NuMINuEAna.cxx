#include "sbnana/SBNAna/Cuts/NuMIXSecCuts.h"
#include "sbnana/SBNAna/Vars/NuMIXSecVars.h"
#include "sbnana/SBNAna/Cuts/NuMINuEAna.h"

namespace ana {

  const TruthVar kTruth_ElectronIndex([](const caf::SRTrueInteractionProxy *nu) -> int {
    double max_E(-999);
    int truth_idx(-1);
    for(std::size_t i(0); i < nu->prim.size(); ++i){
      // primary
      if( nu->prim.at(i).start_process!=0 ) continue;
      // electron
      if( abs(nu->prim.at(i).pdg)!=11 ) continue;
      // non-nan genE
      if(isnan(nu->prim.at(i).genE)) continue;

      double this_E = nu->prim.at(i).genE;
      // if larger E, update
      if(this_E>max_E){
        max_E = this_E;
        truth_idx = i;
      }
    }
    return truth_idx;
  });

  const TruthVar kTruth_ElectronKE([](const caf::SRTrueInteractionProxy *nu) -> double {
    double ret(-5.f);

    int truth_idx = kTruth_ElectronIndex(nu);
    if(truth_idx>=0){
      ret = nu->prim.at(truth_idx).genE - M_ELECTRON;
    }

    return ret;
  });

  bool Is1eNp(const caf::Proxy<caf::SRTrueInteraction>& true_int){

    if ( true_int.index < 0 ) return false;
    if ( abs(true_int.pdg) != 12 ||
         !true_int.iscc ||
         std::isnan(true_int.position.x) || std::isnan(true_int.position.y) || std::isnan(true_int.position.z) ||
         !isInFV(true_int.position.x, true_int.position.y, true_int.position.z) )
      return false; // not signal

    unsigned int nEl(0), nP(0);
    for ( auto const& prim : true_int.prim ) {
      if ( prim.start_process != 0 ) continue;

      double momentum = sqrt( (prim.genp.x*prim.genp.x) + (prim.genp.y*prim.genp.y) + (prim.genp.z*prim.genp.z) );

      if ( abs(prim.pdg) == 11 ) {
        nEl += 1;
      }

      if ( abs(prim.pdg) == 2212 ) {
        if (momentum>0.4){
          nP+=1;
        }
      }

    }

    bool Is1eNp = nEl==1 && nP>0;

    return Is1eNp;

  }

  const TruthCut kTruthCut_Is1eNp([](const caf::SRTrueInteractionProxy* nu) {
    return Is1eNp(*nu);
  });

}
