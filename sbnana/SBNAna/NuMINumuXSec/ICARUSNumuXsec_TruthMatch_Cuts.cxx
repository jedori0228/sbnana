#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_TruthMatch_Cuts.h"

using namespace std;
using namespace ana;

namespace ICARUSNumuXsec{

namespace TruthMatch{

  const Cut TruthProtonStopped([](const caf::SRSliceProxy* slc) {
    int truth_idx = TruthProtonIndex(slc);
    if(truth_idx>=0){
      const auto& prim = slc->truth.prim.at(truth_idx);
      if(prim.end_process==2) return true;
      else return false;
    }
    else return false;
  });

  const Cut TruthProtonInelastic([](const caf::SRSliceProxy* slc) {
    int truth_idx = TruthProtonIndex(slc);
    if(truth_idx>=0){
      const auto& prim = slc->truth.prim.at(truth_idx);
      if(prim.end_process==9 || prim.end_process==10) return true;
      else return false;
    }
    else return false;
  });

  const Cut TruthChargedPionContained([](const caf::SRSliceProxy* slc) {
    int truth_idx = TruthChargedPionIndex(slc);
    if(truth_idx>=0){
      if(slc->truth.prim.at(truth_idx).contained) return true;
      else return false;
    }
    else return false;
  });
  // endprocess
  const Cut TruthChargedPionStopped([](const caf::SRSliceProxy* slc) {
    int truth_idx = TruthChargedPionIndex(slc);
    if(truth_idx>=0){
      const auto& prim = slc->truth.prim.at(truth_idx);
      if(prim.end_process==2) return true;
      else return false;
    }
    else return false;
  });
  const Cut TruthChargedPionDecay([](const caf::SRSliceProxy* slc) {
    int truth_idx = TruthChargedPionIndex(slc);
    if(truth_idx>=0){
      const auto& prim = slc->truth.prim.at(truth_idx);
      if(prim.end_process==3) return true;
      else return false;
    }
    else return false;
  });
  const Cut TruthChargedPionInelastic([](const caf::SRSliceProxy* slc) {
    int truth_idx = TruthChargedPionIndex(slc);
    if(truth_idx>=0){
      const auto& prim = slc->truth.prim.at(truth_idx);
      if(prim.end_process==9 || prim.end_process==10) return true;
      else return false;
    }
    else return false;
  });
  const Cut TruthChargedPionOther([](const caf::SRSliceProxy* slc) {
    int truth_idx = TruthChargedPionIndex(slc);
    if(truth_idx>=0){
      const auto& prim = slc->truth.prim.at(truth_idx);
      if(prim.end_process==2 ||
         prim.end_process==3 ||
         prim.end_process==9 || prim.end_process==10) return true;
      else return false;
    }
    else return false;
  });

  const Cut TruthMuonHasMichel([](const caf::SRSliceProxy* slc) {
    return TruthMuonMichelIndex(slc)>=0;
  });
  const Cut TruthChargedPionHasMichel([](const caf::SRSliceProxy* slc) {
    return TruthChargedPionMichelMatchedShowerIndex(slc)>=0;
  });

} // end namespace TruthMatch

} // end namespace ICARUSNumuXsec
