#pragma once

#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_TwoTrack_Cuts.h"

using namespace ana::PrimaryUtil;

namespace ICARUSNumuXsec{

namespace Efficiency{

  std::map<int, double> GetNuIndexTruthVarPair(
    const caf::SRSpillProxy* sr,
    std::function<bool(const TrueInteraction&)> isSignal,
    std::function<double(const TrueInteraction&) > trueth_var
  );

  std::vector<double> GetVectorForEfficiency(
    const caf::SRSpillProxy* sr,
    std::function<bool(const TrueInteraction&)> isSignal,
    std::function<double(const TrueInteraction&) > trueth_var,
    const Cut& RecoCut=kNoCut,
    const SpillCut& RecoSpillCut=kNoSpillCut
  );

  //extern const SpillMultiVar Eff_MuonP_Den;
  //extern const SpillMultiVar Eff_MuonP_Num;

}

}
