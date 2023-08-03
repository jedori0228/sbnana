#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_Weights.h"

using namespace ana;

namespace ICARUSNumuXsec{

Var NuMIPPFXCVWeight([](const caf::SRSliceProxy* slc) -> double {

  static const NuMIPPFXWeightTool nppfxwt = NuMIPPFXWeightTool::Instance();
  return nppfxwt.GetWeight(slc);

});

SpillVar NuMIPPFXCVFirstNuSpillWeight([](const caf::SRSpillProxy *sr) -> double {

  static const NuMIPPFXWeightTool nppfxwt = NuMIPPFXWeightTool::Instance();
  return nppfxwt.GetFirstNuWeight(sr);

});


}

