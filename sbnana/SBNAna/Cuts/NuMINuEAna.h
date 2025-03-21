#pragma once

#include "sbnana/CAFAna/Core/Cut.h"

namespace ana
{

  // Truth
  extern const TruthVar kTruth_ElectronIndex;
  extern const TruthVar kTruth_ElectronKE;
  bool Is1eNp(const caf::Proxy<caf::SRTrueInteraction>& true_int);
  extern const TruthCut kTruthCut_Is1eNp;


}
