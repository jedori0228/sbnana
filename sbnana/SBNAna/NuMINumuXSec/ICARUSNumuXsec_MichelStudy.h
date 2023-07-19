#pragma once

#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_TruthMatch_Variables.h"

using namespace ana;
using namespace std;

namespace ICARUSNumuXsec{

namespace MichelStudy{

  extern const SpillMultiVar SliceIndiciesWithMichelElectron;
  extern const SpillMultiVar TruthChargedPionWithMichel_TruthChargedPionKEs;
  extern const SpillMultiVar TruthChargedPionWithMichel_TruthMichelElectronKEs;
  extern const SpillMultiVar TruthChargedPionWithMichel_TruthMichelElectronMatchedRecoShowerEnergies;

  extern const Var MichelElectronMatchedRecoShowerIndex;

  extern const Var MichelTest;


} // end namespace TwoTrack

} // end namespace ICARUSNumuXsec
