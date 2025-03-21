#pragma once

#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_Contants.h"
#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_TruthMatch_Variables.h"
#include "sbnana/SBNAna/Vars/Vars.h"
#include "sbnana/CAFAna/Core/Cut.h"
#include "sbnana/CAFAna/Core/MultiVar.h"
#include <iostream>

using namespace std;
using namespace ana;

namespace ICARUSNumuXsec{

namespace TruthMatch{

  extern const Cut TruthProtonStopped;
  extern const Cut TruthProtonInelastic;

  extern const Cut TruthChargedPionContained;
  // endprocess
  extern const Cut TruthChargedPionStopped;
  extern const Cut TruthChargedPionDecay;
  extern const Cut TruthChargedPionInelastic;
  extern const Cut TruthChargedPionOther;

  extern const Cut TruthMuonHasMichel;
  extern const Cut TruthChargedPionHasMichel;
  extern const Cut HasTruthChargedPion;
  extern const Cut HasTruthNeutralPion;

} // end namespace TruthMatch

} // end namespace ICARUSNumuXsec
