#pragma once

#include "sbnana/CAFAna/Core/Tree.h"
#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/SBNAna/Cuts/NuMIXSecRecoEventVars.h"
#include "sbnana/SBNAna/Cuts/NuMIXSecTrueEventVars.h"
#include "sbnana/SBNAna/Cuts/NuMIXSecSysts.h"
#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_Variables.h"

using namespace ana;
using namespace std;

namespace ana{

  std::vector<std::string> GetNuMITrueTreeLabels();
  std::vector<TruthVar> GetNuMITrueTreeVars();

  std::vector<std::string> GetNuMIRecoTreeLabels();
  std::vector<Var> GetNuMIRecoTreeVars();

  std::vector<std::string> GetLowQ2RecStudyLabels();
  std::vector<Var> GetLowQ2RecStudyVars();

}
