#pragma once

#include "sbnana/SBNAna/Cuts/Cuts.h"
#include "sbnana/SBNAna/Vars/Vars.h"
#include "sbnana/SBNAna/Cuts/VolumeDefinitions.h"

#include "sbnana/CAFAna/StandardRecord/Proxy/SRProxy.h"

#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_Contants.h"
#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_Utilities.h"
#include "sbnana/SBNAna/Cuts/NumuCutsIcarus202106.h"
#include "TVector3.h"
#include <iostream>

#define M_MUON 0.1057
#define M_PROTON 0.938272
#define M_NEUTRON 0.939565
#define E_EffNuclB 0.040

using namespace std;
using namespace ana;

namespace ICARUSNumuXsec{

  const TVector3 NuDirection_NuMI(3.94583e-01, 4.26067e-02, 9.17677e-01);
  const FiducialVolumeTool fv = FiducialVolumeTool();

}
