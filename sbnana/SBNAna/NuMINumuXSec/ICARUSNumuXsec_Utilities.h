#pragma once

#include "sbnana/SBNAna/Cuts/Cuts.h"
#include "sbnana/SBNAna/Vars/Vars.h"
#include "sbnana/SBNAna/Cuts/VolumeDefinitions.h"

#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

#include "sbnana/SBNAna/Cuts/NumuCutsIcarus202106.h"
#include "TVector3.h"
#include <iostream>

using namespace std;
using namespace ana;

namespace ICARUSNumuXsec{

  class FiducialVolumeTool{

  public:

    FiducialVolumeTool()
    // : fvCryo0(fvfd_cryo1), fvCryo1(fvfd_cryo2)
    {
      fvCryo0 = {-368.49,  -71.94,  // x
                 -181.86, +134.96,  // y
                 -894.95, +894.85};
      fvCryo1 = {+71.94,  +368.49,  // x
                 -181.86, +134.96,  // y
                 -894.95, +894.85};
      XMargin=25.;
      YMargin=25.;
      ZMarginUp=30.;
      ZMarginDown=50.;
    }

    FidVol fvCryo0, fvCryo1;
    double XMargin, YMargin, ZMarginUp, ZMarginDown;

    bool isContained(double x, double y, double z) const;
    int containedCryo(double x, double y, double z) const;

  };

}
