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
      fvCryo0 = {-358.489, -61.94,  // x
                 -181.86, +134.96,  // y
                 -894.95, +894.85};
      fvCryo1 = { +61.94, +358.489,  // x
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

  }; // END Class FiducialVolumeTool

  int GetMatchedRecoTrackIndex(const caf::SRSliceProxy* slc, int truth_idx);
  int GetMatchedRecoShowerIndex(const caf::SRSliceProxy* slc, int truth_idx);
  int GetMatchedRecoStubIndex(const caf::SRSliceProxy* slc, int truth_idx);


}
