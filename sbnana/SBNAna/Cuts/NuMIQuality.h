#pragma once

#include "sbnana/CAFAna/Core/Var.h"
#include "sbnana/CAFAna/Core/Cut.h"
#include "sbnana/CAFAna/Core/MultiVar.h"
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>


namespace ana
{

  class NuMIQuality
  {
  public:
    NuMIQuality();
    ~NuMIQuality();

    bool IsBadTriggeredSpill(unsigned int run, unsigned int subrun, unsigned int event) const;

  protected:
    std::string filename_bad_triggered_spill;
  };

  // set up to use the flux weight
  static const NuMIQuality kNuMIQuality;
  extern const SpillCut kNuMINotBadTriggeredSpill;

}
