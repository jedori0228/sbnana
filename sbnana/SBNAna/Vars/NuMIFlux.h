// BH - 2023
// HEAVILY based on NuMI flux syst, and thanks to Tony Wood for discussing the right histograms to use

#pragma once

#include "sbnana/CAFAna/Core/Var.h"

#include <string>

class TH1;

namespace ana
{

  class NuMIPpfxFluxWeight
  {
  public:
    NuMIPpfxFluxWeight();
    ~NuMIPpfxFluxWeight();
    mutable TH1* fWeight[2][2][2]; // [fhc/rhc][nue/numu][nu/nubar]
    mutable TH1* fWeight_G3Chase[2][2][2]; // [fhc/rhc][nue/numu][nu/nubar]

  protected:
    std::string fFluxFilePath;
    std::string fFluxFilePath_G3Chase;
  };

  // set up to use the flux weight
  static const NuMIPpfxFluxWeight FluxWeightNuMI;
  extern const Var kGetNuMIFluxWeight;
  extern const TruthVar kGetTruthNuMIFluxWeight;

  extern const Var kGetNuMIFluxWeightWithG3Chase;
  extern const TruthVar kGetTruthNuMIFluxWeightWithG3Chase;

}
