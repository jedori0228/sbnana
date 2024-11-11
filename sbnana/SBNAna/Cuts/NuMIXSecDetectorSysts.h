//ETW May 2020 Fresh version for SBN

#pragma once

#include "sbnana/SBNAna/Cuts/NuMIXSecCuts.h"
#include "sbnana/SBNAna/Vars/NuMIXSecVars.h"
#include "sbnana/CAFAna/Core/ISyst.h"

/*

EXAMPLE)

  const NuMIXSecDetectorSysts kNuMIXSecFrontIndPlaneGainSyst(
    NuMIXSecDetectorSysts::kFrontIndPlaneGain,
    "NuMIXSecFrontIndPlaneGainSyst",
    "Front ind. plane gain #pm10%"
  );
  const NuMIXSecDetectorSysts kNuMIXSecFrontIndPlaneNoiseSyst(
    NuMIXSecDetectorSysts::kFrontIndPlaneNoise,
    "NuMIXSecFrontIndPlaneNoiseSyst",
    "Front ind. plane noise +10%"
  );
  const NuMIXSecDetectorSysts kNuMIXSecFrontIndPlaneSignalShapeSyst(
    NuMIXSecDetectorSysts::kFrontIndPlaneSignalShape,
    "NuMIXSecFrontIndPlaneSignalShapeSyst",
    "Front ind. plane signal shape"
  );
  const NuMIXSecDetectorSysts kNuMIXSecMiddleIndPlaneTransparencySyst(
    NuMIXSecDetectorSysts::kMiddleIndPlaneTransparency,
    "NuMIXSecMiddleIndPlaneTransparencySyst",
    "Middle ind. plane transparency"
  );

*/

namespace ana
{

  class NuMIXSecDetectorSysts: public ISyst
  {
  public:

    enum DetSystType {
      kFrontIndPlaneGain,
      kFrontIndPlaneNoise,
      kFrontIndPlaneSignalShape,
      kFrontIndPlaneSignalShapeFitted,
      kMiddleIndPlaneTransparency,
      kSCE,
    };

    NuMIXSecDetectorSysts(DetSystType detsyst_type, const std::string& name, const std::string& latexName);

    void Shift(double sigma, caf::SRSliceProxy *sr, double& weight) const override;
    void Shift(double sigma, caf::SRTrueInteractionProxy *sr, double& weight) const override;

    double GetSmoothStepFunction(double x, double x_min, double x_max, double offset, double scale) const ;

  private:

    DetSystType kDetSystType;

  };

}
