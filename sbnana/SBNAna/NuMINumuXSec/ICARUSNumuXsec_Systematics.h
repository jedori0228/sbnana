#pragma once

#include "sbnana/CAFAna/Core/ISyst.h"
#include "sbnana/CAFAna/Core/Utilities.h"
#include "sbnana/CAFAna/Analysis/ExpInfo.h"
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_Systematics.h"
#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_Variables.h"

using namespace ana;
using namespace std;
using namespace ICARUSNumuXsec;

namespace ana{

  class MuonMomentumScaleSyst: public ISyst
  {
  public:
    MuonMomentumScaleSyst(double _uncertainty) :
      ISyst("MuonMomentumScaleSyst", "Muon momentum scale systematic"), uncertainty(_uncertainty) {}

    void Shift(double sigma, caf::SRSliceProxy* slc, double& weight) const override
    {
      //==== Only shifting the Muon track
      int muonTrackIndx = varMuonTrackInd(slc);
      if( muonTrackIndx>= 0 ){
        auto const& trk = slc->reco.trk.at(muonTrackIndx);
        bool Contained = ( !isnan(trk.end.x) &&
                    ( trk.end.x < -71.1 - 25 && trk.end.x > -369.33 + 25 ) &&
                    !isnan(trk.end.y) &&
                    ( trk.end.y > -181.7 + 25 && trk.end.y < 134.8 - 25 ) &&
                    !isnan(trk.end.z) &&
                    ( trk.end.z > -895.95 + 30 && trk.end.z < 895.95 - 50 ) );
        if(Contained){
          slc->reco.trk.at(muonTrackIndx).rangeP.p_muon *= 1.0 + uncertainty * sigma;
        }
        else{
          slc->reco.trk.at(muonTrackIndx).mcsP.fwdP_muon *= 1.0 + uncertainty * sigma;
        }
      }
    }

  private:
    double uncertainty;

  }; // END Class


} // namespace ana
