#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_Weight.h"

using namespace std;
using namespace ana;

namespace ICARUSNumuXsec{

  const Var wSterileOscProb([](const caf::SRSliceProxy* slc) -> double {

    if(!isnan(slc->truth.E)){

      double nuprod_NuMI_x = slc->truth.prod_vtx.x/100.;
      double nuprod_NuMI_y = slc->truth.prod_vtx.y/100.;
      double nuprod_NuMI_z = slc->truth.prod_vtx.z/100.;
      TVector3 nuprod_ICARUS = nct.GetICARUSCoord(nuprod_NuMI_x, nuprod_NuMI_y, nuprod_NuMI_z);
      TVector3 nuvtx_ICARUS = TVector3(slc->truth.position.x/100., slc->truth.position.y/100., slc->truth.position.z/100.);

      TVector3 nuProdToVtxVector = nuvtx_ICARUS-nuprod_ICARUS;
      //std::cout << "[wSterileOscProb] L = " << nuProdToVtxVector.Mag() << std::endl;

      double LoE = nuProdToVtxVector.Mag()/( slc->truth.E*1000. );

      return snt.GetOscProb(LoE,1,1);

    }
    else return 1.;

  });

}

