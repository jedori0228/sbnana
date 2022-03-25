#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_Utilities.h"

using namespace ICARUSNumuXsec;

bool FiducialVolumeTool::isContained(double x, double y, double z) const {

  int _containedCryo = containedCryo(x,y,z);
  return (_containedCryo==0 || _containedCryo==1);

}

int FiducialVolumeTool::containedCryo(double x, double y, double z) const {

  int out=-1;
  //==== Cryo0
  if( x <= 0 ){
    if ( !isnan(x) &&
         ( x < fvCryo0.xmax - XMargin && x > fvCryo0.xmin + XMargin ) &&
         !isnan(y) &&
         ( y > fvCryo0.ymin + YMargin && y < fvCryo0.ymax - YMargin ) &&
         !isnan(z) &&
         ( z > fvCryo0.zmin + ZMarginUp && z < fvCryo0.zmax - ZMarginDown )
       ) out = 0;
  }
  //==== Cryo1
  else{
    if ( !isnan(x) &&
         ( x > fvCryo1.xmin + XMargin && x < fvCryo1.xmax - XMargin ) &&
         !isnan(y) &&
         ( y > fvCryo0.ymin + YMargin && y < fvCryo0.ymax - YMargin ) &&
         !isnan(z) &&
         ( z > fvCryo0.zmin + ZMarginUp && z < fvCryo0.zmax - ZMarginDown )
       ) out = 1;
  }
  return out;

}

int ICARUSNumuXsec::GetMatchedRecoTrackIndex(const caf::SRSliceProxy* slc, int truth_idx){

  if(truth_idx>=0){

    int PTrackInd(-999);
    for(std::size_t i(0); i < slc->reco.trk.size(); ++i){
      auto const& trk = slc->reco.trk.at(i);
      if( trk.truth.p.pdg==slc->truth.prim.at(truth_idx).pdg && trk.truth.p.G4ID==slc->truth.prim.at(truth_idx).G4ID ){
        PTrackInd = i;
      }
    }
    return PTrackInd;
  }
  else{
    return -999.;
  }

}

int ICARUSNumuXsec::GetMatchedRecoShowerIndex(const caf::SRSliceProxy* slc, int truth_idx){

  if(truth_idx>=0){

    int PTrackInd(-999);
    for(std::size_t i(0); i < slc->reco.shw.size(); ++i){
      auto const& shw = slc->reco.shw.at(i);
      if( shw.truth.p.pdg==slc->truth.prim.at(truth_idx).pdg && shw.truth.p.G4ID==slc->truth.prim.at(truth_idx).G4ID ){
        PTrackInd = i;
      }
    }
    return PTrackInd;
  }
  else{
    return -999.;
  }

}

int ICARUSNumuXsec::GetMatchedRecoStubIndex(const caf::SRSliceProxy* slc, int truth_idx){

  if(truth_idx>=0){

    int PTrackInd(-999);
    for(std::size_t i(0); i < slc->reco.stub.size(); ++i){
      auto const& stub = slc->reco.stub.at(i);
      if( stub.truth.p.pdg==slc->truth.prim.at(truth_idx).pdg && stub.truth.p.G4ID==slc->truth.prim.at(truth_idx).G4ID ){
        PTrackInd = i;
      }
    }
    return PTrackInd;
  }
  else{
    return -999.;
  }

}
