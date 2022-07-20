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

//==== For a given truth particle, find the reco object
//====   Track : return the longest matched track
int ICARUSNumuXsec::GetMatchedRecoTrackIndex(const caf::SRSliceProxy* slc, int truth_idx){

  if(truth_idx>=0){

    int PTrackInd(-999);
    double LMax(-999.);
    for(std::size_t i(0); i < slc->reco.trk.size(); ++i){
      auto const& trk = slc->reco.trk.at(i);
      if( trk.truth.p.pdg==slc->truth.prim.at(truth_idx).pdg && trk.truth.p.G4ID==slc->truth.prim.at(truth_idx).G4ID ){
        if(trk.len > LMax){
          PTrackInd = i;
          LMax = trk.len;
        }
      }
    }
    return PTrackInd;
  }
  else{
    return -999.;
  }

}

//====   Shower : return (TODO energetic?) shower
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

NuMICoordinateTool::NuMICoordinateTool(){

  TMatrixDRow(rotMatNtoI,0) = {0.921035925, 0.022715103, 0.388814672};
  TMatrixDRow(rotMatNtoI,1) = {0., 0.998297825, -0.058321970};
  TMatrixDRow(rotMatNtoI,2) = {-0.389477631, 0.053716629, 0.919468161};
  rotMatNtoI.Print();

  TMatrixDColumn(tranVecNtoI,0) = {-315.120380, -33.644912, -733.632532};
  tranVecNtoI.Print();

}

TVector3 NuMICoordinateTool::GetICARUSCoord(double x, double y, double z) const {

  TMatrixD coordN(3,1);
  TMatrixDColumn(coordN,0) = {x, y, z};

  TMatrixD ret = (rotMatNtoI*coordN+tranVecNtoI);

  return TVector3(ret(0,0), ret(1,0), ret(2,0));

}
