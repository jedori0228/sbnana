#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_Utilities.h"

using namespace ICARUSNumuXsec;

bool FiducialVolumeTool::isContained(double x, double y, double z) const {

  int _containedCryo = containedCryo(x,y,z);
  return (_containedCryo==0 || _containedCryo==1);

/*
  //==== Cryo0
  if( x <= 0 ){
    return ( !isnan(x) &&
             ( x < fvCryo0.xmax - XMargin && x > fvCryo0.xmin + XMargin ) &&
             !isnan(y) &&
             ( y > fvCryo0.ymin + YMargin && y < fvCryo0.ymax - YMargin ) &&
             !isnan(z) &&
             ( z > fvCryo0.zmin + ZMarginUp && z < fvCryo0.zmax - ZMarginDown ) 
           );
  }
  //==== Cryo1
  else{
    return ( !isnan(x) &&
             ( x > fvCryo1.xmin + XMargin && x < fvCryo1.xmax - XMargin ) &&
             !isnan(y) &&
             ( y > fvCryo0.ymin + YMargin && y < fvCryo0.ymax - YMargin ) &&
             !isnan(z) &&
             ( z > fvCryo0.zmin + ZMarginUp && z < fvCryo0.zmax - ZMarginDown )
           );
  }
*/

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

