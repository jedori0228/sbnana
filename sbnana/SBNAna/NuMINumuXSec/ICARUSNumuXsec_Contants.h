#pragma once

#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_Utilities.h"
#include "sbnana/SBNAna/CRTPMTMatching/ICARUSCRTPMTMatching_Tool.h"

using namespace std;

namespace ICARUSNumuXsec{

  static const TVector3 NuDirection_NuMI(3.94583e-01, 4.26067e-02, 9.17677e-01);

  static const ActiveVolumeTool& av = ActiveVolumeTool::Instance();
  static const VertexContained& fv = VertexContained::Instance();
  static const TrackContained& fv_track = TrackContained::Instance();

  static const NuMICoordinateTool& nct = NuMICoordinateTool::Instance();

  static const dEdXTemplateTool& dedxtempt = dEdXTemplateTool::Instance();

  static const SterileNuTool& snt = SterileNuTool::Instance();

  static const ParticleTool& ptlt = ParticleTool::Instance();

  static const InteractionTool& intt = InteractionTool::Instance();

}

