#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_Utilities.h"

using namespace ICARUSNumuXsec;

VertexContained& VertexContained::Instance(){
  static VertexContained vc;
  return vc;
}

TrackContained& TrackContained::Instance(){
  static TrackContained tc;
  return tc;
}

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

int FiducialVolumeTool::TPCIndex(double x, double y, double z) const {

  double x_cath = x<=0 ? (fvCryo0.xmax+fvCryo0.xmin)/2. : (fvCryo1.xmax+fvCryo1.xmin)/2.;
  double x_from_cath = x-x_cath;
  if(x_from_cath<0 && z<0) return 0;
  else if(x_from_cath<0 && z>0) return 1;
  else if(x_from_cath>0 && z<0) return 2;
  else if(x_from_cath>0 && z>0) return 3;
  else return 4;

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
    return -1;
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
    return -1;
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
    return -1;
  }

}

//==== Get number of daughter track/shower
int ICARUSNumuXsec::GetNDaughterTracks(const caf::SRSliceProxy* slc, int track_idx){

  if(track_idx>=0){

    auto const& trk = slc->reco.trk.at(track_idx);
    TVector3 trk_end(trk.end.x, trk.end.y, trk.end.z);

    int nStTrk(0);
    for(std::size_t i(0); i < slc->reco.trk.size(); ++i){
      if(i==(std::size_t)track_idx) continue;

      auto const& sttrk = slc->reco.trk.at(i);
      TVector3 sttrk_start(sttrk.start.x, sttrk.start.y, sttrk.start.z);
      double this_dist = (trk_end-sttrk_start).Mag();
      if(this_dist<5.){
        nStTrk++;
      }
    }
    return nStTrk;

  }
  else{
    return 0;
  }

}

int ICARUSNumuXsec::GetNDaughterShowers(const caf::SRSliceProxy* slc, int track_idx){

  if(track_idx>=0){
    
    auto const& trk = slc->reco.trk.at(track_idx);
    TVector3 trk_end(trk.end.x, trk.end.y, trk.end.z);
    
    int nStShw(0);
    for(std::size_t i(0); i < slc->reco.shw.size(); ++i){
      auto const& stshw = slc->reco.shw.at(i);
      TVector3 stshw_start(stshw.start.x, stshw.start.y, stshw.start.z);
      double this_dist = (trk_end-stshw_start).Mag();
      if(this_dist<5.){
        nStShw++;
      }
    }
    return nStShw;

  }
  else{
    return 0;
  }

}


//==== Michel tagging
int ICARUSNumuXsec::GetMichelShowerIndex(const caf::SRSliceProxy* slc, int track_idx){

  if(track_idx>=0){

    auto const& trk = slc->reco.trk.at(track_idx);
    TVector3 trk_end(trk.end.x, trk.end.y, trk.end.z);

    int ret(-1);
    double dist(5.);
    for(std::size_t i(0); i < slc->reco.shw.size(); ++i){
      auto const& stshw = slc->reco.shw.at(i);
      TVector3 stshw_start(stshw.start.x, stshw.start.y, stshw.start.z);
      double this_dist = (trk_end-stshw_start).Mag();
      if(this_dist<dist){
        double stshw_e = stshw.bestplane_energy;
        if(stshw_e<0.05){
          ret = i;
          dist = this_dist;
        }
      }
    }
    return ret;

  }
  else{
    return -1;
  }

}
//==== Pion tagging
bool ICARUSNumuXsec::IsPionTagged(const caf::SRSliceProxy* slc, int track_idx){

  if(track_idx>=0){

    auto const& trk = slc->reco.trk.at(track_idx);
    return (trk.chi2pid[trk.bestplane].chi2_pion<20);


    TVector3 trk_end(trk.end.x, trk.end.y, trk.end.z);

    //==== Check if there's a stitched track at the end of the track
    bool HasStitchedProton(false);
    for(std::size_t i(0); i < slc->reco.trk.size(); ++i){
      if(i==(std::size_t)track_idx) continue;

      auto const& sttrk = slc->reco.trk.at(i);
      TVector3 sttrk_start(sttrk.start.x, sttrk.start.y, sttrk.start.z);
      double this_dist = (trk_end-sttrk_start).Mag();
      int sttrk_bp = sttrk.bestplane;
      double sttrk_chi2proton = sttrk.chi2pid[sttrk_bp].chi2_proton;
      if(this_dist<5. &&  sttrk_chi2proton<40.){
        HasStitchedProton = true;
        break;
      }
    }
    return HasStitchedProton;

/*
    //==== Check if there's a stitched shower at the end of the track
    bool HasStitchedShower(false);
    for(std::size_t i(0); i < slc->reco.shw.size(); ++i){
      auto const& stshw = slc->reco.shw.at(i);
      TVector3 stshw_start(stshw.start.x, stshw.start.y, stshw.start.z);
      double this_dist = (trk_end-stshw_start).Mag();
      double stshw_e = stshw.bestplane_energy;
      if(this_dist<5. && stshw_e>0.10){
        HasStitchedShower = true;
        break;
      }
    }
*/

  }
  else{
    return false;
  }

}

double ICARUSNumuXsec::GetEnergyFromStubCharge(double q){
  static double p0 = -0.00236892;
  static double p1 = 0.383538;
  if(q<=0.) return -999.;
  else return (q*23.6e-9-p0)/p1; // to GeV
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

NuMICoordinateTool& NuMICoordinateTool::Instance(){
  static NuMICoordinateTool nct;
  return nct;
}

dEdXTemplateTool::dEdXTemplateTool(){

  std::cout << "[dEdXTemplateTool::dEdXTemplateTool] Setting up.." << std::endl;

  cet::search_path sp("FW_SEARCH_PATH");

  std::string fTemplateFile = "dEdxrestemplates.root";

  sp.find_file(fTemplateFile, fROOTfile);

  TFile *file = TFile::Open(fROOTfile.c_str());
  dedx_range_pro = (TProfile*)file->Get("dedx_range_pro");
  dedx_range_ka  = (TProfile*)file->Get("dedx_range_ka");
  dedx_range_pi  = (TProfile*)file->Get("dedx_range_pi");
  dedx_range_mu  = (TProfile*)file->Get("dedx_range_mu");

}

double dEdXTemplateTool::GetdEdX(double rr, int ptlType) const {

  int bin = dedx_range_pro->FindBin(rr);
  if(bin>dedx_range_pro->GetNbinsX()) bin = dedx_range_pro->GetNbinsX()-1;
  if(bin>=1&&bin<=dedx_range_pro->GetNbinsX()){

    double bincpro = dedx_range_pro->GetBinContent(bin);
    if (bincpro<1e-6){//for 0 bin content, using neighboring bins
      bincpro = (dedx_range_pro->GetBinContent(bin-1)+dedx_range_pro->GetBinContent(bin+1))/2;
    }
    double bincka = dedx_range_ka->GetBinContent(bin);
    if (bincka<1e-6){
      bincka = (dedx_range_ka->GetBinContent(bin-1)+dedx_range_ka->GetBinContent(bin+1))/2;
    }
    double bincpi = dedx_range_pi->GetBinContent(bin);
    if (bincpi<1e-6){
      bincpi = (dedx_range_pi->GetBinContent(bin-1)+dedx_range_pi->GetBinContent(bin+1))/2;
    }
    double bincmu = dedx_range_mu->GetBinContent(bin);
    if (bincmu<1e-6){
      bincmu = (dedx_range_mu->GetBinContent(bin-1)+dedx_range_mu->GetBinContent(bin+1))/2;
    }

    if(ptlType==0) return bincpro;
    else if(ptlType==1) return bincka;
    else if(ptlType==2) return bincpi;
    else if(ptlType==3) return bincmu;
    else{
      std::cout << "[dEdXTemplateTool::GetdEdX] Wrong ptlType : " << ptlType << std::endl;
      abort();
      return -1.;
    }

  }
  else{
    std::cout << "[dEdXTemplateTool::GetdEdX] rr = " << rr << ", bin = " << bin << std::endl;
    return -1.;
  }

}

double dEdXTemplateTool::GetdEdXErr(double rr, int ptlType) const {

  int bin = dedx_range_pro->FindBin(rr);
  if(bin>dedx_range_pro->GetNbinsX()) bin = dedx_range_pro->GetNbinsX()-1;
  if(bin>=1&&bin<=dedx_range_pro->GetNbinsX()){

    double binepro = dedx_range_pro->GetBinError(bin);
    if (binepro<1e-6){
      binepro = (dedx_range_pro->GetBinError(bin-1)+dedx_range_pro->GetBinError(bin+1))/2;
    }
    double bineka = dedx_range_ka->GetBinError(bin);
    if (bineka<1e-6){
      bineka = (dedx_range_ka->GetBinError(bin-1)+dedx_range_ka->GetBinError(bin+1))/2;
    }
    double binepi = dedx_range_pi->GetBinError(bin);
    if (binepi<1e-6){
      binepi = (dedx_range_pi->GetBinError(bin-1)+dedx_range_pi->GetBinError(bin+1))/2;
    }
    double binemu = dedx_range_mu->GetBinError(bin);
    if (binemu<1e-6){
      binemu = (dedx_range_mu->GetBinError(bin-1)+dedx_range_mu->GetBinError(bin+1))/2;
    }

    if(ptlType==0) return binepro;
    else if(ptlType==1) return bineka;
    else if(ptlType==2) return binepi;
    else if(ptlType==3) return binemu;
    else{
      std::cout << "[dEdXTemplateTool::GetdEdXErr] Wrong ptlType : " << ptlType << std::endl;
      abort();
      return -1.;
    }

  }
  else{
    std::cout << "[dEdXTemplateTool::GetdEdX] rr = " << rr << ", bin = " << bin << std::endl;
    return -1.;
  }

}

dEdXTemplateTool& dEdXTemplateTool::Instance(){
  static dEdXTemplateTool dedxtt;
  return dedxtt;
}

SterileNuTool::SterileNuTool(){

  sin2th = 0.10;
  m2 = 7.3;

}

double SterileNuTool::GetOscProb(double LoE, int i, int f) const{

  //==== flav : 0/1/2 = e/m/t

  return 1. - (sin2th*sin2th) * pow( TMath::Sin(1.27 * m2 * LoE), 2 );
  
}

SterileNuTool& SterileNuTool::Instance(){
  static SterileNuTool snt;
  return snt;
}

TrackStitchingTool& TrackStitchingTool::Instance(){
  static TrackStitchingTool tst;
  return tst;
}

TrackStitchingTool::StichOutput TrackStitchingTool::GetStitchedTrack(
  const caf::Proxy<caf::SRTrack>& motherTrack,
  const caf::Proxy<caf::SRSlice>& motherSlice,
  const caf::SRSpillProxy *sr) const
  {

  TrackStitchingTool::StichOutput ret;
  ret.isFound = false;

  TVector3 trk_start(motherTrack.start.x, motherTrack.start.y, motherTrack.start.z);
  TVector3 trk_end(motherTrack.end.x, motherTrack.end.y, motherTrack.end.z);
  TVector3 trk_dir(motherTrack.dir.x, motherTrack.dir.y, motherTrack.dir.z);

  double MinDist = 999999.;

  int sliceCounter = -1;
  for(const auto& slc2: sr->slc){
    sliceCounter++;
    bool currentIsSameSlice = (&motherSlice==&slc2);

    int trackCounter = -1;
    for(const auto& trk2: slc2.reco.trk){
      trackCounter++;

      if(&motherTrack==&trk2) continue;

      TVector3 trk2_start(trk2.start.x, trk2.start.y, trk2.start.z);
      TVector3 trk2_end(trk2.end.x, trk2.end.y, trk2.end.z);
      TVector3 trk2_dir(trk2.dir.x, trk2.dir.y, trk2.dir.z);

      double dist_start_to_start = (trk_start-trk2_start).Mag();
      double dist_start_to_end = (trk_start-trk2_end).Mag();
      double dist_end_to_start = (trk_end-trk2_start).Mag();
      double dist_end_to_end = (trk_end-trk2_end).Mag();

      vector<double> tmp_dists = {dist_start_to_start, dist_start_to_end, dist_end_to_start, dist_end_to_end};
      std::vector<double>::iterator smallest_it = std::min_element(std::begin(tmp_dists), std::end(tmp_dists));
      int smallest_index = std::distance(tmp_dists.begin(), smallest_it);

      double cosDir = trk_dir.Dot(trk2_dir);

      bool DirectionMatched = false;
      double this_dist = +9999999.;
      if(smallest_index==0){
        DirectionMatched = cosDir<0;
        this_dist = dist_start_to_start;
      }
      else if(smallest_index==1){
        DirectionMatched = cosDir>0;
        this_dist = dist_start_to_end;
      }
      else if(smallest_index==2){
        DirectionMatched = cosDir>0;
        this_dist = dist_end_to_start;
      }
      else if(smallest_index==3){
        DirectionMatched = cosDir<0;
        this_dist = dist_end_to_end;
      }

      //if(motherTracktruth.bestmatch.G4ID==trk2.truth.bestmatch.G4ID){

        if(DirectionMatched && this_dist < MinDist){
          MinDist = this_dist;
          ret.isFound = true;
          ret.isSameSlice = currentIsSameSlice;
          ret.minDist = this_dist;
          ret.foundSliceIdx = sliceCounter;
          ret.foundTrackIdx = trackCounter;
          ret.closestMode = smallest_index;
        }

      //}
/*
      std::cout << "[TrackStitchingTool::GetStitchedTrack] Slice index : " << sliceCounter << std::endl;
      std::cout << "[TrackStitchingTool::GetStitchedTrack] - MinDist = " << MinDist << std::endl;
*/
    }
  }

  return ret;

}

CRTPMTMatchingTool::CRTPMTMatchingTool(){
  std::cout << "[CRTPMTMatchingTool::CRTPMTMatchingTool] called" << std::endl;
}

CRTPMTMatchingTool& CRTPMTMatchingTool::Instance(){
  static CRTPMTMatchingTool cpmt;
  return cpmt;
}

void CRTPMTMatchingTool::SetGateType(GateType gt) const {
  GT = gt;
  if(abs(GT)==1){
    timecut_min = 0.;
    timecut_max = 2.2;
  }
  else if(abs(GT)==2){
    timecut_min = 0.;
    timecut_max = 10.1;
  }
  else{
    std::cout << "[CRTPMTMatchingTool] Wrong gate type = " << gt << std::endl;
    abort();
  }

  std::cout << "[CRTPMTMatchingTool::SetGateType] GT = " << GT << ", timecut_min = " << timecut_min << ", timecut_max = " << timecut_max << std::endl;

}

void CRTPMTMatchingTool::SetInTimeRange(double t_min, double t_max) const {
  timecut_min = t_min;
  timecut_max = t_max;
  std::cout << "[CRTPMTMatchingTool::SetInTimeRange] timecut_min = " << timecut_min << ", timecut_max = " << timecut_max << std::endl;
}

bool CRTPMTMatchingTool::IsInTime(double t_gate) const{
  //std::cout << "timecut_min = " << timecut_min << ", t_gate = " << t_gate << ", timecut_max = " << timecut_max << std::endl;
  return ( timecut_min<=t_gate && t_gate<=timecut_max );
}

int CRTPMTMatchingTool::GetMatchedCRTHitIndex(
  double opt,
  const caf::Proxy<std::vector<caf::SRCRTHit> >& crt_hits,
  int mode) const{

  double mindiff = std::numeric_limits<double>::max();
  int ret=-1;
  for(size_t i=0; i<crt_hits.size(); i++){
    const auto& hit = crt_hits.at(i);
    if(hit.plane>=30 && hit.plane<=34){
      double crtt = hit.t1;
      double this_diff = crtt-opt;
      if(fabs(this_diff)<fabs(mindiff)){
        mindiff = this_diff;
        ret = i;
      }
    }
  }

  return ret;

}

int CRTPMTMatchingTool::GetMatchID(
  double opt,
  const caf::Proxy<std::vector<caf::SRCRTHit> >& crt_hits) const{

  int hasCRTHit = 0;
  static double interval = 0.1;
  int topen = 0, topex = 0, sideen = 0, sideex = 0;
  for(size_t i=0; i<crt_hits.size(); i++){
    const auto& crtHit = crt_hits.at(i);
    double tof = crtHit.t1 - opt;

    if(tof<0 && abs(tof)<interval){
      if(crtHit.plane > 36){
        sideen++;
      }
      else{
        topen++;
      }
    }
    else if(tof>=0 && abs(tof)<interval){
      if(crtHit.plane > 36){
        sideex++;
      }
      else{
        topex++;
      }
    }

  }

  // hasCRTHit = 0, no matched CRT
  // hasCRTHit = 1, 1 entering from Top CRT
  // hasCRTHit = 2, 1 entering from Side CRT
  // hasCRTHit = 3, 1 entering from Top and exiting to Side CRT
  // hasCRTHit = 4, No entering; 1 exiting to top
  // hasCRTHit = 5, No entering, 1 exiting to side
  // hasCRTHit = 6, Multiple entering
  // hasCRTHit = 7, Multiple entering and exiting to side
  // hasCRTHit = 8, all other cases

  if (topen == 0 && sideen == 0 && topex == 0 && sideex == 0)
    hasCRTHit = 0;
  else if (topen == 1 && sideen == 0 && topex == 0 && sideex == 0)
    hasCRTHit = 1;
  else if (topen == 0 && sideen == 1 && topex == 0 && sideex == 0)
    hasCRTHit = 2;
  else if (topen == 1 && sideen == 0 && topex == 0 && sideex == 1)
    hasCRTHit = 3;
  else if (topen == 0 && sideen == 0 && topex == 1 && sideex == 0)
    hasCRTHit = 4;
  else if (topen == 0 && sideen == 0 && topex == 0 && sideex == 1)
    hasCRTHit = 5;
  else if (topen >= 1 && sideen >= 1 && topex == 0 && sideex == 0)
    hasCRTHit = 6;
  else if (topen >= 1 && sideen >= 1 && topex == 0 && sideex >= 1)
    hasCRTHit = 7;
  else
    hasCRTHit = 8;

  return hasCRTHit;

}

bool CRTPMTMatchingTool::IsNegativeTOF(double timediff) const{
  return (timediff>-0.1 && timediff<0.);
}

