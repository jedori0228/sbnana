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



