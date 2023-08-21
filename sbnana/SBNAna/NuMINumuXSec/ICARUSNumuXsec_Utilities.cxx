#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_Utilities.h"
#include "sbnana/CAFAna/Core/Utilities.h"

using namespace ICARUSNumuXsec;

//---------------------------------------------------
// ActiveColumeTool

ActiveVolumeTool& ActiveVolumeTool::Instance(){
  static ActiveVolumeTool avt;
  return avt;
}

//---------------------------------------------------
// VertexContained

VertexContained& VertexContained::Instance(){
  static VertexContained vc;
  return vc;
}

//---------------------------------------------------
// TrackContained

TrackContained& TrackContained::Instance(){
  static TrackContained tc;
  return tc;
}

//---------------------------------------------------
// VolumeTool

bool VolumeTool::isContained(double x, double y, double z) const {

  int _containedCryo = containedCryo(x,y,z);
  return (_containedCryo==0 || _containedCryo==1);

}

int VolumeTool::containedCryo(double x, double y, double z) const {

  int out=-1;
  // Cryo0
  if( x <= 0 ){
    if ( !isnan(x) &&
         ( x < fvCryo0.xmax - XMargin && x > fvCryo0.xmin + XMargin ) &&
         !isnan(y) &&
         ( y > fvCryo0.ymin + YMargin && y < fvCryo0.ymax - YMargin ) &&
         !isnan(z) &&
         ( z > fvCryo0.zmin + ZMarginUp && z < fvCryo0.zmax - ZMarginDown )
       ) out = 0;
  }
  // Cryo1
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

int VolumeTool::TPCIndex(double x, double y, double z) const {

  double x_cath = x<=0 ? (fvCryo0.xmax+fvCryo0.xmin)/2. : (fvCryo1.xmax+fvCryo1.xmin)/2.;
  double x_from_cath = x-x_cath;
  if(x_from_cath<0 && z<0) return 0;
  else if(x_from_cath<0 && z>0) return 1;
  else if(x_from_cath>0 && z<0) return 2;
  else if(x_from_cath>0 && z>0) return 3;
  else return 4;

}

double VolumeTool::GetTotalVolume() const {

  double vol0 = (fvCryo0.xmax-fvCryo0.xmin) * (fvCryo0.ymax-fvCryo0.ymin) * (fvCryo0.zmax-fvCryo0.zmin);
  double vol1 = (fvCryo1.xmax-fvCryo1.xmin) * (fvCryo1.ymax-fvCryo1.ymin) * (fvCryo1.zmax-fvCryo1.zmin);

  return vol0+vol1;

}

//---------------------------------------------------
// NuMICoordinateTool

NuMICoordinateTool::NuMICoordinateTool(){

  NuDirection_NuMI.SetXYZ(3.94583e-01, 4.26067e-02, 9.17677e-01);

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


//---------------------------------------------------
// NuMIPPFXWeightTool

NuMIPPFXWeightTool::NuMIPPFXWeightTool(){

  // test

  const std::string fname = "root://fndcadoor.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/persistent/users/jskim/NuMINumuXSec/Flux/2023-07-31_out_450.37_7991.98_79512.66_QEL11.root";
  std::cout << "[NuMIPPFXWeightTool::NuMIPPFXWeightTool] Reading root file: " << fname << std::endl;
  TFile *f_ppfx = TFile::Open(fname.c_str());

  if(f_ppfx->IsZombie()){
    std::cout << "NuMIPPFXWeightTool: Failed to open " << fname << std::endl;
    std::abort();
  }

  std::cout << "[NuMIPPFXWeightTool::NuMIPPFXWeightTool] Filling in fWeight" << std::endl;

  for (int hcIdx : {0, 1}) {
    for (int flavIdx : {0, 1}) {
      for (int signIdx : {0, 1}) {
        std::string hNamePPFX = "ppfx_flux_weights/hweights_";
        if (hcIdx == 0)
          hNamePPFX += "fhc_";
        else
          hNamePPFX += "rhc_";
        if (flavIdx == 0)
          hNamePPFX += "nue";
        else
          hNamePPFX += "numu";
        if (signIdx == 1) hNamePPFX += "bar";

        TH1* h_ppfx = (TH1*)f_ppfx->Get(hNamePPFX.c_str());
        if (!h_ppfx) {
          std::cout << "NuMIPpfxFluxWeight: failed to find " << hNamePPFX << " in " << f_ppfx->GetName()
                    << std::endl;
          std::abort();
        }
        h_ppfx = (TH1*)h_ppfx->Clone(ana::UniqueName().c_str());
        h_ppfx->SetDirectory(0);

        fWeight[hcIdx][flavIdx][signIdx] = h_ppfx;
      }
    }
  }

  std::cout << "[NuMIPPFXWeightTool::NuMIPPFXWeightTool] Done!" << std::endl;

}

NuMIPPFXWeightTool& NuMIPPFXWeightTool::Instance(){
  static NuMIPPFXWeightTool nppfxwt;
  return nppfxwt;
}

double NuMIPPFXWeightTool::GetWeight(const caf::SRSliceProxy* slc) const {

  if ( slc->truth.index < 0 || abs(slc->truth.initpdg) == 16 ) return 1.0;

  if ( !fWeight[0][0][0] ) {
    std::cout << "Trying to access un-available weight array..." << std::endl;
    std::abort();
  }

  unsigned int hcIdx = 0; // assume always FHC for now...
  unsigned int flavIdx = ( abs(slc->truth.initpdg) == 12 ) ? 0 : 1;
  unsigned int signIdx = ( slc->truth.initpdg > 0 ) ? 0 : 1;

  TH1* h = fWeight[hcIdx][flavIdx][signIdx];
  assert(h);

  const int bin = h->FindBin( slc->truth.E );
  if(bin == 0 || bin == h->GetNbinsX()+1) return 1.0;
  return h->GetBinContent(bin);

}

double NuMIPPFXWeightTool::GetFirstNuWeight(const caf::SRSpillProxy* sr) const {

  if(sr->mc.nu.size()==0) return 1.0;
  
  if ( abs(sr->mc.nu[0].initpdg) == 16 ) return 1.0;

  if ( !fWeight[0][0][0] ) {
    std::cout << "Trying to access un-available weight array..." << std::endl;
    std::abort();
  }

  unsigned int hcIdx = 0; // assume always FHC for now...
  unsigned int flavIdx = ( abs(sr->mc.nu[0].initpdg) == 12 ) ? 0 : 1;
  unsigned int signIdx = ( sr->mc.nu[0].initpdg > 0 ) ? 0 : 1;

  TH1* h = fWeight[hcIdx][flavIdx][signIdx];
  assert(h);

  const int bin = h->FindBin( sr->mc.nu[0].E );
  if(bin == 0 || bin == h->GetNbinsX()+1) return 1.0;
  return h->GetBinContent(bin);

}

//---------------------------------------------------
// dEdXTemplateTool

dEdXTemplateTool::dEdXTemplateTool(){

  std::cout << "[dEdXTemplateTool::dEdXTemplateTool] called" << std::endl;

  cet::search_path sp("FW_SEARCH_PATH");

  std::string fTemplateFile = "dEdxrestemplates.root";

  sp.find_file(fTemplateFile, fROOTfile);

  std::cout << "[dEdXTemplateTool::dEdXTemplateTool] Template found from " << fROOTfile << std::endl;

  TFile *file = TFile::Open(fROOTfile.c_str());
  dedx_range_pro = (TProfile*)file->Get("dedx_range_pro");
  dedx_range_ka  = (TProfile*)file->Get("dedx_range_ka");
  dedx_range_pi  = (TProfile*)file->Get("dedx_range_pi");
  dedx_range_mu  = (TProfile*)file->Get("dedx_range_mu");

  std::cout << "[dEdXTemplateTool::dEdXTemplateTool] Done!" << std::endl;

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

double dEdXTemplateTool::CalculateChi2(const caf::Proxy<caf::SRTrackCalo>& calo, int ptlType) const {

  int npt = 0;
  double chi2 = 0;
  for(unsigned i = 0; i < calo.points.size(); ++i) { //hits
    const auto& pt = calo.points[i];
    if (i == 0 || i == calo.points.size() - 1) continue;

    if( pt.rr>=26. ) continue;

    double dedx = GetdEdX(pt.rr, ptlType);
    double dedx_err = GetdEdXErr(pt.rr, ptlType);

    bool UseThisPoint = (dedx>0.5) && (dedx_err>0);
    //bool UseThisPoint = (dedx>0.) && (dedx_err>0);

    if(UseThisPoint){
      double errdedx = 0.04231 + 0.0001783 * pt.dedx * pt.dedx; //resolution on dE/dx
      errdedx *= pt.dedx;
      chi2 += pow( (pt.dedx-dedx)/std::sqrt( pow(dedx_err, 2) + pow(errdedx, 2) ), 2);
      npt++;
    }

  }
  if(npt){
    chi2 /= npt;
  }

  return chi2;

}

double dEdXTemplateTool::CalculateInelasticPionChi2(const caf::Proxy<caf::SRTrackCalo>& calo, int ptlType) const {

  int npt = 0;
  double chi2 = 0;
  for(unsigned i = 0; i < calo.points.size(); ++i) { //hits
    const auto& pt = calo.points[i];
    if (i == 0 || i == calo.points.size() - 1) continue;

    if( pt.rr>=26. ) continue;
    //if( pt.rr<=3. ) continue;

    static double dedx = 1.8814277274698800;
    double dedx_err = 0.3903743442373010;

    bool UseThisPoint = (dedx>0.5) && (dedx_err>0);
    //bool UseThisPoint = (dedx>0.) && (dedx_err>0);

    if(UseThisPoint){
      double errdedx = 0.04231 + 0.0001783 * pt.dedx * pt.dedx; //resolution on dE/dx
      errdedx *= pt.dedx;
      chi2 += pow( (pt.dedx-dedx)/std::sqrt( pow(dedx_err, 2) + pow(errdedx, 2) ), 2);
      npt++;
    }

  }
  if(npt){
    chi2 /= npt;
  }

  return chi2;

}

dEdXTemplateTool& dEdXTemplateTool::Instance(){
  static dEdXTemplateTool dedxtt;
  return dedxtt;
}

//---------------------------------------------------
// SterileNuTool

SterileNuTool::SterileNuTool(){

  sin2th = 0.10;
  m2 = 7.3;

}

double SterileNuTool::GetOscProb(double LoE, int i, int f) const{

  // flav : 0/1/2 = e/m/t

  return 1. - (sin2th*sin2th) * pow( TMath::Sin(1.27 * m2 * LoE), 2 );
  
}

SterileNuTool& SterileNuTool::Instance(){
  static SterileNuTool snt;
  return snt;
}

//---------------------------------------------------
// ParticleTool

ParticleTool::ParticleTool(){
}

ParticleTool& ParticleTool::Instance(){
  static ParticleTool ptlt;
  return ptlt;
}

double ParticleTool::GetMass(int pdg) const {

  const int apdg = std::abs(pdg);
  if(apdg==13) return M_MUON;
  else if(apdg==211) return M_CHARGEDPION;
  else if(apdg==111) return M_PIZERO;
  else if(apdg==2212) return M_PROTON;
  else if(apdg==2112) return M_NEUTRON;
  else return 0.;

}

//---------------------------------------------------
// InteractionTool

InteractionTool::InteractionTool(){

  ClearIndices();
  UseGHepRecord = true;

}
InteractionTool& InteractionTool::Instance(){
  static InteractionTool intt;
  return intt;
}

void InteractionTool::ClearIndices() const {
  MuonIndices.clear();
  ProtonIndices.clear();
  NeutronIndices.clear();
  PipIndices.clear();
  PimIndices.clear();
  Pi0Indices.clear();
}

InteractionTool::NParticles InteractionTool::GetNParticles(const caf::SRSliceProxy* slc) const {

  ClearIndices();

  InteractionTool::NParticles nptls;

  if(UseGHepRecord){

    // using ghep
    int counter = -1;
    for(const auto& ghepptl: slc->truth.ghepptl){
      counter += 1;
      if(ghepptl.gstatus!=1) continue;

      const int pdg = ghepptl.pdg;
      const int apdg = abs(pdg);
      if(apdg==13){
        nptls.NMuon++;
        MuonIndices.push_back(counter);
      }

      if( pdg==2212 ){
        nptls.NProton++;
        ProtonIndices.push_back(counter);
      }
      else if( pdg==2112 ){
        nptls.NNeutron++;
        NeutronIndices.push_back(counter);
      }
      else if( pdg==211 ){
        nptls.NPip++;
        PipIndices.push_back(counter);
      }
      else if( pdg==-211 ){
        nptls.NPim++;
        PimIndices.push_back(counter);
      }
      else if( pdg==111 ){
        Pi0Indices.push_back(counter);
        nptls.NPi0++;
      }

    } // END prim loop

  }
  else{
    // using prim
    int counter = -1;
    for(const auto& prm: slc->truth.prim){
      counter += 1;

      if(prm.start_process!=0){
        continue;
      }

      const int pdg = prm.pdg;
      const int apdg = abs(pdg);
      if(apdg==13){
        nptls.NMuon++;
        MuonIndices.push_back(counter);
      }

      if( pdg==2212 ){
        nptls.NProton++;
        ProtonIndices.push_back(counter);
      }
      else if( pdg==2112 ){
        nptls.NNeutron++;
        NeutronIndices.push_back(counter);
      }
      else if( pdg==211 ){
        nptls.NPip++;
        PipIndices.push_back(counter);
      }
      else if( pdg==-211 ){
        nptls.NPim++;
        PimIndices.push_back(counter);
      }
      else if( pdg==111 ){
        Pi0Indices.push_back(counter);
        nptls.NPi0++;
      }

    } // END prim loop
  }

  return nptls;

}

//---------------------------------------------------
// Printing
void ICARUSNumuXsec::PrintPrimaries(const caf::SRSliceProxy* slc){

  std::cout << "[PrintPrimaries] called" << std::endl;

  for(unsigned int i_prim=0; i_prim<slc->truth.prim.size(); ++i_prim){

    const auto& prim = slc->truth.prim[i_prim];
    printf("  - %d-th prim\n",i_prim);
    printf("    - G4ID = %d\n", prim.G4ID.GetValue());
    printf("    - pdg = %d\n", prim.pdg.GetValue());
    printf("    - Parent ID = %d\n", prim.parent.GetValue());
    printf("    - start = (%1.2f, %1.2f, %1.2f)\n", prim.start.x.GetValue(), prim.start.y.GetValue(), prim.start.z.GetValue());
    printf("    - end = (%1.2f, %1.2f, %1.2f)\n", prim.end.x.GetValue(), prim.end.y.GetValue(), prim.end.z.GetValue());
    const float dist = std::hypot(prim.end.x - prim.start.x, prim.end.y - prim.start.y, prim.end.z - prim.start.z);
    printf("    - |end-start| = %1.2f\n", dist);
    printf("    - start_process = %d\n", prim.start_process.GetValue());
    printf("    - end_process = %d\n", prim.end_process.GetValue());

  }

}

// For a given truth particle, find the reco object
//   Track : return the longest matched track
int ICARUSNumuXsec::GetMatchedRecoTrackIndex(const caf::SRSliceProxy* slc, int truth_idx, double scorecut){

  if(truth_idx>=0){

    int PTrackInd(-999);
    double LMax(-999.);
    for(std::size_t i(0); i < slc->reco.pfp.size(); ++i){
      const auto& pfp = slc->reco.pfp.at(i);
      if(pfp.trackScore<scorecut) continue;
      const auto& trk = pfp.trk;
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

//   Shower : return (TODO energetic?) shower
int ICARUSNumuXsec::GetMatchedRecoShowerIndex(const caf::SRSliceProxy* slc, int truth_idx, double scorecut){

  if(truth_idx>=0){

    int PTrackInd(-999);
    for(std::size_t i(0); i < slc->reco.pfp.size(); ++i){
      const auto& pfp = slc->reco.pfp.at(i);
      if(pfp.trackScore>=scorecut) continue;
      const auto& shw = pfp.shw;
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
std::vector<double> ICARUSNumuXsec::GetMatchedRecoShowerIndices(const caf::SRSliceProxy* slc, int truth_idx, double scorecut){

  std::vector<double> rets;

  if(truth_idx>=0){

    for(std::size_t i(0); i < slc->reco.pfp.size(); ++i){
      const auto& pfp = slc->reco.pfp.at(i);
      if(pfp.trackScore>=scorecut) continue;
      const auto& shw = pfp.shw;
      if( shw.truth.p.pdg==slc->truth.prim.at(truth_idx).pdg && shw.truth.p.G4ID==slc->truth.prim.at(truth_idx).G4ID ){
        rets.push_back( i );
      }
    }
  }

  return rets;

}

int ICARUSNumuXsec::GetMatchedRecoStubIndex(const caf::SRSliceProxy* slc, int truth_idx){

  if(truth_idx>=0){

    int PTrackInd(-999);
    for(std::size_t i(0); i < slc->reco.stub.size(); ++i){
      const auto& stub = slc->reco.stub.at(i);
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


std::vector<std::string> ICARUSNumuXsec::GetGENIEMultisigmaKnobNames(){

  return {
"ZExpA1CCQE",
"ZExpA2CCQE",
"ZExpA3CCQE",
"ZExpA4CCQE",
"RPA_CCQE",
"CoulombCCQE",
"NormCCMEC",
"NormNCMEC",
"MaNCEL",
"EtaNCEL",
"MaCCRES",
"MvCCRES",
"MaNCRES",
"MvNCRES",
"NonRESBGvpCC1pi",
"NonRESBGvpCC2pi",
"NonRESBGvpNC1pi",
"NonRESBGvpNC2pi",
"NonRESBGvnCC1pi",
"NonRESBGvnCC2pi",
"NonRESBGvnNC1pi",
"NonRESBGvnNC2pi",
"NonRESBGvbarpCC1pi",
"NonRESBGvbarpCC2pi",
"NonRESBGvbarpNC1pi",
"NonRESBGvbarpNC2pi",
"NonRESBGvbarnCC1pi",
"NonRESBGvbarnCC2pi",
"NonRESBGvbarnNC1pi",
"NonRESBGvbarnNC2pi",
"RDecBR1gamma",
"RDecBR1eta",
"NormCCCOH",
"NormNCCOH",
"AhtBY",
"BhtBY",
"CV1uBY",
"CV2uBY",
"MFP_pi",
"FrCEx_pi",
"FrInel_pi",
"FrAbs_pi",
"FrPiProd_pi",
"MFP_N",
"FrCEx_N",
"FrInel_N",
"FrAbs_N",
"FrPiProd_N",
  };

}

std::vector<std::string> ICARUSNumuXsec::GetGENIEMorphKnobNames(){
  return {
"VecFFCCQEshape",
"DecayAngMEC",
"Theta_Delta2Npi",
"ThetaDelta2NRad",
  };
}

std::vector<std::string> ICARUSNumuXsec::GetGENIEDependentKnobNames(){
  return {
"ZExpAVariationResponse",
"NCELVariationResponse",
"CCRESVariationResponse",
"NCRESVariationResponse",
"DISBYVariationResponse",
"FSI_pi_VariationResponse",
"FSI_N_VariationResponse",
  };
}

std::vector<std::string> ICARUSNumuXsec::GetGENIEMultisimKnobNames(){

  return {
"ZNormCCQE",
"ZExpAVariationResponse",
"NCELVariationResponse",
"CCRESVariationResponse",
"NCRESVariationResponse",
"NonRESBGvpCC1pi",
"NonRESBGvpCC2pi",
"NonRESBGvpNC1pi",
"NonRESBGvpNC2pi",
"NonRESBGvnCC1pi",
"NonRESBGvnCC2pi",
"NonRESBGvnNC1pi",
"NonRESBGvnNC2pi",
"NonRESBGvbarpCC1pi",
"NonRESBGvbarpCC2pi",
"NonRESBGvbarpNC1pi",
"NonRESBGvbarpNC2pi",
"NonRESBGvbarnCC1pi",
"NonRESBGvbarnCC2pi",
"NonRESBGvbarnNC1pi",
"NonRESBGvbarnNC2pi",
"RDecBR1gamma",
"RDecBR1eta",
"DISBYVariationResponse",
//"FormZone",
"FSI_pi_VariationResponse",
"FSI_N_VariationResponse",
  };

}

//---------------------------------------------------
static double TrackScoreCutValue = 0.45;
bool ICARUSNumuXsec::IsPFPTrack(const caf::SRPFPProxy& pfp){
  return !isnan(pfp.trackScore) && pfp.trackScore>TrackScoreCutValue;
}
bool ICARUSNumuXsec::IsPFPShower(const caf::SRPFPProxy& pfp){
  return !isnan(pfp.trackScore) && pfp.trackScore<=TrackScoreCutValue && pfp.trackScore>0.;
}

