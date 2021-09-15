#ifndef myTempNuMIVariables_h
#define myTempNuMIVariables_h

#include "myConstants.h"

//==== FMScore
const Cut cutTempNuMImyFMScore([](const caf::SRSliceProxy* slc) {
    return ( !isnan(slc->fmatch.score) && slc->fmatch.score < 6.0 );
  });

//==== Muon

const Var varTempNuMIMuonTrackInd([](const caf::SRSliceProxy* slc) -> int {

  //==== The (dis)qualification of a slice is based upon the track level information.
  float Atslc, Longest(0);
  float Chi2Proton, Chi2Muon;
  bool AtSlice, Contained, MaybeMuonExiting, MaybeMuonContained;
  int PTrackInd(-1);

  for(std::size_t i(0); i < slc->reco.trk.size(); ++i){

    auto const& trk = slc->reco.trk.at(i);
    //==== First we calculate the distance of each track to the slice vertex.
    Atslc = sqrt( pow( slc->vertex.x - trk.start.x, 2.0 ) +
            pow( slc->vertex.y - trk.start.y, 2.0 ) +
            pow( slc->vertex.z - trk.start.z, 2.0 ) );
    //===== We require that the distance of the track from the slice is less than
    //==== 10 cm and that the parent of the track has been marked as the primary.
    AtSlice = ( Atslc < 10.0 );//&& trk.parent_is_primary);
    
    if(trk.bestplane == 0){
      Chi2Proton = trk.chi2pid0.chi2_proton;
      Chi2Muon = trk.chi2pid0.chi2_muon;
    }
    else if (trk.bestplane == 1){
      Chi2Proton = trk.chi2pid1.chi2_proton;
      Chi2Muon = trk.chi2pid1.chi2_muon;
    }
    else{
      Chi2Proton = trk.chi2pid2.chi2_proton;
      Chi2Muon = trk.chi2pid2.chi2_muon;
    }

    Contained = ( !isnan(trk.end.x) &&
                ( trk.end.x < -71.1 - 25 && trk.end.x > -369.33 + 25 ) &&
                !isnan(trk.end.y) &&
                ( trk.end.y > -181.7 + 25 && trk.end.y < 134.8 - 25 ) &&
                !isnan(trk.end.z) &&
                ( trk.end.z > -895.95 + 30 && trk.end.z < 895.95 - 50 ) );
    MaybeMuonExiting = ( !Contained && trk.len > 100);

    //==== our current numi sample does not have this chi2
    //MaybeMuonContained = ( Contained && trk.len > 50 );

    MaybeMuonContained = ( Contained && Chi2Proton > 60 && Chi2Muon < 30 && trk.len > 50 ); // MuonSel__MuonChi2LT30_and_ProtonChi2GT60
    //MaybeMuonContained = ( Contained && trk.len > 50 ); // NoChi2

    if( AtSlice && ( MaybeMuonExiting || MaybeMuonContained ) && trk.len > Longest ){
      Longest = trk.len;
      PTrackInd = i;
    }


  }

  return PTrackInd;

});

const Cut cutTempNuMIHasMuon([](const caf::SRSliceProxy* slc) {
  return (varTempNuMIMuonTrackInd(slc) >= 0);
});

const Cut cutTempNuMIMuonContained([](const caf::SRSliceProxy* slc) {

  bool Contained(false);
  if( varTempNuMIMuonTrackInd(slc) >= 0 ){
    auto const& trk = slc->reco.trk.at(varTempNuMIMuonTrackInd(slc));
    Contained = ( !isnan(trk.end.x) &&
                ( trk.end.x < -71.1 - 25 && trk.end.x > -369.33 + 25 ) &&
                !isnan(trk.end.y) && 
                ( trk.end.y > -181.7 + 25 && trk.end.y < 134.8 - 25 ) &&
                !isnan(trk.end.z) &&
                ( trk.end.z > -895.95 + 30 && trk.end.z < 895.95 - 50 ) );
    return Contained;
  }
  return Contained;

});

const Var varTempNuMIMuonRecoP([](const caf::SRSliceProxy* slc) -> float {
  float p(-5.f);
  bool Contained(false);
  
  if( varTempNuMIMuonTrackInd(slc) >= 0 ){
    auto const& trk = slc->reco.trk.at(varTempNuMIMuonTrackInd(slc));
    Contained = ( !isnan(trk.end.x) &&
                ( trk.end.x < -71.1 - 25 && trk.end.x > -369.33 + 25 ) &&
                !isnan(trk.end.y) &&
                ( trk.end.y > -181.7 + 25 && trk.end.y < 134.8 - 25 ) &&
                !isnan(trk.end.z) &&
                ( trk.end.z > -895.95 + 30 && trk.end.z < 895.95 - 50 ) );
    if(Contained){
      p = trk.rangeP.p_muon;
    }
    else{
      p = trk.mcsP.fwdP_muon;
    }
  }
  return p;
});

const Var varTempNuMIMuonTrueP([](const caf::SRSliceProxy* slc) -> float {
  float p(-5.f);
  
  if( varTempNuMIMuonTrackInd(slc) >= 0 ){
    auto const& trk = slc->reco.trk.at(varTempNuMIMuonTrackInd(slc));
    p = sqrt( pow(trk.truth.p.genp.x, 2) + 
              pow(trk.truth.p.genp.y, 2) + 
              pow(trk.truth.p.genp.z, 2) );
  }
  return p;
});

const Var varTempNuMIMuonTruePDG([](const caf::SRSliceProxy* slc) -> int {
  int pdg(-9999999);

  if( varTempNuMIMuonTrackInd(slc) >= 0 ){
    auto const& trk = slc->reco.trk.at(varTempNuMIMuonTrackInd(slc));
    pdg = trk.truth.p.pdg;
  }
  return pdg;
});

const Cut cutTempNuMIMuonMatchedToMuon([](const caf::SRSliceProxy* slc) {
  return ( abs(varTempNuMIMuonTruePDG(slc)) == 13 );
});

const Cut cutTempNuMIMuonMatchedToProton([](const caf::SRSliceProxy* slc) {
  return ( abs(varTempNuMIMuonTruePDG(slc)) == 2212 );
});

const Var varTempNuMIMuonPResidual([](const caf::SRSliceProxy* slc) -> float {
  float pr(-999);

  if( varTempNuMIMuonTrackInd(slc) >= 0 ){
    double RecoP = varTempNuMIMuonRecoP(slc);
    double TrueP = varTempNuMIMuonTrueP(slc);
    pr = RecoP-TrueP;
  }
  return pr;
});

const Var varTempNuMIMuonPResidualFraction([](const caf::SRSliceProxy* slc) -> float {
  float prf(-999);

  if( varTempNuMIMuonTrackInd(slc) >= 0 ){
    double RecoP = varTempNuMIMuonRecoP(slc);
    double TrueP = varTempNuMIMuonTrueP(slc);
    prf = (RecoP-TrueP)/TrueP;
  }
  return prf;
});

const Var varTempNuMIMuonRecoCosineTheta([](const caf::SRSliceProxy* slc) -> float {
  float costh(-5.f);

  if( varTempNuMIMuonTrackInd(slc) >= 0 ){
    auto const& trk = slc->reco.trk.at(varTempNuMIMuonTrackInd(slc));
    costh = trk.costh;
  }
  return costh;
});

const Var varTempNuMIMuonTrueCosineTheta([](const caf::SRSliceProxy* slc) -> float {
  float costh(-5.f);

  if( varTempNuMIMuonTrackInd(slc) >= 0 ){
    auto const& trk = slc->reco.trk.at(varTempNuMIMuonTrackInd(slc));
    TVector3 v3(trk.truth.p.genp.x, trk.truth.p.genp.y, trk.truth.p.genp.z);
    costh = v3.CosTheta();
  }
  return costh;
});

const Var varTempNuMIMuonCosineThetaResidual([](const caf::SRSliceProxy* slc) -> float {
  float costhr(-999);

  if( varTempNuMIMuonTrackInd(slc) >= 0 ){
    double RecoCosineTheta = varTempNuMIMuonRecoCosineTheta(slc);
    double TrueCosineTheta = varTempNuMIMuonTrueCosineTheta(slc);
    costhr = RecoCosineTheta-TrueCosineTheta;
  }
  return costhr;
});

const Var varTempNuMIMuonCosineThetaResidualFraction([](const caf::SRSliceProxy* slc) -> float {
  float costhrf(-999);

  if( varTempNuMIMuonTrackInd(slc) >= 0 ){
    double RecoCosineTheta = varTempNuMIMuonRecoCosineTheta(slc);
    double TrueCosineTheta = varTempNuMIMuonTrueCosineTheta(slc);
    costhrf = (RecoCosineTheta-TrueCosineTheta)/TrueCosineTheta;
  }
  return costhrf;
});

//==== Proton

//====   Track-based

const Var varTempNuMIProtonTrackInd([](const caf::SRSliceProxy* slc) -> int {

  //==== The (dis)qualification of a slice is based upon the track level information.
  float Atslc, Longest(0);
  float Chi2Proton;//, Chi2Muon;
  bool AtSlice, PassPID, Contained;
  int PTrackInd(-1);

  for(std::size_t i(0); i < slc->reco.trk.size(); ++i){

    auto const& trk = slc->reco.trk.at(i);
    //==== First we calculate the distance of each track to the slice vertex.
    Atslc = sqrt( pow( slc->vertex.x - trk.start.x, 2.0 ) +
            pow( slc->vertex.y - trk.start.y, 2.0 ) +
            pow( slc->vertex.z - trk.start.z, 2.0 ) );
    //===== We require that the distance of the track from the slice is less than
    //==== 10 cm and that the parent of the track has been marked as the primary.
    AtSlice = ( Atslc < 10.0 );//&& trk.parent_is_primary);

    Contained = ( !isnan(trk.end.x) &&
                ( trk.end.x < -71.1 - 25 && trk.end.x > -369.33 + 25 ) &&
                !isnan(trk.end.y) &&
                ( trk.end.y > -181.7 + 25 && trk.end.y < 134.8 - 25 ) &&
                !isnan(trk.end.z) &&
                ( trk.end.z > -895.95 + 30 && trk.end.z < 895.95 - 50 ) );

    if(trk.bestplane == 0){
      Chi2Proton = trk.chi2pid0.chi2_proton;
      //Chi2Muon = trk.chi2pid0.chi2_muon;
    }
    else if (trk.bestplane == 1){
      Chi2Proton = trk.chi2pid1.chi2_proton;
      //Chi2Muon = trk.chi2pid1.chi2_muon;
    }
    else{
      Chi2Proton = trk.chi2pid2.chi2_proton;
      //Chi2Muon = trk.chi2pid2.chi2_muon;
    }

    PassPID = Chi2Proton < 60;

    if( AtSlice && Contained && PassPID && trk.len > Longest ){
      Longest = trk.len;
      PTrackInd = i;
    }

  }

  return PTrackInd;

});

const Cut cutTempNuMIHasProton([](const caf::SRSliceProxy* slc) {
  return (varTempNuMIProtonTrackInd(slc) >= 0);
});

const Cut cutTempNuMIProtonContained([](const caf::SRSliceProxy* slc) {

  bool Contained(false);
  if( varTempNuMIProtonTrackInd(slc) >= 0 ){
    auto const& trk = slc->reco.trk.at(varTempNuMIProtonTrackInd(slc));
    Contained = ( !isnan(trk.end.x) &&
                ( trk.end.x < -71.1 - 25 && trk.end.x > -369.33 + 25 ) &&
                !isnan(trk.end.y) && 
                ( trk.end.y > -181.7 + 25 && trk.end.y < 134.8 - 25 ) &&
                !isnan(trk.end.z) &&
                ( trk.end.z > -895.95 + 30 && trk.end.z < 895.95 - 50 ) );
    return Contained;
  }
  return Contained;

});

const Var varTempNuMIProtonCaloP([](const caf::SRSliceProxy* slc) -> float {
  float p(-5.f);

  if( varTempNuMIProtonTrackInd(slc) >= 0 ){
    auto const& trk = slc->reco.trk.at(varTempNuMIProtonTrackInd(slc));

    double best_calo = 0.;
    if(trk.bestplane == 0){
      best_calo = trk.calo0.ke;
    }
    else if(trk.bestplane == 1){
      best_calo = trk.calo1.ke;
    }
    else{
      best_calo = trk.calo2.ke;
    }
    double E_Proton = best_calo+M_PROTON;
    return sqrt( E_Proton*E_Proton - M_PROTON*M_PROTON );
  }
  return p;

});

const Var varTempNuMIProtonRecoP([](const caf::SRSliceProxy* slc) -> float {
  float p(-5.f);
  bool Contained(false);
  
  if( varTempNuMIProtonTrackInd(slc) >= 0 ){
    auto const& trk = slc->reco.trk.at(varTempNuMIProtonTrackInd(slc));

    Contained = ( !isnan(trk.end.x) &&
                ( trk.end.x < -71.1 - 25 && trk.end.x > -369.33 + 25 ) &&
                !isnan(trk.end.y) &&
                ( trk.end.y > -181.7 + 25 && trk.end.y < 134.8 - 25 ) &&
                !isnan(trk.end.z) &&
                ( trk.end.z > -895.95 + 30 && trk.end.z < 895.95 - 50 ) );
    if(Contained){
      p = trk.rangeP.p_proton;
    }
    else{
      p = trk.mcsP.fwdP_proton;
    }

  }
  return p;
});

const Var varTempNuMIProtonTrueP([](const caf::SRSliceProxy* slc) -> float {
  float p(-5.f);
  
  if( varTempNuMIProtonTrackInd(slc) >= 0 ){
    auto const& trk = slc->reco.trk.at(varTempNuMIProtonTrackInd(slc));
    p = sqrt( pow(trk.truth.p.genp.x, 2) + 
              pow(trk.truth.p.genp.y, 2) + 
              pow(trk.truth.p.genp.z, 2) );
  }
  return p;
});

const Var varTempNuMIProtonTruePDG([](const caf::SRSliceProxy* slc) -> int {
  int pdg(-9999999);

  if( varTempNuMIProtonTrackInd(slc) >= 0 ){
    auto const& trk = slc->reco.trk.at(varTempNuMIProtonTrackInd(slc));
    pdg = trk.truth.p.pdg;
  }
  return pdg;
});

const Cut cutTempNuMIProtonMatchedToMuon([](const caf::SRSliceProxy* slc) {
  return ( abs(varTempNuMIProtonTruePDG(slc)) == 13 );
});

const Cut cutTempNuMIProtonMatchedToProton([](const caf::SRSliceProxy* slc) {
  return ( abs(varTempNuMIProtonTruePDG(slc)) == 2212 );
});

const Var varTempNuMIProtonPResidual([](const caf::SRSliceProxy* slc) -> float {
  float pr(-999);

  if( varTempNuMIProtonTrackInd(slc) >= 0 ){
    double RecoP = varTempNuMIProtonRecoP(slc);
    double TrueP = varTempNuMIProtonTrueP(slc);
    pr = RecoP-TrueP;
  }
  return pr;
});

const Var varTempNuMIProtonPResidualFraction([](const caf::SRSliceProxy* slc) -> float {
  float prf(-999);

  if( varTempNuMIProtonTrackInd(slc) >= 0 ){
    double RecoP = varTempNuMIProtonRecoP(slc);
    double TrueP = varTempNuMIProtonTrueP(slc);
    prf = (RecoP-TrueP)/TrueP;
  }
  return prf;
});

const Var varTempNuMIProtonRecoCosineTheta([](const caf::SRSliceProxy* slc) -> float {
  float costh(-5.f);

  if( varTempNuMIProtonTrackInd(slc) >= 0 ){
    auto const& trk = slc->reco.trk.at(varTempNuMIProtonTrackInd(slc));
    costh = trk.costh;
  }
  return costh;
});

const Var varTempNuMIProtonTrueCosineTheta([](const caf::SRSliceProxy* slc) -> float {
  float costh(-5.f);

  if( varTempNuMIProtonTrackInd(slc) >= 0 ){
    auto const& trk = slc->reco.trk.at(varTempNuMIProtonTrackInd(slc));
    TVector3 v3(trk.truth.p.genp.x, trk.truth.p.genp.y, trk.truth.p.genp.z);
    costh = v3.CosTheta();
  }
  return costh;
});

const Var varTempNuMIProtonCosineThetaResidual([](const caf::SRSliceProxy* slc) -> float {

  float costhr(-999);

  if( varTempNuMIProtonTrackInd(slc) >= 0 ){
    double RecoCosineTheta = varTempNuMIProtonRecoCosineTheta(slc);
    double TrueCosineTheta = varTempNuMIProtonTrueCosineTheta(slc);
    costhr = RecoCosineTheta-TrueCosineTheta;
  }
  return costhr;
});

const Var varTempNuMIProtonCosineThetaResidualFraction([](const caf::SRSliceProxy* slc) -> float {

  float costhrf(-999);

  if( varTempNuMIProtonTrackInd(slc) >= 0 ){
    double RecoCosineTheta = varTempNuMIProtonRecoCosineTheta(slc);
    double TrueCosineTheta = varTempNuMIProtonTrueCosineTheta(slc);
    costhrf = (RecoCosineTheta-TrueCosineTheta)/TrueCosineTheta;
  }
  return costhrf;
});

//====   Stub-based

const Var varTempNuMINStub([](const caf::SRSliceProxy* slc) -> int {

  return slc->reco.stub.size();

});

//==== Neutrino

const Var varTempNuMINeutrinoCombinedEnergy([](const caf::SRSliceProxy* slc) -> double {

  double E_Nu(-5.f);
  if( varTempNuMIMuonTrackInd(slc) >= 0 && varTempNuMIProtonTrackInd(slc) >= 0 ){

		double P_Muon = varTempNuMIMuonRecoP(slc); // comdinbed
		double E_Muon = sqrt( P_Muon*P_Muon + M_MUON*M_MUON );

		double P_Proton = varTempNuMIProtonRecoP(slc);
    double E_Proton = sqrt( P_Proton*P_Proton + M_PROTON*M_PROTON );
    double KE_Proton = E_Proton - M_PROTON;

		E_Nu = E_Muon + KE_Proton + E_EffNuclB;

  }
  return E_Nu;

});

const Var varTempNuMINeutrinoCombinedEnergyResidual([](const caf::SRSliceProxy* slc) -> double {

  double E_Nur(-999);
  if( varTempNuMIMuonTrackInd(slc) >= 0 && varTempNuMIProtonTrackInd(slc) >= 0 ){

    double RecoE = varTempNuMINeutrinoCombinedEnergy(slc);
    double TrueE = varNeutrinoTruthE(slc);

    E_Nur = RecoE-TrueE;

  }
  return E_Nur;

});

const Var varTempNuMINeutrinoCombinedEnergyResidualFraction([](const caf::SRSliceProxy* slc) -> double {

  double E_Nurf(-999);
  if( varTempNuMIMuonTrackInd(slc) >= 0 && varTempNuMIProtonTrackInd(slc) >= 0 ){

    double RecoE = varTempNuMINeutrinoCombinedEnergy(slc);
    double TrueE = varNeutrinoTruthE(slc);

    E_Nurf = (RecoE-TrueE)/TrueE;

  }
  return E_Nurf;

});

//==== https://s3.cern.ch/inspire-prod-files-9/93642a13c46438d97680971700e2013c
const Var varTempNuMINeutrinoQE([](const caf::SRSliceProxy* slc) -> double {

  int muonTrackIndex = varMuonTrackIndex(slc);
  if(muonTrackIndex>=0){
    auto const& trk_Muon = slc->reco.trk.at(muonTrackIndex);

    double mu_costh = trk_Muon.costh;
    double mu_p = varTempNuMIMuonRecoP(slc);
    double mu_E = sqrt(mu_p*mu_p+M_MUON*M_MUON);

    double EQE_num = M_PROTON*M_PROTON - (M_NEUTRON-E_EffNuclB)*(M_NEUTRON-E_EffNuclB) - M_MUON*M_MUON + 2.*(M_NEUTRON-E_EffNuclB)*mu_E;
    double EQE_den = 2.*(M_NEUTRON - E_EffNuclB - mu_E + mu_p * mu_costh);

    return EQE_num/EQE_den;

  }
  else{
    return -999.;
  }

});

const Var varTempNuMINeutrinoQEResidual([](const caf::SRSliceProxy* slc) -> double {

  int muonTrackIndex = varMuonTrackIndex(slc);
  if(muonTrackIndex>=0){

    double RecoE = varTempNuMINeutrinoQE(slc);
    double TrueE = varNeutrinoTruthE(slc);

    double E_Nur = RecoE-TrueE;

    return E_Nur;

  }
  else{
    return -999.;
  }

});

const Var varTempNuMINeutrinoQEResidualFraction([](const caf::SRSliceProxy* slc) -> double {

  int muonTrackIndex = varMuonTrackIndex(slc);
  if(muonTrackIndex>=0){

    double RecoE = varTempNuMINeutrinoQE(slc);
    double TrueE = varNeutrinoTruthE(slc);

    double E_Nurf = (RecoE-TrueE)/TrueE;

    return E_Nurf;

  }
  else{
    return -999.;
  }

});



#endif
