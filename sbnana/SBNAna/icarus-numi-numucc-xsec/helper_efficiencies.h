// Variables for efficiency calculations...

#include "sbnana/CAFAna/Core/Spectrum.h"
#include "sbnana/CAFAna/Core/Binning.h"
#include "sbnana/CAFAna/Core/Var.h"

#include "sbnanaobj/StandardRecord/SRTrueInteraction.h"
#include "sbnana/SBNAna/Vars/Vars.h"

#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

// ROOT
#include "TMath.h"
#include "TVector3.h"
#include "TGraphAsymmErrors.h"
#include "TRandom3.h"
#include "TPaveText.h"
#include "TColor.h"
#include "TLegend.h"
#include "TH1D.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TGraph.h"

// C++
#include <vector>
#include <string>
#include <map>
#include <utility>
#include <cassert>

// Other numi xsec includes
#include "sbnana/SBNAna/icarus-numi-numucc-xsec/helper_selvars.h"

using namespace ana;

// Functions to fetch the truth
std::map<int, double> GetTrueNuMuCC_1muNp0pi_IndexMuonPs( const caf::SRSpillProxy& sr, const bool containedLepton,
                                                          const bool containedProton=true, double muThresh=-1., double pThresh=-1., 
                                                          double pMaxThresh=9999., double piThresh=-1., double pi0Thresh=-1. ) {
  std::map<int, double> indexPs;

  for ( auto const& nu : sr.mc.nu ) {
    if ( abs(nu.pdg) != 14 ||
				 !nu.iscc ||
				 std::isnan(nu.position.x) || std::isnan(nu.position.y) || std::isnan(nu.position.z) ||
				 !isInFV(nu.position.x,nu.position.y,nu.position.z) )
      continue; // not signal

    unsigned int nMu(0), nP(0), nPi0(0), nChgPi(0);
    float muonP(0.);
    for ( auto const& prim : nu.prim ) {
      if ( prim.start_process != 0 ) continue; // to handle fact that primaries now include more particles.
      double momentum = sqrt( (prim.startp.x*prim.startp.x) + (prim.startp.y*prim.startp.y) + (prim.startp.z*prim.startp.z) );
      if ( abs(prim.pdg) == 13 ) {
        if ( containedLepton && !prim.contained ) continue; // requesting signal muon to be contained and it's not.
        // Note: 50cm ~ 226 MeV/c for muons
        if ( momentum < muThresh ) continue;
        if ( momentum > muonP ) muonP = momentum;
        nMu+=1;
      }
      if ( abs(prim.pdg) == 2212 ){
        if ( containedProton && !prim.contained ) continue;
        if ( momentum < pThresh || momentum > pMaxThresh ) continue;
        nP+=1;
      }
      if ( abs(prim.pdg) == 111 ){
        if ( momentum < pi0Thresh ) continue;
        nPi0+=1;
      }
      if ( abs(prim.pdg) == 211 ){
        if ( momentum < piThresh ) continue;
        nChgPi+=1;
      }
    }
    if ( nMu==1 && nP>=1 && nChgPi==0 && nPi0==0 ) {
      indexPs[ nu.index ] = muonP;
    }
  }

  return indexPs;
}

std::map<int, double> GetTrueNuMuCC_1muNp0pi_IndexProtonPs( const caf::SRSpillProxy& sr, const bool containedLepton,
                                                            const bool containedProton=true, double muThresh=-1., double pThresh=-1., 
                                                            double pMaxThresh=9999., double piThresh=-1., double pi0Thresh=-1. ) {
  std::map<int, double> indexPs;

  for ( auto const& nu : sr.mc.nu ) {
    if ( abs(nu.pdg) != 14 ||
				 !nu.iscc ||
				 std::isnan(nu.position.x) || std::isnan(nu.position.y) || std::isnan(nu.position.z) ||
				 !isInFV(nu.position.x,nu.position.y,nu.position.z) )
      continue; // not signal

    unsigned int nMu(0), nP(0), nPi0(0), nChgPi(0);
    float protonP(0.);
    for ( auto const& prim : nu.prim ) {
      if ( prim.start_process != 0 ) continue; // to handle fact that primaries now include more particles.
      double momentum = sqrt( (prim.startp.x*prim.startp.x) + (prim.startp.y*prim.startp.y) + (prim.startp.z*prim.startp.z) );
      if ( abs(prim.pdg) == 13 ) {
        if ( containedLepton && !prim.contained ) continue; // requesting signal muon to be contained and it's not.
        // Note: 50cm ~ 226 MeV/c for muons
        if ( momentum < muThresh ) continue;
        nMu+=1;
      }
      if ( abs(prim.pdg) == 2212 ){
        if ( containedProton && !prim.contained ) continue;
        if ( momentum < pThresh || momentum > pMaxThresh ) continue;
        if ( momentum > protonP ) protonP = momentum;
        nP+=1;
      }
      if ( abs(prim.pdg) == 111 ){
        if ( momentum < pi0Thresh ) continue;
        nPi0+=1;
      }
      if ( abs(prim.pdg) == 211 ){
        if ( momentum < piThresh ) continue;
        nChgPi+=1;
      }
    }
    if ( nMu==1 && nP>=1 && nChgPi==0 && nPi0==0 ) {
      indexPs[ nu.index ] = protonP;
    }
  }

  return indexPs;
}

std::map<int, double> GetTrueNuMuCC_1muNp0pi_IndexMuonCosThNuMI( const caf::SRSpillProxy& sr, const bool containedLepton,
                                                                 const bool containedProton=true, double muThresh=-1., double pThresh=-1., 
                                                                 double pMaxThresh=9999., double piThresh=-1., double pi0Thresh=-1. ) {
  std::map<int, double> indexPs;

  for ( auto const& nu : sr.mc.nu ) {
    if ( abs(nu.pdg) != 14 ||
				 !nu.iscc ||
				 std::isnan(nu.position.x) || std::isnan(nu.position.y) || std::isnan(nu.position.z) ||
				 !isInFV(nu.position.x,nu.position.y,nu.position.z) )
      continue; // not signal

    unsigned int nMu(0), nP(0), nPi0(0), nChgPi(0);
    float muonP(0.);
    float muonCosTh(0.);
    for ( auto const& prim : nu.prim ) {
      if ( prim.start_process != 0 ) continue; // to handle fact that primaries now include more particles.
      double momentum = sqrt( (prim.startp.x*prim.startp.x) + (prim.startp.y*prim.startp.y) + (prim.startp.z*prim.startp.z) );
      if ( abs(prim.pdg) == 13 ) {
        if ( containedLepton && !prim.contained ) continue; // requesting signal muon to be contained and it's not.
        // Note: 50cm ~ 226 MeV/c for muons
        if ( momentum < muThresh ) continue;
        if ( momentum > muonP ) {
          TVector3 muDir(prim.end.x - prim.start.x, prim.end.y - prim.start.y, prim.end.z - prim.start.z);
          muDir = muDir.Unit();
          muonCosTh = TMath::Cos( muDir.Angle(rFromNuMI) );

          muonP = momentum;
        }
        nMu+=1;
      }
      if ( abs(prim.pdg) == 2212 ){
        if ( containedProton && !prim.contained ) continue;
        if ( momentum < pThresh || momentum > pMaxThresh ) continue;
        nP+=1;
      }
      if ( abs(prim.pdg) == 111 ){
        if ( momentum < pi0Thresh ) continue;
        nPi0+=1;
      }
      if ( abs(prim.pdg) == 211 ){
        if ( momentum < piThresh ) continue;
        nChgPi+=1;
      }
    }
    if ( nMu==1 && nP>=1 && nChgPi==0 && nPi0==0 ) {
      indexPs[ nu.index ] = muonCosTh;
    }
  }

  return indexPs;
}

std::map<int, double> GetTrueNuMuCC_1muNp0pi_IndexMuonCosThMuP( const caf::SRSpillProxy& sr, const bool containedLepton,
                                                                const bool containedProton=true, double muThresh=-1., double pThresh=-1., 
                                                                double pMaxThresh=9999., double piThresh=-1., double pi0Thresh=-1. ) {
  std::map<int, double> indexPs;

  for ( auto const& nu : sr.mc.nu ) {
    if ( abs(nu.pdg) != 14 ||
				 !nu.iscc ||
				 std::isnan(nu.position.x) || std::isnan(nu.position.y) || std::isnan(nu.position.z) ||
				 !isInFV(nu.position.x,nu.position.y,nu.position.z) )
      continue; // not signal

    unsigned int nMu(0), nP(0), nPi0(0), nChgPi(0);
    float muonP(0.);
    float protonP(0.);
    float muonX(0.), muonY(0.), muonZ(0.);
    float protonX(0.), protonY(0.), protonZ(0.);
    for ( auto const& prim : nu.prim ) {
      if ( prim.start_process != 0 ) continue; // to handle fact that primaries now include more particles.
      double momentum = sqrt( (prim.startp.x*prim.startp.x) + (prim.startp.y*prim.startp.y) + (prim.startp.z*prim.startp.z) );
      if ( abs(prim.pdg) == 13 ) {
        if ( containedLepton && !prim.contained ) continue; // requesting signal muon to be contained and it's not.
        // Note: 50cm ~ 226 MeV/c for muons
        if ( momentum < muThresh ) continue;
        if ( momentum > muonP ) {
          TVector3 muDir(prim.end.x - prim.start.x, prim.end.y - prim.start.y, prim.end.z - prim.start.z);
          muDir = muDir.Unit();
          muonX = muDir.X(); muonY = muDir.Y(); muonZ = muDir.Z();

          muonP = momentum;
        }
        nMu+=1;
      }
      if ( abs(prim.pdg) == 2212 ){
        if ( containedProton && !prim.contained ) continue;
        if ( momentum < pThresh || momentum > pMaxThresh ) continue;
        if ( momentum > protonP ) {
          TVector3 pDir(prim.end.x - prim.start.x, prim.end.y - prim.start.y, prim.end.z - prim.start.z);
          pDir = pDir.Unit();
          protonX = pDir.X(); protonY = pDir.Y(); protonZ = pDir.Z();

          protonP = momentum;
        }
        nP+=1;
      }
      if ( abs(prim.pdg) == 111 ){
        if ( momentum < pi0Thresh ) continue;
        nPi0+=1;
      }
      if ( abs(prim.pdg) == 211 ){
        if ( momentum < piThresh ) continue;
        nChgPi+=1;
      }
    }
    if ( nMu==1 && nP>=1 && nChgPi==0 && nPi0==0 ) {
      TVector3 muDir(muonX, muonY, muonZ);
      TVector3 pDir(protonX, protonY, protonZ);
      indexPs[ nu.index ] = TMath::Cos( muDir.Angle(pDir) );
    }
  }

  return indexPs;
}

// Implementation
const SpillMultiVar k1muNp0pi_TrueMuonPs( [](const caf::SRSpillProxy *sr) {
  std::vector<double> vals;

  std::map<int, double> indexVals = GetTrueNuMuCC_1muNp0pi_IndexMuonPs(*sr, false, true, 0.226, 0.400, 1., 0., 0. );

  if ( indexVals.size() == 0 ) return vals;

	for ( auto const &[index, thisval] : indexVals ) {
    vals.push_back( thisval );
  }

	return vals;
});

const SpillMultiVar k1muNp0pi_TrueMuonPs_Cuts( [](const caf::SRSpillProxy *sr) {
  std::vector<double> vals;

  std::map<int, double> indexVals = GetTrueNuMuCC_1muNp0pi_IndexMuonPs(*sr, false, true, 0.226, 0.400, 1., 0., 0. );

  if ( indexVals.size() == 0 ) return vals;

  std::map<int, double> selectedVals;

  for ( auto const& slc : sr->slc ) {
    if ( slc.truth.index < 0 ) continue;
    else if ( indexVals.find( slc.truth.index ) == indexVals.end() ) continue;

    if ( kNuMISelection_1muNp0pi(&slc) ) selectedVals[ slc.truth.index ] = indexVals[ slc.truth.index ];
  }

	for ( auto const &[index, thisval] : selectedVals ) {
    vals.push_back( thisval );
  }

	return vals;
});

const SpillMultiVar k1muNp0pi_TrueMuonCosThNuMI( [](const caf::SRSpillProxy *sr) {
  std::vector<double> vals;

  std::map<int, double> indexVals = GetTrueNuMuCC_1muNp0pi_IndexMuonCosThNuMI(*sr, false, true, 0.226, 0.400, 1., 0., 0. );

  if ( indexVals.size() == 0 ) return vals;

	for ( auto const &[index, thisval] : indexVals ) {
    vals.push_back( thisval );
  }

	return vals;
});

const SpillMultiVar k1muNp0pi_TrueMuonCosThNuMI_Cuts( [](const caf::SRSpillProxy *sr) {
  std::vector<double> vals;

  std::map<int, double> indexVals = GetTrueNuMuCC_1muNp0pi_IndexMuonCosThNuMI(*sr, false, true, 0.226, 0.400, 1., 0., 0. );

  if ( indexVals.size() == 0 ) return vals;

  std::map<int, double> selectedVals;

  for ( auto const& slc : sr->slc ) {
    if ( slc.truth.index < 0 ) continue;
    else if ( indexVals.find( slc.truth.index ) == indexVals.end() ) continue;

    if ( kNuMISelection_1muNp0pi(&slc) ) selectedVals[ slc.truth.index ] = indexVals[ slc.truth.index ];
  }

	for ( auto const &[index, thisval] : selectedVals ) {
    vals.push_back( thisval );
  }

	return vals;
});

const SpillMultiVar k1muNp0pi_TrueMuonCosThMuP( [](const caf::SRSpillProxy *sr) {
  std::vector<double> vals;

  std::map<int, double> indexVals = GetTrueNuMuCC_1muNp0pi_IndexMuonCosThMuP(*sr, false, true, 0.226, 0.400, 1., 0., 0. );

  if ( indexVals.size() == 0 ) return vals;

	for ( auto const &[index, thisval] : indexVals ) {
    vals.push_back( thisval );
  }

	return vals;
});

const SpillMultiVar k1muNp0pi_TrueMuonCosThMuP_Cuts( [](const caf::SRSpillProxy *sr) {
  std::vector<double> vals;

  std::map<int, double> indexVals = GetTrueNuMuCC_1muNp0pi_IndexMuonCosThMuP(*sr, false, true, 0.226, 0.400, 1., 0., 0. );

  if ( indexVals.size() == 0 ) return vals;

  std::map<int, double> selectedVals;

  for ( auto const& slc : sr->slc ) {
    if ( slc.truth.index < 0 ) continue;
    else if ( indexVals.find( slc.truth.index ) == indexVals.end() ) continue;

    if ( kNuMISelection_1muNp0pi(&slc) ) selectedVals[ slc.truth.index ] = indexVals[ slc.truth.index ];
  }

	for ( auto const &[index, thisval] : selectedVals ) {
    vals.push_back( thisval );
  }

	return vals;
});


// Proton momenta with no kinematic thresholds

const SpillMultiVar k1muNp0pi_TrueProtonPs( [](const caf::SRSpillProxy *sr) {
  std::vector<double> vals;

  std::map<int, double> indexVals = GetTrueNuMuCC_1muNp0pi_IndexProtonPs(*sr, false, true, 0.226, 0., 9999., 0., 0. );

  if ( indexVals.size() == 0 ) return vals;

	for ( auto const &[index, thisval] : indexVals ) {
    vals.push_back( thisval );
  }

	return vals;
});

const SpillMultiVar k1muNp0pi_TrueProtonPs_Cuts( [](const caf::SRSpillProxy *sr) {
  std::vector<double> vals;

  std::map<int, double> indexVals = GetTrueNuMuCC_1muNp0pi_IndexProtonPs(*sr, false, true, 0.226, 0., 9999., 0., 0. );

  if ( indexVals.size() == 0 ) return vals;

  std::map<int, double> selectedVals;

  for ( auto const& slc : sr->slc ) {
    if ( slc.truth.index < 0 ) continue;
    else if ( indexVals.find( slc.truth.index ) == indexVals.end() ) continue;

    if ( kNuMISelection_1muNp0pi(&slc) ) selectedVals[ slc.truth.index ] = indexVals[ slc.truth.index ];
  }

	for ( auto const &[index, thisval] : selectedVals ) {
    vals.push_back( thisval );
  }

	return vals;
});
