// Defining bins, some useful functions and quantities, and signal definitions...

#include "sbnana/CAFAna/Core/Spectrum.h"
#include "sbnana/CAFAna/Core/Binning.h"
#include "sbnana/CAFAna/Core/Var.h"

#include "sbnanaobj/StandardRecord/SRTrueInteraction.h"

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

using namespace ana;

////////////////////////////////////////////////////////////////
// NuMI direction vector (https://sbn-docdb.fnal.gov/cgi-bin/sso/RetrieveFile?docid=22998&filename=main_v2.pdf&version=2)
////////////////////////////////////////////////////////////////
// First attempts at NuMI nu direction vectors
TVector3 dFromNuMI(315.120380, 33.644912, 733.632532);
double magNuMI = sqrt(315.120380*315.120380 + 33.644912*33.644912 + 733.632532*733.632532);
TVector3 rFromNuMI(315.120380/magNuMI, 33.644912/magNuMI, 733.632532/magNuMI);

////////////////////////////////////////////////////////////////
// Binning
////////////////////////////////////////////////////////////////
const Binning kBinsRes = Binning::Simple(22,-1.1,1.1);
const Binning kBinsResZoom = Binning::Simple(40,-1.,1.);

const Binning kBinsDist = Binning::Simple(5,0.,20.);
const Binning kBinsDistZoom = Binning::Simple(10,0.,20.);

const Binning kBinsPosResZoom = Binning::Simple(60,-15.,15.);
const Binning kBinsPosResZoomZoom = Binning::Simple(120,-6.,6.);

const Binning kBinsFrac = Binning::Simple(11,0.,1.1);
const Binning kBinsFracZoom = Binning::Simple(20,0.,1.);

const Binning kBinsPosX = Binning::Simple(80,-400.,400.);
const Binning kBinsPosY = Binning::Simple(35,-200.,150.);
const Binning kBinsPosZ = Binning::Simple(90,-900.,900.);

const Binning kBinsAnodeDist = Binning::Simple(15,0.,150.);

const Binning kBinsCosTh = Binning::Simple(20,-1.,1.);
const Binning kBinsP = Binning::Simple(35,0.,3.5);
const Binning kBinsProtonP = Binning::Simple(16,0,1.6);

////////////////////////////////////////////////////////////////
// Basic checks:
////////////////////////////////////////////////////////////////
bool isInFV (double x, double y, double z)
{
  if ( std::isnan(x) || std::isnan(y) || std::isnan(z) ) return false;

  return (( ( x < -61.94 - 25 && x > -358.49 + 25 ) ||
						( x >  61.94 + 25 && x <  358.49 - 25 )) &&
					( ( y > -181.86 + 25 && y < 134.96 - 25 ) &&
						( z > -894.95 + 30 && z < 894.95 - 50 ) ));
}

bool isContainedVol (double x, double y, double z)
{
  if ( std::isnan(x) || std::isnan(y) || std::isnan(z) ) return false;

  return (( ( x < -61.94 - 10. && x > -358.49 + 10. ) ||
						( x >  61.94 + 10. && x <  358.49 - 10. )) &&
					( ( y > -181.86 + 10. && y < 134.96 - 10. ) &&
						( z > -894.95 + 10. && z < 894.95 - 10. ) ));
}

////////////////////////////////////////////////////////////////
// Signal definitions:
////////////////////////////////////////////////////////////////

const Cut kIsSlcNuNC([](const caf::SRSliceProxy* slc) {
    if ( slc->truth.index < 0 ) return false;

    if ( slc->truth.iscc )
      return false; // IS CC

    return true;
  });

const Cut kIsSlcNotNu([](const caf::SRSliceProxy* slc) {
    return ( slc->truth.index < 0 );
  });

// CC: No kinematic threshold
const Cut kIsSlcSignal_1muNp0pi_NoContainment([](const caf::SRSliceProxy* slc) {
    if ( slc->truth.index < 0 ) return false;

    if ( abs(slc->truth.pdg) != 14 ||
				 !slc->truth.iscc ||
				 std::isnan(slc->truth.position.x) || std::isnan(slc->truth.position.y) || std::isnan(slc->truth.position.z) ||
				 !isInFV(slc->truth.position.x,slc->truth.position.y,slc->truth.position.z) )
      return false; // not signal

    // primaries:
    unsigned int nMu(0), nP(0), nPi(0);
    for ( auto const& prim : slc->truth.prim ) {
      if ( prim.start_process != 0 ) continue;

      double momentum = sqrt( (prim.genp.x*prim.genp.x) + (prim.genp.y*prim.genp.y) + (prim.genp.z*prim.genp.z) );
      if ( abs(prim.pdg) == 13 && momentum > 0.226 ) nMu+=1;
      if ( abs(prim.pdg) == 2212 && isContainedVol(prim.end.x,prim.end.y,prim.end.z) ) nP+=1;
      if ( abs(prim.pdg) == 111 || abs(prim.pdg) == 211 ) nPi+=1;
    }
    if ( nMu!=1 || nP==0 || nPi > 0 ) return false;

    return true;
  });

const Cut kIsSlcOtherNuCC_1muNp0pi_NoContainment([](const caf::SRSliceProxy* slc) {
    if ( slc->truth.index < 0 ) return false; // not Nu
    if ( !slc->truth.iscc ) return false; // not CC

    if ( kIsSlcSignal_1muNp0pi_NoContainment(slc) ) return false; // covered by signal
    return true;
  });

// CC: With kinematic threshold of 400MeV/c to 1 GeV/c on proton
const Cut kIsSlcSignal_1muNp0pi_NoContainment_ProtonKinematicThreshold([](const caf::SRSliceProxy* slc) {
    if ( slc->truth.index < 0 ) return false;

    if ( abs(slc->truth.pdg) != 14 ||
				 !slc->truth.iscc ||
				 std::isnan(slc->truth.position.x) || std::isnan(slc->truth.position.y) || std::isnan(slc->truth.position.z) ||
				 !isInFV(slc->truth.position.x,slc->truth.position.y,slc->truth.position.z) )
      return false; // not signal

    // primaries:
    unsigned int nMu(0), nP(0), nPi(0);
    for ( auto const& prim : slc->truth.prim ) {
      if ( prim.start_process != 0 ) continue;

      double momentum = sqrt( (prim.genp.x*prim.genp.x) + (prim.genp.y*prim.genp.y) + (prim.genp.z*prim.genp.z) );
      if ( abs(prim.pdg) == 13 && momentum > 0.226 ) nMu+=1;
      if ( abs(prim.pdg) == 2212 && isContainedVol(prim.end.x,prim.end.y,prim.end.z) && (momentum > 0.4 && momentum < 1.) ) nP+=1;
      if ( abs(prim.pdg) == 111 || abs(prim.pdg) == 211 ) nPi+=1;
    }
    if ( nMu!=1 || nP==0 || nPi > 0 ) return false;

    return true;
  });

const Cut kIsSlcOtherNuCC_1muNp0pi_NoContainment_ProtonKinematicTreshold([](const caf::SRSliceProxy* slc) {
    if ( slc->truth.index < 0 ) return false; // not Nu
    if ( !slc->truth.iscc ) return false; // not CC

    if ( kIsSlcSignal_1muNp0pi_NoContainment_ProtonKinematicThreshold(slc) ) return false; // covered by signal
    return true;
  });
