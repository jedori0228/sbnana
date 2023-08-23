#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_TwoTrack_Efficiency.h"

using namespace ana;

namespace ICARUSNumuXsec{

namespace Efficiency{

  std::map<int, double> GetNuIndexTruthVarPair(
    const caf::SRSpillProxy* sr,
    std::function<bool(const TrueInteraction&)> isSignal,
    std::function<double(const TrueInteraction&) > trueth_var
  ){

    std::map<int, double> indexPs;

    for ( auto const& nu : sr->mc.nu ) {
      if(!isSignal(nu)) continue; // 1muNp0pi
      indexPs[ nu.index ] = trueth_var(nu);
    }

    return indexPs;

  }

  std::vector<double> GetVectorForEfficiency(
    const caf::SRSpillProxy* sr,
    std::function<bool(const TrueInteraction&)> isSignal,
    std::function<double(const TrueInteraction&) > trueth_var,
    const Cut& RecoCut,
    const SpillCut& RecoSpillCut
  ){

    std::vector<double> vals;

    // Applying SpillCut
    if( !RecoSpillCut(sr) ) return vals;

    // Get (nu index, variable) pair
    std::map<int, double> indexVals = GetNuIndexTruthVarPair(sr, isSignal, trueth_var);
    if ( indexVals.size() == 0 ) return vals;

    // Loop over slices that pass Cut
    std::map<int, double> selectedVals;
    for ( auto const& slc : sr->slc ) {
      if ( slc.truth.index < 0 ) continue;
      else if ( indexVals.find( slc.truth.index ) == indexVals.end() ) continue;

      if ( RecoCut(&slc) ) selectedVals[ slc.truth.index ] = indexVals[ slc.truth.index ];
    }

    for ( auto const &[index, thisval] : selectedVals ) {
      vals.push_back( thisval );
    }

    return vals;

  }

}

}

