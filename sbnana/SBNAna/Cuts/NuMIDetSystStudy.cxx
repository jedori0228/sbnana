#include "sbnana/SBNAna/Cuts/NuMIDetSystStudy.h"
#include "sbnana/SBNAna/Vars/NuMIXSecVars.h"
#include "sbnana/SBNAna/Cuts/NuMIXSecCuts.h"

#include <iostream>

namespace ana {

  const MultiVar Hit_dedx([](const caf::SRSliceProxy* slc) -> std::vector<double> {

    std::vector<double> rets;

    if( !kNuMI_1muNp0piStudy_Signal_WithPhaseSpaceCut(slc) ) return rets;

    int muonIndex = kNuMIMuonCandidateIdx(slc);
    if(muonIndex<0) return rets;

    const auto& mu_pfp = slc->reco.pfp.at(muonIndex);
    const auto& mu_trk = mu_pfp.trk;
    for(const auto& pt: mu_trk.calo[2].points){
      rets.push_back(pt.dedx);
    }
    return rets;
  });
  const MultiVar Hit_dqdx([](const caf::SRSliceProxy* slc) -> std::vector<double> {
    
    std::vector<double> rets;
    
    if( !kNuMI_1muNp0piStudy_Signal_WithPhaseSpaceCut(slc) ) return rets;
    
    int muonIndex = kNuMIMuonCandidateIdx(slc);
    if(muonIndex<0) return rets;
    
    const auto& mu_pfp = slc->reco.pfp.at(muonIndex);
    const auto& mu_trk = mu_pfp.trk;
    for(const auto& pt: mu_trk.calo[2].points){
      rets.push_back(pt.dqdx);
    }
    return rets;
  });
  const MultiVar Hit_integral([](const caf::SRSliceProxy* slc) -> std::vector<double> {

    std::vector<double> rets;

    if( !kNuMI_1muNp0piStudy_Signal_WithPhaseSpaceCut(slc) ) return rets;

    int muonIndex = kNuMIMuonCandidateIdx(slc);
    if(muonIndex<0) return rets;

    const auto& mu_pfp = slc->reco.pfp.at(muonIndex);
    const auto& mu_trk = mu_pfp.trk;
    for(const auto& pt: mu_trk.calo[2].points){
      rets.push_back(pt.integral);
    }
    return rets;
  });
  const MultiVar Hit_sumadc([](const caf::SRSliceProxy* slc) -> std::vector<double> {

    std::vector<double> rets;
  
    if( !kNuMI_1muNp0piStudy_Signal_WithPhaseSpaceCut(slc) ) return rets;

    int muonIndex = kNuMIMuonCandidateIdx(slc);
    if(muonIndex<0) return rets;

    const auto& mu_pfp = slc->reco.pfp.at(muonIndex);
    const auto& mu_trk = mu_pfp.trk;
    for(const auto& pt: mu_trk.calo[2].points){
      rets.push_back(pt.sumadc);
    }
    return rets;
  });
  const MultiVar Hit_t([](const caf::SRSliceProxy* slc) -> std::vector<double> {

    std::vector<double> rets;

    if( !kNuMI_1muNp0piStudy_Signal_WithPhaseSpaceCut(slc) ) return rets;

    int muonIndex = kNuMIMuonCandidateIdx(slc);
    if(muonIndex<0) return rets;

    const auto& mu_pfp = slc->reco.pfp.at(muonIndex);
    const auto& mu_trk = mu_pfp.trk;
    for(const auto& pt: mu_trk.calo[2].points){
      rets.push_back(pt.t);
    }
    return rets;
  });

}
