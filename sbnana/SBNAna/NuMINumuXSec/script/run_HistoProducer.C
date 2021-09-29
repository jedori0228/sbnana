#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_HistoProducer.h"
#include "sbnana/SBNAna/Cuts/NumuCutsIcarus202106.h"
#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_Variables.h"
#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_Cuts.h"
#include "sbnana/SBNAna/Cuts/Cuts.h"

#include "myFilelist.h"
#include "myFilelistReCAF.h"
#include "myFilelistShort.h"
#include "myFilelistBNB.h"
#include "myFilelistNuMI.h"


using namespace ana;
using namespace ICARUSNumuXsec;

void run_HistoProducer(int whichSample, TString outDirBase){

  //==== Define Cuts

  vector<Cut> cuts;
  vector<SpillCut> spillcuts;
  vector<TString> cutNames;

  //====   1) CC

  cuts.push_back( kNoCut && kIsNuMuCC );
  spillcuts.push_back( kNoSpillCut );
  cutNames.push_back( "NuMuCC" );

  cuts.push_back( kNoCut && kIsNuMuCC && cutRFiducial );
  spillcuts.push_back( kNoSpillCut );
  cutNames.push_back( "NuMuCC_RecoFiducial" );

  cuts.push_back( kNoCut && kIsNuMuCC && cutRFiducial && kNotClearCosmic );
  spillcuts.push_back( kNoSpillCut );
  cutNames.push_back( "NuMuCC_RecoFiducial_NotClearCosmic" );

  cuts.push_back( kNoCut && kIsNuMuCC && cutRFiducial && kNotClearCosmic );
  spillcuts.push_back( kCRTHitVetoFD );
  cutNames.push_back( "NuMuCC_RecoFiducial_NotClearCosmic_CRTVeto" );

  cuts.push_back( kNoCut && kIsNuMuCC && cutRFiducial && kNotClearCosmic );
  spillcuts.push_back( kCRTHitVetoFD );
  cutNames.push_back( "NuMuCC_RecoFiducial_NotClearCosmic_CRTVeto" );

  cuts.push_back( kNoCut && kIsNuMuCC && cutRFiducial && kNotClearCosmic && cutFMScore );
  spillcuts.push_back( kCRTHitVetoFD );
  cutNames.push_back( "NuMuCC_RecoFiducial_NotClearCosmic_CRTVeto_FMScoreLT6" );

  cuts.push_back( kNoCut && kIsNuMuCC && cutRFiducial && kNotClearCosmic && cutFMScore && kNuScore );
  spillcuts.push_back( kCRTHitVetoFD );
  cutNames.push_back( "NuMuCC_RecoFiducial_NotClearCosmic_CRTVeto_FMScoreLT6_NuScoreGT0p4" );

  cuts.push_back( kNoCut && kIsNuMuCC && cutRFiducial && kNotClearCosmic && cutFMScore && kNuScore && cutHasMuon );
  spillcuts.push_back( kCRTHitVetoFD );
  cutNames.push_back( "NuMuCC_RecoFiducial_NotClearCosmic_CRTVeto_FMScoreLT6_NuScoreGT0p4_HasMuon" );

  //====     Contained/Exiting

  cuts.push_back( kNoCut && kIsNuMuCC && cutRFiducial && kNotClearCosmic && cutFMScore && kNuScore && cutHasMuon && cutMuonContained );
  spillcuts.push_back( kCRTHitVetoFD );
  cutNames.push_back( "NuMuCC_RecoFiducial_NotClearCosmic_CRTVeto_FMScoreLT6_NuScoreGT0p4_HasMuon_Contained" );

  cuts.push_back( kNoCut && kIsNuMuCC && cutRFiducial && kNotClearCosmic && cutFMScore && kNuScore && cutHasMuon && !cutMuonContained);
  spillcuts.push_back( kCRTHitVetoFD );
  cutNames.push_back( "NuMuCC_RecoFiducial_NotClearCosmic_CRTVeto_FMScoreLT6_NuScoreGT0p4_HasMuon_Exiting" );

  //====       Matched to muon/proton

  cuts.push_back( kNoCut && kIsNuMuCC && cutRFiducial && kNotClearCosmic && cutFMScore && kNuScore && cutHasMuon && cutMuonContained && cutMuonMatchedToMuon );
  spillcuts.push_back( kCRTHitVetoFD );
  cutNames.push_back( "NuMuCC_RecoFiducial_NotClearCosmic_CRTVeto_FMScoreLT6_NuScoreGT0p4_HasMuon_Contained_MatchedToMuon" );

  cuts.push_back( kNoCut && kIsNuMuCC && cutRFiducial && kNotClearCosmic && cutFMScore && kNuScore && cutHasMuon && cutMuonContained && cutMuonMatchedToProton );
  spillcuts.push_back( kCRTHitVetoFD );
  cutNames.push_back( "NuMuCC_RecoFiducial_NotClearCosmic_CRTVeto_FMScoreLT6_NuScoreGT0p4_HasMuon_Contained_MatchedToProton" );

  cuts.push_back( kNoCut && kIsNuMuCC && cutRFiducial && kNotClearCosmic && cutFMScore && kNuScore && cutHasMuon && !cutMuonContained && cutMuonMatchedToMuon );
  spillcuts.push_back( kCRTHitVetoFD );
  cutNames.push_back( "NuMuCC_RecoFiducial_NotClearCosmic_CRTVeto_FMScoreLT6_NuScoreGT0p4_HasMuon_Exiting_MatchedToMuon" );

  cuts.push_back( kNoCut && kIsNuMuCC && cutRFiducial && kNotClearCosmic && cutFMScore && kNuScore && cutHasMuon && !cutMuonContained && cutMuonMatchedToProton );
  spillcuts.push_back( kCRTHitVetoFD );
  cutNames.push_back( "NuMuCC_RecoFiducial_NotClearCosmic_CRTVeto_FMScoreLT6_NuScoreGT0p4_HasMuon_Exiting_MatchedToProton" );

  //====     Additional particles

  cuts.push_back( kNoCut && kIsNuMuCC && cutRFiducial && kNotClearCosmic && cutFMScore && kNuScore && cutHasMuon && !cutTruthNoNeutron );
  spillcuts.push_back( kCRTHitVetoFD );
  cutNames.push_back( "NuMuCC_RecoFiducial_NotClearCosmic_CRTVeto_FMScoreLT6_NuScoreGT0p4_HasMuon_WithNeutron" );

  cuts.push_back( kNoCut && kIsNuMuCC && cutRFiducial && kNotClearCosmic && cutFMScore && kNuScore && cutHasMuon && !cutTruthNoChargedPion );
  spillcuts.push_back( kCRTHitVetoFD );
  cutNames.push_back( "NuMuCC_RecoFiducial_NotClearCosmic_CRTVeto_FMScoreLT6_NuScoreGT0p4_HasMuon_WithChargedPion" );

  cuts.push_back( kNoCut && kIsNuMuCC && cutRFiducial && kNotClearCosmic && cutFMScore && kNuScore && cutHasMuon && !cutTruthNoPiZero );
  spillcuts.push_back( kCRTHitVetoFD );
  cutNames.push_back( "NuMuCC_RecoFiducial_NotClearCosmic_CRTVeto_FMScoreLT6_NuScoreGT0p4_HasMuon_WithNeutralPion" );

  cuts.push_back( kNoCut && kIsNuMuCC && cutRFiducial && kNotClearCosmic && cutFMScore && kNuScore && cutHasMuon && cutTruthNoNeutron && cutTruthNoChargedPion && cutTruthNoPiZero );
  spillcuts.push_back( kCRTHitVetoFD );
  cutNames.push_back( "NuMuCC_RecoFiducial_NotClearCosmic_CRTVeto_FMScoreLT6_NuScoreGT0p4_HasMuon_WithoutNeutronOrPion" );

  //====     HasProton added

  cuts.push_back( kNoCut && kIsNuMuCC && cutRFiducial && kNotClearCosmic && cutFMScore && kNuScore && cutHasMuon && cutHasProton );
  spillcuts.push_back( kCRTHitVetoFD );
  cutNames.push_back( "NuMuCC_RecoFiducial_NotClearCosmic_CRTVeto_FMScoreLT6_NuScoreGT0p4_HasMuon_HasProton" );

  //====   2) NC check

  cuts.push_back( kNoCut && cutIsNuMuNC );
  spillcuts.push_back( kNoSpillCut );
  cutNames.push_back( "NuMuNC" );


  //==== Inputfile

  std::vector<std::string> vec_inputFiles;
  vec_inputFiles.clear();

  //====   BNB, reCAF
  if(whichSample==0){
    cout << "[run_Efficiency] Running BNB, my custom reCAF" << endl;
    vec_inputFiles = inputFilesReCAF;
  }
  //====   BNB, new central sample, IcarusProd2021B_BNB_Nu_Cosmics_v09_28_01_01_01_caf
  else if(whichSample==1){
    cout << "[run_Efficiency] Running BNB, IcarusProd2021B_BNB_Nu_Cosmics_v09_28_01_01_01_caf" << endl;
    vec_inputFiles = inputFilesBNB;
  }
  //====   NuMI, new central sample, IcarusProd2021B_NuMI_Nu_Cosmics_v09_28_01_01_01_caf
  else if(whichSample==2){
    cout << "[run_Efficiency] Running NuMI, IcarusProd2021B_NuMI_Nu_Cosmics_v09_28_01_01_01_caf" << endl;
    vec_inputFiles = inputFilesNuMI;
  }

  SpectrumLoader loader(vec_inputFiles);

  HistoProducer m;
  m.outputDir = outDirBase;
  m.initialize();
  //m.setSystematicWeights();

  for(unsigned int i=0; i<cuts.size(); i++){

    Cut myCut = cuts.at(i);
    SpillCut mySpillCut = spillcuts.at(i);
    TString cutName = cutNames.at(i);

    if( m.setCut(cutName) ){

      m.bookSlice(loader, myCut, mySpillCut);
      m.bookTruth(loader, myCut, mySpillCut);
      if(cutName.Contains("HasMuon")){
        m.bookRecoMuon(loader, myCut, mySpillCut);
      }
      if(cutName.Contains("HasProton")){
        m.bookRecoProton(loader, myCut, mySpillCut);
        m.bookRecoNeutrino(loader, myCut, mySpillCut);
      }

    }

  }

  loader.Go();

  m.saveHistograms();

}
