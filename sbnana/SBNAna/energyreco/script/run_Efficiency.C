#include "sbnana/SBNAna/energyreco/myEfficiency.h"
#include "sbnana/SBNAna/Cuts/NumuCutsIcarus202106.h"
#include "sbnana/SBNAna/energyreco/myTempNuMIVariables.h"
#include "sbnana/SBNAna/energyreco/myTool.h"

using namespace ana;

void run_Efficiency(int whichCut, int whichSample, TString outDirBase){

  vector<Cut> cuts;
  vector<SpillCut> spillcuts;
  vector<TString> outDirNames;

  if(whichCut==0){
    cuts.push_back( kNoCut && kIsNuMuCC );
    spillcuts.push_back( kNoSpillCut );
    outDirNames.push_back( "NuMuCC" );
  }
  else if(whichCut==1){
    cuts.push_back( kNoCut && kIsNuMuCC && cutmyRFiducial );
    spillcuts.push_back( kNoSpillCut );
    outDirNames.push_back( "NuMuCC_RecoFiducial" );
  }
  else if(whichCut==2){
    cuts.push_back( kNoCut && kIsNuMuCC && cutmyRFiducial && kNotClearCosmic );
    spillcuts.push_back( kNoSpillCut );
    outDirNames.push_back( "NuMuCC_RecoFiducial_NotClearCosmic" );
  }
  else if(whichCut==3){
    cuts.push_back( kNoCut && kIsNuMuCC && cutmyRFiducial && kNotClearCosmic && kNuScore );
    spillcuts.push_back( kNoSpillCut );
    outDirNames.push_back( "NuMuCC_RecoFiducial_NotClearCosmic_NuScoreGT0p4" );
  }
  else if(whichCut==4){
    cuts.push_back( kNoCut && kIsNuMuCC && cutmyRFiducial && kNotClearCosmic && kNuScore && cutTempNuMImyFMScore );
    spillcuts.push_back( kNoSpillCut );
    outDirNames.push_back( "NuMuCC_RecoFiducial_NotClearCosmic_NuScoreGT0p4_FMScoreLT6" );
  }
  else if(whichCut==5){
    cuts.push_back( kNoCut && kIsNuMuCC && cutmyRFiducial && kNotClearCosmic && kNuScore && cutTempNuMImyFMScore && cutTempNuMIHasMuon );
    spillcuts.push_back( kNoSpillCut );
    outDirNames.push_back( "NuMuCC_RecoFiducial_NotClearCosmic_NuScoreGT0p4_FMScoreLT6_HasMuon" );
  }
  else if(whichCut==6){
    cuts.push_back( kNoCut && kIsNuMuCC && cutmyRFiducial && kNotClearCosmic && kNuScore && cutTempNuMImyFMScore && cutTempNuMIHasMuon);
    spillcuts.push_back( kCRTHitVetoFD );
    outDirNames.push_back( "NuMuCC_RecoFiducial_NotClearCosmic_NuScoreGT0p4_FMScoreLT6_HasMuon_CRTVeto" );
  }
  else if(whichCut==7){
    cuts.push_back( kNoCut && kIsNuMuCC && cutmyRFiducial && kNotClearCosmic && kNuScore && cutTempNuMImyFMScore && cutTempNuMIHasMuon && cutTempNuMIHasProton);
    spillcuts.push_back( kCRTHitVetoFD );
    outDirNames.push_back( "NuMuCC_RecoFiducial_NotClearCosmic_NuScoreGT0p4_FMScoreLT6_HasMuon_CRTVeto_HasProton" );
  }
  //==== Additional particles
  else if(whichCut==8){
    cuts.push_back( kNoCut && kIsNuMuCC && cutmyRFiducial && kNotClearCosmic && kNuScore && cutTempNuMImyFMScore && cutTempNuMIHasMuon && cutTempNuMIHasProton
      && !cutTruthNoNeutron
    );
    spillcuts.push_back( kCRTHitVetoFD );
    outDirNames.push_back( "NuMuCC_RecoFiducial_NotClearCosmic_NuScoreGT0p4_FMScoreLT6_HasMuon_CRTVeto_HasProton_WithNeutron" );
  }
  else if(whichCut==9){
    cuts.push_back( kNoCut && kIsNuMuCC && cutmyRFiducial && kNotClearCosmic && kNuScore && cutTempNuMImyFMScore && cutTempNuMIHasMuon && cutTempNuMIHasProton
      && !cutTruthNoChargedPion
    );
    spillcuts.push_back( kCRTHitVetoFD );
    outDirNames.push_back( "NuMuCC_RecoFiducial_NotClearCosmic_NuScoreGT0p4_FMScoreLT6_HasMuon_CRTVeto_HasProton_WithChargedPion" );
  }
  else if(whichCut==10){
    cuts.push_back( kNoCut && kIsNuMuCC && cutmyRFiducial && kNotClearCosmic && kNuScore && cutTempNuMImyFMScore && cutTempNuMIHasMuon && cutTempNuMIHasProton
      && !cutTruthNoPiZero
    );
    spillcuts.push_back( kCRTHitVetoFD );
    outDirNames.push_back( "NuMuCC_RecoFiducial_NotClearCosmic_NuScoreGT0p4_FMScoreLT6_HasMuon_CRTVeto_HasProton_WithNeutralPion" );
  }
  else if(whichCut==11){
    cuts.push_back( kNoCut && kIsNuMuCC && cutmyRFiducial && kNotClearCosmic && kNuScore && cutTempNuMImyFMScore && cutTempNuMIHasMuon && cutTempNuMIHasProton
      && cutTruthNoNeutron && cutTruthNoChargedPion && cutTruthNoPiZero
    );
    spillcuts.push_back( kCRTHitVetoFD );
    outDirNames.push_back( "NuMuCC_RecoFiducial_NotClearCosmic_NuScoreGT0p4_FMScoreLT6_HasMuon_CRTVeto_HasProton_WithoutNeutronPion" );
  }
  //==== Contained/Exiting
  else if(whichCut==12){
    cuts.push_back( kNoCut && kIsNuMuCC && cutmyRFiducial && kNotClearCosmic && kNuScore && cutTempNuMImyFMScore && cutTempNuMIHasMuon && cutTempNuMIHasProton && cutTempNuMIMuonContained);
    spillcuts.push_back( kCRTHitVetoFD );
    outDirNames.push_back( "NuMuCC_RecoFiducial_NotClearCosmic_NuScoreGT0p4_FMScoreLT6_HasMuon_CRTVeto_HasProton_Contained" );
  }
  else if(whichCut==13){
    cuts.push_back( kNoCut && kIsNuMuCC && cutmyRFiducial && kNotClearCosmic && kNuScore && cutTempNuMImyFMScore && cutTempNuMIHasMuon && cutTempNuMIHasProton && !cutTempNuMIMuonContained);
    spillcuts.push_back( kCRTHitVetoFD );
    outDirNames.push_back( "NuMuCC_RecoFiducial_NotClearCosmic_NuScoreGT0p4_FMScoreLT6_HasMuon_CRTVeto_HasProton_Exiting" );
  }
  //==== Contained/Exiting just after muon selection

  else if(whichCut==14){
    cuts.push_back( kNoCut && kIsNuMuCC && cutmyRFiducial && kNotClearCosmic && kNuScore && cutTempNuMImyFMScore && cutTempNuMIHasMuon && cutTempNuMIMuonContained);
    spillcuts.push_back( kNoSpillCut );
    outDirNames.push_back( "NuMuCC_RecoFiducial_NotClearCosmic_NuScoreGT0p4_FMScoreLT6_HasMuon_Contained" );
  }
  else if(whichCut==15){
    cuts.push_back( kNoCut && kIsNuMuCC && cutmyRFiducial && kNotClearCosmic && kNuScore && cutTempNuMImyFMScore && cutTempNuMIHasMuon && !cutTempNuMIMuonContained);
    spillcuts.push_back( kNoSpillCut );
    outDirNames.push_back( "NuMuCC_RecoFiducial_NotClearCosmic_NuScoreGT0p4_FMScoreLT6_HasMuon_Exiting" );
  }

  std::vector<std::string> vec_inputFiles;
  vec_inputFiles.clear();
  //==== BNB
  if(whichSample==0){
    //==== ReCAF file
    vec_inputFiles = inputFiles;
    std::cout << "[run_Efficiency] Running on BNB files" << std::endl;
  }
  //==== BNB, reCAF
  else if(whichSample==1){
    vec_inputFiles = inputFilesReCAF;
    std::cout << "[run_Efficiency] Running on BNB, ReCAF files" << std::endl;
  }

  for(unsigned int i=0; i<cuts.size(); i++){

    Cut myCut = cuts.at(i);
    SpillCut mySpillCut = spillcuts.at(i);
    TString outDirName = outDirNames.at(i);

    myEfficiency m;

    //==== BNB
    m.outputDir = outDirBase;
    SpectrumLoader loader(vec_inputFiles);

    m.initialize();
    m.setSystematicWeights();

    m.bookTruth(loader, myCut, mySpillCut);
    if(outDirName.Contains("HasMuon")){
      m.bookRecoMuon(loader, myCut, mySpillCut);
    }
    if(outDirName.Contains("HasProton")){
      m.bookRecoProton(loader, myCut, mySpillCut);
      m.bookRecoNeutrino(loader, myCut, mySpillCut);
    }

    loader.Go();

    m.saveHistograms();

  }

}
