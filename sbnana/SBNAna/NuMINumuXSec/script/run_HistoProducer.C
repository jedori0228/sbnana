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

void run_HistoProducer(int whichInputFile, TString outDirBase){

  //==== Inputfile

  string inputDef = "";
  //====   BNB, new central sample, IcarusProd2021B_BNB_Nu_Cosmics_v09_28_01_01_01_caf
  if(whichInputFile==0){
    cout << "[run_Efficiency] Running BNB, IcarusProd2021B_BNB_Nu_Cosmics_v09_28_01_01_01_caf" << endl;
    inputDef = "IcarusProd2021B_BNB_Nu_Cosmics_v09_28_01_01_01_caf";
  }
  //====   NuMI, new central sample, IcarusProd2021B_NuMI_Nu_Cosmics_v09_28_01_01_01_caf
  else if(whichInputFile==2){
    cout << "[run_Efficiency] Running NuMI, IcarusProd2021B_NuMI_Nu_Cosmics_v09_28_01_01_01_caf" << endl;
    inputDef = "IcarusProd2021B_NuMI_Nu_Cosmics_v09_28_01_01_01_caf";
  }
  //====   Cosmic, IcarusProd2021B_Intime_Cosmic_v09_28_01_01_01_caf
  else if(whichInputFile==3){
    cout << "[run_Efficiency] Running In-time BNB cosmic" << endl;
    inputDef = "IcarusProd2021B_Intime_Cosmic_v09_28_01_01_01_caf";
  }

  cout << "[run_Efficiency] Input definition : " << inputDef << endl;
  SpectrumLoader loader(inputDef);

  //SpectrumLoader loader("/pnfs/sbn/data/sbn_fd/poms_production/ICARUS_NuMI_Nu_Cosmics/mc/reconstructed/icaruscode/v09_28_01_01_01/caf/00/49/prodcorsika_genie_protononly_icarus_numi_volDetEnclosure_tpc_gen_filter_g4_detsim_24037195_1096_reco1_20210913T155758_reco2.caf-6f94e5ce-28b2-4c27-bc19-525bf3c76d7f.root");

  //==== Define samples;
  vector<TString> baseSampleNames;
  vector<Cut> baseSampleCuts;
  vector<SpillCut> baseSampleSpillCuts;

  baseSampleNames.push_back("NuMuCC");
  baseSampleCuts.push_back(kIsNuMuCC);
  baseSampleSpillCuts.push_back(kNoSpillCut);

  baseSampleNames.push_back("NuMuNC");
  baseSampleCuts.push_back(cutIsNuMuNC);
  baseSampleSpillCuts.push_back(kNoSpillCut);

  baseSampleNames.push_back("Cosmic");
  baseSampleCuts.push_back(kIsCosmic);
  baseSampleSpillCuts.push_back(kNoSpillCut);

  //==== Define Cuts

  vector<Cut> cuts;
  vector<SpillCut> spillcuts;
  vector<TString> cutNames;

  //====   1) Numu-cc selection

  cuts.push_back( kNoCut );
  spillcuts.push_back( kNoSpillCut );
  cutNames.push_back( "NoCut" );

  cuts.push_back( kNoCut && cutRFiducial );
  spillcuts.push_back( kNoSpillCut );
  cutNames.push_back( "RecoFiducial" );

  cuts.push_back( kNoCut && cutRFiducial && kNotClearCosmic );
  spillcuts.push_back( kNoSpillCut );
  cutNames.push_back( "NotClearCosmic" );

  cuts.push_back( kNoCut && cutRFiducial && kNotClearCosmic );
  spillcuts.push_back( kCRTHitVetoFD );
  cutNames.push_back( "CRTVeto" );

  cuts.push_back( kNoCut && cutRFiducial && kNotClearCosmic && cutFMScore );
  spillcuts.push_back( kCRTHitVetoFD );
  cutNames.push_back( "FMScoreLT6" );

  cuts.push_back( kNoCut && cutRFiducial && kNotClearCosmic && cutFMScore && kNuScore );
  spillcuts.push_back( kCRTHitVetoFD );
  cutNames.push_back( "NuScoreGT0p4" );

  cuts.push_back( kNoCut && cutRFiducial && kNotClearCosmic && cutFMScore && kNuScore && cutHasMuon );
  spillcuts.push_back( kCRTHitVetoFD );
  cutNames.push_back( "HasMuon" );

  //====     Contained/Exiting

  cuts.push_back( kNoCut && cutRFiducial && kNotClearCosmic && cutFMScore && kNuScore && cutHasMuon && cutMuonContained );
  spillcuts.push_back( kCRTHitVetoFD );
  cutNames.push_back( "HasMuon_Contained" );

  cuts.push_back( kNoCut && cutRFiducial && kNotClearCosmic && cutFMScore && kNuScore && cutHasMuon && !cutMuonContained);
  spillcuts.push_back( kCRTHitVetoFD );
  cutNames.push_back( "HasMuon_Exiting" );

  //====       Matched to muon/proton

  cuts.push_back( kNoCut && cutRFiducial && kNotClearCosmic && cutFMScore && kNuScore && cutHasMuon && cutMuonContained && cutMuonMatchedToMuon );
  spillcuts.push_back( kCRTHitVetoFD );
  cutNames.push_back( "HasMuon_Contained_MatchedToMuon" );

  cuts.push_back( kNoCut && cutRFiducial && kNotClearCosmic && cutFMScore && kNuScore && cutHasMuon && cutMuonContained && cutMuonMatchedToProton );
  spillcuts.push_back( kCRTHitVetoFD );
  cutNames.push_back( "HasMuon_Contained_MatchedToProton" );

  cuts.push_back( kNoCut && cutRFiducial && kNotClearCosmic && cutFMScore && kNuScore && cutHasMuon && !cutMuonContained && cutMuonMatchedToMuon );
  spillcuts.push_back( kCRTHitVetoFD );
  cutNames.push_back( "HasMuon_Exiting_MatchedToMuon" );

  cuts.push_back( kNoCut && cutRFiducial && kNotClearCosmic && cutFMScore && kNuScore && cutHasMuon && !cutMuonContained && cutMuonMatchedToProton );
  spillcuts.push_back( kCRTHitVetoFD );
  cutNames.push_back( "HasMuon_Exiting_MatchedToProton" );

  //====     Additional particles

  cuts.push_back( kNoCut && cutRFiducial && kNotClearCosmic && cutFMScore && kNuScore && cutHasMuon && !cutTruthNoNeutron );
  spillcuts.push_back( kCRTHitVetoFD );
  cutNames.push_back( "HasMuon_WithNeutron" );

  cuts.push_back( kNoCut && cutRFiducial && kNotClearCosmic && cutFMScore && kNuScore && cutHasMuon && !cutTruthNoChargedPion );
  spillcuts.push_back( kCRTHitVetoFD );
  cutNames.push_back( "HasMuon_WithChargedPion" );

  cuts.push_back( kNoCut && cutRFiducial && kNotClearCosmic && cutFMScore && kNuScore && cutHasMuon && !cutTruthNoPiZero );
  spillcuts.push_back( kCRTHitVetoFD );
  cutNames.push_back( "HasMuon_WithNeutralPion" );

  cuts.push_back( kNoCut && cutRFiducial && kNotClearCosmic && cutFMScore && kNuScore && cutHasMuon && cutTruthNoNeutron && cutTruthNoChargedPion && cutTruthNoPiZero );
  spillcuts.push_back( kCRTHitVetoFD );
  cutNames.push_back( "HasMuon_WithoutNeutronOrPion" );

  //====     HasProton added

  cuts.push_back( kNoCut && cutRFiducial && kNotClearCosmic && cutFMScore && kNuScore && cutHasMuon && cutHasProton );
  spillcuts.push_back( kCRTHitVetoFD );
  cutNames.push_back( "HasMuon_HasProton" );

  //==== Declare HistoProducer

  HistoProducer m;
  m.outputDir = outDirBase;
  m.initialize();
  //m.setSystematicWeights(); // TODO Run this to get systematic variations

  for(unsigned int is=0; is<baseSampleNames.size(); is++){

    TString baseSampleName = baseSampleNames.at(is);
    Cut baseSampleCut = baseSampleCuts.at(is);
    SpillCut baseSampleSpillCut = baseSampleSpillCuts.at(is);

    bool isCosmic = (baseSampleName=="Cosmic");

    for(unsigned int i=0; i<cuts.size(); i++){

      Cut myCut = baseSampleCut && cuts.at(i);
      SpillCut mySpillCut = baseSampleSpillCut && spillcuts.at(i);
      TString cutName = baseSampleName+"_"+cutNames.at(i);

      if( m.setCut(cutName) ){

        if(isCosmic){

          m.bookSlice(loader, myCut, mySpillCut);
          if(cutName.Contains("HasMuon")){
            m.bookRecoMuon(loader, myCut, mySpillCut);
          }

          //==== For cosmic, remove truth-related variables to prevent warning messages
          map< TString, vector<Spectrum *> > map_cutName_to_vec_Spectrums;
          map< TString, vector<EnsembleSpectrum *> > map_cutName_to_vec_EnsembleSpectrums;
          for(std::map< TString, vector<Spectrum *> >::iterator mapit=map_cutName_to_vec_Spectrums.begin(); mapit!=map_cutName_to_vec_Spectrums.end(); mapit++){
            bool isTruthVariable = false;
            if(mapit->first.Contains("True")) isTruthVariable = true;
            if(mapit->first.Contains("Residual")) isTruthVariable = true;
            if(mapit->first.Contains("Truth")) isTruthVariable = true;
            map_cutName_to_vec_Spectrums.erase(mapit->first);
          }
          for(std::map< TString, vector<EnsembleSpectrum *> >::iterator mapit=map_cutName_to_vec_EnsembleSpectrums.begin(); mapit!=map_cutName_to_vec_EnsembleSpectrums.end(); mapit++){
            bool isTruthVariable = false;
            if(mapit->first.Contains("True")) isTruthVariable = true;
            if(mapit->first.Contains("Residual")) isTruthVariable = true;
            if(mapit->first.Contains("Truth")) isTruthVariable = true;
            map_cutName_to_vec_EnsembleSpectrums.erase(mapit->first);
          }

        }
        else{

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

      } // END if setCut

    } // END loop cut

  } // END sample cut

  loader.Go();

  m.saveHistograms();

}
