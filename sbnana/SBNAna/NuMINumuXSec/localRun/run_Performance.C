#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_HistoProducer.h"
#include "sbnana/SBNAna/Cuts/NumuCutsIcarus202106.h"
#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_Variables.h"
#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_Cuts.h"
#include "sbnana/SBNAna/Cuts/Cuts.h"

#include "myFilelist_MC_NUMI_Nu_Cosmics.h"
#include "myFilelist_MC_NUMI_InTimeCosmics.h"

using namespace ana;
using namespace ICARUSNumuXsec;

void run_Performance(int whichInputFile, TString outDirBase){

  //==== Inputfile

  string inputDef = "";
  vector<string> vec_inputs;
  bool isDataInput = false;
  bool isCosmicInput = false;
  double TargetPOT = 6e20;
  TString str_TargetPOT = "6e20 POT";

  //==== NuMI, MC, Beam+Cosmic 
  if(whichInputFile==0){
    cout << "[run_Efficiency] Running NuMI+Cosmic" << endl;
    inputDef = "/pnfs/icarus/scratch/users/jskim/mc/NUMI_Nu_Cosmics/flatcaf/v09_41_00/Try1/flatcaf_*.root";
    vec_inputs = inputFiles_MC_NUMI_Nu_Cosmics;

    TargetPOT = 6e20;
    str_TargetPOT = "6e20 POT";
  }
  //==== NuMI, MC, Cosmics in-time
  else if(whichInputFile==1){
    cout << "[run_Efficiency] Running NuMI In-time cosmic" << endl;
    inputDef = "/pnfs/icarus/scratch/users/jskim/mc/NUMI_in-time_Cosmics/flatcaf/v09_41_00/Try1/flatcaf_*.root";
    vec_inputs = inputFiles_MC_NUMI_InTimeCosmics;

    isCosmicInput = true;
  }

  cout << "[run_Efficiency] Input definition : " << inputDef << endl;
  SpectrumLoader loader(vec_inputs);

  //==== Define samples;

  vector<TString> baseSampleNames;
  vector<Cut> baseSampleCuts;
  vector<SpillCut> baseSampleSpillCuts;

  baseSampleNames.push_back("NuMuCC");
  baseSampleCuts.push_back(kIsNuMuCC);
  baseSampleSpillCuts.push_back(kNoSpillCut);

  //==== cut

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


  cuts.push_back( kNoCut && cutRFiducial && kNotClearCosmic && cutHasMuon );
  spillcuts.push_back( kNoSpillCut );
  cutNames.push_back( "HasMuon" );
  cuts.push_back( kNoCut && cutRFiducial && kNotClearCosmic && cutHasMuon && cutMuonMatchedToMuon );
  spillcuts.push_back( kNoSpillCut );
  cutNames.push_back( "HasMuon_MatchedToMuon" );
  cuts.push_back( kNoCut && cutRFiducial && kNotClearCosmic && cutHasMuon && !cutMuonMatchedToMuon );
  spillcuts.push_back( kNoSpillCut );
  cutNames.push_back( "HasMuon_NotMatchedToMuon" );

  cuts.push_back( kNoCut && cutRFiducial && kNotClearCosmic && cutHasMuon && cutMuonContained );
  spillcuts.push_back( kNoSpillCut );
  cutNames.push_back( "HasMuon_Contained" );
  cuts.push_back( kNoCut && cutRFiducial && kNotClearCosmic && cutHasMuon && cutMuonContained && cutMuonMatchedToMuon );
  spillcuts.push_back( kNoSpillCut );
  cutNames.push_back( "HasMuon_Contained_MatchedToMuon" );
  cuts.push_back( kNoCut && cutRFiducial && kNotClearCosmic && cutHasMuon && cutMuonContained && !cutMuonMatchedToMuon );
  spillcuts.push_back( kNoSpillCut );
  cutNames.push_back( "HasMuon_Contained_NotMatchedToMuon" );
  cuts.push_back( kNoCut && cutRFiducial && kNotClearCosmic && cutHasMuon && cutMuonContained && cutRecoMuonTruthContained );
  spillcuts.push_back( kNoSpillCut );
  cutNames.push_back( "HasMuon_Contained_TruthContained" );
  cuts.push_back( kNoCut && cutRFiducial && kNotClearCosmic && cutHasMuon && cutMuonContained && !cutRecoMuonTruthContained );
  spillcuts.push_back( kNoSpillCut );
  cutNames.push_back( "HasMuon_Contained_TruthExiting" );

  cuts.push_back( kNoCut && cutRFiducial && kNotClearCosmic && cutHasMuon && !cutMuonContained );
  spillcuts.push_back( kNoSpillCut );
  cutNames.push_back( "HasMuon_Exiting" );
  cuts.push_back( kNoCut && cutRFiducial && kNotClearCosmic && cutHasMuon && !cutMuonContained && cutMuonMatchedToMuon );
  spillcuts.push_back( kNoSpillCut );
  cutNames.push_back( "HasMuon_Exiting_MatchedToMuon" );
  cuts.push_back( kNoCut && cutRFiducial && kNotClearCosmic && cutHasMuon && !cutMuonContained && !cutMuonMatchedToMuon );
  spillcuts.push_back( kNoSpillCut );
  cutNames.push_back( "HasMuon_Exiting_NotMatchedToMuon" );
  cuts.push_back( kNoCut && cutRFiducial && kNotClearCosmic && cutHasMuon && !cutMuonContained && cutRecoMuonTruthContained );
  spillcuts.push_back( kNoSpillCut );
  cutNames.push_back( "HasMuon_Exiting_TruthContained" );
  cuts.push_back( kNoCut && cutRFiducial && kNotClearCosmic && cutHasMuon && !cutMuonContained && !cutRecoMuonTruthContained );
  spillcuts.push_back( kNoSpillCut );
  cutNames.push_back( "HasMuon_Exiting_TruthExiting" );


  cuts.push_back( kNoCut && cutRFiducial && kNotClearCosmic && cutHasMuon && cutHasProton );
  spillcuts.push_back( kNoSpillCut );
  cutNames.push_back( "HasProton" );

  cuts.push_back( kNoCut && cutRFiducial && kNotClearCosmic && cutHasMuon && cutHasProton && cutProtonMatchedToProton );
  spillcuts.push_back( kNoSpillCut );
  cutNames.push_back( "HasProton_MatchedToProton" );

  cuts.push_back( kNoCut && cutRFiducial && kNotClearCosmic && cutHasMuon && cutHasProton && !cutProtonMatchedToProton );
  spillcuts.push_back( kNoSpillCut );
  cutNames.push_back( "HasProton_NotMatchedToProton" );

  //==== Declare HistoProducer

  HistoProducer m;
  m.outputDir = outDirBase;

  m.initialize();

  m.TargetPOT = TargetPOT;
  m.str_TargetPOT = str_TargetPOT;

  for(unsigned int is=0; is<baseSampleNames.size(); is++){

    TString baseSampleName = baseSampleNames.at(is);
    Cut baseSampleCut = baseSampleCuts.at(is);
    SpillCut baseSampleSpillCut = baseSampleSpillCuts.at(is);

    for(unsigned int i=0; i<cuts.size(); i++){

      Cut myCut = baseSampleCut && cuts.at(i);
      SpillCut mySpillCut = baseSampleSpillCut && spillcuts.at(i);
      TString cutName = baseSampleName+"_"+cutNames.at(i);

      if( m.setCut(cutName) ){

        m.bookRecoPerformance(loader, mySpillCut, myCut);

      } // END if setCut

    } // END loop cut

  } // END sample cut

  loader.Go();

  m.saveHistograms();

}
