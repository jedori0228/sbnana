#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_HistoProducer.h"
#include "sbnana/SBNAna/Cuts/NumuCutsIcarus202106.h"
#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_Variables.h"
#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_Cuts.h"
#include "sbnana/SBNAna/Cuts/Cuts.h"

#include "myFilelist_MC_NUMI_Nu_Cosmics.h"
#include "myFilelist_MC_NUMI_InTimeCosmics.h"
#include "myFilelist_DATA_Run7568_NUMI.h"
#include "myFilelist_DATA_Run7568_OffbeamNUMI.h"

using namespace ana;
using namespace ICARUSNumuXsec;

void run_QuickEfficiency(int whichInputFile, TString outDirBase){

  //==== Inputfile

  string inputDef = "";
  vector<string> vec_inputs;
  bool isDataInput = false;
  bool isCosmicInput = false;
  double TargetPOT = 6e20;
  TString str_TargetPOT = "6e20 POT";

  SpillCut cut_SideCRTVeto = spillcutSideCRTHitVetoFD;
  SpillCut cut_CRTVeto = kCRTHitVetoFD;

  //==== v09_45_00
  //==== NuMI, MC, Beam+Cosmic 
  if(whichInputFile==0){
    cout << "[run_Efficiency] Running NuMI+Cosmic" << endl;
    inputDef = "/pnfs/icarus/persistent/users/jskim/mc/NUMI_Nu_Cosmics/flatcaf/v09_45_00/220315_CRTT0Fixed/flatcaf_*.root";
    vec_inputs = inputFiles_MC_NUMI_Nu_Cosmics;

    TargetPOT = 6e20;
    str_TargetPOT = "6e20 POT";

    cut_SideCRTVeto = spillcutSideCRTHitVetoFDUpdatedT0;
    cut_CRTVeto = spillcutCRTHitVetoFDUpdatedT0;

  }
  //==== v09_45_00
  //==== NuMI, MC, Cosmics in-time
  else if(whichInputFile==1){
    cout << "[run_Efficiency] Running NuMI In-time cosmic" << endl;
    inputDef = "/pnfs/icarus/persistent/users/jskim/mc/NUMI_in-time_Cosmics/flatcaf/v09_45_00/220315_CRTT0Fixed/flatcaf_*.root";
    vec_inputs = inputFiles_MC_NUMI_InTimeCosmics;

    isCosmicInput = true;

    cut_SideCRTVeto = spillcutSideCRTHitVetoFDUpdatedT0;
    cut_CRTVeto = spillcutCRTHitVetoFDUpdatedT0;

  }

  //==== Production
  //==== NuMI, MC, Beam+Cosmic 
  if(whichInputFile==2){
    cout << "[run_Efficiency] Running NuMI+Cosmic" << endl;
    inputDef = "/pnfs/icarus/persistent/users/jskim/mc/NUMI_Nu_Cosmics/flatcaf/v09_37_01_03p01/220315_SameAsProduction/flatcaf_*.root";
    vec_inputs.clear();
    for(unsigned int i=0; i<20; i++){
      TString ts_fp = "/pnfs/icarus/persistent/users/jskim/mc/NUMI_Nu_Cosmics/flatcaf/v09_37_01_03p01/220315_SameAsProduction/flatcaf_"+TString::Itoa(i,10)+".root";
      string s_fp = ts_fp.Data();
      vec_inputs.push_back( s_fp );
    }

    TargetPOT = 6e20;
    str_TargetPOT = "6e20 POT";
  }
  //==== Production
  //==== NuMI, MC, Cosmics in-time
  else if(whichInputFile==3){
    cout << "[run_Efficiency] Running NuMI In-time cosmic" << endl;
    inputDef = "/pnfs/icarus/persistent/users/jskim/mc/NUMI_in-time_Cosmics/flatcaf/v09_37_01_03p01/220315_SameAsProduction/flatcaf_*.root";
    vec_inputs.clear();
    for(unsigned int i=0; i<20; i++){
      TString ts_fp = "/pnfs/icarus/persistent/users/jskim/mc/NUMI_in-time_Cosmics/flatcaf/v09_37_01_03p01/220315_SameAsProduction/flatcaf_"+TString::Itoa(i,10)+".root";
      string s_fp = ts_fp.Data();
      vec_inputs.push_back( s_fp );
    }

    isCosmicInput = true;
  }

  //==== v09_45_00
  //==== NuMI, MC, Beam+Cosmic, NoSCE
  if(whichInputFile==4){
    cout << "[run_Efficiency] Running NuMI+Cosmic" << endl;
    inputDef = "/pnfs/icarus/persistent/users/jskim/mc/NUMI_Nu_Cosmics/flatcaf/v09_45_00/220312_TestMCGeneration_NoSCE/flatcaf_*.root";
    vec_inputs.clear();
    for(unsigned int i=0; i<20; i++){
      TString ts_fp = "/pnfs/icarus/persistent/users/jskim/mc/NUMI_Nu_Cosmics/flatcaf/v09_45_00/220312_TestMCGeneration_NoSCE/flatcaf_"+TString::Itoa(i,10)+".root";
      string s_fp = ts_fp.Data();
      vec_inputs.push_back( s_fp );
    }

    TargetPOT = 6e20;
    str_TargetPOT = "6e20 POT";

    //cut_SideCRTVeto = spillcutSideCRTHitVetoFDUpdatedT0;
    //cut_CRTVeto = spillcutCRTHitVetoFDUpdatedT0;

  }

  //==== v09_45_00
  //==== NuMI, Data, Run7568, NUMI
  else if(whichInputFile==5){
    cout << "[run_Efficiency] Running run7568" << endl;
    inputDef = "/pnfs/icarus/persistent/users/jskim/data/run_7568/flatcaf/v09_45_00/220315_CRTT0Fixed/NUMI/flatcaf_*";
    vec_inputs = inputFiles_DATA_Run7568_NUMI;

    TargetPOT = 6e20;
    str_TargetPOT = "6e20 POT";
    isDataInput = true;

    cut_SideCRTVeto = spillcutSideCRTHitVetoFDUpdatedT0;
    cut_CRTVeto = spillcutCRTHitVetoFDUpdatedT0;

  }
  //==== v09_45_00
  //==== NuMI, Data, Run7568, OffbeamNUMI
  else if(whichInputFile==6){
    cout << "[run_Efficiency] Running run7568" << endl;
    inputDef = "/pnfs/icarus/persistent/users/jskim/data/run_7568/flatcaf/v09_45_00/220315_CRTT0Fixed/OffbeamNUMI/flatcaf_*";
    vec_inputs = inputFiles_DATA_Run7568_OffbeamNUMI;

    isDataInput = true;

    cut_SideCRTVeto = spillcutSideCRTHitVetoFDUpdatedT0;
    cut_CRTVeto = spillcutCRTHitVetoFDUpdatedT0;

  }

  cout << "[run_Efficiency] Input definition : " << inputDef << endl;
  SpectrumLoader loader(vec_inputs);

  //==== Define samples;

  vector<TString> baseSampleNames;
  vector<Cut> baseSampleCuts;
  vector<SpillCut> baseSampleSpillCuts;

  baseSampleNames.push_back("AllSamples");
  baseSampleCuts.push_back(kNoCut);
  baseSampleSpillCuts.push_back(kNoSpillCut);

  baseSampleNames.push_back("NuMuCC");
  baseSampleCuts.push_back(kIsNuMuCC);
  baseSampleSpillCuts.push_back(kNoSpillCut);

  baseSampleNames.push_back("NuMuCCSignalDef");
  baseSampleCuts.push_back(cutIsNuMuCCSignalDef);
  baseSampleSpillCuts.push_back(kNoSpillCut);

  baseSampleNames.push_back("NonNuMuCC");
  baseSampleCuts.push_back(!kIsNuMuCC);
  baseSampleSpillCuts.push_back(kNoSpillCut);

  baseSampleNames.push_back("NuMuNC");
  baseSampleCuts.push_back(cutIsNuMuNC);
  baseSampleSpillCuts.push_back(kNoSpillCut);

  baseSampleNames.push_back("NuECC");
  baseSampleCuts.push_back(cutIsNuECC);
  baseSampleSpillCuts.push_back(kNoSpillCut);

  baseSampleNames.push_back("Cosmic");
  baseSampleCuts.push_back(kIsCosmic);
  baseSampleSpillCuts.push_back(kNoSpillCut);

  //==== cut

  vector<Cut> cuts;
  vector<SpillCut> spillcuts;
  vector<TString> cutNames;

  //==== Base selecitons

  cuts.push_back( kNoCut );
  spillcuts.push_back( kNoSpillCut );
  cutNames.push_back( "NoCut" );

  cuts.push_back( cutRFiducial );
  spillcuts.push_back( kNoSpillCut );
  cutNames.push_back( "RecoFiducial" );

  cuts.push_back( cutRFiducial && kNotClearCosmic );
  spillcuts.push_back( kNoSpillCut );
  cutNames.push_back( "NotClearCosmic" );

  cuts.push_back( cutRFiducial && kNotClearCosmic && cutFMScore);
  spillcuts.push_back( kNoSpillCut );
  cutNames.push_back( "FMScore" );

  //==== CRTVeto on NoCut
  cuts.push_back( kNoCut );
  spillcuts.push_back( cut_SideCRTVeto );
  cutNames.push_back( "NoCut_SideCRTVeto" );
  cuts.push_back( kNoCut );
  spillcuts.push_back( cut_CRTVeto );
  cutNames.push_back( "NoCut_CRTVeto" );

  //==== CRTVeto on RecoFiducial
  cuts.push_back( cutRFiducial );
  spillcuts.push_back( cut_SideCRTVeto );
  cutNames.push_back( "RecoFiducial_SideCRTVeto" );
  cuts.push_back( cutRFiducial );
  spillcuts.push_back( cut_CRTVeto );
  cutNames.push_back( "RecoFiducial_CRTVeto" );

  //==== CRTVeto on NotClearCosmic
  cuts.push_back( cutRFiducial && kNotClearCosmic );
  spillcuts.push_back( cut_SideCRTVeto );
  cutNames.push_back( "NotClearCosmic_SideCRTVeto" );
  cuts.push_back( cutRFiducial && kNotClearCosmic );
  spillcuts.push_back( cut_CRTVeto );
  cutNames.push_back( "NotClearCosmic_CRTVeto" );

  //==== CRTVeto on FMScore
  cuts.push_back( cutRFiducial && kNotClearCosmic && cutFMScore);
  spillcuts.push_back( cut_SideCRTVeto );
  cutNames.push_back( "FMScore_SideCRTVeto" );
  cuts.push_back( cutRFiducial && kNotClearCosmic && cutFMScore);
  spillcuts.push_back( cut_CRTVeto );
  cutNames.push_back( "FMScore_CRTVeto" );

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

        m.bookQuickCutEfficiency(loader, mySpillCut, myCut);

      } // END if setCut

    } // END loop cut

  } // END sample cut

  loader.Go();

  m.saveHistograms();

}
