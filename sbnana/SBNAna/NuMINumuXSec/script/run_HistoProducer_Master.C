#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_HistoProducer.h"
#include "sbnana/SBNAna/Cuts/NumuCutsIcarus202106.h"
#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_Variables.h"
#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_Cuts.h"
#include "sbnana/SBNAna/Cuts/Cuts.h"

#include "myFilelistNuMIBeamPlusCosmicMC.h"
#include "myFilelistNuMIInTimeComsicMC.h"

using namespace ana;
using namespace ICARUSNumuXsec;

void run_HistoProducer_Master(int whichInputFile, int runType, TString outDirBase){

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
    vec_inputs = inputFilesNuMIBeamPlusCosmicMC;

    TargetPOT = 6e20;
    str_TargetPOT = "6e20 POT";
  }
  //==== NuMI, MC, Cosmics in-time
  else if(whichInputFile==1){
    cout << "[run_Efficiency] Running NuMI In-time cosmic" << endl;
    inputDef = "/pnfs/icarus/scratch/users/jskim/mc/NUMI_in-time_Cosmics/flatcaf/v09_41_00/Try1/flatcaf_*.root";
    vec_inputs = inputFilesNuMIInTimeComsicMC;

    isCosmicInput = true;
  }
/*
  //==== NuMI, Data, Run5679
  else if(whichInputFile==2){
    cout << "[run_Efficiency] Running run5679" << endl;
    vec_inputs = inputFilesNuMIDataRun5679;

    TargetPOT = 6e20;
    str_TargetPOT = "6e20 POT";
    isDataInput = true;
  }
  //==== NoBeam, Data, Run6819
  else if(whichInputFile==3){
    cout << "[run_Efficiency] Running run6819" << endl;
    vec_inputs = inputFilesNoBeamDataRun6819;

    TargetPOT = 6e20;
    str_TargetPOT = "6e20 POT";
    isDataInput = true;
  }
*/
  cout << "[run_Efficiency] Input definition : " << inputDef << endl;
  //SpectrumLoader loader("/pnfs/icarus/scratch/users/jskim/mc/NUMI_Nu_Cosmics/flatcaf/v09_41_00/Try1/flatcaf_0.root");
  //SpectrumLoader loader("/pnfs/icarus/scratch/users/jskim/mc/NUMI_in-time_Cosmics/flatcaf/v09_41_00/Try2/flatcaf_0.root");
  //SpectrumLoader loader("test_BNBMC.root");
  SpectrumLoader loader("test_NewGENIESyst.root");
  //SpectrumLoader loader(vec_inputs);

  //==== runType
  //====  0 : Only nominal;     Incl+PerInt.; +Cutflow
  //==== 11 : Only systematics; Incl        ; X
  //==== 21 : Only ststematics;      PerInt.; X

  bool fillNominal = false;
  bool fillSystematic = false;

  bool runInclusiveProcess = false;
  bool runPerInteraction = false;

  bool runCutflow = false;

  if(runType==0){

    fillNominal = true;
    fillSystematic = false;

    runInclusiveProcess = true;
    runPerInteraction = true;

    runCutflow = true;

  }
  if(runType==11){

    fillNominal = false;
    fillSystematic = true;

    runInclusiveProcess = true;
    runPerInteraction = false;

    runCutflow = false;

  }
  if(runType==21){

    fillNominal = false;
    fillSystematic = true;

    runInclusiveProcess = false;
    runPerInteraction = true;

    runCutflow = false;

  }

  //==== Define samples;
  vector<TString> baseSampleNames;
  vector<Cut> baseSampleCuts;
  vector<SpillCut> baseSampleSpillCuts;

  baseSampleNames.push_back("AllSamples");
  baseSampleCuts.push_back(kNoCut);
  baseSampleSpillCuts.push_back(kNoSpillCut);

  vector<Cut> cuts;
  vector<SpillCut> spillcuts;
  vector<TString> cutNames;

  //====   1) Numu-cc selection

  cuts.push_back( kNoCut );
  spillcuts.push_back( kNoSpillCut );
  cutNames.push_back( "NoCut" );

  cuts.push_back( kNoCut );
  spillcuts.push_back( kCRTHitVetoFD );
  cutNames.push_back( "CRTVeto" );

  //==== Declare HistoProducer

  HistoProducer m;
  m.outputDir = outDirBase;
  m.fillNominal = fillNominal;


  m.initialize();
  if(fillSystematic){
    m.setSystematicWeights();
  }

  m.TargetPOT = TargetPOT;
  m.str_TargetPOT = str_TargetPOT;

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

        m.bookSpectrums(loader, mySpillCut, myCut);

      } // END if setCut

    } // END loop cut

  } // END sample cut

  loader.Go();

  m.saveHistograms();

}
