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

void run_HistoProducer_Master(int whichInputFile, int runType, TString outDirBase){

  //==== Inputfile

  string inputDef = "";
  vector<string> vec_inputs;
  bool isDataInput = false;
  bool isCosmicInput = false;
  double TargetPOT = 6e20;
  TString str_TargetPOT = "6e20 POT";

  //==== CRTVeto differs between releases (T1 vs T0)
  SpillCut cut_SideCRTVeto = spillcutSideCRTHitVetoFDUpdatedT0;
  SpillCut cut_CRTVeto = spillcutCRTHitVetoFDUpdatedT0;
  //cut_CRTVeto = spillcutSideCRTHitVetoFDUpdatedT0; // SideCRTOnly

  //==== v09_45_00
  //==== NuMI, MC, Beam+Cosmic 
  if(whichInputFile==0){
    cout << "[run_Efficiency] Running NuMI+Cosmic" << endl;
    inputDef = "/pnfs/icarus/persistent/users/jskim/mc/NUMI_Nu_Cosmics/flatcaf/v09_45_00/220315_CRTT0Fixed/flatcaf_*.root";
    vec_inputs = inputFiles_MC_NUMI_Nu_Cosmics;

    TargetPOT = 6e20;
    str_TargetPOT = "6e20 POT";

  }
  //==== v09_45_00
  //==== NuMI, MC, Cosmics in-time
  else if(whichInputFile==1){
    cout << "[run_Efficiency] Running NuMI In-time cosmic" << endl;
    inputDef = "/pnfs/icarus/persistent/users/jskim/mc/NUMI_in-time_Cosmics/flatcaf/v09_45_00/220315_CRTT0Fixed/flatcaf_*.root";
    vec_inputs = inputFiles_MC_NUMI_InTimeCosmics;

    isCosmicInput = true;

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

    cut_SideCRTVeto = spillcutSideCRTHitVetoFD;
    cut_CRTVeto = kCRTHitVetoFD;

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

    cut_SideCRTVeto = spillcutSideCRTHitVetoFD;
    cut_CRTVeto = kCRTHitVetoFD;

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

    //==== T0 not fixed
    cut_SideCRTVeto = spillcutSideCRTHitVetoFD;
    cut_CRTVeto = kCRTHitVetoFD;

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

  }
  //==== v09_45_00
  //==== NuMI, Data, Run7568, OffbeamNUMI
  else if(whichInputFile==6){
    cout << "[run_Efficiency] Running run7568" << endl;
    inputDef = "/pnfs/icarus/persistent/users/jskim/data/run_7568/flatcaf/v09_45_00/220315_CRTT0Fixed/OffbeamNUMI/flatcaf_*";
    vec_inputs = inputFiles_DATA_Run7568_OffbeamNUMI;

    isDataInput = true;

  }


  cout << "[run_Efficiency] Input definition : " << inputDef << endl;
  SpectrumLoader loader(vec_inputs);

  //==== runType
  //====  0 : Only nominal;     Incl+PerInt.; +Cutflow
  //==== 11 : Only systematics; Incl        ; X
  //==== 21 : Only ststematics;      PerInt.; X

  bool fillNominal = false;
  bool fillSystematic = false;
  bool fillBeamInfo = false;

  bool runInclusiveProcess = false;
  bool runPerInteraction = false;

  bool runCutflow = false;

  if(runType==0){

    fillNominal = true;
    fillSystematic = false;
    fillBeamInfo = true;

    runInclusiveProcess = true;
    runPerInteraction = true;

    runCutflow = true;

  }
  if(runType==11){

    fillNominal = false;
    fillSystematic = true;
    fillBeamInfo = false;

    runInclusiveProcess = true;
    runPerInteraction = false;

    runCutflow = false;

  }
  if(runType==21){

    fillNominal = false;
    fillSystematic = true;
    fillBeamInfo = false;

    runInclusiveProcess = false;
    runPerInteraction = true;

    runCutflow = false;

  }

  //==== TODO
  runPerInteraction = false;

  //==== Define samples;
  vector<TString> baseSampleNames;
  vector<Cut> baseSampleCuts;
  vector<SpillCut> baseSampleSpillCuts;

  baseSampleNames.push_back("AllSamples");
  baseSampleCuts.push_back(kNoCut);
  baseSampleSpillCuts.push_back(kNoSpillCut);

  if(!isDataInput && !isCosmicInput){

    if(runInclusiveProcess){

      baseSampleNames.push_back("NuMuCC");
      baseSampleCuts.push_back(kIsNuMuCC);
      baseSampleSpillCuts.push_back(kNoSpillCut);

      baseSampleNames.push_back("NuMuCCSignalDef");
      baseSampleCuts.push_back(cutIsNuMuCCSignalDef);
      baseSampleSpillCuts.push_back(kNoSpillCut);

      baseSampleNames.push_back("NuMuCCSignalDefContained");
      baseSampleCuts.push_back(cutIsNuMuCCSignalDef&&cutTruthMuonContained);
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

    }
    if(runPerInteraction){

      baseSampleNames.push_back("CC");
      baseSampleCuts.push_back(cutIsCC);
      baseSampleSpillCuts.push_back(kNoSpillCut);

      baseSampleNames.push_back("NC");
      baseSampleCuts.push_back(cutIsNC);
      baseSampleSpillCuts.push_back(kNoSpillCut);

      baseSampleNames.push_back("CCQE");
      baseSampleCuts.push_back(!kIsCosmic && cutIsQE && cutIsCC);
      baseSampleSpillCuts.push_back(kNoSpillCut);

      baseSampleNames.push_back("NCQE");
      baseSampleCuts.push_back(!kIsCosmic && cutIsQE && cutIsNC);
      baseSampleSpillCuts.push_back(kNoSpillCut);


      baseSampleNames.push_back("CCRes");
      baseSampleCuts.push_back(!kIsCosmic && cutIsRes && cutIsCC);
      baseSampleSpillCuts.push_back(kNoSpillCut);

      baseSampleNames.push_back("NCRes");
      baseSampleCuts.push_back(!kIsCosmic && cutIsRes && cutIsNC);
      baseSampleSpillCuts.push_back(kNoSpillCut);


      baseSampleNames.push_back("CCDIS");
      baseSampleCuts.push_back(!kIsCosmic && cutIsDIS && cutIsCC);
      baseSampleSpillCuts.push_back(kNoSpillCut);

      baseSampleNames.push_back("NCDIS");
      baseSampleCuts.push_back(!kIsCosmic && cutIsDIS && cutIsNC);
      baseSampleSpillCuts.push_back(kNoSpillCut);


      baseSampleNames.push_back("CCCoh");
      baseSampleCuts.push_back(!kIsCosmic && cutIsCoh && cutIsCC);
      baseSampleSpillCuts.push_back(kNoSpillCut);

      baseSampleNames.push_back("NCCoh");
      baseSampleCuts.push_back(!kIsCosmic && cutIsCoh && cutIsNC);
      baseSampleSpillCuts.push_back(kNoSpillCut);


      baseSampleNames.push_back("CCCohElastic");
      baseSampleCuts.push_back(!kIsCosmic && cutIsCohElastic && cutIsCC);
      baseSampleSpillCuts.push_back(kNoSpillCut);

      baseSampleNames.push_back("NCCohElastic");
      baseSampleCuts.push_back(!kIsCosmic && cutIsCohElastic && cutIsNC);
      baseSampleSpillCuts.push_back(kNoSpillCut);


      baseSampleNames.push_back("CCElectronScattering");
      baseSampleCuts.push_back(!kIsCosmic && cutIsElectronScattering && cutIsCC);
      baseSampleSpillCuts.push_back(kNoSpillCut);

      baseSampleNames.push_back("NCElectronScattering");
      baseSampleCuts.push_back(!kIsCosmic && cutIsElectronScattering && cutIsNC);
      baseSampleSpillCuts.push_back(kNoSpillCut);


      baseSampleNames.push_back("CCIMDAnnihilation");
      baseSampleCuts.push_back(!kIsCosmic && cutIsIMDAnnihilation && cutIsCC);
      baseSampleSpillCuts.push_back(kNoSpillCut);

      baseSampleNames.push_back("NCIMDAnnihilation");
      baseSampleCuts.push_back(!kIsCosmic && cutIsIMDAnnihilation && cutIsNC);
      baseSampleSpillCuts.push_back(kNoSpillCut);


      baseSampleNames.push_back("CCInverseBetaDecay");
      baseSampleCuts.push_back(!kIsCosmic && cutIsInverseBetaDecay && cutIsCC);
      baseSampleSpillCuts.push_back(kNoSpillCut);

      baseSampleNames.push_back("NCInverseBetaDecay");
      baseSampleCuts.push_back(!kIsCosmic && cutIsInverseBetaDecay && cutIsNC);
      baseSampleSpillCuts.push_back(kNoSpillCut);


      baseSampleNames.push_back("CCGlashowResonance");
      baseSampleCuts.push_back(!kIsCosmic && cutIsGlashowResonance && cutIsCC);
      baseSampleSpillCuts.push_back(kNoSpillCut);

      baseSampleNames.push_back("NCGlashowResonance");
      baseSampleCuts.push_back(!kIsCosmic && cutIsGlashowResonance && cutIsNC);
      baseSampleSpillCuts.push_back(kNoSpillCut);


      baseSampleNames.push_back("CCAMNuGamma");
      baseSampleCuts.push_back(!kIsCosmic && cutIsAMNuGamma && cutIsCC);
      baseSampleSpillCuts.push_back(kNoSpillCut);

      baseSampleNames.push_back("NCAMNuGamma");
      baseSampleCuts.push_back(!kIsCosmic && cutIsAMNuGamma && cutIsNC);
      baseSampleSpillCuts.push_back(kNoSpillCut);


      baseSampleNames.push_back("CCMEC");
      baseSampleCuts.push_back(!kIsCosmic && cutIsMEC && cutIsCC);
      baseSampleSpillCuts.push_back(kNoSpillCut);

      baseSampleNames.push_back("NCMEC");
      baseSampleCuts.push_back(!kIsCosmic && cutIsMEC && cutIsNC);
      baseSampleSpillCuts.push_back(kNoSpillCut);


      baseSampleNames.push_back("CCDiffractive");
      baseSampleCuts.push_back(!kIsCosmic && cutIsDiffractive && cutIsCC);
      baseSampleSpillCuts.push_back(kNoSpillCut);

      baseSampleNames.push_back("NCDiffractive");
      baseSampleCuts.push_back(!kIsCosmic && cutIsDiffractive && cutIsNC);
      baseSampleSpillCuts.push_back(kNoSpillCut);


      baseSampleNames.push_back("CCEM");
      baseSampleCuts.push_back(!kIsCosmic && cutIsEM && cutIsCC);
      baseSampleSpillCuts.push_back(kNoSpillCut);

      baseSampleNames.push_back("NCEM");
      baseSampleCuts.push_back(!kIsCosmic && cutIsEM && cutIsNC);
      baseSampleSpillCuts.push_back(kNoSpillCut);


      baseSampleNames.push_back("CCWeakMix");
      baseSampleCuts.push_back(!kIsCosmic && cutIsWeakMix && cutIsCC);
      baseSampleSpillCuts.push_back(kNoSpillCut);

      baseSampleNames.push_back("NCWeakMix");
      baseSampleCuts.push_back(!kIsCosmic && cutIsWeakMix && cutIsNC);
      baseSampleSpillCuts.push_back(kNoSpillCut);

      baseSampleNames.push_back("UnknownInteractionType1");
      baseSampleCuts.push_back(!kIsCosmic && cutIsUnknownInteractionType1);
      baseSampleSpillCuts.push_back(kNoSpillCut);

      baseSampleNames.push_back("UnknownInteractionType2");
      baseSampleCuts.push_back(!kIsCosmic && cutIsUnknownInteractionType2);
      baseSampleSpillCuts.push_back(kNoSpillCut);

      baseSampleNames.push_back("UnknownInteractionType3");
      baseSampleCuts.push_back(!kIsCosmic && cutIsUnknownInteractionType3);
      baseSampleSpillCuts.push_back(kNoSpillCut);

    }

  }
  else{

  }

  //==== Define Cuts

  vector<Cut> cuts;
  vector<SpillCut> spillcuts;
  vector<TString> cutNames;

  //====   1) Numu-cc selection

  cuts.push_back( kNoCut );
  spillcuts.push_back( kNoSpillCut );
  cutNames.push_back( "NoCut" );

  if(runCutflow){

    cuts.push_back( kNoCut && cutRFiducial );
    spillcuts.push_back( kNoSpillCut );
    cutNames.push_back( "RecoFiducial" );

    cuts.push_back( kNoCut && cutRFiducial && kNotClearCosmic );
    spillcuts.push_back( kNoSpillCut );
    cutNames.push_back( "NotClearCosmic" );

    cuts.push_back( kNoCut && cutRFiducial && kNotClearCosmic && cutNuScore );
    spillcuts.push_back( kNoSpillCut );
    cutNames.push_back( "NuScore" );

    cuts.push_back( kNoCut && cutRFiducial && kNotClearCosmic && cutNuScore && cutFMScore );
    spillcuts.push_back( kNoSpillCut );
    cutNames.push_back( "FMScore" );

    cuts.push_back( kNoCut && cutRFiducial && kNotClearCosmic && cutNuScore && cutFMScore );
    spillcuts.push_back( cut_SideCRTVeto );
    cutNames.push_back( "SideCRTVeto" );

    cuts.push_back( kNoCut && cutRFiducial && kNotClearCosmic && cutNuScore && cutFMScore );
    spillcuts.push_back( cut_CRTVeto );
    cutNames.push_back( "CRTVeto" );

  }

  cuts.push_back( kNoCut && cutRFiducial && kNotClearCosmic && cutNuScore && cutFMScore && cutHasMuon );
  spillcuts.push_back( cut_CRTVeto );
  cutNames.push_back( "HasMuon" );

  cuts.push_back( kNoCut && cutRFiducial && kNotClearCosmic && cutNuScore && cutFMScore && cutHasMuon && !cutZeroMomentum );
  spillcuts.push_back( cut_CRTVeto );
  cutNames.push_back( "HasMuon_NonZeroP" );

  //====     Contained/Exiting

  cuts.push_back( kNoCut && cutRFiducial && kNotClearCosmic && cutNuScore && cutFMScore && cutHasMuon && cutMuonContained );
  spillcuts.push_back( cut_CRTVeto );
  cutNames.push_back( "HasMuon_Contained" );
  cuts.push_back( kNoCut && cutRFiducial && kNotClearCosmic && cutNuScore && cutFMScore && cutHasMuon && !cutMuonContained );
  spillcuts.push_back( cut_CRTVeto );
  cutNames.push_back( "HasMuon_Exiting" );
  cuts.push_back( kNoCut && cutRFiducial && kNotClearCosmic && cutNuScore && cutFMScore && cutHasMuon && !cutMuonContained );
  spillcuts.push_back( cut_CRTVeto );
  cutNames.push_back( "HasMuon_Exiting_NonZeroP" );
  cuts.push_back( kNoCut && cutRFiducial && kNotClearCosmic && cutNuScore && cutFMScore && cutHasMuon && !cutMuonContained && cutZeroMomentum );
  spillcuts.push_back( cut_CRTVeto );
  cutNames.push_back( "HasMuon_Exiting_ZeroP" );

  //====     HasProton added
  cuts.push_back( kNoCut && cutRFiducial && kNotClearCosmic && cutNuScore && cutFMScore && cutHasMuon && cutHasProton );
  spillcuts.push_back( cut_CRTVeto );
  cutNames.push_back( "HasMuon_HasProton" );

  //==== No CRT version for run5679

  cuts.push_back( kNoCut && cutRFiducial && kNotClearCosmic && cutNuScore && cutFMScore && cutHasMuon );
  spillcuts.push_back( kNoSpillCut );
  cutNames.push_back( "NoCRT_HasMuon" );

  //====     Contained/Exiting

  cuts.push_back( kNoCut && cutRFiducial && kNotClearCosmic && cutNuScore && cutFMScore && cutHasMuon && cutMuonContained );
  spillcuts.push_back( kNoSpillCut );
  cutNames.push_back( "NoCRT_HasMuon_Contained" );
  cuts.push_back( kNoCut && cutRFiducial && kNotClearCosmic && cutNuScore && cutFMScore && cutHasMuon && !cutMuonContained );
  spillcuts.push_back( kNoSpillCut );
  cutNames.push_back( "NoCRT_HasMuon_Exiting" );

  //====     HasProton added
  cuts.push_back( kNoCut && cutRFiducial && kNotClearCosmic && cutNuScore && cutFMScore && cutHasMuon && cutHasProton );
  spillcuts.push_back( kNoSpillCut );
  cutNames.push_back( "NoCRT_HasMuon_HasProton" );




  //==== Declare HistoProducer

  HistoProducer m;
  m.outputDir = outDirBase;
  m.fillNominal = fillNominal;
  m.fillBeamInfo = fillBeamInfo;

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
