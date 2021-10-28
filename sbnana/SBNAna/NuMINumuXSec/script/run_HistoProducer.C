#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_HistoProducer.h"
#include "sbnana/SBNAna/Cuts/NumuCutsIcarus202106.h"
#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_Variables.h"
#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_Cuts.h"
#include "sbnana/SBNAna/Cuts/Cuts.h"

using namespace ana;
using namespace ICARUSNumuXsec;

void run_HistoProducer(int whichInputFile, TString outDirBase){

  //==== Inputfile

  string inputDef = "";
  bool isDataInput = false;
  bool isCosmicInput = false;
  double TargetPOT = 6.6e20;
  TString str_TargetPOT = "6.6e20 POT";
  //====   BNB, new central sample, IcarusProd2021B_BNB_Nu_Cosmics_v09_28_01_01_01_caf
  if(whichInputFile==0){
    cout << "[run_Efficiency] Running BNB, IcarusProd2021B_BNB_Nu_Cosmics_v09_28_01_01_01_caf" << endl;
    inputDef = "IcarusProd2021B_BNB_Nu_Cosmics_v09_28_01_01_01_caf";
  }
  //====   NuMI, new central sample, IcarusProd2021B_NuMI_Nu_Cosmics_v09_28_01_01_01_caf
  else if(whichInputFile==1){
    cout << "[run_Efficiency] Running NuMI, IcarusProd2021B_NuMI_Nu_Cosmics_v09_28_01_01_01_caf" << endl;
    inputDef = "IcarusProd2021B_NuMI_Nu_Cosmics_v09_28_01_01_01_caf";
    TargetPOT = 6e20;
    str_TargetPOT = "6e20 POT";
  }
  //====   Cosmic, IcarusProd2021B_Intime_Cosmic_v09_28_01_01_01_caf
  else if(whichInputFile==2){
    isCosmicInput = true;
    cout << "[run_Efficiency] Running In-time BNB cosmic" << endl;
    inputDef = "IcarusProd2021B_Intime_Cosmic_v09_28_01_01_01_caf";
  }
  //====   Cosmic Data, howard_howard_numi_cosmic_run6480_full_2_CAFMaker_20210913_224007_821563
  else if(whichInputFile==3){
    isCosmicInput = true;
    isDataInput = true;
    cout << "[run_Efficiency] Running Cosmic data" << endl;
    inputDef = "howard_howard_numi_cosmic_run6480_full_2_CAFMaker_20210913_224007_821563";
    TargetPOT = 6e20;
    str_TargetPOT = "6e20 POT";
  }
  //====   cosmics in-time (is a little old)
  else if(whichInputFile==4){
    isCosmicInput = true;
    cout << "[run_Efficiency] Running /pnfs/icarus/persistent/users/dmendez/Tests/cosmics_hadded0.flatcaf.root" << endl;
    inputDef = "/pnfs/icarus/persistent/users/dmendez/Tests/cosmics_hadded0.flatcaf.root";
  }
  //====   cosmics out-time (release v09_17_01):
  else if(whichInputFile==5){
    isCosmicInput = true;
    cout << "[run_Efficiency] Running /icarus/data/users/mmr/CAFMaker/v09_20_00/NuMIBeamWindow/icarus_numi_offaxis_cosmic_v09_17_01/icarus_numi_offaxis_cosmic_v09_17_01_CAFMaker_out_flat.root" << endl;
    inputDef = "/icarus/data/users/mmr/CAFMaker/v09_20_00/NuMIBeamWindow/icarus_numi_offaxis_cosmic_v09_17_01/icarus_numi_offaxis_cosmic_v09_17_01_CAFMaker_out_flat.root";
    TargetPOT = 6e20;
    str_TargetPOT = "6e20 POT";
  }
  //====   old NuMI sample
  else if(whichInputFile==6){
    cout << "[run_Efficiency] Running /icarus/data/users/mmr/CAFMaker/v09_20_00/NuMIBeamWindow/icarus_numi_offaxis_cosmic_v09_17_01/icarus_numi_offaxis_cosmic_v09_17_01_CAFMaker_out_flat.root" << endl;
    inputDef = "/icarus/data/users/mmr/CAFMaker/v09_20_00/NuMIBeamWindow/icarus_numi_offaxis_cosmic_v09_17_01/icarus_numi_offaxis_cosmic_v09_17_01_CAFMaker_out_flat.root";
    TargetPOT = 6e20;
    str_TargetPOT = "6e20 POT";
  }
  else if(whichInputFile==7){
    cout << "[run_Efficiency] Running ReCAF NuMI+Cosmic, jsmedley_jsmedley_numi_cafmaker_CAFMaker_20210913_224007_821563" << endl;
    inputDef = "jsmedley_jsmedley_numi_cafmaker_CAFMaker_20210913_224007_821563";
    //inputDef = "/pnfs/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_*.root";
    TargetPOT = 6e20;
    str_TargetPOT = "6e20 POT";
  }

  cout << "[run_Efficiency] Input definition : " << inputDef << endl;
  vector<string> vecf = {
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_0.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_1.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_10.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_100.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_101.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_102.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_11.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_12.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_13.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_14.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_15.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_16.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_17.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_18.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_19.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_2.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_20.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_21.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_22.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_23.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_24.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_25.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_26.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_27.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_28.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_29.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_3.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_30.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_31.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_32.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_33.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_34.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_35.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_36.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_37.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_38.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_39.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_4.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_40.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_41.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_42.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_43.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_44.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_45.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_46.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_47.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_48.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_49.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_5.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_50.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_51.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_52.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_53.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_54.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_55.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_56.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_57.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_58.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_59.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_6.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_60.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_61.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_62.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_63.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_64.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_65.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_66.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_67.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_68.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_69.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_7.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_70.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_71.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_72.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_73.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_74.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_75.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_76.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_77.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_78.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_79.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_8.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_80.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_81.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_82.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_83.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_84.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_85.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_86.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_87.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_88.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_89.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_9.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_90.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_91.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_92.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_93.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_94.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_95.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_96.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_97.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_98.root",
"root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_99.root",
};
  //SpectrumLoader loader(inputDef);
  SpectrumLoader loader(vecf);
  //SpectrumLoader loader("root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_0.root");

  //==== Define samples;
  vector<TString> baseSampleNames;
  vector<Cut> baseSampleCuts;
  vector<SpillCut> baseSampleSpillCuts;

  baseSampleNames.push_back("AllSamples");
  baseSampleCuts.push_back(kNoCut);
  baseSampleSpillCuts.push_back(kNoSpillCut);

  if(!isDataInput && !isCosmicInput){


    baseSampleNames.push_back("NuMuCC");
    baseSampleCuts.push_back(kIsNuMuCC);
    baseSampleSpillCuts.push_back(kNoSpillCut);

    baseSampleNames.push_back("NuMuNC");
    baseSampleCuts.push_back(cutIsNuMuNC);
    baseSampleSpillCuts.push_back(kNoSpillCut);

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

    baseSampleNames.push_back("Cosmic");
    baseSampleCuts.push_back(kIsCosmic);
    baseSampleSpillCuts.push_back(kNoSpillCut);

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

  cuts.push_back( kNoCut && cutRFiducial );
  spillcuts.push_back( kNoSpillCut );
  cutNames.push_back( "RecoFiducial" );

  cuts.push_back( kNoCut && cutRFiducial && kNotClearCosmic );
  spillcuts.push_back( kNoSpillCut );
  cutNames.push_back( "NotClearCosmic" );

  cuts.push_back( kNoCut && cutRFiducial && kNotClearCosmic && kNuScore );
  spillcuts.push_back( kNoSpillCut );
  cutNames.push_back( "NuScoreGT0p4" );

  cuts.push_back( kNoCut && cutRFiducial && kNotClearCosmic && kNuScore && cutFMScore );
  spillcuts.push_back( kNoSpillCut );
  cutNames.push_back( "FMScoreLT6" );

  cuts.push_back( kNoCut && cutRFiducial && kNotClearCosmic && kNuScore && cutFMScore );
  spillcuts.push_back( spillcutSideCRTHitVetoFD );
  cutNames.push_back( "SideCRTVeto" );

  cuts.push_back( kNoCut && cutRFiducial && kNotClearCosmic && kNuScore && cutFMScore );
  spillcuts.push_back( kCRTHitVetoFD );
  cutNames.push_back( "CRTVeto" );

  cuts.push_back( kNoCut && cutRFiducial && kNotClearCosmic && kNuScore && cutFMScore && cutHasMuon );
  spillcuts.push_back( kCRTHitVetoFD );
  cutNames.push_back( "HasMuon" );

  //====     Contained/Exiting

  cuts.push_back( kNoCut && cutRFiducial && kNotClearCosmic && kNuScore && cutFMScore && cutHasMuon && cutMuonContained );
  spillcuts.push_back( kCRTHitVetoFD );
  cutNames.push_back( "HasMuon_Contained" );
  cuts.push_back( kNoCut && cutRFiducial && kNotClearCosmic && kNuScore && cutFMScore && cutHasMuon && !cutMuonContained );
  spillcuts.push_back( kCRTHitVetoFD );
  cutNames.push_back( "HasMuon_Exiting" );

/*
  if(!isDataInput && !isCosmicInput){

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

  }
*/
  //====     HasProton added
  cuts.push_back( kNoCut && cutRFiducial && kNotClearCosmic && kNuScore && cutFMScore && cutHasMuon && cutHasProton );
  spillcuts.push_back( kCRTHitVetoFD );
  cutNames.push_back( "HasMuon_HasProton" );

  //==== Declare HistoProducer

  HistoProducer m;
  m.outputDir = outDirBase;
  m.initialize();
  m.setSystematicWeights(); // TODO Run this to get systematic variations


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

        if(isCosmic){

          m.bookSlice(loader, mySpillCut, myCut);
          if(cutName.Contains("HasMuon")){
            m.bookRecoMuon(loader, mySpillCut, myCut);
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

          m.bookSlice(loader, mySpillCut, myCut);
          m.bookTruth(loader, mySpillCut, myCut);
          if(cutName.Contains("HasMuon")){
            m.bookRecoMuon(loader, mySpillCut, myCut);
          }
          if(cutName.Contains("HasProton")){
            m.bookRecoProton(loader, mySpillCut, myCut);
            m.bookRecoNeutrino(loader, mySpillCut, myCut);
          }

        }

      } // END if setCut

    } // END loop cut

  } // END sample cut

  loader.Go();

  m.saveHistograms();

}
