# Setup grid submission

whichInputFile=$1
echo "@@ whichInputFile : ${whichInputFile}"
runType=$2
echo "@@ runType : ${runType}"

jobName="HistoProducer"

#selectionDir="FMScoreRemoved_NuScoreRemoved"
#selectionDir="FMScoreLT12_NuScoreRemoved"
#selectionDir="FMScoreLT12_NuScoreRemoved_SideCRTOnly"
#selectionDir="FMScoreLT12_NuScoreRemoved_CRTNuMIWindow_SideCRTOnly"
selectionDir="FMScoreLT12_NuScoreRemoved_CRTNuMIWindow"
#selectionDir="FMScoreLT12_NuScoreGT0p4_CRTNuMIWindow"
#selectionDir="FMScoreLT12_NuScoreGT0p5_CRTNuMIWindow"

inputSampleName=""
sampleProdRelease="NULL"
if [[ "$whichInputFile" == 0 ]]; then
  echo "@@ [Input] NUMI_Nu_Cosmics"
  inputSampleName="NUMI_Nu_Cosmics"
  sampleProdRelease="v09_45_00"
fi
if [[ "$whichInputFile" == 1 ]]; then
  echo "@@ [Input] NUMI_in-time_Cosmics"
  inputSampleName="NUMI_in-time_Cosmics"
  sampleProdRelease="v09_45_00"
fi
if [[ "$whichInputFile" == 2 ]]; then
  echo "@@ [Input] NUMI_Nu_Cosmics"
  inputSampleName="NUMI_Nu_Cosmics"
  sampleProdRelease="v09_37_01_03p01"
fi
if [[ "$whichInputFile" == 3 ]]; then
  echo "@@ [Input] NUMI_in-time_Cosmics"
  inputSampleName="NUMI_in-time_Cosmics"
  sampleProdRelease="v09_37_01_03p01"
fi
if [[ "$whichInputFile" == 4 ]]; then
  echo "@@ [Input] NUMI_Nu_Cosmics, NoSEC"
  inputSampleName="NUMI_Nu_Cosmics"
  sampleProdRelease="v09_45_00_NoSCE"
fi
if [[ "$whichInputFile" == 5 ]]; then
  echo "@@ [Input] run7568, NUMI"
  inputSampleName="Run7568_NUMI"
  sampleProdRelease="v09_45_00"
fi
if [[ "$whichInputFile" == 6 ]]; then
  echo "@@ [Input] run7568, OffbeamNUMI"
  inputSampleName="Run7568_OffbeamNUMI"
  sampleProdRelease="v09_45_00"
fi

outputNameToTransfer="output.root"
if [[ "$runType" == 0 ]]; then
  echo "@@ [runType]  0 : Only nominal;     Incl+PerInt.; +Cutflow"
  outputNameToTransfer="output__Nominal__Incl_PerInt__CutFlow.root"
fi
if [[ "$runType" == 11 ]]; then
  echo "@@ [runType] 11 : Only systematics; Incl        ; X"
  outputNameToTransfer="output__Systematic__Incl.root"
fi
if [[ "$runType" == 21 ]]; then
  echo "@@ [runType] 21 : Only systematics;      PerInt.; X"
  outputNameToTransfer="output__Systematic__PerInt.root"
fi

outDir=/pnfs/icarus/persistent/users/jskim/NuMINumuXSec/${jobName}/${selectionDir}/${inputSampleName}/${sampleProdRelease}/

echo "@@ outDir : ${outDir}"
mkdir -p ${outDir}

echo "@@ date :"
date
echo "@@ Running cafe -bq run_HistoProducer_Master.C ${whichInputFile} ${runType} \"${outDir}\""
cafe -bq run_HistoProducer_Master.C ${whichInputFile} ${runType} "${outDir}"
echo "@@ date :"
date
echo "@@ Checking output file"
ls -alh ${outDir}/*

mv ${outDir}/output.root ${outDir}/${outputNameToTransfer}

echo "@@ Done!"
