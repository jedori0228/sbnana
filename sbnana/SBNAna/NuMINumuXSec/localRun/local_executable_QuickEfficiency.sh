# Setup grid submission

whichInputFile=$1
echo "@@ whichInputFile : ${whichInputFile}"

#selectionDir="FMScoreRemoved_NuScoreRemoved_DefaultCRTHitVeto"
#selectionDir="FMScoreRemoved_NuScoreRemoved_NUMIWindowCRTHitVeto"
#selectionDir="FMScoreRemoved_NuScoreRemoved_NUMIWindowPE150CRTHitVeto"
#selectionDir="FMScoreLT12_NuScoreRemoved"
selectionDir="FMScoreLT12_NuScoreRemoved_CRTNuMIWindow"

jobName="QuickEfficiency"

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

outDir=/pnfs/icarus/persistent/users/jskim/NuMINumuXSec/${jobName}/${selectionDir}/${inputSampleName}/${sampleProdRelease}/

echo "@@ outDir : ${outDir}"
mkdir -p ${outDir}

echo "@@ date :"
date
echo "@@ Running cafe -bq run_QuickEfficiency.C ${whichInputFile} \"${outDir}\""
cafe -bq run_QuickEfficiency.C ${whichInputFile} "${outDir}"
echo "@@ date :"
date
echo "@@ Checking output file"
ls -alh ${outDir}/*

echo "@@ Done!"
