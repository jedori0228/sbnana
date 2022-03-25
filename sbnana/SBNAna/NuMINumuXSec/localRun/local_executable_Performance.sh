# Setup grid submission

whichInputFile=$1
echo "@@ whichInputFile : ${whichInputFile}"

jobName="Performance"

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

outDir=/pnfs/icarus/persistent/users/jskim/NuMINumuXSec/${jobName}/${inputSampleName}/${sampleProdRelease}/

echo "@@ outDir : ${outDir}"
mkdir -p ${outDir}

echo "@@ date :"
date
echo "@@ Running cafe -bq run_Performance.C ${whichInputFile} \"${outDir}\""
cafe -bq run_Performance.C ${whichInputFile} "${outDir}"
echo "@@ date :"
date
echo "@@ Checking output file"
ls -alh ${outDir}/*

echo "@@ Done!"
