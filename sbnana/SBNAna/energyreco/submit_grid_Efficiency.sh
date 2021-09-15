# Setup grid submission

echo "@@ pwd"
pwd
echo "@@ mkdir output"
mkdir output
echo "@@ Done!"
thisOutputCreationDir=`pwd`/output/

nProcess=$PROCESS
echo "@@ nProcess : "${nProcess}

userArg=$1
echo "@@ userArg : "${userArg}

echo "@@ setup_icarus.sh"
source /cvmfs/icarus.opensciencegrid.org/products/icarus/setup_icarus.sh

codeVERSION=v09_28_02

cd grid_dir/
export MRB_BUILDDIR=$(pwd)
export MRB_PROJECT="larsoft"
export MRB_PROJECT_VERSION=${codeVERSION}
export MRB_QUALS="e20:prof"
export MRB_TOP=$(pwd)
export MRB_TOP_BUILD=$(pwd)
export MRB_SOURCE=$(pwd)
export MRB_INSTALL=$(pwd)
export PRODUCTS="${MRB_INSTALL}:${PRODUCTS}"
mrbslp

echo "@@ running, setup icaruscode v09_28_02 -q e20:prof, to use samweb"
setup icaruscode v09_28_02 -q e20:prof

outDir=""

Processing test.C...
if [[ "$nProcess" == 0 ]]; then
  outDir="NuMuCC"
fi
if [[ "$nProcess" == 1 ]]; then
  outDir="NuMuCC_RecoFiducial"
fi
if [[ "$nProcess" == 2 ]]; then
  outDir="NuMuCC_RecoFiducial_NotClearCosmic"
fi
if [[ "$nProcess" == 3 ]]; then
  outDir="NuMuCC_RecoFiducial_NotClearCosmic_NuScoreGT0p4"
fi
if [[ "$nProcess" == 4 ]]; then
  outDir="NuMuCC_RecoFiducial_NotClearCosmic_NuScoreGT0p4_FMScoreLT6"
fi
if [[ "$nProcess" == 5 ]]; then
  outDir="NuMuCC_RecoFiducial_NotClearCosmic_NuScoreGT0p4_FMScoreLT6_HasMuon"
fi
if [[ "$nProcess" == 6 ]]; then
  outDir="NuMuCC_RecoFiducial_NotClearCosmic_NuScoreGT0p4_FMScoreLT6_HasMuon_CRTVeto"
fi
if [[ "$nProcess" == 7 ]]; then
  outDir="NuMuCC_RecoFiducial_NotClearCosmic_NuScoreGT0p4_FMScoreLT6_HasMuon_CRTVeto_MatchedToMuon"
fi
if [[ "$nProcess" == 8 ]]; then
  outDir="NuMuCC_RecoFiducial_NotClearCosmic_NuScoreGT0p4_FMScoreLT6_HasMuon_CRTVeto_MatchedToProton"
fi
if [[ "$nProcess" == 9 ]]; then
  outDir="NuMuCC_RecoFiducial_NotClearCosmic_NuScoreGT0p4_FMScoreLT6_HasMuon_CRTVeto_Contained"
fi
if [[ "$nProcess" == 10 ]]; then
  outDir="NuMuCC_RecoFiducial_NotClearCosmic_NuScoreGT0p4_FMScoreLT6_HasMuon_CRTVeto_Exiting"
fi
if [[ "$nProcess" == 11 ]]; then
  outDir="NuMuCC_RecoFiducial_NotClearCosmic_NuScoreGT0p4_FMScoreLT6_HasMuon_CRTVeto_Contained_MatchedToMuon"
fi
if [[ "$nProcess" == 12 ]]; then
  outDir="NuMuCC_RecoFiducial_NotClearCosmic_NuScoreGT0p4_FMScoreLT6_HasMuon_CRTVeto_Exiting_MatchedToMuon"
fi
if [[ "$nProcess" == 13 ]]; then
  outDir="NuMuCC_RecoFiducial_NotClearCosmic_NuScoreGT0p4_FMScoreLT6_HasMuon_CRTVeto_Contained_MatchedToProton"
fi
if [[ "$nProcess" == 14 ]]; then
  outDir="NuMuCC_RecoFiducial_NotClearCosmic_NuScoreGT0p4_FMScoreLT6_HasMuon_CRTVeto_Exiting_MatchedToProton"
fi
if [[ "$nProcess" == 15 ]]; then
  outDir="NuMuCC_RecoFiducial_NotClearCosmic_NuScoreGT0p4_FMScoreLT6_HasMuon_CRTVeto_HasProton"
fi
if [[ "$nProcess" == 16 ]]; then
  outDir="NuMuCC_RecoFiducial_NotClearCosmic_NuScoreGT0p4_FMScoreLT6_HasMuon_CRTVeto_HasProton_WithNeutron"
fi
if [[ "$nProcess" == 17 ]]; then
  outDir="NuMuCC_RecoFiducial_NotClearCosmic_NuScoreGT0p4_FMScoreLT6_HasMuon_CRTVeto_HasProton_WithChargedPion"
fi
if [[ "$nProcess" == 18 ]]; then
  outDir="NuMuCC_RecoFiducial_NotClearCosmic_NuScoreGT0p4_FMScoreLT6_HasMuon_CRTVeto_HasProton_WithNeutralPion"
fi
if [[ "$nProcess" == 19 ]]; then
  outDir="NuMuCC_RecoFiducial_NotClearCosmic_NuScoreGT0p4_FMScoreLT6_HasMuon_CRTVeto_HasProton_WithoutNeutronPion"
fi

selectionDir="MuonSel__MuonChi2LT30_and_ProtonChi2GT60"
#selectionDir="NoChi2"

if [[ "$userArg" == 0 ]]; then
  echo "@@ Running on BNB"
  outDir=/pnfs/icarus/scratch/users/jskim/CAFOutput/myAnalyzer/Efficinecy/${selectionDir}/BNB/${outDir}/
fi
if [[ "$userArg" == 1 ]]; then
  echo "@@ Running on BNB, ReCAF"
  outDir=/pnfs/icarus/scratch/users/jskim/CAFOutput/myAnalyzer/Efficinecy/${selectionDir}/BNBReCAF/${outDir}/
fi

echo "@@ outDir : "${outDir}
echo "@@ ifdh  mkdir_p "${outDir}
ifdh  mkdir_p ${outDir}

echo "@@ Running cafe -bq grid_dir/run_Efficiency.C ${nProcess} ${userArg} "${thisOutputCreationDir}""
cafe -bq run_Efficiency.C ${nProcess} ${userArg} "${thisOutputCreationDir}"
echo "@@ Copying output"
ifdh cp ${thisOutputCreationDir}/output.root ${outDir}/output.root
echo "@@ Done!"
