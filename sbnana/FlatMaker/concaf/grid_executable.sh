# Setup grid submission

outDir=$1
echo "@@ outDir : ${outDir}"

echo "@@ pwd"
pwd
echo "@@ mkdir output"
mkdir output
echo "@@ Done!"
thisOutputCreationDir=`pwd`/output/

nProcess=$PROCESS
echo "@@ nProcess : "${nProcess}

echo "@@ setup_icarus.sh"
source /cvmfs/icarus.opensciencegrid.org/products/icarus/setup_icarus.sh

setup icaruscode v09_37_02_04 -q e20:prof
setup sbnana v09_37_02_03 -q e20:prof

cd grid_dir/

echo "@@ outDir : "${outDir}
echo "@@ ifdh  mkdir_p "${outDir}
ifdh  mkdir_p ${outDir}

echo "@@ source run_"${nProcess}".sh "${thisOutputCreationDir}
source run_${nProcess}.sh ${thisOutputCreationDir}
echo "@@ Check output : "${thisOutputCreationDir}"/flatcaf_"${nProcess}".root"
ls -alh ${thisOutputCreationDir}"/flatcaf_"${nProcess}".root"

outFILE=${thisOutputCreationDir}/flatcaf_${nProcess}.root
if [ -f "$outFILE" ]; then
  echo "ifdh cp "${thisOutputCreationDir}"/flatcaf_"${nProcess}".root "${outDir}"/flatcaf_"${nProcess}".root"
  ifdh cp ${thisOutputCreationDir}/flatcaf_${nProcess}.root ${outDir}/flatcaf_${nProcess}.root
  echo "@@ Done!"
else
  echo "File not exist"
fi

