#!/bin/bash

echo "@@ Setting NuMINumuXSec..."

export NUMIXSECWD=`pwd`
export PostJobWD=${NUMIXSECWD}/PostJob/
export PlotDir=${PostJobWD}/plots/
export gridWD=${NUMIXSECWD}/GridSub/

export PYTHONPATH=${PYTHONPATH}:${PostJobWD}/Libs/:${PostJobWD}/Configs/:${gridWD}/python/:${gridWD}/Configs/



export ANANAME=v09_75_03
export gridJobDir=/pnfs/icarus/scratch/users/${USER}/NuMINumuXSec/GridJob/
export gridBinDir=${gridWD}/bin/
export gridLibDir=${gridWD}/lib/
#export gridLibDirPNFS=/pnfs/icarus/resilient/users/jskim/NuMINumuXSec/${ANANAME}/
export gridLibDirPNFS=/pnfs/icarus/scratch/users/jskim/NuMINumuXSec/tars/${ANANAME}/
export gridDataDir=${gridWD}/data/
export gridOutputDir=/pnfs/icarus/persistent/users/${USER}/NuMINumuXSec/


export PATH=${gridBinDir}:${PYTHONDIR}:${PATH}

echo "ANANAME: "${ANANAME}
echo "gridLibDirPNFS: "${gridLibDirPNFS}

