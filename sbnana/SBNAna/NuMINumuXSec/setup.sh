#!/bin/bash

echo "@@ Setting NuMINumuXSec..."

export NUMIXSECWD=`pwd`
export PostJobWD=${NUMIXSECWD}/PostJob/
export gridWD=${NUMIXSECWD}/gridSubmission/

export PYTHONPATH=${PYTHONPATH}:${PostJobWD}/Libs/:${PostJobWD}/Configs/:${gridWD}/python/:${gridWD}/Configs/



export gridJobDir=/pnfs/icarus/scratch/users/${USER}/NuMINumuXSec/GridJob/
export gridBinDir=${gridWD}/bin/
export gridLibDir=${gridWD}/lib/
export gridDataDir=${gridWD}/data/
export gridOutputDir=/pnfs/icarus/persistent/users/${USER}/NuMINumuXSec/


export PATH=${gridBinDir}:${PYTHONDIR}:${PATH}
