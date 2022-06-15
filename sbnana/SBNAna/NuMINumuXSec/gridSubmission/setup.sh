#!/bin/bash

export gridWD=`pwd`
export gridJobDir=/pnfs/icarus/scratch/users/${USER}/NuMINumuXSec/GridJob
mkdir -p ${gridJobDir}
export gridBinDir=${gridWD}/bin/
export gridLibDir=${gridWD}/lib/
export gridDataDir=${gridWD}/data/
export gridOutputDir=/pnfs/icarus/persistent/users/${USER}/NuMINumuXSec/


export PYTHONPATH=${PYTHONPATH}:${gridWD}/python/:${gridWD}/Configs/
export PATH=${gridBinDir}:${PYTHONDIR}:${PATH}
