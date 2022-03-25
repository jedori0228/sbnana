#!/bin/bash
for is in 0 1 5 6
#for is in 0 1
do
  source local_executable_HistoProducer_Master.sh ${is} 0
done
