#!/bin/bash

# submit jobs
for dir in lig{1..8}; do
  cd $dir
  qsub ../run.pbs # submit job
  cd ..
done
