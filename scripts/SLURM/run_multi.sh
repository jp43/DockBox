#!/bin/bash

# submit jobs
for dir in lig{1..8}; do
  cd $dir
  sbatch ../run.slurm # submit job
  cd ..
done
