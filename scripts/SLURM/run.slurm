#!/bin/bash
#SBATCH --time=2-00:00
#SBATCH --partition=serial
#SBATCH --job-name="docking"
#SBATCH --cpus-per-task=1
#SBATCH --array=1-11

cd $SLURM_SUBMIT_DIR

cd rec$SLURM_ARRAY_TASK_ID
for isodir in isomer*; do
  cd $isodir
  rundock -f config.ini -l lig.mol2 -r rec.pdb
  cd ..
done
