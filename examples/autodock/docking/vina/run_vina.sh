#!/bin/bash
set -e

MGLPATH=`which prepare_ligand4.py`
MGLPATH=`python -c "print '/'.join('$MGLPATH'.split('/')[:-3])"`
export PYTHONPATH=$PYTHONPATH:$MGLPATH

# prepare ligand
prepare_ligand4.py -l /home/ciniero/src/DockBox/examples/autodock/docking/1a30_ligand_dbx.mol2 -o lig.pdbqt
python check_ligand_pdbqt.py lig.pdbqt

# prepare receptor
prepare_receptor4.py -U nphs_lps_waters -r /home/ciniero/src/DockBox/examples/autodock/docking/1a30_protein.pdb -o target.pdbqt &> prepare_receptor4.log
python check_ions.py target.pdbqt prepare_receptor4.log

# run vina
vina --config vina.config 1> vina.out 2> vina.err