import os
import method

required_programs = ['MMPBSA.py', 'cpptraj', 'babel']

class Mmgbsa(method.ScoringMethod):

    def write_rescoring_script(self, filename, file_r, file_l):

        locals().update(self.options)

        with open(filename, 'w') as file:
            script ="""#!/bin/bash
set -e
# create temporary folder
rm -rf tmp; mkdir tmp
# copy ligand file
cp %(file_l)s tmp/lig0.mol2
echo "from MOBPred.tools import mol2
from MOBPred.amber import minimz as mn
mol2.update_mol2file('%(file_l)s', 'tmp/lig0.mol2', ligname='LIG')
mn.prepare_receptor('tmp/rec.pdb', '%(file_r)s', False)" > prepare_files1.py
python prepare_files1.py
#rm -rf prepare_files1.py
cd tmp
# run antechamber
antechamber -fi mol2 -i lig0.mol2 -fo mol2 -o lig.mol2
# run parmchk
parmchk -i lig.mol2 -f mol2 -o lig.frcmod
babel lig0.mol2 lig.pdb
# prepare complex.pdb
echo "with open('complex.pdb', 'w') as cpdb:
    # copy contents of rec.pdb
    with open('rec.pdb', 'r') as rpdb:
        for line in rpdb:
            if line.startswith(('ATOM','TER')):
                cline = line
                cpdb.write(line)
    if not cline.startswith('TER'):
        cpdb.write('TER\\n')
    # copy contents of lig.pdb
    with open('lig.pdb', 'r') as lpdb:
        for line in lpdb:
            if line.startswith(('ATOM','HETATM','TER')):
                cpdb.write(line)
    cpdb.write('END\\n')" > prepare_files2.py
python prepare_files2.py
echo "source leaprc.ff14SB
source leaprc.gaff
LIG = loadmol2 lig.mol2
loadamberparams lig.frcmod
t = loadpdb rec.pdb
p = loadpdb complex.pdb
saveamberparm LIG lig.prmtop lig.inpcrd
saveamberparm t rec.prmtop rec.inpcrd
saveamberparm p complex.prmtop complex.inpcrd
savePdb p complex.pdb
quit" > leap.in
tleap -f leap.in > leap.log
echo "Sample input file for GB calculation
&general
startframe=1, endframe=1, interval=1,
verbose=2, keep_files=0,
/
&gb
igb=2, saltcon=0.150,
/" > mm.in
MMPBSA.py -O -i mm.in -o mm.out -sp complex.prmtop -cp complex.prmtop -rp rec.prmtop -lp lig.prmtop -y complex.pdb > mm.log
cd .."""% locals()
            file.write(script)

    def extract_rescoring_results(self, file_s):

        with open(file_s, 'a') as sf:
            if os.path.isfile('tmp/mm.out'):
                with open('tmp/mm.out', 'r') as mmf:
                    for line in mmf:
                        if line.startswith("DELTA TOTAL"):
                            print >> sf, line.split()[2]
                            break
            else:
                print >> sf, 'NaN'

    def cleanup(self):
        pass
