#!/usr/bin/python
import os
import sys
import shutil
import stat
from glob import glob
import subprocess

from mdtools.utility import reader
from mdtools.utility import mol2
from mdtools.amber import ambertools
from mdtools.amber import clustr

if not os.path.exists('poses'):
    print "No poses directory detected. Quitting..."
    sys.exit()

pwd = os.getcwd()
pdbid = pwd.split('/')[-1]

shutil.rmtree('rmsd', ignore_errors=True)
os.mkdir('rmsd')
os.chdir('rmsd')

# prepare receptor
ambertools.prepare_receptor('protein.pdb', glob('../*_protein_min.pdb')[0])

if not os.path.isfile('../'+pdbid+'_ligand_min_u.mol2'):
    mol2.update_mol2file('../'+pdbid+'_ligand_min.mol2', '../'+pdbid+'_ligand_min_u.mol2', unique=True, ligname='LIG')

shutil.copyfile(glob('../*_ligand_min_u.mol2')[0], 'ligand.mol2')
ambertools.prepare_ligand('protein.pdb', 'ligand.mol2', 'protein-ligand.pdb')
n_mol2files = len(glob('../poses/*.mol2'))
files_l = []
for idx in range(n_mol2files):
   files_l.append(os.path.abspath('../poses/lig-%i.mol2'%(idx+1)))

os.mkdir('PDB')
files_rl = []
for idx, file_l in enumerate(files_l):
    file_rl = 'PDB/protein-ligand-%s.pdb'%(idx+1)
    ambertools.prepare_ligand('protein.pdb', file_l, file_rl)
    files_rl.append(file_rl)

clustr.prepare_leap_config_file('leap.in', ['protein.pdb'], ['ligand.mol2']+files_l, ['protein-ligand.pdb']+files_rl)
subprocess.check_output('tleap -f leap.in', shell=True)

lines_trajin = ""
for file_rl in files_rl:
    lines_trajin += "trajin %s\n"%(file_rl)

# remove last \n
lines_trajin = lines_trajin[:-1]

rmol2 = reader.open('ligand.mol2')
ligname = rmol2.ligname

# write cpptraj config file to cluster frames
with open('cpptraj.in', 'w') as file:
    contents ="""parm protein-ligand.prmtop
trajin protein-ligand.pdb
%(lines_trajin)s
rms first @CA,C,N&!:%(ligname)s
rms nofit :%(ligname)s&!@H= out rmsd.dat
"""% locals()
    file.write(contents)
subprocess.check_output('cpptraj -i cpptraj.in > cpptraj.log', shell=True)
