import sys
import shutil
from tempfile import mkstemp

from mdkit.amber.ambertools import load_atomic_ions

# first all residues are supposed to be recognized
are_unrecognized_residues = False

# check if and which atoms were not recognized
unrecognized_residues = []
with open(sys.argv[2], 'r') as logf:
    for line in logf:
        if line.startswith('Sorry, there are no Gasteiger parameters available for atom'):
            are_unrecognized_residues = True
            resname = line.split()[-1].split(':')[0]
            resname = ''.join([i for i in resname if not i.isdigit()])
            unrecognized_residues.append(resname)

if are_unrecognized_residues:

    ions_amber = load_atomic_ions()
    print "No charges specified for ion(s) " + ', '.join(unrecognized_residues)
    print "Attributing formal charges..."

    # update .pdbqt file for the receptor
    fh, abs_path = mkstemp()

    with open(abs_path, 'w') as tempf:
        with open(sys.argv[1], 'r') as ff:

            for line in ff:
                is_ion = False

                if line.startswith(('ATOM', 'HETATM')):
                    resname = line[17:20].strip()
                    if resname in unrecognized_residues:
                        assert resname in ions_amber
                        charge = "%.3f"%ions_amber[resname]
                        is_ion = True

                if is_ion:
                    tempf.write(line[:70] + ' '*(6-len(charge)) + charge + line[76:])
                else:
                    tempf.write(line)

    shutil.move(abs_path, sys.argv[1])