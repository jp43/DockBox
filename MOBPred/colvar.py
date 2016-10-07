import os
import sys
import method

mandatory_settings = ['type']
known_types = ['distance', 'volume']
default_settings = {'type': 'distance', 'residues': None}

class Colvar(method.ScoringMethod):
    """ScoringMethod class to compute a collective variable from ligand/complex structure"""

    def __init__(self, instance, site, options):

        super(Colvar, self).__init__(instance, site, options)

        if self.options['type'] == 'distance':
            if 'residues' not in self.options:
                raise ValueError('Residues option not specified in instance %s'%self.instance)
            else:
                self.options['residues'] = '[' + self.options['residues'] + ']'

    def write_rescoring_script(self, filename, file_r, file_l):

        locals().update(self.options)

        if self.options['type'] == 'distance':
            with open(filename, 'w') as file:
                script ="""#!/bin/bash
set -e
echo "from MOBPred.tools import mol2, PDB
from MOBPred.tools import reader
from MOBPred.tools import util
import numpy as np

ligrd = reader.open('%(file_l)s')
coords_lig = [map(float,line[2:5]) for line in ligrd.next()['ATOM']]
coords_lig = np.array(coords_lig)
natoms_lig = len(coords_lig)

recrd = reader.open('%(file_r)s')
residue = %(residues)s

coords_rec = [map(float,line[4:7]) for line in recrd.next()['ATOM'] if line[3] in map(str,residue)]
coords_rec = np.array(coords_rec)
natoms_rec = len(coords_rec)

dist = np.zeros((natoms_lig, natoms_rec))
for idx, cl in enumerate(coords_lig):
    for jdx, cr in enumerate(coords_rec):
        dist[idx, jdx] = np.sqrt(np.sum((cl - cr)**2))

# write min distance
with open('cv.out', 'w') as ff:
    ff.write(str(dist.min()))" > get_distance.py 
python get_distance.py"""% locals()
                file.write(script)

        elif self.options['type'] == 'volume':
            with open(filename, 'w') as file:
                script ="""#!/bin/bash
set -e
# use Schrodinger's utility volume_calc.py 
structconvert -imol2 %(file_l)s -omae lig.mae 

$SCHRODINGER/run volume_calc.py -imae lig.mae > cv.out"""% locals()
                file.write(script)

    def extract_rescoring_results(self, file_s):

        with open(file_s, 'a') as sf:
            if os.path.isfile('cv.out'):
                with open('cv.out', 'r') as ff:
                    if self.options['type'] == 'distance':
                        try:
                            print >> sf, ff.next().strip()
                        except StopIteration:
                            print >> sf, 'NaN'
                    elif self.options['type'] == 'volume':
                        try:
                            ff.next()
                            ff.next()
                            print >> sf, ff.next().split(',')[1].replace('\n','') 
                        except StopIteration:
                            print >> sf, 'NaN'
            else:
                print >> sf, 'NaN'
