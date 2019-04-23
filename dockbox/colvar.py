import os
import sys
import method

mandatory_settings = ['type']
known_types = ['distance', 'volume', 'sasa']
default_settings = {'type': 'distance', 'residues': None, 'distance_mode': 'cog'}

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
            if self.options['distance_mode'] == 'cog':
                with open(filename, 'w') as file:
                    script ="""#!/bin/bash
set -e
echo "from mdkit.utility import mol2, PDB
from mdkit.utility import reader
from mdkit.utility import utils
import numpy as np

ligrd = reader.open('%(file_l)s')
coords_lig = [map(float,line[2:5]) for line in ligrd.next()['ATOM']]
coords_lig = np.array(coords_lig)
cog_lig = utils.center_of_geometry(coords_lig)

recrd = reader.open('%(file_r)s')
residue = %(residues)s

dist_min = 1e10
rec = recrd.next()['ATOM']
for rs in residue:

    coords_rec = [map(float,line[4:7]) for line in rec if line[3] == str(rs)]
    coords_rec = np.array(coords_rec)
    cog_rec = utils.center_of_geometry(coords_rec)

    dist = np.sqrt(np.sum((cog_lig - cog_rec)**2))
    if dist < dist_min:
        dist_min = dist

# write min distance
with open('cv.out', 'w') as ff:
    ff.write(str(dist_min))" > get_distance.py 
python get_distance.py"""% locals()
                    file.write(script)

            elif self.options['distance_mode'] == 'min':
                with open(filename, 'w') as file:
                    script ="""#!/bin/bash
set -e
echo "from mdkit.utility import mol2, PDB
from mdkit.utility import reader
from mdkit.utility import utils
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

        elif self.options['type'] == 'sasa':
            files_l_joined = ' '.join(file_l)
            with open(filename, 'w') as file:
                script ="""#!/bin/bash
set -e
cat %(files_l_joined)s > lig.mol2

structconvert -ipdb %(file_r)s -omae rec.mae
structconvert -imol2 lig.mol2 -omae lig.mae

cat rec.mae lig.mae > complex.mae

# use Schrodinger's utility bindind_sasa.py 
$SCHRODINGER/run binding_sasa.py complex.mae -f > cv.out
structconvert -n 2: -imae complex-out_pv.mae -osd lig_out.sdf"""% locals()
                file.write(script)

    def extract_rescoring_results(self, file_s, nligands=None):

        if self.options['type'] in ['distance', 'volume']:
            with open(file_s, 'a') as sf:
                try:
                    with open('cv.out', 'r') as ff:
                        if self.options['type'] == 'distance':
                            print >> sf, ff.next().strip()
                        elif self.options['type'] == 'volume':
                            ff.next()
                            ff.next()
                            print >> sf, ff.next().split(',')[1].replace('\n','') 
                except:
                    print >> sf, 'NaN'

        elif self.options['type'] == 'sasa':
            with open(file_s, 'w') as sf:
                try:
                    with open('lig_out.sdf', 'r') as ff:
                        is_sasa_line = False
                        for line in ff:
                            if line.startswith('> <r_user_sasa_ligand_total_delta>'):
                                is_sasa_line = True
                            elif is_sasa_line:
                                is_sasa_line = False
                                sf.write(line)
                except:
                    for idx in range(nligands):
                        print >> sf, 'NaN'
