import os
import sys
import method
from mdkit.utility import mol2

required_programs = ['dsx']

default_settings = {'pot_dir': None, 'other_flags': None}

class Dsx(method.ScoringMethod):

    def write_rescoring_script(self, filename, file_r, file_l):

        locals().update(self.options)

        if self.options['pot_dir']:
            pot_dir_str = ' -D ' + self.options['pot_dir']
        else:
            pot_dir_str = ''

        if self.options['other_flags']:
            other_flags_str = ' ' + self.options['other_flags']
        else:
            other_flags_str = ''

        # write vina script
        with open(filename, 'w') as file:
            script ="""#!/bin/bash
set -e
# remove pre-existing result file
rm -rf dsx.out

cp %(file_r)s protein.pdb
cp %(file_l)s ligand.mol2

# execute DSX
dsx -P protein.pdb -L ligand.mol2 -F dsx.out%(pot_dir_str)s%(other_flags_str)s
"""% locals()
            file.write(script)

    def extract_rescoring_results(self, filename):

        with open(filename, 'a') as sf:
            is_score = False
            if os.path.isfile('dsx.out'):
                with open('dsx.out', 'r') as outf:
                    for line in outf:
                        if line.startswith(" 0"):
                            sf.write(line.split('|')[3].strip()+'\n')
                            is_score = True
                            break
                if not is_score:
                    sf.write('NaN\n')
            else:
                sf.write('NaN\n')
