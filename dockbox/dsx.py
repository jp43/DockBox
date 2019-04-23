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
rm -rf dsx.txt

cp %(file_r)s protein.pdb
cp %(file_l)s ligand.mol2

# execute DSX
dsx -P protein.pdb -L ligand.mol2 -F dsx.txt%(pot_dir_str)s%(other_flags_str)s
"""% locals()
            file.write(script)

    def extract_rescoring_results(self, file_s):

        with open(file_s, 'a') as sf:
            if os.path.isfile('dsx.txt'):
                with open('dsx.txt', 'r') as txtf:
                    for line in txtf:
                        if line.startswith(" 0"):
                            print >> sf, line.split('|')[3].strip()
                            break
            else:
                print >> sf, 'NaN'
