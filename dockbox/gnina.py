import os
import sys
import method
import subprocess

from mdkit.utility import mol2

required_programs = ['gnina']
default_settings = {'type': 'CNNscore', 'cnn': None}

class Gnina(method.ScoringMethod):

    def write_rescoring_script(self, filename, file_r, file_l):

        if self.options['cnn'] is None or self.options['cnn'].lower() in ["none", "no"]:
            cnn_flag = ""
        else:
            cnn_flag = " --cnn %s"%self.options['cnn']

        # write vina script
        with open(filename, 'w') as file:
            script ="""#!/bin/bash

rm -rf gnina.out

# execute GNINA
gnina -r %(file_r)s -l %(file_l)s%(cnn_flag)s --score_only > gnina.out\n"""% locals()
            file.write(script)

    def extract_rescoring_results(self, filename):

        with open(filename, 'a') as sf:
            is_score = False
            if os.path.isfile('gnina.out'):
                with open('gnina.out', 'r') as outf:
                    for line in outf:
                        if line.startswith(self.options['type']):
                            sf.write(line.split()[1]+'\n')
                            is_score = True
                            break
                if not is_score:
                    sf.write('NaN\n')
            else:
                sf.write('NaN\n')
