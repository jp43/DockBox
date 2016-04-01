import sys
import os
import shutil
import ConfigParser
import time
import glob
import shutil
import fileinput
import subprocess
import licence.check as chkl

ligprep_default_options = {'tautomerizer': False, 'ring_conf': False, 'stereoizer': False, 'tarfiles': False, 'ionization': '1'}

# the first element of each tuple corresponds to the value True or False when the flag applies
ligprep_bool_flags = {'tautomerizer': (False, '-nt'), 'ring_conf': (False, '-nr'), 'stereoizer': (False, '-ns'), 'tarfiles': (False, '-nz')}
ligprep_value_flags = {'ionization': '-i'}

known_formats = ['.sdf', '.smi', '.txt']
known_systems = ['herg', 'herg-cut', 'herg-inactivated']

def prepare_ligand(file_l, config):

    # construct the command to execute ligprep
    ligprep_flags_str = ""
    for key, value in config.ligprep.items():
        if key in ligprep_bool_flags:
            if value == ligprep_bool_flags[key][0]:
                ligprep_flags_str += ligprep_bool_flags[key][1] + " "
        elif key in ligprep_value_flags:
            ligprep_flags_str += ligprep_value_flags[key] + " " + str(value) + " "
        else:
            raise ValueError("Option %s does not look like a ligprep option!"%key)

    sdffile = os.path.basename(file_l)
    shutil.copyfile(file_l, sdffile)

    outputfile = generate_3D_structure(sdffile, ligprep_flags_str, config)

    # generate SDF files
    suffix = os.path.splitext(outputfile)[0]
    outputsdffile = suffix + '_.sdf'
    subprocess.check_call('babel %s %s -m 2>/dev/null'%(outputfile,outputsdffile), shell=True)

    # clean up
    os.remove(outputfile)

def generate_3D_structure(ff, flags, config):

    ext = os.path.splitext(ff)[1]
    if ext == '.sdf':
        inputflag = '-isd'
    elif ext in ['.smi', '.txt']:
        inputflag = '-ismi'

    suffix = (os.path.splitext(ff)[0]).split('/')[-1]
    outputfile = suffix + ".prep.sdf"

    ligprep_cmd = chkl.eval("ligprep -WAIT -W e,-ph,7.0,-pht,2.0 -epik -r 1 -bff 14 -ac %(inputflag)s %(ff)s -osd %(outputfile)s %(flags)s"%locals(), 'glide')

    with open('ligprep.sh', 'w') as file:
        script ="""#!/bin/bash
%(ligprep_cmd)s"""% locals()
        file.write(script)

    subprocess.check_output('bash ligprep.sh', shell=True)

    for logf in glob.glob(suffix +'*.log'):
       os.remove(logf)

    return outputfile

def prepare_receptor(file_r, config):

    # find new file name
    new_file_r = os.path.basename(file_r)
    pref, ext = os.path.splitext(new_file_r)
    new_file_r = pref + '.prep' + ext

    prepwizard_cmd = chkl.eval("prepwizard -WAIT -fix %(file_r)s %(new_file_r)s"%locals(), 'glide')

    with open('prepwizard.sh', 'w') as file:
        script ="""#!/bin/bash
%(prepwizard_cmd)s"""% locals()
        file.write(script)

    subprocess.check_output('bash prepwizard.sh', shell=True)
    return new_file_r
