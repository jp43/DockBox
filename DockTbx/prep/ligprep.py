import sys
import os
import shutil
import stat
import glob
import subprocess

from DockTbx.licence import check as chkl

def prepare_ligand(file_l, flags):

    # copy ligand file in current directory
    input_file = os.path.basename(file_l)
    shutil.copyfile(file_l, input_file)
   
    # Generate 3D structure using ligprep
    output_file = generate_3D_structure(input_file, flags)

    return output_file

def generate_3D_structure(file_l, flags):

    ext = os.path.splitext(file_l)[1]
    if ext == '.sdf':
        input_format_flag = '-isd'
    elif ext in ['.smi', '.txt']:
        input_format_flag = '-ismi'
    else:
        raise IOError("Format %s not recognized!"%(ext[1:]))

    suffix = (os.path.splitext(file_l)[0]).split('/')[-1]
    output_file = suffix + ".prep.sdf"

    # write ligprep command
    cmd = chkl.eval("ligprep -WAIT %(flags)s %(input_format_flag)s %(file_l)s -osd %(output_file)s"%locals(), 'schrodinger')

    script_name = 'run_ligprep.sh'
    with open(script_name, 'w') as file:
        script ="""#!/bin/bash
set -e
%(cmd)s"""% locals()
        file.write(script)
    os.chmod(script_name, stat.S_IRUSR | stat.S_IWUSR | stat.S_IRGRP | stat.S_IROTH | stat.S_IXUSR)
    # execute ligprep
    subprocess.check_output('./' + script_name +" &> ligprep.log", shell=True, executable='/bin/bash')

    # cleanup directory
    for logf in glob.glob(suffix +'*.log'):
       os.remove(logf)

    return output_file

def prepare_receptor(file_r, flags):

    # find new file name
    new_file_r = os.path.basename(file_r)
    pref, ext = os.path.splitext(new_file_r)
    new_file_r = pref + '.prep' + ext

    # write ligprep command
    cmd = chkl.eval("prepwizard -WAIT %(flags)s %(file_r)s %(new_file_r)s"%locals(), 'schrodinger')
    script_name = 'run_prepwizard.sh'
    with open(script_name, 'w') as file:
        script ="""#!/bin/bash
set -e
%(cmd)s"""% locals()
        file.write(script)
    os.chmod(script_name, stat.S_IRUSR | stat.S_IWUSR | stat.S_IRGRP | stat.S_IROTH | stat.S_IXUSR)

    subprocess.check_output('./' + script_name + " &> recprep.log", shell=True, executable='/bin/bash')
    return os.path.abspath(new_file_r)
