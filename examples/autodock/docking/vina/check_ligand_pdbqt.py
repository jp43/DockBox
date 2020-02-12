import os
import sys
import shutil

input_file = sys.argv[1]

filename, ext = os.path.splitext(input_file)
file_tmp = filename + '_tmp.pdbqt'

lines_to_be_removed = []

has_branch_started = False
with open(input_file, 'r') as ff:
    for line in ff:
        if has_branch_started:
            has_branch_started = False
            branch_num = start_branch_line.split()[-1]
            if line.split()[1] != branch_num:
                lines_to_be_removed.append(start_branch_line)
                lines_to_be_removed.append('END' + start_branch_line)
        if line.startswith('BRANCH'):
            start_branch_line = line
            has_branch_started = True

if lines_to_be_removed:
    with open(input_file, 'r') as ff:
        with open(file_tmp, 'w') as of:
            for line in ff:
                if line.startswith(('BRANCH', 'ENDBRANCH')) and line in lines_to_be_removed:
                    pass
                else:
                    of.write(line)
    shutil.move(file_tmp, input_file)