import sys
import shutil

def give_unique_atom_names(file_l):

    tmpfile = 'tmp.mol2'
    with open(file_l, 'r') as oldf:
        newf = open(tmpfile, 'w')
        known_atom_types = []
        atom_numbers = []
        for line in oldf:
            if line.startswith('@<TRIPOS>ATOM'):
                is_structure = True
                newf.write(line)
            elif line.startswith('@<TRIPOS>'):
                is_structure = False
                newf.write(line)
            elif is_structure:
                line_s = line.split()
                atom = line_s[1]
                atom_type = ''.join([ch for ch in atom if not ch.isdigit()])
                if atom_type not in known_atom_types:
                    known_atom_types.append(atom_type)
                    atom_number = '1'
                    atom_numbers.append(1)
                else:
                    idx = known_atom_types.index(atom_type)
                    atom_number = str(atom_numbers[idx]+1)
                    atom_numbers[idx] += 1
                new_atom_type = line_s[1]+atom_number
                newline = ' '*(7-len(line_s[0])) + line_s[0] + 2 * ' ' + new_atom_type + \
' '*(10-len(new_atom_type)) + line[19:]
                newf.write(newline)
            else:
                newf.write(line)

    shutil.move(tmpfile, file_l)
