from PyBabel.addh import *

class ArrangeHydrogens(AddHydrogens):

    def addHydrogens(self, atoms, atomsNoH):
        """ """
        idxs_ha = [-1 for atom in atoms]
        for idx, a in enumerate(atoms):
            for jdx, anoh in enumerate(atomsNoH):
                if a.coords == anoh.coords:
                    idxs_ha[idx] = jdx
                    break
        Hat = []
        nonHat = []
        for idx, a in enumerate(atoms):
            nbonds = len(a.bonds)
            if a.element != 'H':
                nHbonds = 0
                for bond in a.bonds:
                    if bond.atom1.element == 'H' or bond.atom2.element == 'H':
                        nHbonds += 1
                anoh = atomsNoH[idxs_ha[idx]]
                if nHbonds > 0:
                    nonHat.extend([idx]*nHbonds)
                nHat_p = len(Hat)
                if a.element == 'C': # treat the atom as C3
                    if nbonds == 4:
                        if nHbonds == 1:
                            Hat = Hat + self.add_tertiary_hydrogen(anoh, SP3_C_H_DIST)
                        elif nHbonds == 2:
                            Hat = Hat + self.add_methylene_hydrogens(anoh ,SP3_C_H_DIST)
                        elif nHbonds == 3:
                            Hat = Hat + self.add_methyl_hydrogen(anoh, SP3_C_H_DIST)
                            Hat = Hat + self.add_methylene_hydrogens(anoh, SP3_C_H_DIST,
                                                             Hat[-1])
                    elif nbonds == 3: # treat the atom as C2
                        if nHbonds == 1:
                            Hat = Hat + self.add_sp2_hydrogen(anoh ,SP2_C_H_DIST)
                        elif nHbonds == 2:
                            Hat = Hat + self.add_vinyl_hydrogens(anoh ,SP2_C_H_DIST)
                    elif nbonds == 2: # treat the atom as C1
                        if nHbonds == 1:
                            Hat = Hat + self.add_sp_hydrogen(anoh ,SP_C_H_DIST)
                elif a.element == 'N' and nbonds == 4: # treat the atom as N3+
                    if nHbonds == 1:
                        Hat = Hat + self.add_tertiary_hydrogen(anoh, SP3_N_H_DIST)
                    elif nHbonds == 2:
                        Hat = Hat + self.add_methylene_hydrogens(anoh ,SP3_N_H_DIST)
                    elif nHbonds == 3:
                        Hat = Hat + self.add_methyl_hydrogen(anoh, SP3_N_H_DIST)
                        Hat = Hat + self.add_methylene_hydrogens(anoh ,SP3_N_H_DIST,
                                                             Hat[-1])
                elif a.element == 'O' and nbonds == 2: # treat the atom as O3
                    if nHbonds == 1:
                        Hat = Hat + self.add_methyl_hydrogen(anoh, SP3_O_H_DIST)
                elif nbonds == 3:
                    if nHbonds == 1:
                        Hat = Hat + self.add_sp2_hydrogen(anoh ,SP2_N_H_DIST)
                    if nHbonds == 2:
                        Hat = Hat + self.add_vinyl_hydrogens(anoh ,SP2_N_H_DIST)
                elif nbonds == 2:
                    if nHbonds == 1:
                        Hat = Hat + self.add_sp_hydrogen(anoh ,SP2_N_H_DIST)
                #print a.element, a.babel_type, nbonds, nHbonds, len(Hat)-nHat_p

        if len(Hat) != len(nonHat):
            sys.exit("Number of hydrogens arranged do not match the initial number in the input file!")

        return Hat, nonHat
