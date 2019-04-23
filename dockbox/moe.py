import os
import sys
import glob
import shutil
import subprocess
import method
import license

from mdkit.utility import reader
from mdkit.utility import mol2

required_programs = ['moebatch']

default_settings = {'placement': 'Triangle Matcher', 'placement_nsample': '10', 'placement_maxpose': '250',  'scoring': 'London dG',
'maxpose': '30', 'remaxpose': '5', 'gtest': '0.01', 'rescoring': 'GBVI/WSA dG'}

known_scorings = ['ASE', 'Affinity dG', 'Alpha HB', 'GBVI/WSA dG', 'London dG', 'None']
known_placements = ['Alpha PMI', 'Alpha Triangle', 'Proxy Triangle', 'Triangle Matcher'] 

known_settings = {'placement': known_placements, 'scoring': known_scorings, 'rescoring': known_scorings}

class Moe(method.DockingMethod):

    def __init__(self, instance, site, options):

        super(Moe, self).__init__(instance, site, options)

        # set box center
        self.options['center_bs'] = '[' + ', '.join(map(str.strip, site[1].split(','))) + ']'

        # set box size
        self.options['boxsize_bs'] = '[' + ', '.join(map(str.strip, site[2].split(','))) + ']'

    def write_docking_script(self, filename, file_r, file_l):
   
        self.write_moe_docking_script('moe_dock.svl')
    
        convertmol2_cmd = license.wrap_command("moebatch -exec \"mdb_key = db_Open ['lig.mdb','create']; db_Close mdb_key;\
db_ImportMOL2 ['%(file_l)s','lig.mdb', 'molecule']\""%locals(), 'moe') # create mdb for ligand
    
        dock_cmd = license.wrap_command("moebatch -run moe_dock.svl -rec %(file_r)s -lig lig.mdb"%locals(), 'moe') # cmd for docking

        # write script
        with open(filename, 'w') as file:
            script ="""#!/bin/bash
# convert .mol2 file to mdb
%(convertmol2_cmd)s

# run docking
%(dock_cmd)s
"""% locals()
            file.write(script)
    
    def write_moe_docking_script(self, filename):
    
        locals().update(self.options)
    
        # write vina script
        with open(filename, 'w') as file:
            script ="""#svl
function DockAtoms, DockFile;
function DockMDBwAtoms, DockMDBwFile;

global argv;
function ArgvPull;

local function main []

    // Set potential and setup parameters
    pot_Load '$MOE/lib/Amber10EHT.ff';

    pot_Setup [
        strEnable: 1,
        angEnable: 1,
        stbEnable: 1,
        oopEnable: 1,
        torEnable: 1,
        vdwEnable: 1,
        eleEnable: 1,
        solEnable: 0,
        resEnable: 1,
        strWeight: 1,
        angWeight: 1,
        stbWeight: 1,
        oopWeight: 1,
        torWeight: 1,
        vdwWeight: 1,
        eleWeight: 1,
        solWeight: 1,
        resWeight: 1,
        cutoffEnable: 1,
        cutoffOn: 8,
        cutoffOff: 10,
        eleDist: 2,
        vdwScale14: 0.5,
        vdwBuffer1: 0,
        vdwBuffer2: 0,
        eleScale14: 0.833333,
        eleDielectric: 1,
        eleBuffer: 0,
        solDielectric: 80,
        solDielectricOffset: 0,
        state0: 1,
        state1: 0,
        state2: 1,
        threadCount: 0
    ];

ArgvReset ArgvExpand argv;
    local [recmdb, ligmdb, ph4file, outf] = ArgvPull [
        ['-rec','-lig','-ph4','-o'],
        1
    ];

    // If no receptor given as argument use default rec.moe
    if isnull recmdb then
        recmdb = 'rec.moe';
    endif

    local basename = fbase recmdb;
    local extension = fext recmdb;

    // output docking database file
    outf = 'dock.mdb';

    // Receptor file or database
    // Assume that the file is a moe or pdb file extract chains atoms

    local chains = ReadAuto [recmdb, []];
    local rec = cat cAtoms chains; // extract atom info from atom

    // get residues involved in the binding site
    local center_bs = %(center_bs)s; // center for the binding site
    local boxsize_bs = %(boxsize_bs)s; // size of the box for the binding site
    local residues_bs = []; // residues involved in binding site

    local idx, jdx;
    local com, dist;
    local isinbox;

    local rec_bs = cat cResidues chains; // extract residues info
    for idx = 1, length rec_bs loop
        com = oCenterOfMass rec_bs(idx);
        dist = sqrt add pow[sub[center_bs, com], 2];
        isinbox = 1;
        for jdx = 1, 3 loop
            if abs(center_bs(jdx)-com(jdx)) > 0.5*boxsize_bs(jdx) then
                isinbox = 0;
            endif
        endloop
        if isinbox == 1 then
            residues_bs = append [residues_bs, rec_bs(idx)];
        endif
    endloop

    rec_bs = cat rAtoms residues_bs;
    View (Atoms[]);

    local alpha_sites = run['sitefind.svl', [rec_bs, []], 'AlphaSites'];

    // Take first/highest scoring pocket alpha_sites(1)
    // Take fpos data alpha_sites(1)(1)
    // Take only coords of fpos data alpha_sites(1)(1)(2)
    local a_sites = apt cat alpha_sites(1)(1)(2); // x, y, z coords

    // Make dummy He atoms for alpha site
    // local dummy, x, y, z;
    // for x = 1, length a_sites(1) loop
    //    dummy(x) = sm_Build ['[He]'];
    //    aSetPos [dummy(x), [a_sites(1)(x), a_sites(2)(x), a_sites(3)(x)]];
    //endloop

    // Make dummy He atoms for alpha site
    local dummy, x, y, z;
    for x = 1, length a_sites loop
        dummy(x) = sm_Build ['[He]'];
        aSetPos [dummy(x), a_sites(x)];
    endloop

    // Make a collection of site atoms to send to docking
    // from the alpha site
    oSetCollection ['Site', dummy];
    local site = oGetCollection 'Site';

    // Ligand database
    local lmdb = _db_Open [ligmdb, 'read'];
    if lmdb == 0 then
        exit twrite ['Cannot read ligand mdb file {}', ligmdb];
    endif

    local ent = 0; // must have this set to zero
    while ent = db_NextEntry[lmdb, ent] loop; //loop through ligand database
        local ligdata = db_Read[lmdb, ent]; //read data for each entry
        local ligmoldata = ligdata.mol; // extract into moldata
        local ligchains = mol_Create ligmoldata; //create molecule in window
        local lig = cat cAtoms ligchains; // extract atom info from atom
    endloop

    // Set options for docking and refinement
    // maxpose is set to accept 50 poses, change as required
    local opt = [
                outrmsd: 1,
                sel_ent_only_rec: 0,
                sel_ent_only: 0,
                wall: [ '', 0, [ 0, 0, 0 ], [ 1000000, 1000000, 1000000 ], 0 ],
                csearch: 1,
                placement: '%(placement)s',
                placement_opt: [nsample : %(placement_nsample)s, maxpose : %(placement_maxpose)s ],
                scoring: '%(scoring)s',
                scoring_opt: [ train : 0 ],
                dup_placement: 1,
                maxpose: %(maxpose)s,
                refine: 'Forcefield',
                refine_opt: [ cutoff : 6, wholeres : 1, mmgbvi : 1, fixrec : 'Fix', tether : 10, gtest : %(gtest)s,
                maxit : 500, OverrideSetup : 1, k_potl : 100, offset : 0.4 ],
                rescoring: '%(rescoring)s',
                rescoring_opt: [ train : 0 ],
                dup_refine: 1,
                remaxpose: %(remaxpose)s,
                descexpr: '',
                receptor_mfield: '',
                ligand_mfield: '',
                tplate: [  ],
                tplateSel: [  ],
                //ph4: ph4file,
                ligmdbname: ligmdb,
                recmdbname: recmdb
    ];

    //Perform the docking
    DockFile [rec, site, ligmdb, outf, opt];

    oDestroy ligchains;
    db_Close lmdb;
    write ['Docking finished at {}.\\n', asctime []];

endfunction;"""% locals()
            file.write(script)
    
    def extract_docking_results(self, file_s, input_file_r, input_file_l):

        subprocess.check_output(license.wrap_command("moebatch -exec \"db_ExportTriposMOL2 ['dock.mdb', 'lig.mol2', 'mol', []]\"", 'moe'), shell=True, executable='/bin/bash')

        if os.path.exists('lig.mol2'):
            ligname = reader.open(input_file_l).ligname
            mol2.update_mol2file('lig.mol2', 'lig-.mol2', ligname=ligname, multi=True)
            os.remove('lig.mol2')

            # get SDF to extract scores
            sdffile = 'lig.sdf'
            subprocess.check_output(license.wrap_command("moebatch -exec \"db_ExportSD ['dock.mdb', '%s', ['mol','S'], []]\""%sdffile, 'moe'), shell=True, executable='/bin/bash')
            with open(sdffile, 'r') as sdff:
                with open(file_s, 'w') as sf:
                    for line in sdff:
                        if line.startswith("> <S>"):
                            print  >> sf, sdff.next().strip()
            os.remove(sdffile)
    
    def write_rescoring_script(self, filename, file_r, file_l):

        locals().update(self.options)

        if self.options['rescoring'] == 'prolig':
            rescoring_cmd = license.wrap_command("moebatch -run moe_rescoring.svl -rec %(file_r)s -lig %(file_l)s"%locals(), 'moe') # cmd for docking

            with open(filename, 'w') as file:
                script ="""#!/bin/bash
echo "#svl
function prolig_Calculate;

global argv;
function ArgvPull;

local function main[]

    ArgvReset ArgvExpand argv;
    local [recmdb, ligmdb, outf] = ArgvPull [
        ['-rec','-lig','-o'],
        1
    ];
    local lk = ReadTriposMOL2 [ligmdb, []];

    // Load pdb
    local rk = ReadAuto [recmdb, []];

    local itypes = ['hbond', 'metal', 'ionic', 'covalent', 'arene', 'distance'];
    local iract = prolig_Calculate [itypes, lk, rk, []];
    //local iract_v = Formulate2DInteractions [lk, rk, []];

    local idx;
    local interaction_energy = 0.;
    for idx = 1, length iract(1) loop
        if iract(1)(idx) == 'distance' then
            break;
        else
            interaction_energy = interaction_energy + iract(4)(idx);
        endif
    endloop

    write ['Interaction energy: {f.2} kCal/mol \\n', interaction_energy];

endfunction;" > moe_rescoring.svl

%(rescoring_cmd)s
""" %locals()
                file.write(script)

        else:
            convertmol2_cmd = license.wrap_command("moebatch -exec \"mdb_key = db_Open ['lig.mdb','create']; db_Close mdb_key;\
db_ImportMOL2 ['%(file_l)s','lig.mdb', 'molecule']\""%locals(), 'moe') # create mdb for ligand
            rescoring_cmd = license.wrap_command("moebatch -run moe_rescoring.svl -rec %(file_r)s -lig lig.mdb"%locals(), 'moe') # cmd for docking

            # write vina script
            with open(filename, 'w') as file:
                script ="""#!/bin/bash

%(convertmol2_cmd)s

echo "#svl
function DockAtoms, DockFile;
function DockMDBwAtoms, DockMDBwFile;

global argv;
function ArgvPull;

local function main []

    // Set potential and setup parameters
    pot_Load '$MOE/lib/Amber10EHT.ff';

    pot_Setup [
        strEnable: 1,
        angEnable: 1,
        stbEnable: 1,
        oopEnable: 1,
        torEnable: 1,
        vdwEnable: 1,
        eleEnable: 1,
        solEnable: 0,
        resEnable: 1,
        strWeight: 1,
        angWeight: 1,
        stbWeight: 1,
        oopWeight: 1,
        torWeight: 1,
        vdwWeight: 1,
        eleWeight: 1,
        solWeight: 1,
        resWeight: 1,
        cutoffEnable: 1,
        cutoffOn: 8,
        cutoffOff: 10,
        eleDist: 2,
        vdwScale14: 0.5,
        vdwBuffer1: 0,
        vdwBuffer2: 0,
        eleScale14: 0.833333,
        eleDielectric: 1,
        eleBuffer: 0,
        solDielectric: 80,
        solDielectricOffset: 0,
        state0: 1,
        state1: 0,
        state2: 1,
        threadCount: 0
    ];

ArgvReset ArgvExpand argv;
    local [recmdb, ligmdb, ph4file, outf] = ArgvPull [
        ['-rec','-lig','-ph4','-o'],
        1
    ];

    // If no receptor given as argument use default rec.moe
    if isnull recmdb then
        recmdb = 'rec.moe';
    endif

    local basename = fbase recmdb;
    local extension = fext recmdb;

    // output docking database file
    outf = 'dock.mdb';

    // Receptor file or database
    // Assume that the file is a moe or pdb file extract chains atoms

    local chains = ReadAuto [recmdb, []];
    local rec = cat cAtoms chains; // extract atom info from atom

    local alpha_sites = run['sitefind.svl', [rec, []], 'AlphaSites'];

    // Take first/highest scoring pocket alpha_sites(1)
    // Take fpos data alpha_sites(1)(1)
    // Take only coords of fpos data alpha_sites(1)(1)(2)
    local a_sites = apt cat alpha_sites(1)(1)(2); // x, y, z coords

    // Make dummy He atoms for alpha site
    local dummy, x, y, z;
    for x = 1, length a_sites loop
        dummy(x) = sm_Build ['[He]'];
        aSetPos [dummy(x), a_sites(x)];
    endloop

    // Make a collection of site atoms to send to docking
    // from the alpha site
    oSetCollection ['Site', dummy];
    local site = oGetCollection 'Site';

    // Ligand database
    local lmdb = _db_Open [ligmdb, 'read'];
    if lmdb == 0 then
        exit twrite ['Cannot read ligand mdb file {}', ligmdb];
    endif

    local ent = 0; // must have this set to zero
    while ent = db_NextEntry[lmdb, ent] loop; //loop through ligand database
        local ligdata = db_Read[lmdb, ent]; //read data for each entry
        local ligmoldata = ligdata.mol; // extract into moldata
        local ligchains = mol_Create ligmoldata; //create molecule in window
        local lig = cat cAtoms ligchains; // extract atom info from atom
    endloop

    // Set options for docking and refinement
    // maxpose is set to accept 50 poses, change as required
    local opt = [
                outrmsd: 1,
                sel_ent_only_rec: 0,
                sel_ent_only: 0,
                wall: [ '', 0, [ 0, 0, 0 ], [ 1000000, 1000000, 1000000 ], 0 ],
                csearch: 1,
                placement: 'None',
                scoring: 'None',
                dup_placement: 1,
                rescoring: '%(rescoring)s',
                rescoring_opt: [ train : 0 ],
                dup_refine: 1,
                remaxpose: 1,
                descexpr: '',
                receptor_mfield: '',
                ligand_mfield: '',
                tplate: [  ],
                tplateSel: [  ],
                ligmdbname: ligmdb,
                recmdbname: recmdb
    ];

    //Perform the docking
    DockFile [rec, site, ligmdb, outf, opt];

    oDestroy ligchains;
    db_Close lmdb;
    write ['Docking finished at {}.\\n', asctime []];

endfunction;" > moe_rescoring.svl

%(rescoring_cmd)s"""% locals()
                file.write(script)

    def extract_rescoring_results(self, file_s):

        locals().update(self.options)

        if self.options['rescoring'] == 'prolig': 
            with open(file_s, 'a') as sf:
                if os.path.exists('moebatch.log'):
                    with open('moebatch.log', 'r') as logf:
                        is_interaction_energy = False
                        for line in logf:
                            if line.startswith("Interaction energy:"):
                                print >> sf, line.split()[-2]
                                is_interaction_energy = True
                                break
                        if not is_interaction_energy:
                            print >> sf, 'NaN'
                else:
                     print >> sf, 'NaN'
        else:
            # get SDF to extract scores
            sdffile = 'lig.sdf'
            subprocess.check_output(license.wrap_command("moebatch -exec \"db_ExportSD ['dock.mdb', '%s', ['mol','S'], []]\""%sdffile, 'moe'), shell=True, executable='/bin/bash')
            with open(file_s, 'a') as sf:
                if os.path.exists(sdffile):
                    with open(sdffile, 'r') as sdff:
                        for line in sdff:
                            if line.startswith("> <S>"):
                                print >> sf, sdff.next().strip()
                                break
                        os.remove(sdffile)
                else:
                    print >> sf, 'NaN'

def write_sitefinder_script(filename, file_r, args):
    
    write_moe_sitefinder_script('sitefinder.svl', file_r, args)
    sitefinder_cmd = license.wrap_command("moebatch -run sitefinder.svl", 'moe') # cmd for docking

    # write script
    with open(filename, 'w') as file:
        script ="""#!/bin/bash
# run docking
%(sitefinder_cmd)s
"""% locals()
        file.write(script)

def write_moe_sitefinder_script(filename, file_r, args):

    if args.nsitesmax == 0:
        nsitesmax = 'length alpha_sites'
    else:
        nsitesmax = str(args.nsitesmax)
    minplb = args.minplb

    # write svl script
    with open(filename, 'w') as file:
        script ="""#svl

local function main []
    local chains = ReadAuto ['%(file_r)s', []];
    local rec = cat cAtoms chains; // extract atom info from atom

    // locate alpha sites
    local alpha_sites = run['sitefind.svl', [rec, []], 'AlphaSites'];

    local dummy, x, dist;
    local a_sites, plb;
    local minplb = %(minplb)s, maxdist;
    local idx;
    local nsites;
    local cog; // center of geometry

    write ['#ID PLB  x  y  z  radius\\n'];

    for idx = 1, length alpha_sites loop
        plb = alpha_sites(idx)(4)(2);

        if (plb > minplb or idx == 1) and idx <= %(nsitesmax)s then
            a_sites = alpha_sites(idx)(1)(2);
            nsites = length a_sites(1);

            // get center of geometry of the alpha sites
            cog = [0.0, 0.0, 0.0];
            for x = 1, nsites loop
                cog = add[[a_sites(1)(x), a_sites(2)(x), a_sites(3)(x)], cog];
            endloop
            cog = div[cog, nsites];
            maxdist = 0;

            // get distance to the farthest atom
            for x = 1, nsites loop
                dist = sqrt add pow[sub[[a_sites(1)(x), a_sites(2)(x), a_sites(3)(x)], cog], 2];
                if dist > maxdist then
                    maxdist = dist;
                endif
            endloop
            write ['{f.0} {f.2} {f.3} {f.3}\\n', idx, plb, cog, maxdist];
        endif
    endloop
endfunction;""" %locals()
        file.write(script)
