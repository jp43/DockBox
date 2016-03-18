import os
import sys
import glob
import shutil
import subprocess

import util.check_license as chkl

required_programs = ['moebatch']

default_settings = {'placement': 'Alpha PMI', 'placement_nsample': '10', 'placement_maxpose': '250', 
'scoring': 'London dG', 'maxpose': '20', 'gtest': '0.01', 'rescoring': 'London dG', 'remaxpose': '20', 'binding_radius': '10000'}

known_settings = {'herg': {'binding_site': '3.966, 8.683, 11.093'}, \
'herg-cut': {'binding_site': '3.966, 8.683, 11.093'}, \
'herg-inactivated': {'binding_site': '0.000, 0.000, -5.000'}}

required_settings_names = ['binding_site']

def write_docking_script(filename, input_file_r, input_file_l, config):

    write_moe_docking_script('moe_dock.svl', config)

    convertsdf_cmd = chkl.eval("moebatch -exec \"mdb_key = db_Open ['lig.mdb','create']; db_Close mdb_key;\
        db_ImportSD ['lig.mdb','%(input_file_l)s','mol']\""%locals(), 'moe') # create mdb for ligand 

    dock_cmd = chkl.eval("moebatch -run moe_dock.svl -rec %(input_file_r)s -lig lig.mdb"%locals(), 'moe') # cmd for docking

    # write autodock script
    with open(filename, 'w') as file:
        script ="""#!/bin/bash
set -e
# convert sdf file to mdb
%(convertsdf_cmd)s

# run docking
%(dock_cmd)s
"""% locals()
        file.write(script)

def write_moe_docking_script(filename, config):

    binding_site_moe = '[' + config.options['moe']['binding_site'] + ']'

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
    local binding_site = %(binding_site_moe)s; // center for the binding site
    local binding_radius = %(binding_radius)s; // radius of the binding site
    local binding_res = []; // residues involved in binding site

    local idx;
    local com, dist;

    local rrec;
    rrec = cat cResidues chains; // extract residues info
    for idx = 1, length rrec loop
        com = oCenterOfMass rrec(idx);
        dist = sqrt add pow[sub[binding_site, com],2];
        if dist < binding_radius then
            binding_res = append [binding_res, rrec(idx)];
        endif
    endloop
    rrec = cat rAtoms binding_res;

    View (Atoms[]);

    local alpha_sites =  run ['sitefind.svl', [rrec, []], 'AlphaSites'];

    // Take first/highest scoring pocket alpha_sites(1)
    // Take fpos data alpha_sites(1)(1)
    // Take only coords of fpos data alpha_sites(1)(1)(2)
    local a_sites = apt cat alpha_sites(1)(1)(2); // x,y,z coords


    // Make dummy He atoms for alpha site
    local dummy, x,y,z;
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

endfunction;"""%dict(dict(locals()).items()+config.options['moe'].items())
        file.write(script)

def extract_docking_results(file_r, file_l, file_s, config):

    # copy receptor structure
    shutil.copyfile(config.input_file_r, file_r)

    sdffile = os.path.splitext(file_l)[0] + '.sdf'
    subprocess.check_call(chkl.eval("moebatch -exec \"db_ExportSD ['dock.mdb', '%s', ['mol','S'], []]\""%sdffile, 'moe'), shell=True)

    # save ligand structures in PDB file
    subprocess.check_call("babel %s %s"%(sdffile, file_l), shell=True)

    # save scores
    with open(sdffile, 'r') as sdff:
        with open(file_s, 'w') as sf:
            for line in sdff:
                if line.startswith("> <S>"):
                    print  >> sf, sdff.next().strip()

def cleanup(config):
    pass
