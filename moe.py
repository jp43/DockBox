import os
import sys
import glob
import shutil
import subprocess

import licence.check as chkl

required_programs = ['moebatch']

default_settings = {'placement': 'Alpha PMI', 'placement_nsample': '10', 'placement_maxpose': '250', 
'scoring': 'London dG', 'maxpose': '10', 'gtest': '0.01', 'rescoring': 'London dG', 'binding_radius': '10000'}

default_sitefind_settings = {'minplb': '1.0'}

def set_site_options(site, options):

    # set box center
    center = site[1]
    options['center_bs'] = '[' + ', '.join(map(str.strip, center.split(','))) + ']'

    # set box size
    boxsize = site[2]
    boxsize = map(float, map(str.strip, boxsize.split(',')))
    options['radius_bs'] = max(boxsize)*1./2

def write_sitefinder_script(filename, file_r, config):
    
    write_moe_sitefinder_script('sitefinder.svl', file_r, config)

    sitefinder_cmd = chkl.eval("moebatch -run sitefinder.svl"%locals(), 'moe') # cmd for docking

    # write script
    with open(filename, 'w') as file:
        script ="""#!/bin/bash

# run docking
%(sitefinder_cmd)s
"""% locals()
        file.write(script)

def write_moe_sitefinder_script(filename, file_r, config):

    # write vina script
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

    write ['#ID    PLB     x     y    z   radius\\n'];

    for idx = 1, length alpha_sites loop
        plb = alpha_sites(idx)(4)(2);

        if plb > minplb or idx == 1 then
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
endfunction;""" %dict(dict(locals()).items()+config.site.items())
        file.write(script)


def write_docking_script(filename, input_file_r, input_file_l, site, options):

    # set box center
    center = site[1]
    center_bs = '[' + ', '.join(map(str.strip, center.split(','))) + ']'

    # set box size
    boxsize = site[2]
    boxsize = map(float, map(str.strip, boxsize.split(',')))
    radius_bs = max(boxsize)*1./2

    write_moe_docking_script('moe_dock.svl', options, center_bs, radius_bs)

    convertsdf_cmd = chkl.eval("moebatch -exec \"mdb_key = db_Open ['lig.mdb','create']; db_Close mdb_key;\
        db_ImportSD ['lig.mdb','%(input_file_l)s','mol']\""%locals(), 'moe') # create mdb for ligand 

    dock_cmd = chkl.eval("moebatch -run moe_dock.svl -rec %(input_file_r)s -lig lig.mdb"%locals(), 'moe') # cmd for docking

    # write script
    with open(filename, 'w') as file:
        script ="""#!/bin/bash
set -e
# convert sdf file to mdb
%(convertsdf_cmd)s

# run docking
%(dock_cmd)s
"""% locals()
        file.write(script)

def write_moe_docking_script(filename, options, center_bs, radius_bs):

    locals().update(options)

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
    local radius_bs = %(radius_bs)s; // radius of the binding site
    local residues_bs = []; // residues involved in binding site

    local idx;
    local com, dist;

    local rec_bs = cat cResidues chains; // extract residues info
    for idx = 1, length rec_bs loop
        com = oCenterOfMass rec_bs(idx);
        dist = sqrt add pow[sub[center_bs, com], 2];
        if dist < radius_bs then
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
                remaxpose: %(maxpose)s,
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

def extract_docking_results(file_r, file_l, file_s, input_file_r, extract):

    # copy receptor structure
    shutil.copyfile(input_file_r, file_r)

    sdffile = os.path.splitext(file_l)[0] + '.sdf'
    subprocess.check_call(chkl.eval("moebatch -exec \"db_ExportSD ['dock.mdb', '%s', ['mol','S'], []]\""%sdffile, 'moe'), shell=True)

    # save ligand structures in PDB file
    subprocess.check_call("babel %s %s &>/dev/null"%(sdffile, file_l), shell=True)

    # save scores
    with open(sdffile, 'r') as sdff:
        with open(file_s, 'w') as sf:
            for line in sdff:
                if line.startswith("> <S>"):
                    print  >> sf, sdff.next().strip()

def cleanup(config):
    pass
