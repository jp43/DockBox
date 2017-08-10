Mode-of-binding Predictor
=========================

Mode-of-binding Predictor (MOBPred) is a python package used to facilitate the use of popular docking softwares (including structure preparation, docking and rescoring). The package is particularly suitable to compare docking results obtained from different softwares or combine them in a consensus docking or consensus scoring strategy.

Note that the softwares used for structure preparation or docking are not part of the MOBPred package. Whatever software needs to be used, it should be installed separately on the same machine MOBPred is set up.

Below is a list of all the programs which can be used together with MOBPred. 

* Structure preparation/optimization:

  * ligprep (Schrodinger 2015, https://www.schrodinger.com/ligprep) used to prepare compounds (generate protonation states, isomers conformers...)

  * antechamber, prmchk, tleap, sander (AMBER12 or later) used to assign partial charges, minimize structures...

  * moebatch (MOE2015) used to identify probable binding sites

  * Open Babel (http://openbabel.org/wiki/Main_Page)

* Docking:

  * Autodock (http://autodock.scripps.edu)

  * Autodock Vina (http://autodock.scripps.edu)

  * Glide (https://www.schrodinger.com/glide)

  * MOE2015 (https://www.chemcomp.com/MOE-Molecular_Operating_Environment.htm)

  * DOCK 6 (http://dock.compbio.ucsf.edu/DOCK_6/index.htm)

Prerequisites
-------------

Before installing the MOBPred package, make sure that you have the following packages installed:

* NumPy; version 1.4.1 or later

* pandas; version 0.18.1 or later

* AmberTools; version 12 or later

Any software intended to be used in conjunction with MOBPred should be installed separetely and should work as a standalone program. In addition, make sure the applications mentioned below are in your PATH, depending on which docking softwares is used:

* Autodock: Babel and autodock4 set aside, all the executables below can be found in the AutoDockTools package (http://autodock.scripps.edu/downloads/resources/adt/index_html): **prepare_ligand4.py**, **prepare_receptor4.py**, **prepare_dpf4.py**, **prepare_gpf4.py**, **autogrid4**, **autodock4**, **babel**.

* Autodock Vina: **prepare_ligand4.py**, **prepare_receptor4.py**, **vina**, **babel**.

* Glide: All the executables can be found in the Schrodinger package: **prepwizard**, **glide**, **ligprep**, **glide_sort**, **pdbconvert**.

* MOE: **moebatch**.

* DOCK 6: **chimera**, **dms**, **sphgen_cpp**, **sphere_selector**, **showbox**, **grid**, **dock6**.

*Pharmamatrix users*: on the pharmamatrix cluster, the majority of the docking softwares mentioned above have been already installed. Here is an example on how the PATH environment variable can be updated to include Autodock, Vina, Glide, MOE, DOCK 6:

    # Amber 14
    export AMBERHOME=/pmshare/amber/ambertools14_ibm_gnu-20140820
    PATH=$AMBERHOME/bin:$PATH

    # AD and ADV
    PATH=/opt/mgltools/1.5.4:/opt/mgltools/1.5.4/MGLToolsPckgs/AutoDockTools/Utilities24:$PATH
    PATH=/pmshare/vina/autodock_vina_1_1_2_linux_x86/bin:$PATH
    PATH=/opt/autodock/4.2.3/bin:$PATH

    # Glide
    export SCHRODINGER=/nfs/r720-1/preto/tools/schrodinger2015-4
    PATH=/nfs/r720-1/preto/tools/schrodinger2015-4:$PATH
    PATH=/nfs/r720-1/preto/tools/schrodinger2015-4/utilities:$PATH

    # MOE 2015
    export MOE=/nfs/r720-1/preto/tools/moe2015
    PATH=/nfs/r720-1/preto/tools/moe2015/bin:$PATH

    # DOCK 6
    PATH=/nfs/r720-1/preto/tools/UCSF-Chimera64-1.10.2/bin:$PATH
    PATH=/nfs/r720-1/preto/tools/dms:$PATH
    PATH=/nfs/r720-1/preto/tools/src/dock6/bin:$PATH

    export PATH

Installation
------------

The Python Distutils are used to build and install MOBPred, so it is fairly simple to get things ready to go. Following are very simple instructions on how to proceed:

1. First, make sure that you have the NumPy and pandas modules. If not, get them from http://numpy.scipy.org/, http://pandas.pydata.org. Compile/install them.

2. Make sure AmberTools is installed and that standard executables (e.g., sander, tleap,...)  are accessible through your PATH variable. For pharmamatrix users, see section *prerequisites*.

3. From the main MOBPred distribution directory run this command (plus any extra flags, e.g., --prefix or --user to specify the installation directory):

        python setup.py install

After installation, make sure that the bin folder within your installation directory is included to your PATH. It contains the executables *prepvs*, *rundock* and *runanlz* that are used for virtual screening preparation, docking and docking analysis, respectively.


Commands
--------

After including to your PATH the bin folder from your installation, you can use the following commands from any location: 


1. prepvs

    prepvs is used to prepare the ligand structures (carried out with Schrodinger's ligprep utility) and to create folders (one per receptor and per ligand) intended to facilitate docking or virtual screening runs performed with rundock. When typing "prepvs -h" on the command line, the following help message will pop up:

        prepvs -h
        usage: prepvs [-h] [-l INPUT_FILES_L [INPUT_FILES_L ...]]
                      [-r INPUT_FILES_R [INPUT_FILES_R ...]] [-f CONFIG_FILE]
                      [-lpflags LPFLAGS] [-ligprep] [-site SITE] [-noprep]
        
        Prepare files for Docking or Virtual Screening
        
        optional arguments:
          -h, --help            show this help message and exit
          -l INPUT_FILES_L [INPUT_FILES_L ...]
                                Ligand coordinate file(s): .sdf, .smi
          -r INPUT_FILES_R [INPUT_FILES_R ...]
                                Receptor coordinate file(s): .pdb
          -f CONFIG_FILE        config file: .ini
          -lpflags LPFLAGS      Ligprep (Schrodinger) flags for ligand preparation.
                                Default: "-ph 7.0 -pht 2.0 -i 2 -s 8 -t 4"
          -ligprep              Prepare compounds using ligprep
          -site SITE            Update binding sites info in config file from file
          -noprep               No structure preparation, update directories and files
                                only (used debbuging)

2. rundock

    rundock is used to dock and eventually minimize and rescore the poses. When typing "rundock -h" on the command line, the following help message will pop up:

        usage: rundock [-h] -l INPUT_FILE_L -r INPUT_FILE_R -f CONFIG_FILE
                       [-q CHARGE_FILE] [-rescore_only] [-extract_only] [-d POSEDIR]
                       [-norun]
        
        rundock : dock with multiple softwares -------- Requires one file for the
        ligand (1 struct.) and one file for the receptor (1 struct.)
        
        optional arguments:
          -h, --help       show this help message and exit
          -l INPUT_FILE_L  Ligand coordinate file(s): .mol2
          -r INPUT_FILE_R  Receptor coordinate file(s): .pdb
          -f CONFIG_FILE   config file containing docking parameters
          -q CHARGE_FILE   File with partial charges of non-standard residues
          -rescore_only    Run rescoring only
          -extract_only    Extract structures only (usually used for debugging)
          -d POSEDIR       Directory containing poses to rescore (should be used with
                           rescore_only option)
          -norun           Do not run the scripts for docking (simply generate the
                           files)

    * Mandatory arguments

        Besides one .mol2 file containing the ligand structure (-l flag) and one .pdb file containing the receptor structure (-r flag), rundock requires another mandatory input file, namely, a configuration file (-f flag) where all the parameters needed for the docking procedure are specified.

        **N.B.**: *rundock* can only be used to run docking and scoring procedures with single protein and ligand structures. If multiple protein or/and ligand structures need to be used, the *prepvs* command can be used to create folders for each protein-ligand pair (see the above section *prepvs*). 

        **Configuration to dock with multiple softwares on a single binding site and eventually minimize the poses**

        Below is an example of configuration file that can be used as an input of *rundock*. The docking procedure is carried out on a single binding site specied as a box with dimensions 30.0 x 30.0 x 30.0 centered at the position (x, y, z) = 8.446, 25.365, 4.394.

            [DOCKING]
            program = autodock, vina, dock, glide, moe, moe1, moe2
            rescoring = no
            minimize = yes
            cleanup = no
            
            [AUTODOCK]
            ga_run = 50
            spacing = 0.3
            
            [VINA]
            num_modes = 20
            
            [DOCK]
            nposes = 200
            
            [GLIDE]
            poses_per_lig = 200
            pose_rmsd = 2.0
            precision = SP
            use_prepwizard = False
            
            [MOE]
            scoring = London dG
            maxpose = 100
            remaxpose = 50
            
            [MOE1]
            scoring = GBVI/WSA dG
            maxpose = 100
            remaxpose = 50
            
            [MOE2]
            scoring = Affinity dG
            maxpose = 100
            remaxpose = 50
            
            [SITE]
            center = 8.446, 25.365, 4.394
            boxsize = 30.0, 30.0, 30.0

        * The DOCKING section includes the softwares that should be used for docking, and if minimization, rescoring and/or cleanup should be performed. The docking softwares should be specified with coma separation through the keyword *programs*.

            Keywords

            * **programs**: specifies the softwares which are used for docking (autodock, vina, dock6, glide and/or moe). Options relative to each program (or instance) are specfied within the section of the same name. For example, if autodock is in the list of programs, options associated with autodock should be specified in the AUTODOCK section of the configuration file.

    * Optional arguments

    Preferably do not use flags other than -l, -r and -f

3. runanlz


    runanlz is used to anlyze the docking poses obtained with the rundock command

LigPrep
-------

Used to prepare the ligand structure

default flags: ligprep -WAIT -W e,-ph,7.0,-pht,2.0 -s 8 -t 4
These flags aim at generating a few low-risk variations on the input structures (p.40 of ligprep manual)

Steps:

    sdconvert
        -- Converts the input sdf or smi to the schrodinger format

    applyhtreat

        -- Adds (or deletes) hydrogen atoms following treatment
        -- Chemical structures often are specified with implicit hydrogens
        -- The default treatment should be fine "All-atom with No-Lp" (lone pair)
        -- Note that for AutoDock you need to remove non-polar hydrogens, but this will be taken care of later by the ligand preparation script for AutoDock
        -- Also if you are preparing the ligands for a particular force field you may want to select a different treatment, or again you can post-process it

    desalter
        -- Normally you should just leave this on
        -- This will remove the counter-ions that you sometimes find in chemical database structures
        -- Also rarely there might be multiple unbonded molecules stored as a single "structure", this will just pick the single largest molecule (for example, this happens in drugbank with some "drugs" that are mixtures)

    neutralizer
        -- The default is to neutralize, that is normally what you want
        -- It will do this by adding/removing protons
        -- Can check the manual for the exact list of changes that it may make

    ionizer
        -- This doesn't run by default
        -- For docking normally a neutral state only is what you want... at least that's what we've done in the past

    tautomerizer
        -- This will generate multiple isomers from the input structures by moving protons & double bonds
        -- The default is up to 8 tautomers
        -- The default is to exclude tautomers with probability < 0.01

    stereoizer
        -- This will generate multiple stereoisomers (e.g. at carbon stereo centers or double-bonds)
        -- It will keep the chirality from the input structures where it is specified, but where it is not specified it will generate most possible stereoisomers (up to the max stereoisomers allowed)
        -- There are some restrictions it will apply by default, i.e. it will exclude some states are not achievable for geometric reasons or are atypical for some types of natural products (e.g. peptides and steroids).
        -- The default is up to 32 stereoisomers

    ring_conf
        -- For non-flexible rings it will always use the input conformation
        -- By default this will only generate a single (most likely) ring conformation
        -- Might be worth trying to increase the max number of ring conformations, e.g. add the ligprep option "-r 3"

    premin & bmin
        -- Uses a forcefield to generate a 3D conformation
        -- One reasonable conformation should be fine, the docking program will explore other conformations
        -- A few input structures may be filtered by premin, these are problematic structures that it couldn't generated a conformation for, should be ok to exclude these

Glide
----- 

parameters
* outerbox: box within which the grids are calculated. This is also the box within which all the ligand atoms must be contained. The maximum size of the enclosing box is 50Å.
* innerbox: box explored by the ligand center (restricted to a cube whose sides cannot be longer than 40Å)

* DOCKING_METHOD = confgen ensure flexible docking
