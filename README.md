Mode-of-binding Predictor
=========================

Mode-of-binding Predictor (MOBPred) is a python package used to facilitate the use of popular docking software (including structure preparation, docking and rescoring). The package is particularly suitable to compare docking results obtained from different software or combine them in a consensus docking or consensus scoring strategy.

Note that the software that can be used by MOBPred are not part of the current repository. Whatever program needs to be utilized, it should be downloaded from its official website and installed separately on the same machine MOBPred is set up.

Below is a list of all the programs which can be used by MOBPred. 

* **Structure preparation/optimization**:

  * antechamber, prmchk, tleap, sander (AMBER12 or later, http://ambermd.org) used to assign partial charges, minimize structures...
  * ligprep (Schrodinger 2015, https://www.schrodinger.com/ligprep) used to prepare compounds (generate protonation states, isomers, conformers,...)
  * moebatch (MOE2015) used to identify probable binding sites
  * Open Babel (http://openbabel.org/wiki/Main_Page)

* **Docking**:

  * Autodock (http://autodock.scripps.edu)
  * Autodock Vina (http://autodock.scripps.edu)
  * DOCK 6 (http://dock.compbio.ucsf.edu/DOCK_6/index.htm)
  * Glide (https://www.schrodinger.com/glide)
  * GOLD (https://www.ccdc.cam.ac.uk/solutions/csd-discovery/components/gold/)
  * Induced Fit (https://www.schrodinger.com/induced-fit)
  * MOE2015 (https://www.chemcomp.com/MOE-Molecular_Operating_Environment.htm)

* **Scoring**:

  * Autodock (http://autodock.scripps.edu) 
  * Autodock Vina (http://autodock.scripps.edu)
  * DSX (http://pc1664.pharmazie.uni-marburg.de/drugscore/)
  * Glide (https://www.schrodinger.com/glide)
  * MOE2015 (https://www.chemcomp.com/MOE-Molecular_Operating_Environment.htm)

Table of contents
=================

  * [Prerequisites](#prerequisites)
  * [Installation](#installation)
  * [Commands](#commands)
    * [prepvs](#prepvs)
    * [rundock](#rundock)
    * [runanlz](#runanlz)
  * [Preparing the rundock config file](#preparing-the-rundock-config-file)
    * [General sections](#general-sections)
    * [Sections relative to each software](#sections-relative-to-each-software)
    * [Examples]

Prerequisites
=============

Before installing the MOBPred package, make sure that you have the following packages installed:

* NumPy; version 1.4.1 or later

* pandas; version 0.18.1 or later

* AmberTools; version 12 or later

Any software intended to be used in conjunction with MOBPred should be installed separetely and should work as a standalone program. In addition, make sure the applications mentioned below are in your PATH, depending on which docking/scoring software will be used:

* **Autodock**: except for babel and autodock4, all the executables below can be found in the AutoDockTools package (http://autodock.scripps.edu/downloads/resources/adt/index_html):
  * autodock4
  * autogrid4
  * babel
  * prepare_dpf4.py
  * prepare_gpf4.py
  * prepare_ligand4.py
  * prepare_receptor4.py

* **Autodock Vina**: 
  * babel
  * prepare_ligand4.py
  * prepare_receptor4.py
  * vina

* **DOCK 6**:
  * chimera (http://www.cgl.ucsf.edu/chimera)
  * dms (http://www.cgl.ucsf.edu/chimera/docs/UsersGuide/midas/dms1.html)
  * dock6
  * grid
  * showbox
  * sphere_selector
  * sphgen_cpp

* **DSX**:
  * dsx (symbolic link to dsx_linux_64.lnx or similar executables)

* **Glide**: All the executables can be found in the Schrodinger package (https://www.schrodinger.com):
  * glide
  * glide_sort
  * ligprep
  * pdbconvert
  * prepwizard

* **Gold**:
  * gold_auto

* **Induced Fit**: All the executables can be found in the Schrodinger package (https://www.schrodinger.com): 
  * ifd
  * prepwizard

* **MOE2015**:
  * moebatch

***Pharmamatrix users***: on the pharmamatrix cluster, the majority of the docking programs mentioned above have already been installed. Below is an example on how to set the PATH environment variable in order to use Autodock, Vina, Glide, MOE2015 and DOCK 6:

    # Amber 14
    export AMBERHOME=/pmshare/amber/ambertools14_ibm_gnu-20140820
    PATH=$PATH:$AMBERHOME/bin

    # AD and ADV
    PATH=$PATH:/opt/mgltools/1.5.4:/opt/mgltools/1.5.4/MGLToolsPckgs/AutoDockTools/Utilities24
    PATH=$PATH:/pmshare/vina/autodock_vina_1_1_2_linux_x86/bin
    PATH=$PATH:/opt/autodock/4.2.3/bin

    # Glide
    export SCHRODINGER=/nfs/r720-1/preto/tools/schrodinger2015-4
    PATH=$PATH:/nfs/r720-1/preto/tools/schrodinger2015-4
    PATH=$PATH:/nfs/r720-1/preto/tools/schrodinger2015-4/utilities

    # MOE2015
    export MOE=/nfs/r720-1/preto/tools/moe2015
    PATH=$PATH:/nfs/r720-1/preto/tools/moe2015/bin

    # DOCK 6
    PATH=$PATH:/nfs/r720-1/preto/tools/UCSF-Chimera64-1.10.2/bin
    PATH=$PATH:/nfs/r720-1/preto/tools/dms
    PATH=$PATH:/nfs/r720-1/preto/src/dock6/bin

    export PATH

Installation
============

The Python Distutils are used to build and install MOBPred, so it is fairly simple to get things ready to go. Following are very simple instructions on how to proceed:

1. First, make sure that you have the NumPy and pandas modules. If not, get them from http://numpy.scipy.org/, http://pandas.pydata.org. Compile/install them.

2. Make sure AmberTools is installed and that standard executables (e.g., sander, tleap,...) are accessible through your PATH variable. For pharmamatrix users, see section **prerequisites**.

3. From the main MOBPred distribution directory run this command (plus any extra flags, e.g., --prefix or --user to specify the installation directory):

        python setup.py install

After installation, make sure that the bin folder within your installation directory is accessible via your PATH variable. It contains the executables **prepvs**, **rundock** and **runanlz** that are used for virtual screening preparation, docking and docking analysis, respectively.


Commands
========

After adding the bin folder of your installation directory to your PATH variable, the following commands can be run: 


prepvs
------

prepvs is used to prepare the ligand structures (carried out with Schrodinger's ligprep utility) and to create folders (one folder per receptor and per ligand) intended to facilitate docking or virtual screening runs performed with rundock. When typing "prepvs -h" on the command line, the following help message will pop up:

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

rundock
-------

rundock is used to dock a ligand to a protein structure and eventually minimize and rescore the output poses. When typing "rundock -h" on the command line, the following help message will pop up:

    usage: rundock [-h] -l INPUT_FILE_L -r INPUT_FILE_R -f CONFIG_FILE
                   [-q CHARGE_FILE] [-rescore_only] [-extract_only] [-d POSEDIR]
                   [-norun]
    
    rundock : dock with multiple software -------- Requires one file for the
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

    * -l INPUT_FILE_L: **.mol2** file containing the coordinates of the ligand (only one structure allowed)

    * -r INPUT_FILE_R: **.pdb** file containing the receptor coordinates (only one structure allowed)

    * -f CONFIG_FILE: **.ini** configuration file containing the docking parameters (see the section **preparing the rundock configuration file**)

* Optional arguments

    Preferably do not use any flags other than -l, -r and -f

Thus, a typical use of **rundock** is done through the following command:

    rundock -f config.ini -r receptor.pdb -l ligand.mol2


Preparing the rundock configuration file
========================================

Besides one **.mol2** file containing the ligand structure (-l flag) and one **.pdb** file containing the receptor structure (-r flag), running **rundock** requires a configuration file (-f flag) that specifies all the parameters needed for the docking procedure.

**Note**: **rundock** can only be used to run docking and scoring procedures with a single protein and ligand structure. If multiple protein or/and ligand structures need to be used, the **prepvs** command can be used to create folders for each protein-ligand pair (see the above section **prepvs**). 

The rundock configuration file should be a .ini file (https://en.wikipedia.org/wiki/INI_file), i.e., the file should be split in sections, each section name appearing on a line by itself, in square brackets ("[" and "]"). Each section contains a certain number of keys which refer to specific options used; all keys after the section declaration are associated with that section. Finally, every key should have a name (option name) and a value (option value), delimited by an equals sign (=).

Below is an example of configuration file used to dock on two binding sites and rescore with DrugScoreX (dsx), Autodock and Autodock Vina.

    [DOCKING]
    site = site1, site2
    program = autodock, vina, dock, glide
    rescoring = yes
    minimize = yes
    cleanup = yes
    
    [RESCORING]
    program = dsx, autodock, vina
    
    [DSX]
    pot_dir = /pmshare/jordane/CSD_potentials/DSX_CSD_Potentials_v0511/csd_pot_0511/
    other_flags = -T0 1.0 -T1 1.0 -T3 1.0 -j
    
    [AUTODOCK]
    ga_run = 20
    spacing = 0.4
    
    [VINA]
    num_modes = 20
    
    [DOCK]
    nposes = 20
    
    [GLIDE]
    poses_per_lig = 20
    
    [SITE1]
    center = 75.5, 80.0, 31.0
    boxsize = 40.0, 40.0, 40.0
    
    [SITE2]
    center = 75.5, 40.0, 50.0
    boxsize = 40.0, 40.0, 40.0


General sections
----------------

* The **DOCKING** section includes the software that should be used for docking, and if minimization, rescoring and/or cleanup should be performed. The docking software should be specified with coma separation through the key **programs**. The keys relative to the **DOCKING** section are:

    * **programs**: specifies the software which are used for docking (autodock, dock6, glide, gold, moe and/or vina). Options relative to each program (or instance) are specfied within the section of the same name. For example, if autodock is in the list of programs, options associated with autodock should be specified in the **AUTODOCK** section. In case the same software needs to be used multiple times, numbering can be appended to the name of the program (e.g., in the first example below, multiple runs of MOE are performed using different scoring methods: moe, moe1, moe2).

    * **minimization**: performs minimization on the generated poses (yes or no).

    * **rescoring**: performs rescoring on the generated poses (yes or no). I strongly recommend to enable minimization in case rescoring is done. This will avoid a lot clashes, especially when the software used for rescoring are different from those used for docking. If the rescoring option is enabled, a section RESCORING should be created that contains all the options relative to that step (see below).

    * **cleanup**: specifies if big intermediate files should be removed (yes or no).

    * **site**: specifies the labels for the binding sites in case multiple binding sites are considered (site1, site2,...). See the example configuration to dock on multiple binding site, minimize and rescore the poses with multiple software.


    Below is a list of all the programs that can be used by MOBPred specifying if they can be used for docking or/and rescoring.


    |    Software   |    Docking    |   Rescoring   |
    | :-----------  |:-------------:|:-------------:|
    |   Autodock    |      Yes      |      Yes      |
    |    Dock 6     |      Yes      |      Yes      |
    |     DSX       |      No       |      Yes      |
    |     Glide     |      Yes      |      Yes      |
    |     Gold      |      Yes      |      No       |
    |  Induced Fit  |      Yes      |      No       |
    |     MOE       |      Yes      |      Yes      |
    |     Vina      |      Yes      |      Yes      |

    Docking and rescoring options relative to each program are detailed in the section **Docking/scoring options relative to each software**

* The **SITE** section includes the information about the box to spot the binding site. The keys are the following:

    *  **center**: x, y, z coordinates of the center of the binding box (in Å).

    *  **boxsize**: size of the box along each dimension x, y, z. The dimensions of the box should be no more than 50.0, 50.0, 50.0 (in Å).


* The **RESCORING** section has only one key specifying the programs used to rescore:

    *  **program**: specifies the software which are used for docking (autodock, dock6, glide, gold, moe and/or vina). Options relative to each program (or instance) are specfied within the section of the same name. For example, if autodock is in the list of programs, options associated with autodock should be specified in the **AUTODOCK** section. In case the same software needs to be used multiple times, numbering can be appended to the name of the program (e.g., in the example below, multiple runs of MOE are performed using different scoring methods: moe, moe1, moe2).

Docking/scoring options relative to each software
-------------------------------------------------

Each section relative to a docking/scoring program should be named the way it appears through the keys **program** of the **DOCKING** and/or **RESCORING** section. Below is a list of all the options per software that can be specified in the configuration file.

* **Autodock** (docking/scoring method)

    * ga_run (default: 100): number of autodock runs = targeted number of final poses
    * spacing (default: 0.3): grid spacing

    **Note 1**: the partial charges of the ligand are obtained from the Gasteiger method using the AutodockTools command *prepare_ligand4.py*

    **Note 2**: the number of energy evalutations *ga_num_evals* is automatically calculated from the number of torsions angles in the ligand structure via the formula:

        ga_num_evals = min(25000000, 987500 * n_torsion_angles + 125000)

    **Note 3**: As is usually the case for Autodock, non polar hydrogens in the ligand structure are removed prior to docking in order to properly use the Autodock force field. Once the docking has been performed, nonpolar hydrogens are reattributed in a way consistent with the input structure. Unless the *minimize* option in the configuration file is set to *yes*, no minimization is performed on those hydrogens.

    **Note 4** Final poses are extracted from the .dlg file using Open Babel via the following command:

        babel -ad -ipdbqt dock.dlg -omol2 lig-.mol2 -m

* **Autodock Vina** (docking/scoring method)

    * cpu (default: 1)
    * energy_range (default: 3)
    * num_modes (default: 9): targeted number of final poses

    **Note 1**: the partial charges of the ligand are obtained from the Gasteiger method using the AutodockTools command *prepare_ligand4.py*

    **Note 2**: As is usually the case for Autodock Vina, non polar hydrogens in the ligand structure are removed prior to docking in order to properly use the Autodock force field. Once the docking has been performed, nonpolar hydrogens are reattributed in a way consistent with the input structure. Unless the *minimize* option in the configuration file is set to *yes*, no minimization is performed on those hydrogens.


* **DOCK 6** (docking method)

    * attractive_exponent (default: 6)
    * extra_margin (default: 2.0)
    * grid_spacing (default: 0.3)
    * maximum_sphere_radius (default: 4.0)
    * max_orientations (default: 10000)
    * minimum_sphere_radius (default: 1.4)
    * nposes (default: 20): targeted number of final poses
    * num_scored_conformers (default 5000)
    * probe_radius (default: 1.4)
    * repulsive_exponent (default: 12)

* **DSX** (scoring method)

* **Glide** (docking/scoring)

    * pose_rmsd (default: 0.5):
    * poses_per_lig (default: 10): targeted number of final poses
    * precision (default: SP):
    * use_prepwizard (default: True):

* **GOLD**

    * nposes (default: 20)

* **MOE**

    * gtest (default: 0.01)
    * maxpose (default: 5)
    * placement (default: Triangle Matcher)
    * placement_maxpose (default: 250)
    * placement_nsample (default: 10)
    * remaxpose (default: 1)
    * rescoring (default: GBVI/WSA dG)
    * scoring (default: London dG)



Examples
--------

Docking with multiple software on a single binding site and minimize the poses
-------------------------------------------------------------------------------

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



Docking on multiple binding site, minimize and rescore the poses with multiple software
----------------------------------------------------------------------------------------

Below is another example of configuration file for *rundock* used to dock on two binding sites and rescore with DrugScoreX (dsx), Autodock and Autodock Vina.

    [DOCKING]
    site = site1, site2
    program = autodock, vina, dock, glide
    rescoring = yes
    minimize = yes
    cleanup = yes
    
    [RESCORING]
    program = dsx, autodock, vina
    
    [DSX]
    pot_dir = /pmshare/jordane/CSD_potentials/DSX_CSD_Potentials_v0511/csd_pot_0511/
    other_flags = -T0 1.0 -T1 1.0 -T3 1.0 -j
    
    [AUTODOCK]
    ga_run = 20
    spacing = 0.4
    
    [VINA]
    num_modes = 20
    
    [DOCK]
    nposes = 20
    
    [GLIDE]
    poses_per_lig = 20
    
    [SITE1]
    center = 75.5, 80.0, 31.0
    boxsize = 40.0, 40.0, 40.0
    
    [SITE2]
    center = 75.5, 40.0, 50.0
    boxsize = 40.0, 40.0, 40.0

* Note that the DOCKING section includes the label of the binding sites through the keyword *site*, here, site1 and site2. Each label refers to the section of the same name SITE1 and SITE2, respectively. 


Docking/scoring options relative to each software
-------------------------------------------------

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

