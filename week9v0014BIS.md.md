
# Internship - Week 9 - Summary

At this point, and considering that we do not count with CPU time at Genius, it seems like a good idea trying to gather input and output files corresponding to all the calculations carried out so far with the Ab-Initio/MD software. On the one hand, the list of tar files uploaded to the "computationalchemistry" repository is completed. On the other, job submission scripts and input files are exhibited in this document. For cases of sequential calculations or parallel calculations on the laptop, no submission script is shown (only the main input file).  

 On the other hand, installation work (on DSI guest2 laptop), is summarised below, as well as work on the first runs with the software below:

 - GROMACS
 - PLUMED
 - AlphaFold (not installed by the intern, simply run on Genius, when/if possible)
 - OPENMM
 - OPENPLATHSAMPLING
 
The last two codes are not separately installed on Genius. Actually, on the HPC version, they take part in the distribution of other software. Firstly, PLUMED is a free-energy calculation code that makes use of OPENPATHSAMPLING. Secondly, AlphaFold is a protein structure determination sofware that uses OPENMM as its force-engine, and therefore, OpenMM is part of its distribution. Regarding the laptop installations of OPENMM and OPENPATHSAMPLING, explicit installation have been carried out. After looking at the AlphaFold documentation, and seeing that its local installation would require a huge amount of storage space (about 3 TB for the database), not available on the laptop, no attempt of installing AlphaFold is made.

## HPC calculations: tar files uploaded and inputs/job submission scripts - Development of supercomputing training material regarding Ab-Initio/MD calculations

Following up the same order appearing in the computationalchemistry repository (README.md), job submission scripts and input files are shown next:

### ABINIT

 - **Sequential** ABINIT calculation of the as provided example "tbasepar_1", located within abinit-test/tutorial ("Lead crystal - parallelism over k-points"): 
   *sequentialABINITtbasepar_1.tar.gz* (Genius),
```
#
# Lead crystal
#

# Simulation parameters
ecut  30.0
acell 10.0 10.0 10.0
rprim
    0.0  0.5  0.5
    0.5  0.0  0.5
    0.5  0.5  0.0

# K-points
ngkpt 8 8 8
nshiftk 4
shiftk
   0.5 0.5 0.5
   0.5 0.0 0.0
   0.0 0.5 0.0
   0.0 0.0 0.5
occopt 7
tsmear 0.01

# System description
ntypat 1
znucl  82
natom  1
typat  1 
xred
  0.000   0.000   0.000
nband  4

# SCF procedure
nstep   3
tolvrs  1.0d-10


## After modifying the following section, one might need to regenerate the pickle database with runtests.py -r
#%%<BEGIN TEST_INFO>
#%% [setup]
#%% executable = abinit
#%% [files]
#%% files_to_test = tbasepar_1.out, tolnlines=0, tolabs=0.0, tolrel=0.0
#%% psp_files= HGH/82pb.4.hgh
#%% [paral_info]
#%% max_nprocs = 4
#%% [extra_info]
#%% authors = Unknown
#%% keywords = NC
#%% description = Lead crystal. Parallelism over k-points
#%%<END TEST_INFO>
```
 - **Sequential** ABINIT calculation of the as provided example "tbasepar_2", located within abinit-test/tutorial ("FCC 4-atom A1-phonon deformed ferromagnetic Fe 
   crystal" - parallelism over spin): *sequentialABINITtbasepar_2.tar.gz* (Genius)
```
# FCC Fe (ferromagnetic for fun) with four atoms per cell
# Distorted with a A1 phonon, so as to keep the symmetry ...
# Only one k point in the IBZ
# Test the parallelism over the spins 

  tolvrs  1.0d-13
   ngkpt  2 2 2    

 ecut 30
 acell 3*7.00  

 ixc    1
 natom  4 
 nband 40
 nline 5
 nsppol 2
 spinat 0.0 0.0 3.0
        0.0 0.0 3.0
        0.0 0.0 3.0
        0.0 0.0 3.0
 nstep 5
 ntypat  1
 occopt 7
 shiftk 0.5 0.5 0.5
 typat  4*1
 xred  0.01 0.01 0.01
       0.49 0.49 0.01
       0.49 0.01 0.49
       0.01 0.49 0.49
 znucl 26.0



## After modifying the following section, one might need to regenerate the pickle database with runtests.py -r
#%%<BEGIN TEST_INFO>
#%% [setup]
#%% executable = abinit
#%% [files]
#%% files_to_test = tbasepar_2.out, tolnlines=0, tolabs=0.0, tolrel=0.0
#%% psp_files= HGH/26fe.8.hgh
#%% [paral_info]
#%% max_nprocs = 2
#%% [extra_info]
#%% keywords = NC
#%% authors = Unknown
#%% description = 
#%%   FCC Fe (ferromagnetic for fun) with four atoms per cell
#%%   Distorted with a A1 phonon, so as to keep the symmetry ...
#%%   Only one k point in the IBZ
#%%   Test the parallelism over the spins , 
#%%<END TEST_INFO>
```
 - **MPI** ABINIT calculation of the as provided example "tdfpt", located within abinit-test/tutoparal ("DFPT phonon calculation for Aluminium crystal - step 1: 
   ground-state electronic structure calcuation - step 2: phonon calculation": *mpiABINITtdfpt1and2.tar.gz* (DSI laptop guest2; no batch job script used)
```bash
#!/bin/bash

mpirun abinit < tdfpt_01.files > tdfpt_01.log
cp tdfpt_01.o_WFK tdfpt_02.i_WFK
cp tdfpt_01.o_WFK tdfpt_02.i_WFQ
mpirun abinit < tdfpt_02.files > tdfpt_02.log
```
```
#   FCC Al; 10 special points

#timopt -1

 acell 3*7.56
 densty 1.2
 ecut 10

 enunit 2  

 localrdwf 1
 ngkpt 8 8 8 
 nshiftk 4
 shiftk 0.5 0.5 0.5
        0.5 0.0 0.0
        0.0 0.5 0.0
        0.0 0.0 0.5

 natom  1 nband 5
 nline 3  nstep 20
 ntypat  1
 occopt  4  prtden 1   prtvol 10
 rprim   0 .5 .5  .5 0 .5  .5 .5 0
 timopt 2
 tnons   72*0.0d0
 tolvrs 1.0d-18
 typat  1  
 xred  0.0 0.0 0.0
 znucl 13.0
```
```
#   FCC Al; phonon at 1/4 1/8 1/8

 istatr 1000
 irdwfk 1
 irdwfq 1

 nbdbuf 3
 rfphon 1
 rfatpol 1 1
 rfdir  1 1 1
 nqpt   1
 qpt    0.25 -0.125 0.125

 acell 3*7.56
 densty 1.2

 ecut 10

 enunit 2  

 kptopt 3
 localrdwf 1
 ngkpt 8 8 8  
 nshiftk 4
 shiftk 0.5 0.5 0.5
        0.5 0.0 0.0
        0.0 0.5 0.0
        0.0 0.0 0.5 

 natom  1 nband 5

 nstep 20

 ntypat  1
 amu    26.96
 occopt  4 
 rprim   0 .5 .5  .5 0 .5  .5 .5 0
 timopt 2
 tolvrs 1.0d-10
 typat  1  
 xred  0.0 0.0 0.0
 znucl 13.0
```
### ASE

 - **Sequential** ASE convex hull determination for a CuxPt1-x(111) surface slab by means of a genetic algorithm: 
   *sequentialASEgenalgCuPt111slab.tar.gz* (DSI guest2 laptop)
```python
from pathlib import Path
import random

from ase.build import fcc111
from ase.calculators.emt import EMT
from ase.data import atomic_numbers, reference_states
from ase.ga.data import PrepareDB
from ase.ga import set_raw_score


def get_avg_lattice_constant(syms):
    a = 0.
    for m in set(syms):
        a += syms.count(m) * lattice_constants[m]
    return a / len(syms)


metals = ['Cu', 'Pt']
# Use experimental lattice constants
lattice_constants = dict((m, reference_states[atomic_numbers[m]]['a'])
                         for m in metals)

# Create the references (pure slabs) manually
pure_slabs = []
refs = {}
print('Reference energies:')
for m in metals:
    slab = fcc111(m, size=(2, 4, 3), a=lattice_constants[m],
                  vacuum=5, orthogonal=True)
    slab.calc = EMT()

    # We save the reference energy as E_A / N
    e = slab.get_potential_energy()
    e_per_atom = e / len(slab)
    refs[m] = e_per_atom
    print('{0} = {1:.3f} eV/atom'.format(m, e_per_atom))

    # The mixing energy for the pure slab is 0 by definition
    set_raw_score(slab, 0.0)
    pure_slabs.append(slab)

# The population size should be at least the number of different compositions
pop_size = 2 * len(slab)

# We prepare the db and write a few constants that we are going to use later
target = Path('hull.db')
if target.exists():
    target.unlink()
db = PrepareDB(target, population_size=pop_size,
               reference_energies=refs, metals=metals,
               lattice_constants=lattice_constants)

# We add the pure slabs to the database as relaxed because we have already
# set the raw_score
for slab in pure_slabs:
    db.add_relaxed_candidate(slab,
                             atoms_string=''.join(slab.get_chemical_symbols()))


# Now we create the rest of the candidates for the initial population
for i in range(pop_size - 2):
    # How many of each metal is picked at random, making sure that
    # we do not pick pure slabs
    nA = random.randint(0, len(slab) - 2)
    nB = len(slab) - 2 - nA
    symbols = [metals[0]] * nA + [metals[1]] * nB + metals

    # Making a generic slab with the correct lattice constant
    slab = fcc111('X', size=(2, 4, 3),
                  a=get_avg_lattice_constant(symbols),
                  vacuum=5, orthogonal=True)

    # Setting the symbols and randomizing the order
    slab.set_chemical_symbols(symbols)
    random.shuffle(slab.numbers)

    # Add these candidates as unrelaxed, we will relax them later
    atoms_string = ''.join(slab.get_chemical_symbols())
    db.add_unrelaxed_candidate(slab, atoms_string=atoms_string)
```
```python
import numpy as np
from ase.ga.population import RankFitnessPopulation
from ase.ga.data import DataConnection
from ase.ga.offspring_creator import OperationSelector
from ase.ga.slab_operators import (CutSpliceSlabCrossover,
                                   RandomSlabPermutation,
                                   RandomCompositionMutation)
from ase.ga import set_raw_score

from ase.calculators.emt import EMT

# Connect to the database containing all candidates
db = DataConnection('hull.db')

# Retrieve saved parameters
pop_size = db.get_param('population_size')
refs = db.get_param('reference_energies')
metals = db.get_param('metals')
lattice_constants = db.get_param('lattice_constants')


def get_mixing_energy(atoms):
    # Set the correct cell size from the lattice constant
    new_a = get_avg_lattice_constant(atoms.get_chemical_symbols())
    # Use the orthogonal fcc cell to find the current lattice constant
    current_a = atoms.cell[0][0] / np.sqrt(2)
    atoms.set_cell(atoms.cell * new_a / current_a, scale_atoms=True)

    # Calculate the energy
    atoms.calc = EMT()
    e = atoms.get_potential_energy()

    # Subtract contributions from the pure element references
    # to get the mixing energy
    syms = atoms.get_chemical_symbols()
    for m in set(syms):
        e -= syms.count(m) * refs[m]
    return e


def get_avg_lattice_constant(syms):
    a = 0.
    for m in set(syms):
        a += syms.count(m) * lattice_constants[m]
    return a / len(syms)


def get_comp(atoms):
    return atoms.get_chemical_formula()


# Specify the number of generations this script will run
num_gens = 10

# Specify the procreation operators for the algorithm and
# how often each is picked on average
# The probability for an operator is the prepended integer divided by the sum
# of integers
oclist = [(3, CutSpliceSlabCrossover()),
          (1, RandomSlabPermutation()),
          (1, RandomCompositionMutation())
          ]
operation_selector = OperationSelector(*zip(*oclist))

# Pass parameters to the population instance
# A variable_function is required to divide candidates into groups here we use
# the chemical composition
pop = RankFitnessPopulation(data_connection=db,
                            population_size=pop_size,
                            variable_function=get_comp)

# Evaluate the starting population
# The only requirement of the evaluation is to set the raw_score
# Negative mixing energy means more stable than the pure slabs
# The optimization always progress towards larger raw score,
# so we take the negative mixing energy as the raw score
print('Evaluating initial candidates')
while db.get_number_of_unrelaxed_candidates() > 0:
    a = db.get_an_unrelaxed_candidate()
    set_raw_score(a, -get_mixing_energy(a))
    db.add_relaxed_step(a)
pop.update()

# Below is the iterative part of the algorithm
gen_num = db.get_generation_number()
for i in range(num_gens):
    print('Creating and evaluating generation {0}'.format(gen_num + i))
    new_generation = []
    for _ in range(pop_size):
        # Select parents for a new candidate
        parents = pop.get_two_candidates()

        # Select an operator and use it
        op = operation_selector.get_operator()
        offspring, desc = op.get_new_individual(parents)
        # An operator could return None if an offspring cannot be formed
        # by the chosen parents
        if offspring is None:
            continue

        set_raw_score(offspring, -get_mixing_energy(offspring))
        new_generation.append(offspring)

    # We add a full relaxed generation at once, this is faster than adding
    # one at a time
    db.add_more_relaxed_candidates(new_generation)

    # update the population to allow new candidates to enter
    pop.update()
```
```python
import numpy as np
import matplotlib.pyplot as plt
from ase.phasediagram import PhaseDiagram
from ase.db import connect
from ase.io import write

db = connect('hull.db')

# Select the evaluated candidates and retrieve the chemical formula and mixing
# energy for the phase diagram
refs = []
dcts = list(db.select('relaxed=1'))
for dct in dcts:
    refs.append((dct.formula, -dct.raw_score))

pd = PhaseDiagram(refs)
ax = pd.plot(show=not True,  # set to True to show plot
             only_label_simplices=True)
plt.savefig('hull.png')

# View the simplices of the convex hull
simplices = []
toview = sorted(np.array(dcts)[pd.hull], key=lambda x: x.mass)
for dct in toview:
    simplices.append(dct.toatoms())

write('hull.traj', simplices)
```

   ![](hull.png)

 - **MPI** ASE water-box equilibration as a strong-scaling **MPI** test: *mpiASE_STRONGSCALING_waterboxequi.tar.gz* (Genius)

### LAMMPS

 - **MPI** LAMMPS calculation of c-HfO2: *mpiLAMMPScHfO2.tar.gz* (Genius)
 - **OpenMP** LAMMPS ethanol-box equlibration as a strong-scaling **OpenMP** test: *openmpLAMMPSstrongscalingETHANOLBOX.tar.gz* (Genius)
 - **KOKKOS-HYBRID** LAMMPS ethanol-box equilibration as a KOKKOS-hybrid (36/2) parallel calculation case: KOKKOSPARALLhybrid36and2ethanol.tar.gz (Genius)
 - **Sequential** LAMMPS Monte-Carlo relaxation of a two-dimensional deformed hexagonal lattice: *sequentialLAMMPSmc.tar.gz* (Genius)
 - **MPI** LAMMPS calculation of quartz amorphisation via melting/quenching temperature ramps, as an **MPI** strong-scaling case: *mpiLAMMPSstrongscalsilicaamorph.tar.gz* (Genius)

### QUANTUM-ESPRESSO

 - **Sequential** GaAs PWscf calculation under PAW: *sequentialQEgaas.tar.gz* (Genius),
 ```
 &control
   calculation   = 'scf'
   restart_mode  = 'from_scratch'
   prefix        = 'GaAs'
   pseudo_dir    = '.'
   outdir        = '.'
   verbosity     = 'high'
   etot_conv_thr = 1.d-5
   forc_conv_thr = 1.d-4
/
&system
   ibrav = 2,
   celldm(1) = 10.6867,
   nat = 2,
   ntyp = 2,
   ecutwfc = 50,
   ecutrho = 600,
   nspin = 1
   nbnd  = 16
   la2F  = .true.
/
&electrons
   mixing_beta      = 0.3
   conv_thr         = 1.0d-8
   electron_maxstep = 200,
/
&ions
/
&cell
  press = 0.0
/
ATOMIC_SPECIES
   Ga  69.723    Ga.pbe-dn-kjpaw_psl.0.2.upf
   As  74.921595 As.pbe-n-kjpaw_psl.0.2.upf
ATOMIC_POSITIONS
 Ga  0.00 0.00 0.00
 As  0.25 0.25 0.25
K_POINTS {automatic}
   32 32 32 0 0 0
```
 - **MPI** SiH4 -molecule in a box- PWscf calculation under Hamann's ONCV pseudopotentials: *mpiQEsih4.tar.gz* (Genius),
```bash
 #!/bin/bash

#SBATCH --job-name=mpiesp
#SBATCH --nodes=1
#SBATCH --time=00:05:00
#SBATCH --account=lp_h_vsc35663
#SBATCH --ntasks=8           # TOTAL Number of cores
#SBATCH --cpus-per-task=1    # Number of OpenMP threads per task
#SBATCH --cluster=genius

module purge
module use /apps/leuven/skylake/2021a/modules/all
module load QuantumESPRESSO/6.8-intel-2021a

mpirun pw.x < silane.in > scf.out
```
```
&control
calculation  = 'scf'
restart_mode = 'from_scratch'
pseudo_dir   = './'
outdir       = './'
prefix       = 'silane'
wf_collect   = .TRUE.
/
&system
ibrav           = 1
celldm(1)       = 20
nat             = 5
ntyp            = 2
ecutwfc         = 25.0
nbnd            = 10
assume_isolated ='mp'
/
&electrons
diago_full_acc = .TRUE.
/
ATOMIC_SPECIES
Si 28.0855  Si_ONCV_PBE-1.2.upf
H  1.00794   H_ONCV_PBE-1.2.upf
ATOMIC_POSITIONS bohr
Si      10.000000   10.000000  10.000000
H       11.614581   11.614581  11.614581
H        8.385418    8.385418  11.614581
H        8.385418   11.614581   8.385418
H       11.614581    8.385418   8.385418
K_POINTS {gamma}
``` 
 - **MPI** *Immm* Ag2PdO2 PWscf calculation under ultrasoft pseudopotentials: *mpiQE_Ag2PdO2.tar.gz* (Genius),
```bash
#!/bin/bash

#SBATCH --job-name=mpiesp
#SBATCH --nodes=1
#SBATCH --time=00:59:00
#SBATCH --account=lp_h_vsc35663
#SBATCH --ntasks=8           # TOTAL Number of cores
#SBATCH --cpus-per-task=1    # Number of OpenMP threads per task
#SBATCH --cluster=genius

module purge
module use /apps/leuven/skylake/2021a/modules/all
module load QuantumESPRESSO/6.8-intel-2021a

mpirun pw.x < scf.in > scf.out
```
```
&CONTROL
  calculation = 'scf'
  etot_conv_thr =   1.0000000000d-04
  forc_conv_thr =   1.0000000000d-04
  outdir = './out/'
  prefix = ''
  pseudo_dir = './'
  tprnfor = .true.
  tstress = .true.
  verbosity = 'high'
/
&SYSTEM
  degauss =   1.0000000000d-03
  ecutwfc         =       59.00000000
  ecutrho         =      545.00000000
  ibrav = 0
  nat = 10
  nosym = .false.
  ntyp = 3
  occupations = 'smearing'
  smearing        =  'f-d'
  noncolin        =  .true.
  lspinorb        =  .true.
/
&ELECTRONS
  conv_thr =   2.0000000000d-09
  electron_maxstep = 80
  mixing_beta =   4.0000000000d-01
/
ATOMIC_SPECIES
Ag     107.8682 Ag.rel-pbe-n-rrkjus_psl.1.0.0.UPF
O      15.9994 O.rel-pbe-n-rrkjus_psl.1.0.0.UPF
Pd     106.42 Pd.rel-pbe-spn-rrkjus_psl.1.0.0.UPF
ATOMIC_POSITIONS crystal
Pd           0.0000000000       0.0000000000       0.0000000000
Pd           0.5000000000       0.5000000000       0.5000000000
O            0.5000000000       0.0000000000       0.3627000000
O            0.5000000000       0.0000000000       0.6373000000
O            0.0000000000       0.5000000000       0.8627000000
O            0.0000000000       0.5000000000       0.1373000000
Ag           0.0000000000       0.0000000000       0.3575800000
Ag           0.0000000000       0.0000000000       0.6424200000
Ag           0.5000000000       0.5000000000       0.8575800000
Ag           0.5000000000       0.5000000000       0.1424200000
K_POINTS automatic
7 11 4 0 0 0
CELL_PARAMETERS angstrom
      4.5552300000       0.0000000000       0.0000000000
      0.0000000000       3.0080300000       0.0000000000
      0.0000000000       0.0000000000       9.8977000000
```
 - **MPI** PWscf calculation of 112 atom-cell Gold surface, from the PRACE benchmark (UEABS): *mpiQEgoldsurface.tar.gz* (Genius),
```bash 
#!/bin/bash

#SBATCH --job-name=qempi
#SBATCH --nodes=1
#SBATCH --time=00:59:00
#SBATCH --account=lp_h_vsc35663
#SBATCH --ntasks=36          # TOTAL Number of cores
#SBATCH --cpus-per-task=1    # Number of OpenMP threads per task
#SBATCH --cluster=genius

module purge
module use /apps/leuven/skylake/2021a/modules/all
module load QuantumESPRESSO/6.8-intel-2021a

mpirun pw.x < ausurf.in > ausurf.out
```
```
&CONTROL
  title = ' DEISA pw benchmark ',
  calculation = 'scf',
  restart_mode = 'from_scratch', ! 'restart',
  tprnfor = .TRUE.,
  etot_conv_thr = 1.d-8,
  prefix = 'ausurf'
  pseudo_dir = './'
  outdir = './out/'
/

&SYSTEM
  ibrav = 8, 
  celldm(1) = 38.7583, 
  celldm(2) = 0.494393, 
  celldm(3) = 1.569966, 
  nat = 112,
  ntyp = 1,
  nbnd = 800, 
  ecutwfc = 25.0,
  ecutrho = 200.0,
  occupations='smearing', smearing='marzari-vanderbilt', degauss=0.05
/

&ELECTRONS
    diagonalization='david'
    mixing_beta = 0.7
/

&IONS
  ion_dynamics = 'none',
/

&CELL
  cell_dynamics = 'none',
/

ATOMIC_SPECIES
 AU  196.96  Au.pbe-nd-van.UPF  

K_POINTS (automatic)
2 2 1 1 1 0

ATOMIC_POSITIONS (angstrom)
AU       29.285000       40.578999        7.173000
AU       29.285000       35.511002        7.173000
AU       32.214001       40.578999        7.173000
AU       32.214001       35.511002        7.173000
AU       35.141998       40.578999        7.173000
AU       35.141998       35.511002        7.173000
AU       38.070999       40.578999        7.173000
AU       38.070999       35.511002        7.173000
AU       40.999001       40.578999        7.173000
AU       40.999001       35.511002        7.173000
AU       43.928001       40.578999        7.173000
AU       43.928001       35.511002        7.173000
AU       46.855999       40.578999        7.173000
AU       46.855999       35.511002        7.173000
AU       29.285000       42.270000        4.782000
AU       29.285000       37.202000        4.782000
AU       32.214001       42.270000        4.782000
AU       32.214001       37.202000        4.782000
AU       35.141998       42.270000        4.782000
AU       35.141998       37.202000        4.782000
AU       38.070999       42.270000        4.782000
AU       38.070999       37.202000        4.782000
AU       40.999001       42.270000        4.782000
AU       40.999001       37.202000        4.782000
AU       43.928001       42.270000        4.782000
AU       43.928001       37.202000        4.782000
AU       46.855999       42.270000        4.782000
AU       46.855999       37.202000        4.782000
AU       30.749001       43.115002        7.173000
AU       30.749001       38.047001        7.173000
AU       33.678001       43.115002        7.173000
AU       33.678001       38.047001        7.173000
AU       36.605999       43.115002        7.173000
AU       36.605999       38.047001        7.173000
AU       39.535000       43.115002        7.173000
AU       39.535000       38.047001        7.173000
AU       42.464001       43.115002        7.173000
AU       42.464001       38.047001        7.173000
AU       45.391998       43.115002        7.173000
AU       45.391998       38.047001        7.173000
AU       48.320999       43.115002        7.173000
AU       48.320999       38.047001        7.173000
AU       30.749001       34.666000        4.782000
AU       30.749001       39.737999        4.782000
AU       33.678001       34.666000        4.782000
AU       33.678001       39.737999        4.782000
AU       36.605999       34.666000        4.782000
AU       36.605999       39.737999        4.782000
AU       39.535000       34.666000        4.782000
AU       39.535000       39.737999        4.782000
AU       42.464001       34.666000        4.782000
AU       42.464001       39.737999        4.782000
AU       45.391998       34.666000        4.782000
AU       45.391998       39.737999        4.782000
AU       48.320999       34.666000        4.782000
AU       48.320999       39.737999        4.782000
AU       29.285000       40.578999        0.000000
AU       29.285000       35.511002        0.000000
AU       32.214001       40.578999        0.000000
AU       32.214001       35.511002        0.000000
AU       35.141998       40.578999        0.000000
AU       35.141998       35.511002        0.000000
AU       38.070999       40.578999        0.000000
AU       38.070999       35.511002        0.000000
AU       40.999001       40.578999        0.000000
AU       40.999001       35.511002        0.000000
AU       43.928001       40.578999        0.000000
AU       43.928001       35.511002        0.000000
AU       46.855999       40.578999        0.000000
AU       46.855999       35.511002        0.000000
AU       30.749001       41.424000        2.391000
AU       30.749001       36.355999        2.391000
AU       33.678001       41.424000        2.391000
AU       33.678001       36.355999        2.391000
AU       36.605999       41.424000        2.391000
AU       36.605999       36.355999        2.391000
AU       39.535000       41.424000        2.391000
AU       39.535000       36.355999        2.391000
AU       42.464001       41.424000        2.391000
AU       42.464001       36.355999        2.391000
AU       45.391998       41.424000        2.391000
AU       45.391998       36.355999        2.391000
AU       48.320999       41.424000        2.391000
AU       48.320999       36.355999        2.391000
AU       29.285000       43.959999        2.391000
AU       29.285000       38.893002        2.391000
AU       32.214001       43.959999        2.391000
AU       32.214001       38.893002        2.391000
AU       35.141998       43.959999        2.391000
AU       35.141998       38.893002        2.391000
AU       38.070999       43.959999        2.391000
AU       38.070999       38.893002        2.391000
AU       40.999001       43.959999        2.391000
AU       40.999001       38.893002        2.391000
AU       43.928001       43.959999        2.391000
AU       43.928001       38.893002        2.391000
AU       46.855999       43.959999        2.391000
AU       46.855999       38.893002        2.391000
AU       30.749001       43.115002        0.000000
AU       30.749001       38.047001        0.000000
AU       33.678001       43.115002        0.000000
AU       33.678001       38.047001        0.000000
AU       36.605999       43.115002        0.000000
AU       36.605999       38.047001        0.000000
AU       39.535000       43.115002        0.000000
AU       39.535000       38.047001        0.000000
AU       42.464001       43.115002        0.000000
AU       42.464001       38.047001        0.000000
AU       45.391998       43.115002        0.000000
AU       45.391998       38.047001        0.000000
AU       48.320999       43.115002        0.000000
AU       48.320999       38.047001        0.000000
```
 - **MPI** PWscf calculation of 443-atom cell Iridium carbide under PAW, from the PRACE benchmark: *mpiQEiridiumcarbide.tar.gz* (Genius),
```bash
#!/bin/bash

#SBATCH --job-name=mpi1152
#SBATCH --nodes=32
#SBATCH --time=02:00:00
#SBATCH --account=lp_h_vsc35663
#SBATCH --ntasks=1152        # TOTAL Number of cores
#SBATCH --cpus-per-task=1    # Number of OpenMP threads per task
#SBATCH --cluster=genius

module purge
module use /apps/leuven/skylake/2021a/modules/all
module load QuantumESPRESSO/6.8-intel-2021a

mpirun pw.x -pd .true. -npool 4 -ndiag 16 -input pw.in > med.out
```
```
 &control
    calculation = 'scf'
    prefix='GRIR'
    max_seconds=20000 !42000
    restart_mode='from_scratch'
    pseudo_dir='./',
    outdir='./',
    disk_io = 'none'
 /
 &system
    ibrav=  4
    celldm(1) = 46.5334237988185d0
    celldm(3) =  1.15276608084821
    nat=  443
    ntyp= 2,
    ecutwfc=30
    occupations = 'smearing'
        smearing='mv'
        degauss=0.025d0
    nspin = 2 
    starting_magnetization(1) = +.00
    starting_magnetization(2) = +.00
 /
 &electrons
    conv_thr =  1.0d-5
    mixing_beta=0.3d0
    mixing_mode='local-TF'
    startingwfc='atomic'
    diagonalization='david'
    electron_maxstep = 1
  /
ATOMIC_SPECIES
 C    12.010   C.pbe-paw_kj-x.UPF
 Ir  192.22   Ir.pbe-paw_kj.UPF
K_POINTS {automatic}
2 2 2 0 0 0
ATOMIC_POSITIONS (crystal)
C   0.00000000000    0.00000000000    0.29309767430
 C  0.03333333333    0.06666666667    0.29309767430
C   0.10000000000    0.00000000000    0.29309767430
 C  0.13333333333    0.06666666667    0.29309767430
C   0.20000000000    0.00000000000    0.29309767430
 C  0.23333333333    0.06666666667    0.29309767430
C   0.30000000000    0.00000000000    0.29309767430
 C  0.33333333333    0.06666666667    0.29309767430
C   0.40000000000    0.00000000000    0.29309767430
 C  0.43333333333    0.06666666667    0.29309767430
C   0.50000000000    0.00000000000    0.29309767430
 C  0.53333333333    0.06666666667    0.29309767430
C   0.60000000000    0.00000000000    0.29309767430
 C  0.63333333333    0.06666666667    0.29309767430
C   0.70000000000    0.00000000000    0.29309767430
 C  0.73333333333    0.06666666667    0.29309767430
C   0.80000000000    0.00000000000    0.29309767430
 C  0.83333333333    0.06666666667    0.29309767430
C   0.90000000000    0.00000000000    0.29309767430
 C  0.93333333333    0.06666666667    0.29309767430
C   0.00000000000    0.10000000000    0.29309767430
 C  0.03333333333    0.16666666667    0.29309767430
C   0.10000000000    0.10000000000    0.29309767430
 C  0.13333333333    0.16666666667    0.29309767430
C   0.20000000000    0.10000000000    0.29309767430
 C  0.23333333333    0.16666666667    0.29309767430
C   0.30000000000    0.10000000000    0.29309767430
 C  0.33333333333    0.16666666667    0.29309767430
C   0.40000000000    0.10000000000    0.29309767430
 C  0.43333333333    0.16666666667    0.29309767430
C   0.50000000000    0.10000000000    0.29309767430
 C  0.53333333333    0.16666666667    0.29309767430
C   0.60000000000    0.10000000000    0.29309767430
 C  0.63333333333    0.16666666667    0.29309767430
C   0.70000000000    0.10000000000    0.29309767430
 C  0.73333333333    0.16666666667    0.29309767430
C   0.80000000000    0.10000000000    0.29309767430
 C  0.83333333333    0.16666666667    0.29309767430
C   0.90000000000    0.10000000000    0.29309767430
 C  0.93333333333    0.16666666667    0.29309767430
C   0.00000000000    0.20000000000    0.29309767430
 C  0.03333333333    0.26666666667    0.29309767430
C   0.10000000000    0.20000000000    0.29309767430
 C  0.13333333333    0.26666666667    0.29309767430
C   0.20000000000    0.20000000000    0.29309767430
 C  0.23333333333    0.26666666667    0.29309767430
C   0.30000000000    0.20000000000    0.29309767430
 C  0.33333333333    0.26666666667    0.29309767430
C   0.40000000000    0.20000000000    0.29309767430
 C  0.43333333333    0.26666666667    0.29309767430
C   0.50000000000    0.20000000000    0.29309767430
 C  0.53333333333    0.26666666667    0.29309767430
C   0.60000000000    0.20000000000    0.29309767430
 C  0.63333333333    0.26666666667    0.29309767430
C   0.70000000000    0.20000000000    0.29309767430
 C  0.73333333333    0.26666666667    0.29309767430
C   0.80000000000    0.20000000000    0.29309767430
 C  0.83333333333    0.26666666667    0.29309767430
C   0.90000000000    0.20000000000    0.29309767430
 C  0.93333333333    0.26666666667    0.29309767430
C   0.00000000000    0.30000000000    0.29309767430
 C  0.03333333333    0.36666666667    0.29309767430
C   0.10000000000    0.30000000000    0.29309767430
 C  0.13333333333    0.36666666667    0.29309767430
C   0.20000000000    0.30000000000    0.29309767430
 C  0.23333333333    0.36666666667    0.29309767430
C   0.30000000000    0.30000000000    0.29309767430
 C  0.33333333333    0.36666666667    0.29309767430
C   0.40000000000    0.30000000000    0.29309767430
 C  0.43333333333    0.36666666667    0.29309767430
C   0.50000000000    0.30000000000    0.29309767430
 C  0.53333333333    0.36666666667    0.29309767430
C   0.60000000000    0.30000000000    0.29309767430
 C  0.63333333333    0.36666666667    0.29309767430
C   0.70000000000    0.30000000000    0.29309767430
 C  0.73333333333    0.36666666667    0.29309767430
C   0.80000000000    0.30000000000    0.29309767430
 C  0.83333333333    0.36666666667    0.29309767430
C   0.90000000000    0.30000000000    0.29309767430
 C  0.93333333333    0.36666666667    0.29309767430
C   0.00000000000    0.40000000000    0.29309767430
 C  0.03333333333    0.46666666667    0.29309767430
C   0.10000000000    0.40000000000    0.29309767430
 C  0.13333333333    0.46666666667    0.29309767430
C   0.20000000000    0.40000000000    0.29309767430
 C  0.23333333333    0.46666666667    0.29309767430
C   0.30000000000    0.40000000000    0.29309767430
 C  0.33333333333    0.46666666667    0.29309767430
C   0.40000000000    0.40000000000    0.29309767430
 C  0.43333333333    0.46666666667    0.29309767430
C   0.50000000000    0.40000000000    0.29309767430
 C  0.53333333333    0.46666666667    0.29309767430
C   0.60000000000    0.40000000000    0.29309767430
 C  0.63333333333    0.46666666667    0.29309767430
C   0.70000000000    0.40000000000    0.29309767430
 C  0.73333333333    0.46666666667    0.29309767430
C   0.80000000000    0.40000000000    0.29309767430
 C  0.83333333333    0.46666666667    0.29309767430
C   0.90000000000    0.40000000000    0.29309767430
 C  0.93333333333    0.46666666667    0.29309767430
C   0.00000000000    0.50000000000    0.29309767430
 C  0.03333333333    0.56666666667    0.29309767430
C   0.10000000000    0.50000000000    0.29309767430
 C  0.13333333333    0.56666666667    0.29309767430
C   0.20000000000    0.50000000000    0.29309767430
 C  0.23333333333    0.56666666667    0.29309767430
C   0.30000000000    0.50000000000    0.29309767430
 C  0.33333333333    0.56666666667    0.29309767430
C   0.40000000000    0.50000000000    0.29309767430
 C  0.43333333333    0.56666666667    0.29309767430
C   0.50000000000    0.50000000000    0.29309767430
 C  0.53333333333    0.56666666667    0.29309767430
C   0.60000000000    0.50000000000    0.29309767430
 C  0.63333333333    0.56666666667    0.29309767430
C   0.70000000000    0.50000000000    0.29309767430
 C  0.73333333333    0.56666666667    0.29309767430
C   0.80000000000    0.50000000000    0.29309767430
 C  0.83333333333    0.56666666667    0.29309767430
C   0.90000000000    0.50000000000    0.29309767430
 C  0.93333333333    0.56666666667    0.29309767430
C   0.00000000000    0.60000000000    0.29309767430
 C  0.03333333333    0.66666666667    0.29309767430
C   0.10000000000    0.60000000000    0.29309767430
 C  0.13333333333    0.66666666667    0.29309767430
C   0.20000000000    0.60000000000    0.29309767430
 C  0.23333333333    0.66666666667    0.29309767430
C   0.30000000000    0.60000000000    0.29309767430
 C  0.33333333333    0.66666666667    0.29309767430
C   0.40000000000    0.60000000000    0.29309767430
 C  0.43333333333    0.66666666667    0.29309767430
C   0.50000000000    0.60000000000    0.29309767430
 C  0.53333333333    0.66666666667    0.29309767430
C   0.60000000000    0.60000000000    0.29309767430
 C  0.63333333333    0.66666666667    0.29309767430
C   0.70000000000    0.60000000000    0.29309767430
 C  0.73333333333    0.66666666667    0.29309767430
C   0.80000000000    0.60000000000    0.29309767430
 C  0.83333333333    0.66666666667    0.29309767430
C   0.90000000000    0.60000000000    0.29309767430
 C  0.93333333333    0.66666666667    0.29309767430
C   0.00000000000    0.70000000000    0.29309767430
 C  0.03333333333    0.76666666667    0.29309767430
C   0.10000000000    0.70000000000    0.29309767430
 C  0.13333333333    0.76666666667    0.29309767430
C   0.20000000000    0.70000000000    0.29309767430
 C  0.23333333333    0.76666666667    0.29309767430
C   0.30000000000    0.70000000000    0.29309767430
 C  0.33333333333    0.76666666667    0.29309767430
C   0.40000000000    0.70000000000    0.29309767430
 C  0.43333333333    0.76666666667    0.29309767430
C   0.50000000000    0.70000000000    0.29309767430
 C  0.53333333333    0.76666666667    0.29309767430
C   0.60000000000    0.70000000000    0.29309767430
 C  0.63333333333    0.76666666667    0.29309767430
C   0.70000000000    0.70000000000    0.29309767430
 C  0.73333333333    0.76666666667    0.29309767430
C   0.80000000000    0.70000000000    0.29309767430
 C  0.83333333333    0.76666666667    0.29309767430
C   0.90000000000    0.70000000000    0.29309767430
 C  0.93333333333    0.76666666667    0.29309767430
C   0.00000000000    0.80000000000    0.29309767430
 C  0.03333333333    0.86666666667    0.29309767430
C   0.10000000000    0.80000000000    0.29309767430
 C  0.13333333333    0.86666666667    0.29309767430
C   0.20000000000    0.80000000000    0.29309767430
 C  0.23333333333    0.86666666667    0.29309767430
C   0.30000000000    0.80000000000    0.29309767430
 C  0.33333333333    0.86666666667    0.29309767430
C   0.40000000000    0.80000000000    0.29309767430
 C  0.43333333333    0.86666666667    0.29309767430
C   0.50000000000    0.80000000000    0.29309767430
 C  0.53333333333    0.86666666667    0.29309767430
C   0.60000000000    0.80000000000    0.29309767430
 C  0.63333333333    0.86666666667    0.29309767430
C   0.70000000000    0.80000000000    0.29309767430
 C  0.73333333333    0.86666666667    0.29309767430
C   0.80000000000    0.80000000000    0.29309767430
 C  0.83333333333    0.86666666667    0.29309767430
C   0.90000000000    0.80000000000    0.29309767430
 C  0.93333333333    0.86666666667    0.29309767430
C   0.00000000000    0.90000000000    0.29309767430
 C  0.03333333333    0.96666666667    0.29309767430
C   0.10000000000    0.90000000000    0.29309767430
 C  0.13333333333    0.96666666667    0.29309767430
C   0.20000000000    0.90000000000    0.29309767430
 C  0.23333333333    0.96666666667    0.29309767430
C   0.30000000000    0.90000000000    0.29309767430
 C  0.33333333333    0.96666666667    0.29309767430
C   0.40000000000    0.90000000000    0.29309767430
 C  0.43333333333    0.96666666667    0.29309767430
C   0.50000000000    0.90000000000    0.29309767430
 C  0.53333333333    0.96666666667    0.29309767430
C   0.60000000000    0.90000000000    0.29309767430
 C  0.63333333333    0.96666666667    0.29309767430
C   0.70000000000    0.90000000000    0.29309767430
 C  0.73333333333    0.96666666667    0.29309767430
C   0.80000000000    0.90000000000    0.29309767430
 C  0.83333333333    0.96666666667    0.29309767430
C   0.90000000000    0.90000000000    0.29309767430
 C  0.93333333333    0.96666666667    0.29309767430
 Ir  0.00000000000   0.00000000000   0.00000000000
 Ir  0.00000000000   0.11111111111   0.00000000000
 Ir  0.00000000000   0.22222222222   0.00000000000
 Ir  0.00000000000   0.33333333333   0.00000000000
 Ir  0.00000000000   0.44444444444   0.00000000000
 Ir  0.00000000000   0.55555555556   0.00000000000
 Ir  0.00000000000   0.66666666667   0.00000000000
 Ir  0.00000000000   0.77777777778   0.00000000000
 Ir  0.00000000000   0.88888888889   0.00000000000
 Ir  0.11111111111   0.00000000000   0.00000000000
 Ir  0.11111111111   0.11111111111   0.00000000000
 Ir  0.11111111111   0.22222222222   0.00000000000
 Ir  0.11111111111   0.33333333333   0.00000000000
 Ir  0.11111111111   0.44444444444   0.00000000000
 Ir  0.11111111111   0.55555555556   0.00000000000
 Ir  0.11111111111   0.66666666667   0.00000000000
 Ir  0.11111111111   0.77777777778   0.00000000000
 Ir  0.11111111111   0.88888888889   0.00000000000
 Ir  0.22222222222   0.00000000000   0.00000000000
 Ir  0.22222222222   0.11111111111   0.00000000000
 Ir  0.22222222222   0.22222222222   0.00000000000
 Ir  0.22222222222   0.33333333333   0.00000000000
 Ir  0.22222222222   0.44444444444   0.00000000000
 Ir  0.22222222222   0.55555555556   0.00000000000
 Ir  0.22222222222   0.66666666667   0.00000000000
 Ir  0.22222222222   0.77777777778   0.00000000000
 Ir  0.22222222222   0.88888888889   0.00000000000
 Ir  0.33333333333   0.00000000000   0.00000000000
 Ir  0.33333333333   0.11111111111   0.00000000000
 Ir  0.33333333333   0.22222222222   0.00000000000
 Ir  0.33333333333   0.33333333333   0.00000000000
 Ir  0.33333333333   0.44444444444   0.00000000000
 Ir  0.33333333333   0.55555555556   0.00000000000
 Ir  0.33333333333   0.66666666667   0.00000000000
 Ir  0.33333333333   0.77777777778   0.00000000000
 Ir  0.33333333333   0.88888888889   0.00000000000
 Ir  0.44444444444   0.00000000000   0.00000000000
 Ir  0.44444444444   0.11111111111   0.00000000000
 Ir  0.44444444444   0.22222222222   0.00000000000
 Ir  0.44444444444   0.33333333333   0.00000000000
 Ir  0.44444444444   0.44444444444   0.00000000000
 Ir  0.44444444444   0.55555555556   0.00000000000
 Ir  0.44444444444   0.66666666667   0.00000000000
 Ir  0.44444444444   0.77777777778   0.00000000000
 Ir  0.44444444444   0.88888888889   0.00000000000
 Ir  0.55555555556   0.00000000000   0.00000000000
 Ir  0.55555555556   0.11111111111   0.00000000000
 Ir  0.55555555556   0.22222222222   0.00000000000
 Ir  0.55555555556   0.33333333333   0.00000000000
 Ir  0.55555555556   0.44444444444   0.00000000000
 Ir  0.55555555556   0.55555555556   0.00000000000
 Ir  0.55555555556   0.66666666667   0.00000000000
 Ir  0.55555555556   0.77777777778   0.00000000000
 Ir  0.55555555556   0.88888888889   0.00000000000
 Ir  0.66666666667   0.00000000000   0.00000000000
 Ir  0.66666666667   0.11111111111   0.00000000000
 Ir  0.66666666667   0.22222222222   0.00000000000
 Ir  0.66666666667   0.33333333333   0.00000000000
 Ir  0.66666666667   0.44444444444   0.00000000000
 Ir  0.66666666667   0.55555555556   0.00000000000
 Ir  0.66666666667   0.66666666667   0.00000000000
 Ir  0.66666666667   0.77777777778   0.00000000000
 Ir  0.66666666667   0.88888888889   0.00000000000
 Ir  0.77777777778   0.00000000000   0.00000000000
 Ir  0.77777777778   0.11111111111   0.00000000000
 Ir  0.77777777778   0.22222222222   0.00000000000
 Ir  0.77777777778   0.33333333333   0.00000000000
 Ir  0.77777777778   0.44444444444   0.00000000000
 Ir  0.77777777778   0.55555555556   0.00000000000
 Ir  0.77777777778   0.66666666667   0.00000000000
 Ir  0.77777777778   0.77777777778   0.00000000000
 Ir  0.77777777778   0.88888888889   0.00000000000
 Ir  0.88888888889   0.00000000000   0.00000000000
 Ir  0.88888888889   0.11111111111   0.00000000000
 Ir  0.88888888889   0.22222222222   0.00000000000
 Ir  0.88888888889   0.33333333333   0.00000000000
 Ir  0.88888888889   0.44444444444   0.00000000000
 Ir  0.88888888889   0.55555555556   0.00000000000
 Ir  0.88888888889   0.66666666667   0.00000000000
 Ir  0.88888888889   0.77777777778   0.00000000000
 Ir  0.88888888889   0.88888888889   0.00000000000
Ir   0.03703703704  -0.03703703704   0.07854140049
Ir   0.03703703704   0.07407407407   0.07854140049
Ir   0.03703703704   0.18518518519   0.07854140049
Ir   0.03703703704   0.29629629630   0.07854140049
Ir   0.03703703704   0.40740740741   0.07854140049
Ir   0.03703703704   0.51851851852   0.07854140049
Ir   0.03703703704   0.62962962963   0.07854140049
Ir   0.03703703704   0.74074074074   0.07854140049
Ir   0.03703703704   0.85185185185   0.07854140049
Ir   0.14814814815  -0.03703703704   0.07854140049
Ir   0.14814814815   0.07407407407   0.07854140049
Ir   0.14814814815   0.18518518519   0.07854140049
Ir   0.14814814815   0.29629629630   0.07854140049
Ir   0.14814814815   0.40740740741   0.07854140049
Ir   0.14814814815   0.51851851852   0.07854140049
Ir   0.14814814815   0.62962962963   0.07854140049
Ir   0.14814814815   0.74074074074   0.07854140049
Ir   0.14814814815   0.85185185185   0.07854140049
Ir   0.25925925926  -0.03703703704   0.07854140049
Ir   0.25925925926   0.07407407407   0.07854140049
Ir   0.25925925926   0.18518518519   0.07854140049
Ir   0.25925925926   0.29629629630   0.07854140049
Ir   0.25925925926   0.40740740741   0.07854140049
Ir   0.25925925926   0.51851851852   0.07854140049
Ir   0.25925925926   0.62962962963   0.07854140049
Ir   0.25925925926   0.74074074074   0.07854140049
Ir   0.25925925926   0.85185185185   0.07854140049
Ir   0.37037037037  -0.03703703704   0.07854140049
Ir   0.37037037037   0.07407407407   0.07854140049
Ir   0.37037037037   0.18518518519   0.07854140049
Ir   0.37037037037   0.29629629630   0.07854140049
Ir   0.37037037037   0.40740740741   0.07854140049
Ir   0.37037037037   0.51851851852   0.07854140049
Ir   0.37037037037   0.62962962963   0.07854140049
Ir   0.37037037037   0.74074074074   0.07854140049
Ir   0.37037037037   0.85185185185   0.07854140049
Ir   0.48148148148  -0.03703703704   0.07854140049
Ir   0.48148148148   0.07407407407   0.07854140049
Ir   0.48148148148   0.18518518519   0.07854140049
Ir   0.48148148148   0.29629629630   0.07854140049
Ir   0.48148148148   0.40740740741   0.07854140049
Ir   0.48148148148   0.51851851852   0.07854140049
Ir   0.48148148148   0.62962962963   0.07854140049
Ir   0.48148148148   0.74074074074   0.07854140049
Ir   0.48148148148   0.85185185185   0.07854140049
Ir   0.59259259259  -0.03703703704   0.07854140049
Ir   0.59259259259   0.07407407407   0.07854140049
Ir   0.59259259259   0.18518518519   0.07854140049
Ir   0.59259259259   0.29629629630   0.07854140049
Ir   0.59259259259   0.40740740741   0.07854140049
Ir   0.59259259259   0.51851851852   0.07854140049
Ir   0.59259259259   0.62962962963   0.07854140049
Ir   0.59259259259   0.74074074074   0.07854140049
Ir   0.59259259259   0.85185185185   0.07854140049
Ir   0.70370370370  -0.03703703704   0.07854140049
Ir   0.70370370370   0.07407407407   0.07854140049
Ir   0.70370370370   0.18518518519   0.07854140049
Ir   0.70370370370   0.29629629630   0.07854140049
Ir   0.70370370370   0.40740740741   0.07854140049
Ir   0.70370370370   0.51851851852   0.07854140049
Ir   0.70370370370   0.62962962963   0.07854140049
Ir   0.70370370370   0.74074074074   0.07854140049
Ir   0.70370370370   0.85185185185   0.07854140049
Ir   0.81481481481  -0.03703703704   0.07854140049
Ir   0.81481481481   0.07407407407   0.07854140049
Ir   0.81481481481   0.18518518519   0.07854140049
Ir   0.81481481481   0.29629629630   0.07854140049
Ir   0.81481481481   0.40740740741   0.07854140049
Ir   0.81481481481   0.51851851852   0.07854140049
Ir   0.81481481481   0.62962962963   0.07854140049
Ir   0.81481481481   0.74074074074   0.07854140049
Ir   0.81481481481   0.85185185185   0.07854140049
Ir   0.92592592593  -0.03703703704   0.07854140049
Ir   0.92592592593   0.07407407407   0.07854140049
Ir   0.92592592593   0.18518518519   0.07854140049
Ir   0.92592592593   0.29629629630   0.07854140049
Ir   0.92592592593   0.40740740741   0.07854140049
Ir   0.92592592593   0.51851851852   0.07854140049
Ir   0.92592592593   0.62962962963   0.07854140049
Ir   0.92592592593   0.74074074074   0.07854140049
Ir   0.92592592593   0.85185185185   0.07854140049
 Ir  0.00000000000  -0.11111111111   0.15570700293   0  0  0
 Ir  0.00000000000   0.00000000000   0.15570700293   0  0  0
 Ir  0.00000000000   0.11111111111   0.15570700293   0  0  0
 Ir  0.00000000000   0.22222222222   0.15570700293   0  0  0
 Ir  0.00000000000   0.33333333333   0.15570700293   0  0  0
 Ir  0.00000000000   0.44444444444   0.15570700293   0  0  0
 Ir  0.00000000000   0.55555555556   0.15570700293   0  0  0
 Ir  0.00000000000   0.66666666667   0.15570700293   0  0  0
 Ir  0.00000000000   0.77777777778   0.15570700293   0  0  0
 Ir  0.11111111111  -0.11111111111   0.15570700293   0  0  0
 Ir  0.11111111111   0.00000000000   0.15570700293   0  0  0
 Ir  0.11111111111   0.11111111111   0.15570700293   0  0  0
 Ir  0.11111111111   0.22222222222   0.15570700293   0  0  0
 Ir  0.11111111111   0.33333333333   0.15570700293   0  0  0
 Ir  0.11111111111   0.44444444444   0.15570700293   0  0  0
 Ir  0.11111111111   0.55555555556   0.15570700293   0  0  0
 Ir  0.11111111111   0.66666666667   0.15570700293   0  0  0
 Ir  0.11111111111   0.77777777778   0.15570700293   0  0  0
 Ir  0.22222222222  -0.11111111111   0.15570700293   0  0  0
 Ir  0.22222222222   0.00000000000   0.15570700293   0  0  0
 Ir  0.22222222222   0.11111111111   0.15570700293   0  0  0
 Ir  0.22222222222   0.22222222222   0.15570700293   0  0  0
 Ir  0.22222222222   0.33333333333   0.15570700293   0  0  0
 Ir  0.22222222222   0.44444444444   0.15570700293   0  0  0
 Ir  0.22222222222   0.55555555556   0.15570700293   0  0  0
 Ir  0.22222222222   0.66666666667   0.15570700293   0  0  0
 Ir  0.22222222222   0.77777777778   0.15570700293   0  0  0
 Ir  0.33333333333  -0.11111111111   0.15570700293   0  0  0
 Ir  0.33333333333   0.00000000000   0.15570700293   0  0  0
 Ir  0.33333333333   0.11111111111   0.15570700293   0  0  0
 Ir  0.33333333333   0.22222222222   0.15570700293   0  0  0
 Ir  0.33333333333   0.33333333333   0.15570700293   0  0  0
 Ir  0.33333333333   0.44444444444   0.15570700293   0  0  0
 Ir  0.33333333333   0.55555555556   0.15570700293   0  0  0
 Ir  0.33333333333   0.66666666667   0.15570700293   0  0  0
 Ir  0.33333333333   0.77777777778   0.15570700293   0  0  0
 Ir  0.44444444444  -0.11111111111   0.15570700293   0  0  0
 Ir  0.44444444444   0.00000000000   0.15570700293   0  0  0
 Ir  0.44444444444   0.11111111111   0.15570700293   0  0  0
 Ir  0.44444444444   0.22222222222   0.15570700293   0  0  0
 Ir  0.44444444444   0.33333333333   0.15570700293   0  0  0
 Ir  0.44444444444   0.44444444444   0.15570700293   0  0  0
 Ir  0.44444444444   0.55555555556   0.15570700293   0  0  0
 Ir  0.44444444444   0.66666666667   0.15570700293   0  0  0
 Ir  0.44444444444   0.77777777778   0.15570700293   0  0  0
 Ir  0.55555555556  -0.11111111111   0.15570700293   0  0  0
 Ir  0.55555555556   0.00000000000   0.15570700293   0  0  0
 Ir  0.55555555556   0.11111111111   0.15570700293   0  0  0
 Ir  0.55555555556   0.22222222222   0.15570700293   0  0  0
 Ir  0.55555555556   0.33333333333   0.15570700293   0  0  0
 Ir  0.55555555556   0.44444444444   0.15570700293   0  0  0
 Ir  0.55555555556   0.55555555556   0.15570700293   0  0  0
 Ir  0.55555555556   0.66666666667   0.15570700293   0  0  0
 Ir  0.55555555556   0.77777777778   0.15570700293   0  0  0
 Ir  0.66666666667  -0.11111111111   0.15570700293   0  0  0
 Ir  0.66666666667   0.00000000000   0.15570700293   0  0  0
 Ir  0.66666666667   0.11111111111   0.15570700293   0  0  0
 Ir  0.66666666667   0.22222222222   0.15570700293   0  0  0
 Ir  0.66666666667   0.33333333333   0.15570700293   0  0  0
 Ir  0.66666666667   0.44444444444   0.15570700293   0  0  0
 Ir  0.66666666667   0.55555555556   0.15570700293   0  0  0
 Ir  0.66666666667   0.66666666667   0.15570700293   0  0  0
 Ir  0.66666666667   0.77777777778   0.15570700293   0  0  0
 Ir  0.77777777778  -0.11111111111   0.15570700293   0  0  0
 Ir  0.77777777778   0.00000000000   0.15570700293   0  0  0
 Ir  0.77777777778   0.11111111111   0.15570700293   0  0  0
 Ir  0.77777777778   0.22222222222   0.15570700293   0  0  0
 Ir  0.77777777778   0.33333333333   0.15570700293   0  0  0
 Ir  0.77777777778   0.44444444444   0.15570700293   0  0  0
 Ir  0.77777777778   0.55555555556   0.15570700293   0  0  0
 Ir  0.77777777778   0.66666666667   0.15570700293   0  0  0
 Ir  0.77777777778   0.77777777778   0.15570700293   0  0  0
 Ir  0.88888888889  -0.11111111111   0.15570700293   0  0  0
 Ir  0.88888888889   0.00000000000   0.15570700293   0  0  0
 Ir  0.88888888889   0.11111111111   0.15570700293   0  0  0
 Ir  0.88888888889   0.22222222222   0.15570700293   0  0  0
 Ir  0.88888888889   0.33333333333   0.15570700293   0  0  0
 Ir  0.88888888889   0.44444444444   0.15570700293   0  0  0
 Ir  0.88888888889   0.55555555556   0.15570700293   0  0  0
 Ir  0.88888888889   0.66666666667   0.15570700293   0  0  0
 Ir  0.88888888889   0.77777777778   0.15570700293   0  0  0
```

## HPC calculations: Installation and first runs with GROMACS, PLUMED, ALPHAFOLD, OPENMM and OPENPATHSAMPLING

This section constitutes a work plan in itself, since none of those codes have been run as yet. We are currently gathering information, particulary regarding input file syntax and generic input files to be used as sample cases. The PDB file ("geometryCIF_PDB_XSF_FASTA" ("k-cl-sol.pdb") repository) corresponding to the structure below (two molecules forming a complex in a KCL (saline) water solution) provides a good sample to be run with most codes, except for AlphaFold (not a protein). The example has been taken from the PAPRIKA distribution [^1]. See below:

![](complexationinwater.png)

Moreover, the list of software would not be complete if we did not include all the packages installed under conda (WSL Ubuntu 22.04 LTS) in a series of virtual environments, on the DSI guest2 laptop:

- "geniushpc": ASE, GPAW, LAMMPS, QE, CP2K, NWCHEM
- "abinitbigdft": ABINIT, BIGDFT and PYTEST;
- "nbsjarvistools": GRIP, JupyterLab, JARVIS-TOOLS, PYIRON, LAMMPS;
- "pkinetrareev": PyEMMA, OPENMM, OPENPATHSAMPLING, GROMACS, PLUMED

The ABINIT installation has one purpose: To see whether version 9.8 has already the new test suite implementation. One final point to be noted is that the presence of generic geometry files for solvated molecules (format PDB, CIF, etc) will be useful for running the ENVIRON solvation module within quantum-espresso [^2].

[^1]: https://github.com/GilsonLabUCSD/pAPRika/blob/master/docs/index.rst
[^2]: http://www.quantum-environ.org/