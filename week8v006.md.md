# Internship - Week 8 - Summary

## Computational Chemistry

Work carried out with Quantum Espresso [^1] is described next, and the tables depticting our progress with the Computational Chemistry packages are updated accordingly (repository "computationalchemistry"). A set of independent runs has been carried out both sequentially and under the MPI parallel mode. On the other hand, the as-provided test-suite is currently being tried.

Incidentally, it must be noted that all our CPU time on GENIUS has been consumed, and a new application for five million extra credits submitted [^2].

### Quantum Espresso calculations

Two as-provided independent systems as well as a system with input files built from scratch, have been run: 

 - GaAs PWSCF calculation under PAW (as provided),
 - SiH4 (molecule in a box) PWSCF calculation under norm-conserving pseudopotentials (as provided),
 - *Immm* Ag2PdO2 PWSCF calculation under ultrasoft pseudopotentials: The materialscloud.org website (Quantum ESPRESSO input generator and structure visualizer) lets us transform the corresponding CIF file (downloaded from the topological materials database) on QE inputs. The full set of input and output files from this quantum espresso run are uploaded to the "computationalchemistry" repository, as well as to the "geometryfilesCIF_PDB_XSF_FASTA" repository.

For the sake of completeness, the QE files corresponding to Ag2PdO2 are shown below:

 - job submission script
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
 - main input file
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
### Specific CIF files and geometry visualisation 

The Ag(16 10 9) surface structure [^3] is depicted below.

 ![](Ag16109.png)
 
On the other hand, the crystal structure for the above topological material Ag2PdO2 is depicted next:

 ![](Ag2PdO2.png)

Finally, the modulated structure of the plagioclase feldspar CaxNa1-xAl1+xSi3-xO8 is shown below [^4]. This represents the textbook example of a solid solution.

 ![](plagioclasefeldspar.png)

[^1]: www.quantum-espresso.org
[^2]: https://admin.kuleuven.be/icts/onderzoek/hpc/extra-project-credits
[^3]: Both CIF file and atomic structure image of Ag(16 10 9) have been provided by Dr Stephen J Jenkins.
[^4]: It has been downloaded from The Bilbao Incommensurate Crystal Structure Database at https://www.cryst.ehu.es/
