
# Internship - Week 7 - Summary

LAMMPS related work has been carried out, although some information about ALPHAFOLD has been gathered by the intern during this time period as well.

## Solving the LAMMPS parallel calculation problem in the hybrid mode

Since I identified last week what I thought could be a problem when running LAMMPS in the hybrid mode (OpenMP/MPI), comprehensive discussions were held with GJB. The original indication of a problem consisted in a meaningful difference in wall time for specific LAMMPS calculations (official distribution example directory) in the hybrid mode when the parallel calculation is carried out under the `mpirun lmp -sf omp -pk omp ${OMP_NUM_THREADS}` options, and in comparison with the case in which the specific options have been removed: `mpirun lmp`.

During this week, GJB has provided an extensive description of the well known "racing problem", which is typical of OpenMP. He insisted in this possibility as the most likely source of the above wall time discrepancies found for hybrid mode calculations (see week 6 report). A direct inspection of the outputs has let us see that different C++ routines are called/run in the hybrid mode (when including all the options described in the documentation), as opposed to not using such options. The LAMMPS manual already explains that an explicit subroutine (omp 0) that sets up explicit variables with an extension "omp" in their name is called whenever the `-sf omp` option is utilised, in common with pure OpenMP calculations. On the other hand, the `-pk` option calls an explicit LAMMPS package. It is fair to say that, in principle, running parallel LAMMPS jobs in the hybrid mode, should strictly follow the LAMMPS manual, and, therefore, the results obtained in the hybrid mode, without the appropiate option syntax, must simply be discarded as wrongful!

Still, we decided to seek an explanation. Following GJB's suggestions on the nature of the "OpenMP racing problem", which would involve the appearance of uncertainty in the absolute values of the computed properties as well as on the total wall time, we made the decision of comparing for one specific system ('sic') the output of:
 1. Sequential calculation.
 1. MPI calculation.
 1. OpenMP calculation.
 1. Hybrid calculation with all the correct manual options (-sf omp/-pk omp N).
 1. Hybrid calculation with no options included.
Key to this exercise, provided our hypothesis had been correct, was to include some statistics in the data. Therefore, we computed ten times each of the above calculations. On the other hand, we decided to pick up a specific physical property well known for the dependence of its absolute value on the convergence criteria, ie a second energy derivative, phonons/elastic constants. For that reason, the system we have investigated for this purpose is the Cu2O elestic constants problem, which is located within the examples directory of the LAMMPS distribution ("COMB").

To our surprise, the outcome indicated that the "OpenMP racing problem" had to be discarded as the source of the problem: No statistical variation of the absolute values of the computed elastic constants has been identified at all, all ten calculations calculations for each of the modes produced the same values. Still, we seem to have partly isolated the problem, since, after comparing the elastic constants for all five calculation modes, we find two separate groups of values:

 - Sequential:
 ```
   Elastic constants (GPa):
C11 =    2937.501
C12 =    1552.401
C13 =    1462.783
C14 =     -28.321
C15 =     -34.381
C16 =      -2.030
C33 =    2980.169
C44 =    1509.223
C66 =    1485.259
  B =    1979.012
  G =    1192.555 
 ```
 - MPI:
 ```
 Elastic constants (GPa):
C11 =    2937.501
C12 =    1552.401
C13 =    1462.783
C14 =     -28.321
C15 =     -34.381
C16 =      -2.030
C33 =    2980.169
C44 =    1509.223
C66 =    1485.259
  B =    1979.012
  G =    1192.555
 ```
 - OpenMP
 ```
 Elastic constants (GPa):
C11 =    1483.422
C12 =     805.400
C13 =     767.869
C14 =     -13.036
C15 =     -14.736
C16 =      -1.144
C33 =    1501.464
C44 =     740.086
C66 =     730.685
  B =    1016.731
  G =     583.983
 ```
 - Hybrid under `mpirun lmp -sf omp -pk omp ${OMP_NUM_THREADS}`
 ``` 
  Elastic constants (GPa):
C11 =    1483.422
C12 =     805.400
C13 =     767.869
C14 =     -13.036
C15 =     -14.736
C16 =      -1.144
C33 =    1501.464
C44 =     740.086
C66 =     730.685
  B =    1016.731
  G =     583.983
 ```
- Hybrid with no options: `mpirun lmp`
 ```
 Elastic constants (GPa):
C11 =    2937.501
C12 =    1552.401
C13 =    1462.783
C14 =     -28.321
C15 =     -34.381
C16 =      -2.030
C33 =    2980.169
C44 =    1509.223
C66 =    1485.259
  B =    1979.012
  G =    1192.555
 ```
## Running an MPI weak scaling problem: OH adsorbed on graphene deposited on Cu2O(110)

This is one of the LAMMPS provided examples, which has two-dimensional periodicity. The goal is to carry out weak scaling runs similarly to what I had described in the report of week 6 (LAMMPS benchmarks), but this time extending only the surface dimensions: The original 689-atom geometry is extended according to the 4x4, 8x8, 12x12, ... rule, system size (number of atoms) increasing exactly according to the number of MPI cores in the parallel calculation (16, 64, 144, ...). Results follow:

|number MPI cores |     1     |    16     |    64     |   144     |    400    |
|---------------- |-----------|-----------|-----------|-----------|-----------|
|   atom number   |    682    |   10912   |   43648   |   98208   |   272800  |
|surface extension|     -     |   (4x4)   |   (8x8)   |  (12x12)  |  (20x20)  |
|      time/s     |   1284    |   1288    |   1746    |   1781    |    1654   |

## Information gathering and tool testing

### Geometry tools

Both "CIF2CELL" and "atomsk" have been installed and are currently being tested for CIF/XSF/PDB geometry generation, as well as for producing geometries appropriate to the LAMMPS syntax. At the moment, simple geometries and supercells are being built. The latter is very helpful as starting point for point defect geometry generation, for which a specific tool associated to ASE is to be installed and tested shortly. 

Moreover, the inverse process of converting LAMMPS geometry files provided within the examples/benchmark directories of the distribution, to CIF/PDB/XSF geometries for direct visualisation, is being tried in order to identify some of the provided complex systems (the rhodopsin protein example for instance).

PYMOL for visualising PDB files, as well as XCRYSDEN for XSF files have been installed. On the other hand, the crystal toolkit within the materials project website (materialsproject.org) is used to visualise CIF files. Regarding the latter, it must be noted that some crystal databases contain their own visualisation facility.

### Geometry databases and specific geometry file acquisition

A list of specific websites with structural databases (and CIF files) for crystalline solids is given next:
 1. The "materials explorer" within materialsproject.org.
 1. The "Crystallography Open Database" (COD): crystallography.net.
 1. BILBAO CRYSTALLOGRAPHY SERVER: https://www.cryst.ehu.es/.
 1. The "Materials Cloud three-dimensional crystals database (MC3D)"
 within materialscloud.org.
 1. nomad-lab.eu

Although no specific structural databases have been located for amorphous solids/glasses, the website of the SIMONS COLLABORATION FOR CRACKING THE GLASS PROBLEM (https://scglass.uchicago.edu/) provides an interesting entry point to the subject. A useful aspect of this source is given by sets of LAMMPS files for specific calculations within glass physics.

### AlphaFold

### Thermodynamic integration

### Other LAMMPS packages: phonons, MC

 
