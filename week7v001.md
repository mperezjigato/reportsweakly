
# Internship - Week 7 - Summary

LAMMPS related work has been carried out, although some information about ALPHAFOLD has been gathered by the intern during this time period as well.

## Solving the LAMMPS parallel calculation problem in the hybrid mode

Since I identified last week what I thought could be a problem when running LAMMPS in the hybrid mode (OpenMP/MPI), comprehensive discussions were held with GJB. The original indication of a problem consisted in a meaningful difference in wall time for specific LAMMPS calculations (official distribution example directory) in the hybrid mode when the parallel calculation is carried out under the `mpirun lmp -sf omp -pk omp ${OMP_NUM_THREADS}` options, and in comparison with the case in which the specific options have been removed: `mpirun lmp`.

During this week, GJB has provided an extensive description of the well known "racing problem", which is typical of OpenMP. He insisted in this possibility as the most likely source of the above wall time discrepancies found for hybrid mode calculations (see week 6 report). A direct inspection of the outputs has let us see that different C++ routines are called/run in the hybrid mode (when including all the options described in the documentation), as opposed to not using such options. The LAMMPS manual already explains that an explicit subroutine (omp 0) that sets up explicit variables with an extension "omp" in their name is called whenever the `-sf omp` option is utilised, in common with pure OpenMP calculations. On the other hand, the `-pk` option calls an explicit LAMMPS package. It is fair to say that, in principle, running parallel LAMMPS jobs in the hybrid mode, should strictly follow the LAMMPS manual, and, therefore, the results obtained in the hybrid mode, without the appropiate option syntax, must simply be discarded as wrongful!

## Running an MPI weak scaling problem: OH adsorbed on graphene deposited on Cu2O(110)

This is one of the LAMMPS provided examples, which has two-dimensional periodicity. The goal is to carry out weak scaling runs similarly to what I had described in the report of week 6 (LAMMPS benchmarks), but this time extending only the surface dimensions: The original 689-atom geometry is extended according to the 4x4, 8x8, 12x12, ... rule, system size (number of atoms) increasing exactly according to the number of MPI cores in the parallel calculation (16, 64, 144, ...).

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

 
