
# Internship - Week 9 - Summary

At this point, and considering that we do not count with CPU time at Genius, it seems as a good idea trying to gather input and output files corresponding to all the calculations carried out so far with the Ab-Initio/MD software. The list coincides with the tar files uploaded to the "computationalchemistry" repository, although without large pseudopotential (nor wavefunction) files. 

 On the other hand, installation work (on DSI guest2 laptop), is summarised below, as well as work on the first runs with the software below:

 - GROMACS
 - PLUMED
 - AlphaFold (not implemented by the intern, simply run on Genius, if possible)
 - OPENMM
 - OPENPLATHSAMPLING
 
 The last two codes are not separately implemented. Rather, they take part in the distribution of other software. Firstly, PLUMED is a free-energy calculation code that makes use of OPENPATHSAMPLING. Secondly, AlpaFold is a protein structure determination sofware that uses OPENMM as force-engine, and therefore, OpenMM is part of its distribution.

## HPC calculations: Inputs and outputs - Development of supercomputing training material regarding Ab-Initio/MD calculations

Following up the same order appearing in the computationalchemistry repository (README.md), input and output files are shown below

### ABINIT

 - Sequential ABINIT calculation of the as provided example "tbasepar_1", located within abinit-test/tutorial ("Lead crystal - parallelism over k-points"): 
   *sequentialABINITtbasepar_1.tar.gz* (Genius)
 - Sequential ABINIT calculation of the as provided example "tbasepar_2", located within abinit-test/tutorial ("FCC 4-atom A1-phonon deformed ferromagnetic Fe 
   crystal" - parallelism over spin): *sequentialABINITtbasepar_2.tar.gz* (Genius)
 - **MPI** ABINIT calculation of the as provided example "tdfpt", located within abinit-test/tutoparal ("DFPT phonon calculation for Aluminium crystal - step 1: 
   ground-state electronic structure calcuation - step 2: phonon calculation": *mpiABINITtdfpt1and2.tar.gz* (DSI laptop guest2; no batch job script used)

### ASE

 - Sequential ASE convex hull determination of a CuxPt1-x(111) surface slab by means of a genetic algorithm calculation (global optimization): 
   *sequentialASEgenalgCuPt111slab.tar.gz* (Genius)
 - ASE water-box equilibration calculation as a strong-scaling **MPI** test: *mpiASE_STRONGSCALING_waterboxequi.tar.gz* (Genius)

### LAMMPS

 - **MPI** LAMMPS calculation of c-HfO2: *mpiLAMMPScHfO2.tar.gz* (Genius)
 - LAMMPS ethanol-box equlibration as a strong-scaling **OpenMP** test: *openmpLAMMPSstrongscalingETHANOLBOX.tar.gz* (Genius)
 - LAMMPS ethanol-box equilibration as a *KOKKOS* parallelism *hybrid (36/2)* calculation case: KOKKOSPARALLhybrid36and2ethanol.tar.gz (Genius)
 - Sequential LAMMPS Monte-Carlo relaxation of a two-dimensional deformed hexaganal lattice: *sequentialLAMMPSmc.tar.gz* (Genius)
 - **MPI** LAMMPS calculation of quartz amorphisation via melting/quenching temperature ramps, as an **MPI** strong-scaling case: *mpiLAMMPSstrongscalsilicaamorph.tar.gz* (Genius)

### QUANTUM-ESPRESSO

## HPC calculations: Installation and first runs with GROMACS, PLUMED, ALPHAFOLD, OPENMM and OPENPATHSAMPLING

This section constitutes a work plan in itself, since none of the codes have been run as yet. We are currently gathering information, particulary regarding input file syntax and example input files. 
