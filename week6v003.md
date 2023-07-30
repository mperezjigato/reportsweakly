
# Internship - Week 6 - Summary

## Progress on HPC runs of MD and Ab-Initio codes within the *Development of Supercomputing Training Material for Computational Chemistry* task

Simple materials science problems provided within the examples directory of the LAMMPS distribution have been run (the force-field type is given in brackets).
Wall time tables for the first case are given at the bottom of the document. No data are reported for the other three cases (only the first case takes long enough as to record wall times). The systems that have been run include:
 - Amorphous Silica geometry generation (Vashishta).
 - Cu2O elasticity (COMB).
 - OH adsorbed on Cu2O(110) deposited graphene monolayer (COMB).
 - mHfO2 minimisation (COMB).
On the other hand, the five as provided LAMMPS benchmarks have been run: 1000 time-steps for 32000 atom systems (starting geometry). Namely,
 1. NVE integration of a Lennard-Jones fluid.
 2. NVE integration of a bead-spring polymer melt under the FENE approximation.
 3. NVE integration of fcc Cu under the EAM potential.
 4. NVE integration of granular Chutte flow.
 5. NPT integration of the rhodopsin protein in solvated lipid bilayer under the CHARMM force-field and the PPPM Coulomb solution method.
Three types of benchmark have been tried, ie sequential, fix-scale (MPI) and scaled-size (MPI) calculations. Data are reported only for the MPI scaled-size calculations [^1] (see below). 

### Amorphous silica geometry generation (Vashishta)

#### OpenMP

|     |  4  |  8  |  12  |  16  |  20  |  24  |  28  |  32  |  36  |sequent|
|-----|-----|-----|------|------|------|------|------|------|------|-------|
|t/sec|19748|11509| 8834 | 7308 | 6469 | 6272 | 7875 | 7058 | 6070 |       |

#### MPI

|     |  4  |  8  | 12 | 16 | 20 | 24 | 28 | 32 | 36 | 72 | 108| 144|
|-----|-----|-----|----|----|----|----|----|----|----|----|----|----|
|t/sec|20610|11920|8817|6842|5414|5047|4777|4449|4225|2965|2645|2607|

### MPI scaled-size benchmarks






[^1]: On July 25th, 2023, GJB organised a training session regarding parallel computing. He then discussed the OpenMP, the MPI and the hybrid (OpenMP/MPI) parallel modes, as well as the concepts of lattency, bandwidth, strong and weak scaling, speedup and parallel efficiency. Within the weak scaling realm, it is key to investigate size scaling. The herein reported set of benchmarks goes along those lines, as opposed to the fixed-size problem (varying the number of processors whilst keeping the same size, ie strong scaling). In a previous report I had shown OpenMP and MPI scaling results for two separate systems, although they corresponded only to the strong scaling concept (weak scaling had been ignored at that point by the intern).
