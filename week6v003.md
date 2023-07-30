
# Internship - Week 6 - Summary

## Progress on HPC runs of MD and Ab-Initio codes within the *Development of Supercomputing Training Material for Computational Chemistry* task

Following previous work on the Computational Chemistry task, I have decided not to report any hybrid calculations data for the time being. Once the consistency problem we identified last week for hybrid parallel calculations is fully resolved, calculations using that mode will resume. Therefore, only OpenMP and MPI data are discussed now.

On the one hand, simple materials science problems provided within the examples directory of the LAMMPS distribution have been run (the force-field type is given in brackets). Wall time tables for the first case are given at the bottom of the document. No data are reported for the other three cases (only the first case takes long enough as to record wall times). The systems that have been run include:
 - Amorphous Silica geometry generation (Vashishta).
 - Cu2O elasticity (COMB).
 - OH adsorbed on Cu2O(110) deposited graphene monolayer (COMB).
 - mHfO2 minimisation (COMB).

On the other, the five as provided LAMMPS benchmarks have been run: 1000 time-steps for 32000 atom systems (starting geometry). Namely,
 1. NVE integration of a Lennard-Jones fluid.
 2. NVE integration of a bead-spring polymer melt under the FENE approximation.
 3. NVE integration of fcc Cu under the EAM potential.
 4. NVE integration of granular Chutte flow.
 5. NPT integration of the rhodopsin protein in solvated lipid bilayer under the CHARMM force-field and the PPPM Coulomb solution method.
Three types of calculation have been tried for each, ie sequential, fix-scale (MPI) and scaled-size (MPI) calculations. Data are reported only for the MPI scaled-size calculations [^1] (see below). 

Moreover, the `pyiron` HTCMS software has been locally implemented by the intern on the guest2 laptop under Ubuntu 22.04 LTS, in order to shortly start generating geometries for single-phase crystalline solid solutions as required for High-Entropy alloys [^2] calculations with LAMMPS. Likewise, one of the traditional problems
of equilibrium molecular dynamics, the calculation of heat transport, seems to have been resolved by Stefano Baroni and co-workers via a new theory based on the cepstral analysis of current time series derived from equilibrium molecular dynamics. We are in the process of installing the SporTran [^3] software for heat transport calculations in combination with LAMMPS.  

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
[^2]: www.pyiron.org ; a tool for geometry generation of crystalline solid solutions is contained within pyiron.
[^3]: https://github.com/sissaschool/sportran.
