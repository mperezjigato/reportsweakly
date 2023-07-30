
# Internship - Week 6 - Summary

## Progress on HPC runs of MD and Ab-Initio codes within the *Development of Supercomputing Training Material for Computational Chemistry* task

Simple materials science problems provided within the examples directory of the LAMMPS distribution have been run (the force-field type is given in brackets):
 - Amorphous Silica geometry generation (Vashishta).
 - Cu2O elasticity (COMB).
 - OH adsorbed on Cu2O(110) deposited graphene monolayer (COMB).
 - mHfO2 minimisation (COMB).

Wall time tables for the Silica amorphisation case are given at the bottom of the document. No data are reported for the other three cases.

On the other hand, the as provided LAMMPS benchmarks have been run:





Three types of benchmark have been tried, ie sequential, fix-scale (MPI) and scaled-size (MPI) calculations. Data are reported only for the scaled-size calculations [^1] (see below). 

# Amorphous silica geometry generation (Vashishta)

## OpenMP

|     |  4  |  8  |  12  |  16  |  20  |  24  |  28  |  32  |  36  |sequent|
|-----|-----|-----|------|------|------|------|------|------|------|-------|
|t/sec|19748|11509| 8834 | 7308 | 6469 | 6272 | 7875 | 7058 | 6070 |       |

## MPI

|     |  4  |  8  | 12 | 16 | 20 | 24 | 28 | 32 | 36 | 72 | 108| 144|
|-----|-----|-----|----|----|----|----|----|----|----|----|----|----|
|t/sec|20610|11920|8817|6842|5414|5047|4777|4449|4225|2965|2645|2607|

# Scaled-size benchmarks






[^1]: GJB organised a training session on July 25th, 2023 regarding parallel computing. He covered the OpenMP, the MPI and the hybrid (OpenMP/MPI) parallel modes, as well as the concepts of lattency, bandwidth, strong and weak scaling, speedup and parallel efficiency. Within the weak scaling realm, it is key to investigate  size scaling. The reported set of benchmarks goes along those lines, as opposed to the fixed-size problem (by varying the number of processors, ie strong scaling). 
