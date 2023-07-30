
# Internship - Week 6 - Summary

## Progress on HPC runs of MD and Ab-Initio codes within the *Development of Supercomputing Training Material for Computational Chemistry* task

Simple materials science problems provided within the examples directory of the LAMMPS distribution have been computed under different submission conditions (the force field type is mentioned in the parenthesis):
 
 - Amorphous Silica geometry generation (Vashishta)
 - Cu2O elasticity (COMB)
 - OH adsorbed on Cu2O(110) deposited graphene monolayer (COMB)
 - mHfO2 minimisation (COMB)

Wall time tables for the respective system are given below.


# Amorphous silica geometry generation (Vashishta)

## OpenMP

|     |  4  |  8  |  12  |  16  |  20  |  24  |  28  |  32  |  36  |sequent|
|-----|-----|-----|------|------|------|------|------|------|------|-------|
|t/sec|19748|11509| 8834 | 7308 | 6469 | 6272 | 7875 | 7058 | 6070 |       |

## MPI

|     |  4  |  8  | 12 | 16 | 20 | 24 | 28 | 32 | 36 | 72 | 108| 144|
|-----|-----|-----|----|----|----|----|----|----|----|----|----|----|
|t/sec|20610|11920|8817|6842|5414|5047|4777|4449|4225|2965|2645|2607|

