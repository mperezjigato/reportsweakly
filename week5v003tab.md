# Internship - Week 5 - Summary

## Progress on HPC runs of MD and Ab-Initio codes within the *Development of Supercomputing Training Material for Computational Chemistry* task

LAMMPS (HPC) calculations of the ethanol box (week 4) are used as a testbed to explore the different parallelisation options of the program.
Pure OpenMP as well as pure MPI and hybrid runs are reported under specific submission parameters. Moreover, HPC runs including the Kokkos acceleration package of LAMMPS are considered. A set of job submission scripts have been uploaded in the "computationalchemistry" GitHub repository and total wall time tables are enclosed below:

# Ethanol box - LAMMPS

## Total wall time - Pure OpenMP vs Pure MPI vs sequential calculation

|          | one core | pOMP:36c | pOMP:20c | pMPI:36c | pMPI:72c | pMPI:144c |
|----------|----------|----------|----------|----------|----------|-----------|
|t/min:sec |  0:08:19 |  0:01:06 |  0:01:15 |  0:02:58 | 0:01:22  |  0:01:02  |

## Total wall time - Hybrid calculations

|          | h:18r/2t | h:9r/4t  | h:6r/6t | h:12r/6t | h:18r/4t | h:36r/2t |
|----------|----------|----------|---------|----------|----------|----------|
|t/min:sec | 0:04:31  | 0:04:12  | 0:06:05 | 0:05:22  | 0:04:13  | 0:02:48  |

Likewise, wall time tables are shown below for *atomic simulation environment* (ASE) water box (week 4) relaxation MD runs, under different conditions. On the other hand, total wall time values for large scale (single-point) ASE calculations for a slab [^1] of a AgAuRuIrPt high-entropy alloy are reported. Embedded atom method force fields have been used, and the calculation has been carried out using LAMMPS in a completely separate run.

Furthermore, and considering ASE as a management system running a collection of both MD and first-principles quantum mechanical codes (as opposed to the cases discussed above), the software LAMMPS, GPAW and QE has been (conda) installed locally (Ubuntu laptop) and alongside ASE. A separate virtual environment for both ABINIT and ASE has been set up as well. Total wall time tables of the single-point calculations for the above alloy slab are shown below.

## Progress on `BASH shell` scripting and templating

M4, renderest and full BASH shell templating software is being discussed in relation to the problem posed in previous reports and that allows the manipulation of LAMMPS input files.

### Templating

The following script transforms a (LAMMPS input file) template with a BASH shell script, ie without any external software:
```bash
#!/usr/bin/env bash

workh='tobe.txt'
cp template.txt $workh
outputf='tb.txt'; touch $outputf

tsys='ethanol'; delimiter=( "molecule" "create_atoms" ) 
printf -v datf "data.%s" $tsys
xFORM="C6H6"; xDATAP="$(pwd)/$datf"

printf -v lmole "grep %s %s" ${delimiter[0]} $workh
printf -v lcrat "grep %s %s" ${delimiter[1]} $workh

lm=`$lmole`
lc=`$lcrat`

while read input_text
do
  case $input_text in
        ${lm} | ${lc}) 
		eval echo $input_text
		;;
        *)              
		echo $input_text
                ;;
   esac
done < $workh > $outputf
```
In this case, the search/replace operation is carried out linewise, via a nested while-case loop. Certainly the `sed` option is still preferable: Carrying out several search/replace operations in a single step as opposed to generating an intermediate template.