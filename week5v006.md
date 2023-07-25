# Internship - Week 5 - Summary

## Progress on HPC runs of MD and Ab-Initio codes within the *Development of Supercomputing Training Material for Computational Chemistry* task

LAMMPS (HPC) calculations of the ethanol box (week 4) are used as a testbed to explore the different parallelisation options of the program.
Pure OpenMP as well as pure MPI and hybrid runs are reported under specific submission parameters. Moreover, HPC runs including the Kokkos acceleration package of LAMMPS are considered. A set of job submission scripts have been uploaded in the "computationalchemistry" GitHub repository and total wall time tables are enclosed below:

### Ethanol box - LAMMPS (total wall time)

|          |      pure OpenMP     |           pure MPI            | one core |
|----------|----------------------|-------------------------------|----------|
|t/min:sec | <table>  <thead>  <tr>  <th>pOMP:36c</th>  <th>pOMP:20c</th>  </tr>  </thead>  <tbody>  <tr>  <td>0:01:06</td>  <td>0:01:15</td>  </tr>  </tbody>  </table>      | <table>  <thead>  <tr>  <th>pMPI:36c</th>  <th>pMPI:72c</th>  <th>pMPI:144c</th>  </tr>  </thead>  <tbody>  <tr>  <td>0:02:58</td>  <td>0:01:22</td>  <td>0:01:02</td>  </tr>  </tbody>  </table>      | 0:08:19 |

|          |                            hybrid                              |
|----------|----------------------------------------------------------------|
|t/min:sec | <table>  <thead>  <tr>  <th>h:18r/2t</th>  <th>h:9r/4t</th>  <th>h:6r/6t</th>  <th>h:12r/6t</th>  <th>h:18r/4t</th>  <th>h:36r/2t</th>  </tr>  </thead>  <tbody>  <tr>  <td>0:04:31</th>  <th>0:04:12</th>  <th>0:06:05</th>  <th>0:05:22</th>  <th>0:04:13</th>  <th>0:02:48</td>  </tr>  </tbody>  </table>      |

|          |               Kokkos                     |
|----------|------------------------------------------|
|t/min:sec |  <table>  <thead>  <tr>  <th>pMPI:36c</th>  <th>pMPI:20c</th>  <th>pMPI:72c</th>  <th>h:36r/2t</th>  </tr>  </thead>  <tbody>  <tr>  <td>0:03:06 </th>  <th>0:03:37</th>  <th>0:01:38</th>  <th>0:03:48</td>  </tr>  </tbody>  </table>      |                                     |

Likewise, wall time tables are shown below for both silica LAMMPS amorphisation (just the step to generate the disordered geometry) [^1] and *atomic simulation environment* (ASE) water box (week 4) relaxation MD runs, under different conditions. Moreover, geometries (in the form of a CIF file) are being gathered for already known metal solid solutions [^2], particularly those constituting equiatomic high-entropy alloys. My intention is to carry out large scale materials optimisation calculations with both LAMMPS and ASE, simply using embedded-atom method calculations.

Furthermore, and considering ASE as a management system running a collection of both MD and first-principles quantum mechanical codes (as opposed to the cases consideres so far), the software LAMMPS, GPAW and QE hasve all been (conda) installed locally (Ubuntu laptop) and alongside ASE. A separate virtual environment for both ABINIT and ASE has been set up as well. Total wall time tables of the single-point calculations for the above alloy slab are shown below.

## Templating

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