# Internship - Week 5 - Summary

Two separate systems are used as testbeds for all available parallelisation options of LAMMPS: The ethanol box calculation (week 4) and a Lennard-Jones fluid calculation [^1].

Pure OpenMP as well as pure MPI and hybrid runs are reported under specific submission parameters. Moreover, HPC runs including the Kokkos acceleration package of LAMMPS are considered. A set of job submission scripts have been uploaded in the "computationalchemistry" GitHub repository and total wall time tables are enclosed below for both cases:

## Ethanol box - LAMMPS (total wall time)

|          |      pure OpenMP     |           pure MPI            | one core |
|----------|----------------------|-------------------------------|----------|
|t/min:sec | <table>  <thead>  <tr>  <th>pOMP:36c</th>  <th>pOMP:20c</th>  </tr>  </thead>  <tbody>  <tr>  <td>0:01:06</td>  <td>0:01:15</td>  </tr>  </tbody>  </table>      | <table>  <thead>  <tr>  <th>pMPI:36c</th>  <th>pMPI:72c</th>  <th>pMPI:144c</th>  </tr>  </thead>  <tbody>  <tr>  <td>0:02:58</td>  <td>0:01:22</td>  <td>0:01:02</td>  </tr>  </tbody>  </table>      | 0:08:19 |

|          |                            hybrid                              |
|----------|----------------------------------------------------------------|
|t/min:sec | <table>  <thead>  <tr>  <th>h:18r/2t</th>  <th>h:9r/4t</th>  <th>h:6r/6t</th>  <th>h:12r/6t</th>  <th>h:18r/4t</th>  <th>h:36r/2t</th>  </tr>  </thead>  <tbody>  <tr>  <td>0:04:31</th>  <th>0:04:12</th>  <th>0:06:05</th>  <th>0:05:22</th>  <th>0:04:13</th>  <th>0:02:48</td>  </tr>  </tbody>  </table>      |

|          |               Kokkos                     |
|----------|------------------------------------------|
|t/min:sec |  <table>  <thead>  <tr>  <th>pMPI:36c</th>  <th>pMPI:20c</th>  <th>pMPI:72c</th>  <th>h:36r/2t</th>  </tr>  </thead>  <tbody>  <tr>  <td>0:03:06 </th>  <th>0:03:37</th>  <th>0:01:38</th>  <th>0:03:48</td>  </tr>  </tbody>  </table>      |                                     |

## Lennard-Jones fluid - LAMMPS (total wall time)

|          |      pure OpenMP     |           pure MPI            | one core |
|----------|----------------------|-------------------------------|----------|
|t/min:sec | <table>  <thead>  <tr>  <th>pOMP:36c</th>  <th>pOMP:16c</th>  </tr>  </thead>  <tbody>  <tr>  <td>0:03:45</td>  <td>0:02:35</td>  </tr>  </tbody>  </table>      | <table>  <thead>  <tr>  <th>pMPI:36c</th>  <th>pMPI:20c</th>  <th>pMPI108</th>  </tr>  </thead>  <tbody>  <tr>  <td>0:01:28</td>  <td>0:01:22</td>  <td>0:00:53</td>  </tr>  </tbody>  </table>      | 0:10:20 |

|           |       hybrid (`mpirun` with "-sf" and "-pk")       |
|-----------|----------------------------------------------------|
| t/min:sec | <table>  <thead>  <tr>  <th>h4r4t</th>  <th>h5r4t</th>  <th>h6r6t</th>  <th>h9r4t</th>  <th>h12r3t</th>  <th>h18r2t</th>  <th>h12r6t</th> <th>h18r4t</th>  <th>h36r2t</th>  <th>h12r12t</th>  <th>h72r2t</th>  </tr>  </thead>  <tbody>  <tr>  <td>0:03:59</td>  <td>0:03:59</td>  <td>0:03:28</td>  <td>0:02:56</td>  <td>0:03:21</td>  <td>0:02:42</td>  <td>0:02:34</td>  <td>0:02:16</td> <td>0:01:29</td>  <td>0:04:02</td>  <td>0:01:34</td>  </tr>  </tbody>  </table>      |

|           |               hybrid (bare `mpirun`)               |
|-----------|----------------------------------------------------|
| t/min:sec | <table>  <thead>  <tr>  <th>h4r4t</th>  <th>h5r4t</th>  <th>h6r6t</th>  <th>h9r4t</th>  <th>h12r3t</th>  <th>h18r2t</th>  <th>h12r6t</th> <th>h18r4t</th>  <th>h36r2t</th>  <th>h12r12t</th>  <th>h72r2t</th>  </tr>  </thead>  <tbody>  <tr>  <td>0:03:40</td>  <td>0:03:58</td>  <td>0:02:23</td>  <td>0:02:10</td>  <td>0:02:03</td>  <td>0:02:23</td>  <td>0:01:49</td>  <td>0:01:47</td> <td>0:01:18</td>  <td>0:01:53</td>  <td>0:01:11</td>  </tr>  </tbody>  </table>      |

It must be noted that two separate wall time entries appear for the LJ hybrid calculations. One of them uses all the `mpirun` options appearing in the manual as appropriate to hybrid (MPI+OpenMP) calculations - including both `-sk omp` and the `-pk omp ${OMP_NUM_THREADS}` -. The other does not utilise them, and as can be observed, it has an affect on performance. After discussions with GJB, we are asserting the accuracy of both submission modes, and comparing their outputs to existing reference output files. In case no difference in the actual outputs there existed, it must be communicated to the authors of LAMMPS. 

In order to move towards more involved materials science applications, information is currently being gathered regarding metal solid solutions equiatomic high-entropy alloys [^2]. Geometry files in either CIF form or in simple coordinate "xyz" form are being sought. The intention is to get started with geometry optimisations on large scale systems using MD calculations (LAMMPS and/or ASE) in combination with empirical/semiempirical potentials, before considering quantum effects (ML potentials of SNAP type; Tight-binding; fully-fledged density-functional theory). Certainly, it must be noted that the embedded-atom method leads to the wrong phonons on common metals, although geometries seem to be correct. Most materials properties, including free energy calculations of defects and of dynamical stabilisation effects (metastable phases at finite temperature) require accurate anharmonic effects, therefore accurate phonon predictions become unavoidable at some point. Furthermore, existing alloy structures seem appropriate for testing purposes, although tools based on the special quasiramdom structure algorithm as well as on the cluster expansion method are to be eventually explored in order to identify unknown phases.

Furthermore, local Ubuntu implementations of several codes that fall within the ASE umbrella are being carried out on the DSI laptop guest2, in order to explore the possibility of identifying the "Genius" test suite problems that appear on both ABINIT and LAMMPS (they seem to have in common a similar Python problem, hopefully disappearing when running the codes from within ASE). The ABINIT test suite is installed on the share directory next to "bin", whilst the LAMMPS test suite is automatically generated on my .`dot-local` as soon as I start running LAMMPS.

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

[^1]: ("2-post-analysis"): A Lennard-Jones fluid simple MD calculation. See the material provided with the tutorial "LAMMPS Course for Intermediate Users" https://epcced.github.io/archer2-advanced-use-of-lammps/

[^2]: Both fcc alloys made out of elemental fcc (Ag, Au, Rh, Pd, Ir, Pt) metals as well as bcc alloys made out of elemental bcc metals (Nb, Mo, Ta, W) are of interest. Whilst the former are well known for their ductility, the latter seem to exhibit high hardness.
