# Internship - Week 3 - Summary: *Development of training material for Computational Chemistry supercomputing applications at VSC*

Efforts have been devoted to the completion of task 3. As a long-term project, a strategy for task 3 has to be devised alongside giving strides towards the actual practical HPC computational runs of the involved codes.

Progress on task 3 comprising HPC runs of the existing computational chemistry codes has met the added benefit of the intern acquiring `BASH shell` scripting skills as well as of learning the concept of templating.

It must be noted that on Friday, July 7th 2023, a meeting was held between the intern and GJB in order to discuss progress and possible further support required by the intern. This document contains the main points of the discussion. 

In order to foster future discussions, seven separate sections have been included in this document: (i) Practical code (Ab-Initio and MD) HPC running, (ii) `Bash shell` scripting, (iii) templating, (iv) information on the LAMMPS templating library, (v) pre- and post-processing tools in MD and ab-initio  codes, (vi) aspects of the strategy to complete task 3, and (vii) information on existing HPC benchmarks for MD and ab-initio codes.  

## Practical code (Ab-Initio and MD) HPC running

Work on ABINIT and LAMMPS is described next. Most of the intern's time has been devoted to LAMMPS, including information gathering and selected independent HPC runs.

Regarding ABINIT, the problem we had met (mentioned already on week's 1 report) and regarding the failure of the Test Suite for versions 8.8.3 and 8.4.4 on *Genius* seems to persist on version 9.6.0: The intern did this time a `sudo` ABINIT installation on the guest2 laptop under Ubuntu 22.04 LTS, and the 9.6.0 version test suite exhibited the same (Python) error message we had met on the HPC with the more obsolete versions. 

According to the official ABINIT documentation, new developments have been carried out by the authors on the Test Suite for the most advanced versions (9.10.0 is the latest version). The next step on ABINIT will be to try the Test Suite of version 9.8.0 on the *guest2* laptop under Ubuntu 22.04 LTS, to be installed via Anaconda.  

On the other hand, a specific LAMMPS in/data input file pair [^1] has been run using different types of job on Genius: (i) single core/sequential, (ii) pure openMP on a single node -36 core- (iii) pure MPI. At this point, it is worth mentioning that the scaling for the 36 core openMP job was below ten times, which is much lower than expected. As advised by GJB, this is a case thats merits further scaling analysis, starting from a smaller number of cores. Regarding the pure MPI jobs, after repeatedly meeting job failure with the conventional `srun` submission command, GJB has informed me that,  at *Genius*, SLURM batch job sumission, only works with the `mpirun` command instead. Further work is currently being carried out with pure MPI job submissions under mpirun.On the other hand, the updated markdown file containing an updated set of tables that describe the current status of this part of the task has been enclosed on the github repository "computationalchemistry".

## Bash shell scripting

`BASH shell` scripting training material from KU Leuven [^2] was recommended by GJB, and is being studied by the intern. The reason for looking into `BASH` is that in order to run LAMMPS, there are two input files, and one of them (in.*) requires manual modification of data.* full path. I understood this to be a good opportunity to refresh my `BASH shell` scripting skills. 

## Templating

Both templating tools *M4* [^3] and *renderest* [^4] are being tested and incorporated into `BASH shell` script examples as well. During discussions with GJB regarding the above BASH shell scripting exercise for LAMMPS, he mentioned the templating concept and particularly recommended M4.The advantage of the second option (*renderest*) is that it constitutes a complete BASH shell script, although the M4 choice is clear.

Incidentally, both templating Python packages DJANGO and JINJA (valid for both HTML and Python files) are worth mentioning as well.

## Information on the *LAMMPS template library*

## Pre- and post-processing tools in MD and Ab-Initio codes

## Aspects of the strategy to complete task 3

## Information on existing HPC benchmarks for MD and Ab-Initio codes


[^1]: in.ethanol/data.ethanol ("1-performance-exercise"): See the material provided with the tutorial "LAMMPS Course for Intermediate Users" https://epcced.github.io/archer2-advanced-use-of-lammps/
[^2]: "Linux scripting", by Mag Selwa, ICTS-Leuven, https://github.com/orgs/hpcleuven/repositorie
[^3]: https://www.gnu.org/software/m4
[^4]: https://github.com/relaxdiego/renderest
