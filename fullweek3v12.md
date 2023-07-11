# Internship - Week 3 - Summary: *Development of training material for Computational Chemistry supercomputing applications at VSC*

Efforts have been devoted to the completion of task 3. As a long-term project, a strategy for task 3 has to be devised alongside giving strides towards the actual practical HPC computational runs of the involved codes.

Progress on task 3 comprising HPC runs of the existing computational chemistry codes has met the added benefit of the intern acquiring `BASH shell` scripting skills as well as of learning the concept of templating.

It must be noted that on Friday, July 7th 2023, a meeting was held between the intern and GJB in order to discuss progress and possible further support required by the intern. This document contains the main points of the discussion. 

In order to foster future discussions, seven separate sections have been included in this document: (i) Practical code (Ab-Initio and MD) HPC running, (ii) `Bash shell` scripting, (iii) templating, (iv) information on the LAMMPS templating library, (v) pre- and post-processing tools in MD and ab-initio  codes, (vi) aspects of the strategy to complete task 3, and (vii) information on existing HPC benchmarks for MD and ab-initio codes.  

## Practical code (Ab-Initio and MD) HPC running

Work on ABINIT and LAMMPS is described next. Most of the intern's time has been devoted to LAMMPS, including information gathering and selected independent HPC runs.

Regarding ABINIT, the problem we had met (mentioned already on week's 1 report) and regarding the failure of the Test Suite for versions 8.8.3 and 8.4.4 on *Genius* seems to persist on version 9.6.0: The intern did this time a `sudo` ABINIT installation on the guest2 laptop under Ubuntu 22.04 LTS, and the 9.6.0 version test suite exhibited the same (Python) error message we had met on the HPC with the more obsolete versions. 

According to the official ABINIT documentation, new developments have been carried out by the authors on the Test Suite for the most advanced versions (9.10.0 is the latest version). The next step on ABINIT will be to try the Test Suite of version 9.8.0 on the *guest2* laptop under Ubuntu 22.04 LTS, to be installed via Anaconda.  

On the other hand, a specific LAMMPS in/data input file pair [^1] has been run using different types of job on Genius: (i) single core/sequential, (ii) pure openMP on a single node -36 core- (iii) pure MPI. At this point, it is worth mentioning that the scaling for the 36 core openMP job was below ten times, which is much lower than expected. As advised by GJB, this is a case thats merits further scaling analysis, starting from a smaller number of cores. Regarding the pure MPI jobs, after repeatedly meeting job failure with the conventional `srun` submission command, GJB has informed me that,  at *Genius*, SLURM batch job sumission, only works with the `mpirun` command instead. Further work is currently being carried out with pure MPI job submissions under mpirun.

The updated markdown file containing a set of tables that describe the current status of this part of the task has been enclosed on the github repository "computationalchemistry". Moreover, a single job script for SLURM batch submission for LAMMPS under pure openMP (36 core), has been enclosed as well on the same repository. 

## Bash shell scripting

`BASH shell` scripting training material from KU Leuven [^2] was recommended by GJB, and is being studied by the intern. The reason for looking into `BASH` is that in order to run AMMPS, there are two input files, and one of them (in.ethanol) requires manual modification of the corresponding data.ethanol full path. I understand this to be a good opportunity to refresh my `BASH shell` scripting skills. 

A very simple script is shown next:
```bash
#/usr/bin/env bash

infile=$(ls in.*)
datafile=$(ls data.*)

mmessage="you need a few changes on the in.* file, and off you go!"
dmm="both files in/data are required (bye)"

if [[ -f $infile && -f $datafile ]]; then
	echo $mmessage
else
	echo $dmm
	break
fi

pathaddress=$(pwd)
echo $pathaddress

pathanddata=$pathaddress/$datafile
echo $pathanddata

# all the entries to be templated within the input in.* file
## flag-variable pair number 1 and 2
m4flag1=xFORM
m4flag2=xDATAP
rendflag1="{{formula}}"
rendflag2="{{fullpathdat}}"
var1=A1B2C3D4E5F6
var2=$pathanddata
# now create the template file
##############keep here all terms to match
match=molecule
printf "the word to match is %s\n" $match
##############keep a log of the full lines to be modified (new file)
lf=linefile.txt
##############
if [ -f $lf ]; then
	printf "a file named %s does already exist" $lf
else
	touch $lf
fi

olx=$(grep -n molecule $infile)
ref=$(grep molecule $infile)
echo $olx >> $lf

while IFS=":" read -r f1 f2
do
        printf "line number: %d/n" $f1
        p1=ol; p2=".txt"; printf -v fbp "%s%d%s" $p1 $f1 $p2
	touch $fbp
	echo $ref > $fbp
done < $lf
```

## Templating

Both templating tools *M4* [^3] and *renderest* [^4] are being tested and incorporated into `BASH shell` script examples as well. During discussions with GJB regarding the above BASH shell scripting exercise for LAMMPS, he mentioned the templating concept and particularly recommended M4.The advantage of the second option (*renderest*) is that it constitutes a complete BASH shell script, although the M4 choice is clear.

Incidentally, both templating Python packages DJANGO and JINJA (valid for both HTML and Python files) are worth mentioning as well.

## Information on the *LAMMPS template library*

The intern has come across information regarding LAMMPS according to which in previous versions there existed a collection of input variable templates [^5] that cover the full functionality: One special script per type of calculation.On the other hand, a recent Python templating system for LAMMPS input variables has been identified as well [^6].   

## Pre- and post-processing tools in MD and Ab-Initio codes

A comprehensive list of pre- and post-processing software is to be enclosed for all the involved  software. At the moment, just the gpaw-tools toolkit is mentioned, as well as the fact that GPAW is developed on the basis of the atomic simulation environment (ASE).

## Aspects of the strategy to complete task 3

In order to make a decision regarding how comprehensively the different codes should be tested and profiled, the model provided by the Sweden
 supercomputing framework [^7] provides a good starting point for discussion.
On the other hand, some of the codes have a wide functionality with different degree of parallelisation, and a decision needs to be made on how far we should go regarding functionality.

## Information on existing HPC benchmarks for MD and Ab-Initio codes

Details of existing benchmarks for large scale calculations will be given shortly regarding the PRACE consortium.It covers a collection of large scale examples for several MD and ab-initio codes, including LAMMPS and Quantum Espresso. On the other hand, the material provided by the MAX centre is to be explored as well.

[^1]: in.ethanol/data.ethanol ("1-performance-exercise"): See the material provided with the tutorial "LAMMPS Course for Intermediate Users" https://epcced.github.io/archer2-advanced-use-of-lammps/
[^2]: "Linux scripting", by Mag Selwa, ICTS-Leuven, https://github.com/orgs/hpcleuven/repositorie
[^3]: https://www.gnu.org/software/m4
[^4]: https://github.com/relaxdiego/renderest
[^5]: "Building a reusable LAMMPS script library", Craig Tenney and Edward Maginn, August 2011
[^6]: https://github.com/anyuzx/LAMMPS-Template-Tool
[^7]: https://enccs.se/

