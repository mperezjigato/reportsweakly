# Internship - Week 1 : Summary 

This document describes the proposed tasks and the work carried out during the first week of internship at the University of Hasselt, as well as the training sessions attended by the intern.

## Tasks and work

1. Carry out a full application for an HPC cluster account at VSC-KU Leuven-ICTS. **COMPLETED**
1. In order to make it fully reproducible for the user, the complete application process (1) is to be documented in the form of computer screenshots. **COMPLETED**
1. As part of the Computational Chemistry work, the intern is to locate the relevant software already installed on GENIUS, as well as to run batch jobs from the set of selected packages. The full task is long term and involves producing training material in the form of documentation specific to the GENIUS HPC system. It includes input files, memory management information, scalability, etc. **WORK IN PROGRESS**
1. The intern is to get familiar with the markdown typesetting system. **WORK IN PROGRESS**
1. The intern is to build a GITHUB repository containing all the computer screenshots from task (2). This task is part of the practicals to be carried out by the intern in order to get familiar with GITHUB/GIT. Following the training session described below and the training material provided by GJB, the intern is to spend some time trying simple cases from scratch. The first step is to open a new GITHUB account under the Hasselt University e-mail. Furthermore, the intern is to fork one of the repositories set up by GJB during the training session, in order to run simple cases. **WORK IN PROGRESS**  

### Task 1: completed

### Task 2: completed

### Task 3: work in progress

Two separate categories of installed sotware within the Computational Chemistry realm can be established, viz first-principles quantum mechanical (otherwise denominated Ab-Initio calculations) codes and molecular dynamics codes. The relevant software currenly located by the intern within the first category on GENIUS includes VASP, ABINIT, QUANTUM EXPRESSO and WANNIER90 (CP2K). Within the second category, the LAMMPS code is to be highlighted. Other MD codes of interest include GROMACS and NAMD.

Further software packages of future interest are described next.
In relation to the second category, the AlphaFold software is currenly installed in GENIUS in combination with the openMM molecular dynamics code. On the other hand, the ELPA eigenproblem solver has been located: A key potential component to future implementations of some/most of the Ab-Initio software packages. No doubt, the efficiency benefits are well known. It is worth mentioning the SLEPC eigenproblem solver implementation as a separate option to ELPA. Incidentally and covering both categories, the atomic simulation environment (ASE) cannot be overlooked as a high-througput computational materials science platform. Kinetic simulations could be performed via PLUMED (an implementation exists on GENIUS) and with further rare event sampling software that will require full installation (see "nested sampling" in materials science -Pymatnest-). As a matter of fact, Prof Peter Bolhuis (the author of openMM) has developed openpathsampling as well, a transition path sampling code for rare event simulations. The uninstalled package P4 (unstable/metastable systems) would also be relevant to kinetic simulations in materials science.

The ABINIT code, on the one hand has been run using the slurm scheduler for small test systems on 18 cores, and on a single node (thin node - CPU). These jobs have been completed successfully in less than a minute of wall time. Furthermore, the intern has tried to run the full test bunch provided by ABINIT, and accompanying the full easybuild implementation (2018a). The combined test system has failed, giving a python error. Further work remains to be carried out on ABINIT.

The VASP code cannot be accessed. This has been actually tested by installing the appropriate 2018a module: Certainly, being a commercial code, VASP requires a licence, and GJB is to set up access to one of the existing licences. Regarding Quantum Espresso and the Wannier90 codes, there has been no progress up to now. One important point related to the latter is that the Wannier90 existing implementation is an outdated one. The intern requires version 3 in order to carry out current developments. The implementation of the various first principles codes in combination with the ELPA software, constitutest another request to be submitted to VSC.

The LAMMPS code is to be tested shortly, including both the CPU and the GPU versions. The ASE and AlphaFold codes remain to be tested as well. 

### Task 4: work in progress

The intern is trying very basic markdown typesetting as evidenced by this document.

### Task 5: work in progress

The intern has created a GITHUB/GIT account and uploaded an SSH public key. He has also forked a folder provided by GJB.

## Training sessions

On June 22, 2023 the intern attended a session provided by GJB to students of the Data Science Institute describing the HPC system GENIUS at VSC-KU Leuven-ICTS. On the other hand, a separate training session was provided on the same day by GJB regarding the GITHUB/GIT version control system. Specific examples were built and provided by GJB for the intern to experiment on. 
