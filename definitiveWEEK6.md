
# Internship - Week 6 - Summary

## Progress on HPC runs of MD and Ab-Initio codes within the *Development of Supercomputing Training Material for Computational Chemistry* task

Following previous work on the Computational Chemistry task, I have decided not to report any hybrid calculations data for the time being. Once the consistency problem we identified last week for hybrid parallel calculations is fully resolved, calculations using that mode will resume. Therefore, only OpenMP and MPI data are discussed below.

On the one hand, simple materials science problems provided within the examples directory of the LAMMPS distribution have been run (the force-field type is given in brackets). Wall time tables for the first case are given at the bottom of the document. No data are reported for the other three cases. The systems that have been run include:
 - Amorphous Silica geometry generation (Vashishta).
 - Cu2O elasticity (COMB).
 - OH adsorbed on Cu2O(110) deposited graphene monolayer (COMB).
 - mHfO2 minimisation (COMB).

On the other, the five as provided LAMMPS benchmarks have been run: 32000 atom systems as starting geometry. Namely,
 1. NVE integration of a Lennard-Jones fluid.
 2. NVE integration of a bead-spring polymer melt under the FENE approximation.
 3. NVE integration of fcc Cu under the EAM potential.
 4. NVE integration of granular Chute flow.
 5. NPT integration of the rhodopsin protein in a solvated lipid bilayer under the CHARMM force-field and the PPPM Coulomb solution method.
    
Three types of calculation have been tried for each, ie sequential, fix-scale (MPI) and scaled-size (MPI) calculations. Data are reported only for the MPI scaled-size calculations [^1] (see below). It should be noted that the number of time-steps has been modified in order to ensure a squential wall-time of ca. 2 minutes in all five cases.  

Moreover, the `pyiron` HTCMS software [^2] has been locally implemented by the intern on the guest2 laptop under Ubuntu 22.04 LTS, in order to shortly start generating geometries for single-phase crystalline solid solutions as seems optimal for High-Entropy alloys. Likewise, one of the traditional problems of equilibrium molecular dynamics, the calculation of heat transport, seems to have been resolved by Stefano Baroni and co-workers via a new theory based on the cepstral analysis of current time series derived from equilibrium molecular dynamics. We are in the process of installing Baroni's SporTran [^3] software for heat transport calculations in combination with LAMMPS.

Finally, it should be mentioned that the authors of VASP have agreed in granting a VASP 6 licence to GJB, and that he is in the process of signing the contract documents.  

### Amorphous Silica geometry generation (Vashishta)

#### OpenMP

|     |  4  |  8  |  12  |  16  |  20  |  24  |  28  |  32  |  36  |seq. |
|-----|-----|-----|------|------|------|------|------|------|------|-----|
|t/sec|19748|11509| 8834 | 7308 | 6469 | 6272 | 7875 | 7058 | 6070 |77556|

#### MPI

|     |  4  |  8  | 12 | 16 | 20 | 24 | 28 | 32 | 36 | 72 | 108| 144|
|-----|-----|-----|----|----|----|----|----|----|----|----|----|----|
|t/sec|20610|11920|8817|6842|5414|5047|4777|4449|4225|2965|2645|2607|

### MPI scaled-size benchmarks

#### LJ fluid (sequential:126.3sec/7534steps/32000 atoms)

|core|16   |36   |64   |81   |100  |144  |196  |225  |256  |324  |400  |441  |
|----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|
|t/s |154.4|189.2|185.3|196.2|190.2|196.8|193.4|195.0|194.8|201.6|195.1|183.3|

|core|484  |576  |625  |676  |729  |784  |900  |1024 |1089 |1156 |1225 |1296 |
|----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|
|t/s |198.8|205.1|203.3|204.9|209.0|205.6|182.4|186.2|188.8|209.2|187.8|188.7|

|core|1444 |1521 |1600 |1764 |1936 |2025 |2116 |2304 |2500 |2601 |2704 |2916 |
|----|---- |-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|
|t/s |211.1|198.5|190.2|217.1|191.6|194.3|189.6|     |     |194.4|196.1|195.5|

#### Bead-Spring Polymer Melt (sequential:124.4sec/12945steps/32000 atoms)

|core|16   |36   |64   |81   |100  |144  |196  |225  |256  |324  |400  |441  |
|----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|
|t/s |174.5|229.9|231.8|249.6|268.1|260.0|258.9|289.8|265.8|273.3|291.8|252.1|

|core|484  |576  |625  |676  |729  |784  |900  |1024 |1089 |1156 |1225 |1296 |
|----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|
|t/s |270.5|281.2|310.7|283.0|285.2|289.5|279.3|265.7|263.3|297.3|271.1|267.5|

|core|1444 |1521 |1600 |1764 |1936 |2025 |2116 |2304 |2500 |2601 |2704 |2916 |
|----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|
|t/s |300.6|277.9|298.0|309.5|279.2|277.3|274.1|     |     |286.4|294.3|291.7|

#### Fcc Copper (sequential:122.5sec/2831steps/32000atoms)

|core|16   |36   |64   |81   |100  |144  |196  |225  |256  |324  |400  |441  |
|----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|
|t/s |150.6|175.2|172.9|180.6|174.9|183.6|184.3|179.5|180.6|186.5|181.9|169.9|

|core|484  |576  |625  |676  |729  |784  |900  |1024 |1089 |1156 |1225 |1296 |
|----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|
|t/s |184.5|190.3|190.5|190.0|191.9|193.4|175.6|174.3|174.9|197.9|175.3|174.3|

|core|1444 |1521 |1600 |1764 |1936 |2025 |2116 |2304 |2500 |2601 |2704 |2916 |
|----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|
|t/s |199.0|179.0|176.8|202.4|179.4|179.0|178.2|     |     |181.2|185.0|183.1|

#### Granular Chute flow (sequential:128.9sec/33670steps/32000atoms) 

|core|16   |36   |64   |81   |100  |144  |196  |225  |256  |324  |400  |441  |
|----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|
|t/s |253.4|379.8|358.2|375.0|372.6|388.7|381.1|400.3|384.7|405.5|397.8|368.1|

|core|484  |576  |625  |676  |729  |784  |900  |1024 |1089 |1156 |1225 |1296 |
|----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|
|t/s |396.8|416.4|409.4|409.7|412.6|440.3|391.2|387.5|390.7|433.9|403.0|399.6|

|core|1444 |1521 |1600 |1764 |1936 |2025 |2116 |2304 |2500 |2601 |2704 |2916 |
|----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|
|t/s |441.9|405.9|409.9|467.7|408.8|432.4|412.1|     |     |418.1|427.1|453.1|

#### Rhodopsin in a solvated lipid bilayer (sequential:117.5sec/500steps/32000atoms)

|core|16   |36   |64   |81   |100  |144  |196  |225  |256  |324  |400  |441  |
|----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|
|t/s |150.4|173.9|176.5|189.8|181.8|186.7|184.9|191.4|193.7|193.6|189.2|179.0|

|core|484  |576  |625  |676  |729  |784  |900  |1024 |1089 |1156 |1225 |1296 |
|----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|
|t/s |193.0|197.8|202.0|200.3|197.6|209.4|179.1|190.4|186.0|210.0|190.7|184.4|

|core|1444 |1521 |1600 |1764 |1936 |2025 |2116 |2304 |2500 |2601 |2704 |2916 |
|----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|
|t/s |212.1|194.1|188.4|218.0|203.3|190.1|207.1|     |     |210.7|215.1|212.6|

[^1]: On July 25th, 2023, GJB organised a training session regarding parallel computing. He then discussed the OpenMP, the MPI and the hybrid (OpenMP/MPI) parallel modes, as well as the concepts of lattency, bandwidth, strong and weak scaling, speedup and parallel efficiency. Within the weak scaling realm, it is key to investigate size scaling. The herein reported set of benchmarks goes along those lines, as opposed to the fixed-size problem (varying the number of processors whilst keeping the same size, as in strong scaling). In a previous report I had shown OpenMP and MPI scaling results for two separate systems, although they actually corresponded only to the strong scaling concept (weak scaling had been ignored at that point by the intern).
[^2]: www.pyiron.org ; a tool for geometry generation of crystalline solid solutions is contained within pyiron.
[^3]: https://github.com/sissaschool/sportran.
