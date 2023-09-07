
# Internship - week 11 - summary

 ASE (under empirical EMT potentials) and GPAW total energy calculations are described next:
 - Ag nanoparticles, including an MPI weak scaling study for GPAW calculations of both octahedral and truncated octahedra Ag (fcc) nanoparticles; 
 - force fields from the "openkim" collection for bcc metals (W) are tested for bulk, surface and nanoparticle cases;
 - interfaces with liquid water (TIP potentials) of both fcc and bcc metals;
 - a comprehensive test suite for GPAW is run and reported.
 
 On the other hand, a couple of examples of atomistic model build-up using `atomsk` are shown below for extended defects in materials science,
 
## Ag Nanostructures - ASE (EMT) geometry optimisation, GPAW (LCAO) self-consistent field electronic structure and GPAW electronic density of states calculations (*count.tar.gz* and *nanoparticlessilver.tar.gz* uploaded to the GitHub repository "computationalchemistry")

One of the valuable aspects of ASE corresponds to its capability to set up atomistic geometry models from scratch. The "ase.cluster.octahedron" module allows the creation of atomistic geometries for this type of nanostructures, including both octahedra and truncated octahedra. The FCC Ag example from the ASE tutorial material has been chosen to be run (based on `ase.cluster.octahedra`). It includes (a) the generation of atomistic geometries and the ASE (EMT) optimisation; (b) the GPAW computation of the self-consistent electronic structure for the optimal cluster geometry, and (c) a GPAW electronic density of states calculation starting from the GPAW electronically converged system. All three steps are (MPI) run on Genius, for nanostructures of increasing size and a full weak-scaling study is presented. Incidentally, GPAW electronic structure calculations can be run using three different modes: (LCAO) atom-centered basis functions, real-space grids/multigrids and planewaves. This test has been run exclusively under (LCAO) atom-centered basis functions.

Although they have been downloaded directly from the ASE and GPAW webpages, the three python/ASE/GPAW input scripts used in this calculation are shown below for the sake of completeness (another subsection of this document will suggest modifications in order to achieve specific HPC procedures):
```python









```
```python








```
```python











```
Two separate tables are shown: One with the number of atoms in the nanoparticle, the other with the MPI wall-time for the second step of the calculation (GPAW self-consistent field electronic structure calculation). Neither the first (EMT geometry optimisation) nor the third step (electronic density of states) are accounted for at this point. As a weak scaling problem, the number of MPI cores increases proportionally to the number of cluster atoms (the number of cores is ca. atom number divided by 10). On the other hand, two types of truncation is considered, ie a `cutoff=2` truncation and a cutoff 40 % truncation.

|NPAR|c=40%|NATpOH|NATtOH2|NATtOH40|Nc(pOH)|Nc(c=2)|Nc(c=40)|
|----|-----|------|-------|--------|-------|-------|--------|
|5   |  2  |85    |    55 |   ---  |     9 |     6 |     -- |
|6   |  -  |146   |   116 |   ---  |    15 |    12 |     -- |
|7   |  -  |231   |   201 |   ---  |    23 |    20 |     -- |
|8   |  -  |344   |   314 |   ---  |    34 |    31 |     -- |
|9   |  -  |489   |   459 |   ---  |    49 |    46 |     -- |
|10  |  4  |670   |   640 |   490  |    67 |    64 |     49 |
|11  |  -  |891   |   861 |   ---  |    89 |    86 |     -- |
|12  |  -  |1156  |  1126 |   ---  |   116 |   113 |     -- |
|13  |  -  |1469  |  1439 |   ---  |   147 |   144 |     -- |
|14  |  -  |1834  |  1804 |   ---  |   183 |   180 |     -- |
|15  |  6  |2255  |  2225 |  1709  |   226 |   223 |    171 |
|20  |  8  |5340  |  5310 |  4116  |   534 |   531 |    412 |

The first column contains the values of the main parameter defining the octahedron size (NPAR). The second one is the value of the truncation parameter c for 40 % of the original octahedron size parameter (NPAR).
The third, fourth and fifth columns contain the main piece of information of this table, ie the number of atoms of each specific cluster (pOH refers to "pure-octahedron"; tOH2 is the truncated octahedron originally built under size parameter NPAR and subsequent truncation variable c=2; tOH40 is the truncated octahedron originally built under size parameter NPAR and truncation variable c=40% of NPAR). The last three columns contain the number of MPI cores for each calculation.

 Regarding the wall-time table:
 
|NPAR|NCpOH|NCc=2|NC c=40|WTpOH [s]|WTtOH2[s]|WTtOH40  |
|----|-----|-----|-------|---------|---------|---------|
|5   |9    |6    |   --- |75.258   |17.607   |  ---    |
|6   |15   |12   |   --- |92.670   |36.106   |  ---    |
|7   |23   |20   |   --- |176.731  |69.871   |  ---    |
|8   |34   |31   |   --- |220.460  |140.867  |  ---    |
|9   |49   |46   |   --- |795.626  |495.511  |  ---    |
|10  |67   |64   |    49 |1574.049 |1046.277 |269.801  |
|11  |89   |86   |   --- |FAIL[^1] |2421.106 |  ---    |

[^1]: The GPAW total energy calculation for this 891 atom (fcc Ag) octahedral nanoparticle (NPAR=11) has simply reached the time limit, without failure. Moreover, it should be noted that *the first step of the calculation, ie the EMT classical geometry optimisation (ASE), has converged in 43 cycles*. Regarding the GPAW total energy calculation, it must be mentioned that 89 electronic iterations had been completed before the calculation reached the time limit, without electronic convergence.

The GPAW scf calculations fail for NPAR above 11. After simple inspection of the total energies, it is apparent that they exhibit huge fluctuations, which, in turn is characteristic of the charge sloshing behaviour of delocalised electronic systems (metals), when LCAO basis functions are used.

### *HPC parameterisation*: An exercise regarding automatisation for ASE-GAPW calculations including BASH scripting and ATOOLS job submission

A more professional approach to the above ASE-GPAW calculation is described next. It actually constitutes a unique chance to try not only your BASH capabilities and learn the ATOOLS scheduler control system, but the following will be brashed off:

  - Python: ASE and GPAW both require python scripts
  - M4: In order to HPC-parameterise the calculation, several python variables need search/replacement
  
As a matter of fact, a new python script is proviby by GJB in order to carry out a particular subtask (search/replacement of a particular variable with special syntax, on a python script) and used directly by the intern, therefore, not much merit has to be attributed to the intern himself. Regarding the M4 tool, it has been used for simple cases.

A couple of comments are in order in relation to M4, before we carry on. Firstly, and since the intern did not explicitly mention it during the LAMMPS input file exercise in previous weekly reports, it should be noted that the case at hand (LAMMPS ethanol input file search/replacement) is easily solved:

>
> $ m4 -D xFORM="CH3COOH" -D xDATAP=$(pwd) templateM4.txt > in.acetic
>
  
where the templateM4.txt file (modified "in.ethanol" file) contains two lines (36 and 37) with the actual variables to be searched/replaced:

```
molecule        xFORM xDATAP
create_atoms    0 region start_box mol xFORM 6871
```
The rest of the document coincides with "in.ethanol".

The second comment is about more involved M4 syntax: It turns out a string `part1_part2` requires a specific macros in order to have "part2" searched/replaced, whilst leaving "part1_" still. If you create a file (preserve.m4) with the contents:
```m4
format(SOMETHING_%s, ELSE)
```
the search/replace operation is now:

>
> $ m4 --define ELSE='123456789' preserve.m4 > output.txt
>

the file "output.txt" successfully containing:
```
SOMETHING_123456789
```
A good example is now in order. Imagine both the "atoms_xEXT" variable and the "atoms_xEXT.xyz" filename appearing within the ASE python script below (s1geo.m4), needs to have the part "xEXT" searched/replaced:

```python
from ase.cluster import Octahedron
from ase.io import write, read
atoms_xEXT = Octahedron('Ag', xNP)
write('atoms_xEXT.xyz', atoms_xEXT)
```
In this case,

```bash
#!/usr/bin/env bash
OLD='atoms_xEXT'
NEW='format(atoms_%s, xEXT)'
sed -i "s#${OLD}#${NEW}#g" s1geo.m4
```
the simple BASH script (name it "file.sh") will create the appropriate M4 macros (bash ./file.sh):

```m4
from ase.cluster import Octahedron
from ase.io import write, read
format(atoms_%s, xEXT) = Octahedron('Ag', xNP)
write(format(atoms_%s, xEXT).xyz', format(atoms_%s, xEXT))
```
which, after applying M4:

>
> $ m4 --define xEXT="123456789" s1geo.m4 > result.py
>

gives the resulting python file:

```python
from ase.cluster import Octahedron
from ase.io import write, read
atoms_123456789 = Octahedron('Ag', xNP)
write(atoms_123456789.xyz', atoms_123456789)
```