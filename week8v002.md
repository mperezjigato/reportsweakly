# Internship - Week 8 - Summary

## Computational Chemistry

Work carried out with Quantum Espresso [^1] is described next, and the tables depticting our progress with the Computational Chemistry packages are updated accordingly (repository "computationalchemistry"). A set of independent runs has been carried out both sequentially and under the MPI parallel mode. On the other hand, the as-provided test-suite is currently being tried.

Incidentally, it must be noted that all our CPU time on GENIUS has been consumed, and a new application for five million extra credits duly submitted [^2].

### Quantum Espresso calculations

Two QE as-provided independent systems as well as a system with input files built from scratch, have been run: 

 - GaAs PWSCF calculation under PAW,
 - SiH4 (molecule in a box) PWSCF calculation under norm-conserving pseudopotentials,
 - Immm Ag2PdO2 PWSCF calculation under ultrasoft pseudopotentials: The materialscloud.org website (Quantum ESPRESSO input generator and structure visualizer) lets us transform the corresponding CIF file on QE inputs, which has been obtained from the topological materials database.

### Specific CIF files and geometry visualisation 

The Ag(16 10 9) surface structure is depicted below as well as the corresponding CIF file.

 ![](Ag16109.png)
 
 Moreover, the incommensurate structure of the metallic solid solution (ss) AuCd and  its  CIF file is shown below [^3]. Incidentally, the AuCd ss are well known for their ferroelastic properties under the 50/50 composition [^4],
 
 Finally, the crystal structure of the above QE example Ag2PdO2 (topological material) is depicted, as well as its CIF file:

[^1]: www.quantum-espresso.org
[^2]: https://admin.kuleuven.be/icts/onderzoek/hpc/extra-project-credits
[^3]: It has been downloaded from The Bilbao Incommensurate Crystal Structure Database at https://www.cryst.ehu.es/
[^4]: EKH Salje, "Phase transitions in ferroelastic and co-elastic crystals", Student edition, 1993, Cambridge University Press