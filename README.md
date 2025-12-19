This repository contains the code used for the MD simulations in:

#### Spermine regulates condensation of Tau and α-Synuclein and promotes the degradation of condensates in the C. elegans model of PD

Xun Sun, Debasis Saha, Cecilia Mörman, Rebecca Sternke-Hoffmann, Juan Atilio Gerez, Fátima Herranz, Roland Riek, Wenwei Zheng, and Jinghui Luo

#### File description:

**asyn/**: Contains the HOOMD script for running the simulations with 96 chains of alpha-Synuclein and different number of SPM in the system. 
The order for scripts' usage is: 
- initial.py for initializing the system
- resize1.py and resize2.py for resizing the box
- extend.py for extending the z-axis
- final.py for productive simulations

**asyn/all-configurations/**: Contains the gsd files for inital and final steps at different conditions. The naming follows the pattern: (prot-spm ratio)-(initial or final state).gsd.

**tau/**: Contains the HOOMD scripts for running the Tau simulations with different SPM ratios using dual epsilon values. The usage of the scripts is same as the ones in **asyn/**.

**tau/all-configurations**: Contains the gsd files for inital and final steps for each system. The naming follows the pattern: (prot-spm ratio)-(initial or final state).gsd.

**all-rg.txt**: Contains the radius of gyration ($R_g$) values in nm unit from all the single chain simulations of aSyn and Tau without SPM and with 1:20 SPM.

**all-densities.txt**: Contains the density in mg/mL unit of condensed phase region in the phase separation simulations for all the systems.  

**karanicolas_dihe_parm.dat**: dihedral potential parameters from Karanicolas and Brooks, Protein Sci. 11:2351 (2002).




