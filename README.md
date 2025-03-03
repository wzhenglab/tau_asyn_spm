This repository contains the code used for the MD simulations in:

#### Spermine regulates condensation of Tau and α-Synuclein and promotes the degradation of condensates in the C. elegans model of PD

Xun Sun, Debasis Saha, Cecilia Mörman, Rebecca Sternke-Hoffmann, Juan Atilio Gerez, Fátima Herranz, Roland Riek, Wenwei Zheng, and Jinghui Luo

#### File description:

**asyn**: Contains the HOOMD script for running the simulations with 96 chains of protein and different number of SPM in the system. The order for scripts' usage is: initial.py > resize1.py > resize2.py > extend.py > final.py. <br>
**asyn/all-configurations**: Contains the gsd files for inital and final steps for each system. The naming follows the pattern: (prot-spm ratio)-(initial or final state).gsd. <br>
**tau**: Contains the HOOMD scripts for running the tau simulations with different SPM ratios using a dual epsilon value. The usage of the scripts is same as asyn. <br>
**tau/all-configurations**: Contains the gsd files for inital and final steps for each system. The naming follows the pattern: (prot-spm ratio)-(initial or final state).gsd. <br>
**all-rg.txt**: Contains the radius of gyration (rg) values in nm unit from all the single chain simulations of asyn and tau without SPM and with 1:20 SPM. <br>
**all-densities.txt**: Contains the density in mg/mL unit of condensed phase region in the phase separation simulations for all the systems.  





