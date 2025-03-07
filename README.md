# Heterogeneous-Bistable-Kirigami-Square-to-Disc
Programming bistable kirigami metamaterials with target heterogeneous shape change (Square-to-Disc)

Please run 'Heter_Disc.m' to generate 'Square-to-Disc' kirigami designs. Figure 1 and Figure 2 present the mapping between the lattice vectors in the reference and deformed plane. Figure 3 and Figure 4 present the 'Disc' designs in reference and deformed plane respectively.

**Heter_Disc.m** contains MATLAB code for 'Disc' model simulation.

**para.m** contains Matlab code for parameterizing all the vectors in kirigami design through the design vectors and prescribed lattice vectors.

**obj.m** shows the objective function used for the very first step during the optimization.

**obj_2.m** shows the objective function used for cells except for the first one. It minimizes the difference of all parameters between the neighboring cells.

**constraint** represents the inequality constraints used for the optimization at the first step (corresponding to obj.m).

**energy.m** and **EnergyBarrier.m** construct the functions for the elastic energy and energy barrier respectively. 

**Rot.m** represents the rotation tensor. **ccc.m** is used for the clearance.
