# Programming bistability in geometrically perturbed mechanical metamaterials (Square-to-Disc shape change)


This repository provides a versatile design framework for designing and optimizing bistable kirigami metamaterials with target shapes and mechanical properties. A Square-to-Disc transformation example is illustrated.



<img src = "https://github.com/Hougeichao/Heterogeneous-Bistable-Kirigami-Square-to-Disc/blob/main/Figure/Github_square-to-disc.png" width="1000"/>


## Paper
This repository contains the code developed for the paper below. Feel free to try this code to create your own planar Kirigami pattern. 

>Peng, Y., Niloy, I., Kam, M., Celli, P., & Plucinsky, P. (2024). "[Programming bistability in geometrically perturbed mechanical metamaterials](https://doi-org.libproxy2.usc.edu/10.1103/PhysRevApplied.22.014073)". Physical Review Applied, 22(1), 014073.



## Getting Started

### Prerequisites
* MATLAB (tested on R2023a and later versions)
* Optimization Toolbox (for `FMINCON`)

### Usage 

* Main Script:

- Run `Heter_Disc.m` to generate 'Square-to-Disc' kirigami designs.
  - **Figure 1 & 2**: Present the mapping between the lattice vectors in the reference and deformed plane.
  - **Figure 3 & 4**: Illustrate the 'Disc' designs in reference and deformed plane, respectively.



### Key Functions

* `Heter_Disc.m`: Contains MATLAB code for 'Disc' model simulation.
* `Para.m`: Contains Matlab code for parameterizing all the vectors in kirigami design through the design vectors and prescribed lattice vectors.
* `Obj.m`: Shows the objective function used for the very first step during the optimization.
* `obj_2.m`: Shows the objective function used for cells except for the first one. It minimizes the difference of all parameters between the neighboring cells.
* `Constraint.m`: Represents the inequality constraints used for the optimization at the first step (corresponding to `Obj.m`).
* `Energy_spr.m` and `EnergyBarrier.m` construct the functions for the elastic energy and energy barrier respectively.
* `Rot.m`: Represents the 2D rotation tensor.
* `ccc.m`: Is used for the clearance.

### Customization and Optimization

* Feel free to experiment with different numbers of quad mesh (default number is 40-by-40 quad mesh)
* If the optimization does not yield a desirable pattern, consider the following:
  * Adjusting fmincon parameters (e.g., MaxFunctionEvaluations, MaxIterations, ConstraintTolerance) to refine the optimization process.
  * Modifying the prescribed Bravais lattice vectors to explore different configurations.
  * Tuning the objective function to better match your design goals.
  * Changing the initial point.

### How to cite
If you use this code in your work, please cite the following [paper](https://doi-org.libproxy2.usc.edu/10.1103/PhysRevApplied.22.014073)

```bibtex
@article{peng2024programming,
  title={Programming bistability in geometrically perturbed mechanical metamaterials},
  author={Peng, Yingchao and Niloy, Imtiar and Kam, Megan and Celli, Paolo and Plucinsky, Paul},
  journal={Physical Review Applied},
  volume={22},
  number={1},
  pages={014073},
  year={2024},
  publisher={APS}
}
```
