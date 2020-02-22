# Exploiting Dynamics Parameter Linearity for Design Optimization in Combined Structural and Dimensional Robot Synthesis


Authors:  
Moritz Schappler, Svenja Tappe, Prof. Dr.-Ing. Tobias Ortmaier

This code is supplemental material to the following submission to the IFToMM World Congress 2019:

```bib
@Conference{SchapplerTapOrt2019,
  author    = {Schappler, M. and Tappe, S. and Ortmaier, T.},
  title     = {Resolution of Functional Redundancy for 3T2R Robot Tasks using Two Sets of Reciprocal Euler Angles},
  booktitle = {Proc. of the 15th IFToMM World Congress},
  year      = {2019},
  month     = {October},
  doi       = {10.1007/978-3-030-20131-9_168},
}
```

## Abstract

Robotic tasks like welding or drilling with three translational and only two rotational degrees of freedom ("3T2R") are of high industrial relevance but are ra\-ther scarcely addressed in scientific publications.
Existing solutions for the resolution of the functional redundancy of robotic manipulators with more than five axes performing these tasks either expand the full kinematic formulation or reduce it in intermediate steps.
This paper presents an approach to reduce the kinematic formulation from the start to solve the problem in a simpler way.
This is done by using a set of reciprocal Euler angles to describe the end-effector orientation and the orientation error in inverse kinematics.

## Contents

This repository contains Maple worksheets to reproduce the equations and validate the theory of the paper:  

* `calculations/ikfr_paper_equations.mw`
  * validates the statement around Equ. 18. (core hypothesis of the paper)
  * creates Matlab functions for Equ. 29
  * creates terms for the implementation of Equ. 31
  * Results are used in [the function constr2grad.m](https://github.com/SchapplM/robsynth-modelgen/blob/master/robot_codegen_scripts/templates_num/robot_constr2grad.m.template), corresponding to Equ. 31
* `calculations/helper/proc_rmat2vec` corresponds to Equ. 27
* `calculations/helper/proc_rvecprod` corresponds to Equ. 28
* `IK_test/IK_3T2R_test.m` performs tests for 4 robots (4DoF to 7DoF)
  * calculates the gradient-based inverse kinematics according to Sec. 3
  * tests the frame transformations from Fig. 1
  * tests the gradient matrix from Equ. 22 with Equ. 18

## Other Implementations and Further Tests

The method described in the paper is employed in the [robotics toolbox](https://github.com/SchapplM/robotics-toolbox) and is used in other related research projects as well, e.g. [the combined structural and dimensional robot synthesis](https://github.com/SchapplM/robsynth-paper_iftommdach2020) or the inverse kinematics of parallel robots.
The functionality of the gradient calculation and the inverse kinematics algorithm is tested to more extent in the test scripts `SerRob_constraints_test.m` and `SerRob_IK_test.m` in the [robotics toolbox](https://github.com/SchapplM/robotics-toolbox).

## Prerequisites

* The Maple worksheet has been tested with Maple version 2018
* For the worksheet `ikfr_paper_equations.mw`, download Maple dependencies containing procedures for Euler angle transformations:
  * https://github.com/SchapplM/robotics-toolbox (tested with rev. bcb2e77; 2020-02-15)
  * https://github.com/SchapplM/matlab_toolbox (tested with rev. 366f149; 2020-02-15)
* Create one path-file for each dependencies as described in the files `calculations/matlab_toolbox_path.template` and `calculations/robotics_toolbox_path.template`
* For the testscript `IK_3T2R_test.m`, download additional dependencies
  * https://github.com/SchapplM/robsynth-serroblib (tested with rev. dff760b; 2020-02-16)
  * https://github.com/SchapplM/robsynth-modelgen (tested with rev. 0f11fd8; 2020-02-15)
  * https://github.com/SchapplM/robotics-dep-ext (tested with rev. bc57a42; 2020-02-17)
