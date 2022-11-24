# data_driven_pos
Data-driven stabilization of positive linear systems. 
A collection of state-input-transition data is collected from a continuous-time or discrete-time LTI dynamical system, under a known L-infinity noise bound. A full-state feedback control u=Kx is computed in order to stabilize all systems that are consistent with the observed data.

Stabilization is certified by the Extended Farkas Lemma (polytope inclusion): the polytope of systems that can be stabilized by the solved controller contains the polytope of plants that are consistent with the data.

Additional features include sign-patterns/structured control and peak-to-peak gain minimization.

## Dependencies

- YALMIP: https://yalmip.github.io/
- Mosek: https://www.mosek.com/ (or any solver compatible with YALMIP)

All code is written and tested on Matlab R2021a.


## Instructions

The stabilization classes using the Extended Farkas Lemma are `posstab_f` (discrete-time) and `posstab_cont_f` (continuous-time).

The simulation classes are `possim` (discrete-time) and `possim_cont` (continuous-time). These methods will optionally generate a random system (`rand_sys`), and will simulate a trajectory (`sim` for discrete-time) or a set of transitions (`sample_slope` for continuous-time) as corrupted by the given L-infinity noise bound. 

The resultant trajectory/collection of data is passed to the stabilization classes. An additional argument to the stabilization classes is the `data_opts` class, with the following options:
- `nontrivial`: eliminate trivial faces from the data-consistency polytope
- `pos_A`: impose prior knowledge that A is a positive system (A Metzler/Nonnegative)
- `pos_B`: impose prior knowledge that B is an internally  positive system (B Nonnegative)
- `gez`: A 0/1 matrix where 1-entries of the controller are greater than or equal to zero
- `lez`: A 0/1 matrix where 1-entries of the controller are less  than or equal to zero

### Robust Optimization

The classes `posstab` and `posstab_cont` use the Yalmip robust optimization package https://yalmip.github.io/tutorial/robustoptimization/ to eliminate the plant parameters (A, B) and derive the robust counterpart (utilizing the duality option). The derived programs are equivalent, but utilizing the `robustoptimization` package adds a preprocessing cost.

The README in the `experiments` folder describes all test scripts.

## Reference
To be uploaded to arXiv.

## Contact
For comments and questions please email [Jared Miller](mailto:miller.jare@northeastern.edu?Subject=data_driven_pos).

