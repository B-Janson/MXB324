# MXB324 Modelling Group 001

This is the MATLAB solution to solve the aquifer in the Dalby region.

## How To Run
The code is all stored within the Verification folder. Please open this to see all scripts and functions. The Newton Solver and Newton-Krylov solver scripts will run the simulation under the parameters specified in VERIF_PARAMS.m. There are extra scripts that have been created to store individual results seen in the report.
<!-- 
There are two parts to this code base:
- The main solver for all parts of the problem
- A verification script to display a simplfied version of the problem, 
that can be shown to be correct to a degree of accuracy.
 -->
## Main Solver
1. Adjust any variables as needed in [INIT_PARAMS](INIT_PARAMS.m).
2. Run the script NEWTON_SOLVER.m or NEWTON_KRYLOV_SOVLER.m.

<!-- ## Verification
All code for verification is contained in the [Verification](Verification/) directory.
1. Adjust any variables as needed in [VERIF_PARAMS](Verification/VERIF_PARAMS.m).
2. Run the script [VERIF_SCRIPT](Verification/VERIF_SCRIPT.m). -->

### Vertices
This section outlines which vertex is which based on a numbering scheme.

---------------- INCLUDE AN IMAGE HERE ---------------------------

Number | Description | Coordinates
------ | ----------- | -----------
1 | Bottom left node | x = 0 && z = 0
2 | Bottom boundary | 0 < x < 500 && z = 0
3 | Bottom right node | x = 500 && z = 0
4 | Left boundary | x = 0 && 0 < z < 80
5 | Interior node | 0 < x < 500 && 0 < z < 80
6 | Right boundary | x = 500 && 0 < z < 80
7 | Top left node | x = 0 && z = 80
8 | Top boundary | 0 < x < 500 && z = 80
9 | Top right node | x = 500 && z = 80

