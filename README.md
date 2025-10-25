# Finite Element Models for Multiscale Theory of Porous Media

PoroelasticityFEM is a porous media Finite Element code developed during the PhD project of Bruna Campos at the University of Waterloo (2021-2024), under the supervision of Prof. Robert Gracie.
The PhD thesis related to this implementation is available [here](https://hdl.handle.net/10012/21213).

The code contains 12 distinct main solvers, encompassing models based on Biot (Bt) annd de la Cruz and Spanos (dCS) porous media theories, dynamic and transient formulations, and five different combinations of main variables.

## Available modules

- `RunSim`: main function to run. This is where the configuration file is read and `main` is called. The existent configuration files are listed within this file.
- `RunTests`: used to test the code functionalities. Contains patch tests and manufactured solution tests for the uncoupled equations.
- `StudyWaveVelAtt` and `StudyWaveVelAttContours`: functions used to plot wave velocity and attenuation patterns; contain various sets of material parameters to be chosen from. 


## Keywords 
Computational Geomechanics, Porous Media, Finite Element Method, Biot Theory, de la Cruz and Spanos Theory, Wave Propagation, Consolidation, Water-saturated Rocks
