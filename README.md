# simple_signac_workflow

This repository provides a simple example of using MoSDeF and Signac/Signac-flow with GROMACS. 

The [Basic release](https://github.com/chrisiacovella/simple_signac_workflow/tree/Basic) does not provide any code to directly support queuing systems on HPC clusters (e.g., PBS, SLURM, etc.). Support for SLURM on Rahman is included in the main branch.

[init.py](init.py): Defines the thermodynamic conditions and initializes the signac workspace.  In this example, each job has the same thermodynamic conditions, but the system is composed of different length alkanes.

[project.py](project.py): A simple demonstration of using signac-flow to encode a workflow for setting up and performing GROMACS simulations using mBuild and Foyer. 

[analyze.py](analyze.py): A very simple framework to demonstrate how to load in your signac project, read in statepoint information, and  to locate the data for each job that would allow you to perform analysis over the dataspace.
