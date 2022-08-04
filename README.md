# simple_signac_workflow

This repository provides a simple example of using MoSDeF and Signac/Signac-flow with GROMACS


[init.py](blob/main/init.py): Defines the thermodynamic conditions and initializes the signac workspace.  In this example, each job has the same thermodynamic conditions, but the system is composed of different length alkanes.

[project.py](blob/main/project.py): A simple demonstration of using signac-flow to encode a workflow for setting up and performing GROMACS simulations using mBuild and Foyer. 

[analysis.py](blob/main/analysis.py): A very simple framework to demonstrate how to load in your signac project, read in statepoint information, and  to locate the data for each job that would allow you to perform analysis over the dataspace.
