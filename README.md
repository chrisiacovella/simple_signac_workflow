# simple_signac_workflow

This repository provides a simple example of using MoSDeF and Signac/Signac-flow with GROMACS. 

Note, if you want a version that does not include support for the SLURM queuing system  on Rahman, check out the [Basic release](https://github.com/chrisiacovella/simple_signac_workflow/tree/Basic).  

The main branch includes support for slurm on Rahman. To use this code on Rahman, you'll need to load the Anaconda and Gromacs modules. 

```
module load anaconda
module load gr


[init.py](init.py): Defines the thermodynamic conditions and initializes the signac workspace.  In this example, each job has the same thermodynamic conditions, but the system is composed of different length alkanes. Current project generates systems of methane (C), ethane (CC), and propane (CCC):

```
python init.py
```

[project.py](project.py): A simple demonstration of using signac-flow to encode a workflow for setting up and performing GROMACS simulations using mBuild and Foyer. Usage (execute each in series)

Initialize and parameterize the system, generating  files need to run gromacs:
```
python project.py exec init 
```

Submit the grommp and mdrun commands to the compute nodes, selecting the partition and hardware:
``` 
python project.py submit -o run --gres=gpu:V100:2 --partition=short-tesla --ntasks=8
```
Check on the status of the jobs.
```
python project.py exec check
```

[analyze.py](analyze.py): A very simple framework to demonstrate how to load in your signac project, read in statepoint information, and  to locate the data for each job that would allow you to perform analysis over the dataspace. No analysis is actually performed here, this is just a simple framework to demonstrate how to access the project info

```
python analyze.py
```
