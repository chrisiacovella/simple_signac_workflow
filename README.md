# simple signac workflow

This repository provides a simple example of using [MoSDeF](https://github.com/mosdef-hub) and [Signac/Signac-flow](https://docs.signac.io/en/latest/index.html) with GROMACS.  Note this was updated to work with the syntax changes for decorators that came with [0.23 in signac-flow](https://docs.signac.io/projects/flow/en/latest/changes.html#id12).

The main branch includes support for slurm on Rahman. To use this code on Rahman, you'll need to load the Anaconda module first.  Note, selection of the appropriate GROMACS module (in this case 2020.6) is set in the run function in the project.py file. 


```
module load anaconda
```

The first time we run this, let us create a conda environment.  We can do this two ways...manually:

```
conda create --name simple_signac_gmx

conda activate simple_signac_gmx

conda install -c conda-forge mbuild foyer rdkit signac signac-flow
```

or create using the environment.yml file (note I include specific versions in the environment file that were used to test the framework):

```
conda env create -f environment.yml

conda activate simple_signac_gmx
```

General usage below:

[init.py](init.py): Defines the thermodynamic conditions and initializes the signac workspace.  In this example, each job has the same thermodynamic conditions, but the system is composed of different length alkanes. The current project generates systems of methane (C), ethane (CC), and propane (CCC):

```
python init.py
```

[project.py](project.py): A simple demonstration of using signac-flow to encode a workflow for setting up and performing GROMACS simulations using mBuild and Foyer. Usage (execute each in series):

Initialize and parameterize the system, generating files needed to run gromacs:
```
python project.py exec init 
```

Submit the grommp and mdrun commands to the compute nodes, selecting the partition and hardware:
``` 
python project.py submit -o run --gres=gpu:V100:2 --partition=short-tesla --ntasks=8
```

The options for each flag are shown below.  Note these are included in the Rahman class in project.py.  

```
--gres: [ 'gpu:GTX980:2', 'gpu:A100:2', 'gpu:V100:2', 
          'gpu:GTX980:1', 'gpu:A100:1', 'gpu:V100:1']
--partition: [ 'short-std', 'short-tesla', 
               'day-long-std', 'day-long-tesla', 
               'week-long-std', 'week-long-tesla', 
               'month-long-std', 'month-long-tesla']
--ntasks: ['1', '2', '4', '8', '16']

```

Check on the status of the jobs.
```
python project.py exec check
```

[analyze.py](analyze.py): A very simple framework to demonstrate how to load in your signac project, read in statepoint information, and locate the data for each job that would allow you to perform analysis over the dataspace. No analysis is actually performed here, this is just a simple framework to demonstrate how to access the project info

```
python analyze.py
```
