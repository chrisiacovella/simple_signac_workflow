""" init.py :: Initialize the signac dataspace.
    ------
    This is a simple example to demonstrate the use of
    MoSDeF and Signac to perform a simulation using GROMACS.
    This file defines the dataspace that will be considered.
    Note, each job in the dataspace will have the same
    thermodynamic inputs, but will be initialized
    with different length alkanes.
    -------
"""
import itertools
import os

import numpy as np
import signac
import unyt as u

# Define a new project
pr = signac.init_project('alkane_screen')

# Define the design space:
# In this case, we will set our design space to be
# the chemical structure, where we increase alkane length.
# For convenience, since are only considering
# short, linear alkaneswe will define the
# molecule dataspace using the SMILES strings.

molecule_strings = ['C', 'CC', 'CCC']

# Define statepoint information:
# We could move any of these variables into the
# design space so we could consider
# e.g., the behavior as a function of temperature
# in addition to molecule structure.

n_molecules = 100
temperature = 300* u.K
velocity_seed = 8675309
run_time = 10000
box_length = 10*u.nm


# Loop over the design space to create an array of statepoints
total_statepoints = []
for molecule_string in molecule_strings:
    statepoint = {
        "molecule_string": molecule_string,
        "temperature": temperature.to_value("K"),
        "velocity_seed": velocity_seed,
        "run_time": run_time,
        "n_molecules": n_molecules,
        "box_length": box_length.to_value("nm"),

    }
    total_statepoints.append(statepoint)

# Initialize the statepoints.
# Since we only are defining small handful,
# we will also print them to the screen.
print('statepoints initialized:')
for sp in total_statepoints:
    print(sp)
    pr.open_job(
            statepoint=sp,
            ).init()
