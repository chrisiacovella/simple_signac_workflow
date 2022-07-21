""" init.py :: Initialize the signac dataspace.
    ------
    This is a simple example to demonstrate the use of
    Signac to perform a simulation using GROMACS.
    This example will perform a simulations at various
    thermodynamics statepoints for a system of bulk hexane.
    We will not use mbuild/foyer directly in this example,
    instead the associated .gro and .top files have
    already been created to simplify the scripts.
    -------
"""
import itertools
import os

import numpy as np
import signac
import unyt as u

# Define a new project and give it a name.
pr = signac.init_project('alkane_screen')

# Define the design space:
# In this case, we will set our design space to be
# the chemical structure, where we increase alkane length

molecule_strings = ['C', 'CC']

# Define statepoint information:
# We of course could move any of these variables into the
# design space so we can loop over, e.g., temperature

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
