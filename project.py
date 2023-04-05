import signac
import os
import pathlib
import sys
from flow import FlowProject
import unyt as u

import foyer
from foyer import Forcefield
import mbuild as mb

# Template mdp files are stored in engine_input/gromacs/mdp.
from engine_input.gromacs import mdp

def santize_path(temp_string):
    # on my machine, dropbox integration ends up with my path have spaces and parenthesis
    #make sure we were passed a string
    temp_string = str(temp_string)
    temp_string = temp_string.replace(' ', '\ ')
    temp_string = temp_string.replace('(', '\(')
    temp_string = temp_string.replace(')', '\)')

class Project(FlowProject):
    pass

# To run on a cluster, you may need to define a template for the scheduler
from flow.environment import DefaultSlurmEnvironment
class Rahman(DefaultSlurmEnvironment):
    # Subclass of DefaultSlurmEnvironment for VU's Rahman cluster.
    # The Slurm template are stored in a "templates" folder in the project
    # directory.
    template = "rahman_gmx.sh"
    hostname_pattern = "head.cl.vanderbilt.edu"

    @classmethod
    def add_args(cls, parser):
        # Add command line arguments to the submit call.
        parser.add_argument(
                            "--gres",
                            choices=[
                                'gpu:GTX980:2',
                                'gpu:A100:2',
                                'gpu:V100:2',
                                'gpu:GTX980:1',
                                'gpu:A100:1',
                                'gpu:V100:1',
                            ],
                            default='gpu:GTX980:2',
                            help="which type of gpu",
                            )
        parser.add_argument(
                            "--partition",
                            choices=[
                                'short-std',
                                'short-tesla',
                                'day-long-std',
                                'day-long-tesla',
                                'week-long-std',
                                'week-long-tesla',
                                'month-long-std',
                                'month-long-tesla',
                            ],
                            default='short-std',
                            help="which queue to run",
                            )
        parser.add_argument(
                            "--ntasks",
                            choices=[
                                '1',
                                '2',
                                '4',
                                '8',
                                '16',
                            ],
                            default='16',
                            help="number of cores",
                            )


# This function will read in jinja template, replace variables, and write out the new file
# This is not called directly in the signac project, but called by the init function
# defined below.
def _setup_mdp(fname, template, data, overwrite=False):
    """Create mdp files based on a template and provided data.
        Parameters
        ----------
        fname: str
        Name of the file to be saved out
        template: str, or jinja2.Template
        Either a jinja2.Template or path to a jinja template
        data: dict
        Dictionary storing data matched with the fields available in the template
        overwrite: bool, optional, default=False
        Options to overwrite (or not) existing mdp file of the
        Returns
        -------
        File saved with names defined by fname
        """
    from jinja2 import Template
    
    if isinstance(template, str):
        with open(template, "r") as f:
            template = Template(f.read())

    if not overwrite:
        if os.path.isfile(fname):
            raise FileExistsError(
                                  f"{fname} already exists. Set overwrite=True to write out."
                                  )
        
    rendered = template.render(data)
    with open(fname, "w") as f:
        f.write(rendered)

    return None


# init function to set up the simulation
# This will call mbuild to construct the system, foyer to atom-type,
# save to the appropriate .top and .gro format,
# and then generate the propopogate an .mdp file for GROMACS using
# the thermodynamic variables defined in init.py
# This operation is considered successful if we have generated the .top, .gro, and .mdp files.
@Project.post(lambda j: j.isfile("system_input.top"))
@Project.post(lambda j: j.isfile("system_input.gro"))
@Project.post(lambda j: j.isfile("system_input.mdp"))
@Project.operation(f'init', with_job=True)
def init_job(job):

    # get the root directory so that we can read in the appropriate force field file later
    # and fetch the appropriate .mdp templates
    project_root = job.project

    # fetch the key information related to system structure parameterization
    molecule_string = job.sp.molecule_string
    box_length = job.sp.box_length
    n_molecules = job.sp.n_molecules
    
    # use mbuild to constract a compound and fill a box
    compound = mb.load(molecule_string, smiles=True)
    box = mb.Box(lengths=[box_length, box_length, box_length])
    compound_system = mb.fill_box(compound, n_compounds=n_molecules, box=box)
    
    # atomtype and save the input files to GROMACS format
    compound_system.save(f"system_input.top", forcefield_files=f"{project_root}/../xml_files/oplsaa_alkanes.xml", overwrite=True)
    compound_system.save(f"system_input.gro", overwrite=True)

    
    # fetch run time variables that will be set in the .mdp file
    run_time = job.sp.run_time
    temperature = job.sp.temperature
    velocity_seed = job.sp.velocity_seed
    
    # aggregate info into a simple dictionary
    mdp_abs_path = f'{project_root}/../engine_input/gromacs/mdp'
    mdp = {
        "fname": "system_input.mdp",
        "template": f"{mdp_abs_path}/system.mdp.jinja",
        "data": {
            "run_time": run_time,
            "velocity_seed": velocity_seed,
            "temperature": temperature,
        }
    }
    # call the function that will read in template, perform replacement, and generate the .mdp file
    _setup_mdp(
               fname=mdp["fname"],
               template=mdp["template"],
               data=mdp["data"],
               overwrite=True,
               )


# This function defines the gmx grompp and gmx mdrun commands
# used to pre-process and run the simulation, respectively.
# By using the flow.cmd decorator, the string in the return statement will be executed
# in the same way we would call it at the command line or in a shell script.
# This avoids the need to import, e.g., os and use popen.

# Note the grompp and mdrun calls could certainly be in a separate signac functions,
# and could be desirable for some workflows.
# One advantage to packaging in a single command is that it allows chaining together a sequence
# of simulations, e.g., a simulation workflow with 3 distinct stages,
# where each stage depends on the input from the prior stage.
# This can be done by simply concatenating together the separate msg statements, before returning.
# Although, caution should be taken when making a single really long string,
# as it may overrun the shell can handle (e.g., getting an "Argument list too long" error)

@Project.post(lambda j: j.isfile(f"system.gro"))
@Project.operation(f'run', cmd=True, with_job=True)
def run_job(job):
    
    module_to_load =f"module load gromacs/2020.6"

    grompp = "gmx grompp -f system_input.mdp -o system_input.tpr -c system_input.gro -p system_input.top --maxwarn 2"
    mdrun ="gmx mdrun -v -deffnm system -s system_input.tpr -cpi system.cpt -nt 16"
    
    msg = f"{module_to_load} && {grompp} && {mdrun}"
    print(msg)
    return(msg)
    
# This is a simple function to check to see if the job has completed, writing to the job.doc.
# This will be used in the analysis.py file to ensure that we are only performing analysis
# on simulations that have completed.
@Project.post(lambda j: j.isfile("system.gro"))
@Project.operation(f'check')
def check_job(job):
    molecule_string = job.sp.molecule_string
    if job.isfile(f"system.gro") == True:
        print(f"{molecule_string} ::  completed")
        if not molecule_string in job.doc.get('completed', []):
            job.doc.setdefault('completed', []).append(molecule_string)
    


if __name__ == "__main__":
    pr = Project(path='alkane_screen')
    pr.main()
