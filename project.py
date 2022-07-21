import signac
import os
import pathlib
import sys
import flow
import unyt as u

import foyer
from foyer import Forcefield
import mbuild as mb

# Template mdp files are stored in engine_input/gromacs/mdp.
from engine_input.gromacs import mdp


class Project(flow.FlowProject):
    """Subclass of FlowProject to provide custom methods and attributes."""
    
    def __init__(self):
        super().__init__()
        current_path = pathlib.Path(os.getcwd()).absolute()

# To run on a cluster, you may need to define a template for the scheduler
# For example, below is what I would use on our local group cluster.
# It is currently commented out as it won't be used in the simple
# example here.
"""
from flow.environment import DefaultTorqueEnvironment
class Rahman(DefaultTorqueEnvironment):
    # Subclass of DefaultPBSEnvironment for VU's Rahman cluster.
    # PBS templates are stored in a "templates" folder in your project
    # directory.
    template = "rahman_gmx.sh"
    
    @classmethod
    def add_args(cls, parser):
        # Add command line arguments to the submit call.
        parser.add_argument(
                            "--walltime",
                            type=float,
                            default=96,
                            help="Walltime for this submission",
                            )

"""
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


@Project.operation(f'init')
@Project.post(lambda j: j.isfile("system_input.top"))
@Project.post(lambda j: j.isfile("system_input.gro"))
@Project.post(lambda j: j.isfile("system_input.mdp"))

# The with_job decorator basically states that this function accepts
# a single job as a parameter
@flow.with_job
def init_job(job):

    #fetch the key information related to system structure parameterization
    molecule_string = job.sp.molecule_string
    project_root = Project().root_directory()
    box_length = job.sp.box_length
    n_molecules = job.sp.n_molecules
    
    #use mbuild to constract a compound and fill a box
    compound = mb.load(molecule_string, smiles=True)
    box = mb.Box(lengths=[box_length, box_length, box_length])
    compound_system = mb.fill_box(compound, n_compounds=n_molecules, box=box)
    
    #atomtype and save the input files to GROMACS format
    compound_system.save(f"system_input.top", forcefield_files=f"{project_root}/xml_files/oplsaa_alkanes.xml", overwrite=True)
    compound_system.save(f"system_input.gro", overwrite=True)

    
    #fetch run time variables
    run_time = job.sp.run_time
    temperature = job.sp.temperature
    velocity_seed = job.sp.velocity_seed
    
    #set up mdp files
    mdp_abs_path = Project().root_directory() + '/engine_input/gromacs/mdp'
    mdp = {
        "fname": "system_input.mdp",
        "template": f"{mdp_abs_path}/system.mdp.jinja",
        "data": {
            "run_time": run_time,
            "velocity_seed": velocity_seed,
            "temperature": temperature,
        }
    }
    _setup_mdp(
               fname=mdp["fname"],
               template=mdp["template"],
               data=mdp["data"],
               overwrite=True,
               )

# This bit will define  the gmx grompp command and gmx mdrun command
# by the flow.cmd decorator, the string in the return statement will be executed
@Project.operation(f'run')
@Project.post(lambda j: j.isfile(f"system.gro"))
@flow.with_job
@flow.cmd
def run_job(job):

    grompp = f"gmx grompp -f system_input.mdp -o system_input.tpr -c system_input.gro -p system_input.top --maxwarn 2"
    mdrun ="gmx mdrun -v -deffnm system -s system_input.tpr -cpi system.cpt -nt 16"
    
    msg = f"{grompp} && {mdrun}"
    print(msg)
    return(msg)
    
@Project.operation(f'check')
@Project.post(lambda j: j.isfile("system.gro"))
@flow.with_job
def check_job(job):
    molecule_string = job.sp.molecule_string
    if job.isfile(f"system.gro") == True:
        print(f"{molecule_string} ::  completed")
        if not molecule_string in job.doc.get('completed', []):
            job.doc.setdefault('completed', []).append(molecule_string)
    


if __name__ == "__main__":
    pr = Project()
    pr.main()
