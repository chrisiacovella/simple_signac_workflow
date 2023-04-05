""" analyze.py :: Initialize the signac dataspace.
    ------
    This is a simple example to demonstrate how to read in
    a Signac project.
    It doesn't actually really do any analysis, but shows how
    to use the job.doc to make sure you are only loading completed
    simulations and shows how to locate your data for analysis.
    -------
"""
import signac

project_local = signac.get_project()

for job in project_local:
    molecule_string = job.sp.molecule_string
    if molecule_string in job.doc.get('completed', []):
        print(f"{molecule_string} : completed")
        # can easily fetch the job_path to be able to access
        # and analyze your data
        job_path = job.path
        print(job_path)
