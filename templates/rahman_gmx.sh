{% extends "base_script.sh" %}
{% block header %}
#!/bin/bash -l
#PBS -j oe
#PBS -l nodes=1:ppn=16
#PBS -l walltime={{ walltime|format_timedelta }}
#PBS -q standard
#PBS -m abe
#PBS -M chris.iacovella@gmail.com

#module load anaconda/v3.7
#source activate signac
module load gromacs/2020.6

{% endblock header %}
