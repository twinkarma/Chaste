#!/bin/bash
#SBATCH --array=1-{{ num_tasks }}
#SBATCH --chdir {{ output_dir }}
{% for bp in batch_params %}
#SBATCH {{ bp }}
{% endfor %}

cd {{ output_dir }}

python runsimulation.py $SLURM_ARRAY_TASK_ID