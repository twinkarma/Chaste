#!/bin/bash
#$ -cwd {{ output_dir }}
#$ -t 1-{{ num_tasks }}
{% for bp in batch_params %}
#$ {{ bp }}
{% endfor %}

cd {{ output_dir }}

python runsimulation.py $SGE_TASK_ID