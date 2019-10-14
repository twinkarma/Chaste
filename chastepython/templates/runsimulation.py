#!/usr/bin/env python
"""Copyright (c) 2005-2019, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""
from __future__ import print_function
import os
import sys
import json
import subprocess

if len(sys.argv) < 2:
    print("A simulation ID must be provided as an argument, e.g. : \n runsimulation.py $SGE_TASK_ID")
    sys.exit(1)

try:
    task_id = int(sys.argv[1])
    with open("params.json", "r") as params_file:
        # Load and collect simulation parameters
        exec_params = json.load(params_file)
        exec_cmd = exec_params["exec_cmd"]
        output_dir = exec_params["output_dir"]
        iteration_param = exec_params["params"][task_id - 1]

        # Create a folder to store simulation results
        simulation_instance_output_dir = os.path.join(output_dir, str(task_id))
        if os.path.exists(simulation_instance_output_dir):
            print("Output directory for simulation id {} already exists, aborting".format(task_id))
            sys.exit(1)
        else:
            os.mkdir(simulation_instance_output_dir)

        # Builds the exec command with parameters and run the simulation
        iteration_param_string = ""
        for key, value in iteration_param.items():
            iteration_param_string = iteration_param_string + " {}={}".format(key, value)
        final_cmd = "{} output_dir={}{}".format(exec_cmd, simulation_instance_output_dir, iteration_param_string)

        print("Running simulation ID {}, outputting to {}".format(task_id, simulation_instance_output_dir))
        sys.exit(subprocess.call(final_cmd, shell=True))


except TypeError as e:
    print("Simulation ID must be an integer")
    sys.exit(1)
except IOError as e:
    print("Could not open parameters file")
    sys.exit(1)
except OSError as e:
    print("OS error")
    sys.exit(1)