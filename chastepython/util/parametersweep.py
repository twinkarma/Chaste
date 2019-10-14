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
import json
import subprocess
import numpy as np
from jinja2 import Environment, PackageLoader, select_autoescape
from .pscan import Scan


class ParamSweeper:
    """
    Class to help execute the parameter sweeping process locally and on the cluster (SGE and SLURM support).
    """

    SGE = 0
    SLURM = 1

    def __init__(self):

        self.params_file_name = "params.json"
        self.sge_batch_file_name = "batch.sge.sh"
        self.slurm_batch_file_name = "batch.slurm.sh"
        self.python_sim_runner_file_name = "runsimulation.py"



    def expand_parameters(self, parameters, joint_lists=[], default_repeats=1, count_funcs=[]):
        """
        Generate an expanded list of parameters from the instance variable.
        :param parameters:
        :param joint_lists:
        :param default_repeats:
        :param count_funcs:
        :return: Array of expanded parameters
        """


        scan = Scan(parameters, joint_lists,
                    default_repeats, count_funcs)

        expanded_output = []

        def add_paramdict_to_output(**kwargs):
            iteration_params = {}
            for kw in kwargs:
                iteration_params[kw] = kwargs[kw]

            expanded_output.append(iteration_params)

        scan.run_scan(add_paramdict_to_output)

        return expanded_output

    def generate_batch_output(self, output_dir, exec_cmd, parameters, scheduler=SGE, joint_lists=[], default_repeats=1, count_funcs=[], batch_params=[]):
        """
        Generate output files needed to run parameter sweep in batch mode
        :param output_dir:
        :param exec_cmd:
        :param parameters:
        :param scheduler:
        :param joint_lists:
        :param default_repeats:
        :param count_funcs:
        :param batch_params:
        :return:
        """


        if output_dir is None:
            raise ValueError("Output directory not specified")

        if exec_cmd is None:
            raise ValueError("Must specify an executable command")

        expanded_output = self.expand_parameters(parameters, joint_lists, default_repeats, count_funcs)
        env = Environment(loader=PackageLoader('chastepython', 'templates'))

        if not os.path.exists(output_dir):
            os.mkdir(output_dir)

        params_output = {"params": expanded_output,
                         "exec_cmd": os.path.abspath(exec_cmd),
                         "output_dir": os.path.abspath(output_dir)}


        json_output_path = os.path.join(output_dir, self.params_file_name)
        sge_batch_output_path = os.path.join(output_dir, self.sge_batch_file_name)
        slurm_batch_output_path = os.path.join(output_dir, self.slurm_batch_file_name)
        python_sim_runner_output_path = os.path.join(output_dir, self.python_sim_runner_file_name)


        # Output json
        with open(json_output_path, 'w') as json_out_file:
            json.dump(params_output, json_out_file)

        context = {
            "num_tasks": len(expanded_output),
            "exec_cmd": os.path.abspath(exec_cmd),
            "output_dir": os.path.abspath(output_dir),
            "batch_params": batch_params
        }

        if scheduler == ParamSweeper.SGE :
            with open(sge_batch_output_path, "w") as sge_batch_file:
                sge_batch_file.write(env.get_template("batch.sge.sh").render(context))
        elif scheduler == ParamSweeper.SLURM:
            with open(slurm_batch_output_path, "w") as slurm_batch_file:
                slurm_batch_file.write(env.get_template("batch.slurm.sh").render(context))
        else:
            raise Exception("Unsupported scheduler {}".format(scheduler))

        with open(python_sim_runner_output_path, "w") as simrunner_file:
            simrunner_file.write(env.get_template("runsimulation.py").render(context))

    def perform_serial_sweep(self, output_dir, exec_cmd, parameters, joint_lists=[], default_repeats=1, count_funcs=[]):
        """
        Runs the sweep serially
        :param output_dir:
        :param exec_cmd:
        :param parameters:
        :param joint_lists:
        :param default_repeats:
        :param count_funcs:
        :return:
        """

        if output_dir is None:
            raise ValueError("Output directory not specified")

        if not os.path.exists(output_dir):
            os.mkdir(output_dir)

        expanded_output = self.expand_parameters(parameters, joint_lists, default_repeats, count_funcs)
        num_iterations = len(expanded_output)
        for i in range(num_iterations):
            iteration_param = expanded_output[i]

            # Create a folder to store simulation results
            simulation_instance_output_dir = os.path.join(output_dir, str(i))
            if os.path.exists(simulation_instance_output_dir):
                print("Output directory for simulation id {} already exists, aborting".format(i))
            else:
                os.mkdir(simulation_instance_output_dir)

                # Builds the exec command with parameters and run the simulation
                iteration_param_string = ""
                for key, value in iteration_param.items():
                    iteration_param_string = iteration_param_string + " {}={}".format(key, value)
                final_cmd = "{} output_dir={}{}".format(exec_cmd, simulation_instance_output_dir,
                                                        iteration_param_string)

                print("Running simulation ID {}, outputting to {}".format(i, simulation_instance_output_dir))
                subprocess.call(final_cmd, shell=True)



    def generate_main_cpp(self, param_name_list, output_path):
        """
        Generates a main.cpp file from the list of parameter names
        :param param_name_list: A list of parameter names
        :param output_path: Path of the output file
        :return:
        """

        env = Environment(loader=PackageLoader('chastepython', 'templates'))
        context = {
            "param_names": param_name_list
        }

        with open(output_path, "w") as main_cpp_file:
            main_cpp_file.write(env.get_template("main.cpp").render(context))


        print("*****************************\n")
        print("**** Code output below ******\n\n")
        print(env.get_template("main.cpp").render(context))