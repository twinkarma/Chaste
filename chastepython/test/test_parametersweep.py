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
import unittest
import numpy as np
import json
import os
import shutil
from chastepython.util.parametersweep import ParamSweeper
from chastepython.util.pscan import JointParameterListSizeError


class TestParameterSweeper(unittest.TestCase):

    def test_expansion(self):

        p = {}
        p['a'] = np.linspace(0, 10, 5)
        p['b'] = np.linspace(0.1, -0.5, 5)
        p['c'] = [10, 20, 30]

        sweeper = ParamSweeper()

        # First example a is combinatorial with b
        expanded_params = sweeper.expand_parameters(p)

        self.assertEqual(len(expanded_params), 75)

        # Where a and b are joined i.e. no combinatorial expansion
        joint_lists = [['a', 'b']]
        expanded_params = sweeper.expand_parameters(parameters=p, joint_lists=joint_lists)
        self.assertEqual(len(expanded_params), 15)

        # Joins must have the same dimension
        joint_lists_fail = [['a', 'b', 'c']]
        with self.assertRaises(JointParameterListSizeError):
            sweeper.expand_parameters(parameters=p, joint_lists=joint_lists_fail)

    def test_batch_generation(self):

        p = {}
        p['a'] = np.linspace(0, 10, 5)
        p['b'] = np.linspace(0.1, -0.5, 5)

        # Feed to param sweeper
        sweeper = ParamSweeper()

        exec_cmd = "chastepython/test/test_params.sh"
        output_dir = "/tmp/myoutdir_sge"
        scheduler = ParamSweeper.SGE

        if os.path.exists(output_dir):
            shutil.rmtree(output_dir)

        # Test output for SGE
        sweeper.generate_batch_output(output_dir=output_dir,
                                      exec_cmd=exec_cmd,
                                      parameters=p,
                                      batch_params=["-M myemail@mydomain.com", "-m bes"])
        self.assertTrue(os.path.exists(os.path.join(output_dir, sweeper.params_file_name)))
        self.assertTrue(os.path.exists(os.path.join(output_dir, sweeper.sge_batch_file_name)))
        self.assertTrue(os.path.exists(os.path.join(output_dir, sweeper.python_sim_runner_file_name)))

        # Test output json
        expanded_params = sweeper.expand_parameters(p)
        with open(os.path.join(output_dir, sweeper.params_file_name), "r") as out_file:
            out_json = json.load(out_file)
            self.assertEqual(len(expanded_params), len(out_json["params"]))

        # Test output for SLURM
        output_dir = "/tmp/myoutdir_slurm"
        if os.path.exists(output_dir):
            shutil.rmtree(output_dir)

        sweeper.generate_batch_output(output_dir=output_dir,
                                      exec_cmd=exec_cmd,
                                      parameters=p,
                                      scheduler=ParamSweeper.SLURM,
                                      batch_params=["-M myemail@mydomain.com", "-m bes"])

        self.assertTrue(os.path.exists(os.path.join(output_dir, sweeper.params_file_name)))
        self.assertTrue(os.path.exists(os.path.join(output_dir, sweeper.slurm_batch_file_name)))
        self.assertTrue(os.path.exists(os.path.join(output_dir, sweeper.python_sim_runner_file_name)))

    def test_serial_sweep(self):

        p = {}
        p['a'] = np.linspace(0, 10, 5)
        p['b'] = np.linspace(0.1, -0.5, 5)

        # Feed to param sweeper
        sweeper = ParamSweeper()

        exec_cmd = "chastepython/test/test_params.sh"
        output_dir = "/tmp/serial_sweep"

        if os.path.exists(output_dir):
            shutil.rmtree(output_dir)

        sweeper.perform_serial_sweep(output_dir, exec_cmd, p)

        for i in range(25):
            self.assertTrue(os.path.exists("{}/{}/testout.txt".format(output_dir, i)))

    def test_main_generation(self):

        out_file = "/tmp/test_main.cpp"
        if os.path.exists(out_file):
            os.remove(out_file)

        sweeper = ParamSweeper()
        sweeper.generate_main_cpp(["param1", "param2", "param3"], out_file)

        self.assertTrue(os.path.exists(out_file))
