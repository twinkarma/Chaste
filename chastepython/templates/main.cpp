/*

Copyright (c) 2005-2019, University of Oxford.
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

*/

/*
 * Note: Do not put any VTK-specific functionality in this file, as we
 * don't ever test it with VTK support turned off!
 */


#include <string>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <iostream>
#include <iterator>
using namespace std;


#include "ExecutableSupport.hpp"

int main(int argc, char *argv[])
{

    try {

        ExecutableSupport::StandardStartup(&argc, &argv);
        int exit_code = ExecutableSupport::EXIT_OK;

        po::options_description desc("Allowed options");
        desc.add_options()
                ("help", "produce help message")
                ("output_dir", "directory to write simulation outputs")
                {% for name in param_names %}
                ("{{ name }}", po::value<double>(), "Set the value of {{ name }}")
                {% endfor %}
                ;

        po::variables_map vm;
        po::store(po::parse_command_line(argc, argv, desc), vm);
        po::notify(vm);

        if (vm.count("help")) {
            cout << desc << "\n";
            return 0;
        }

        // Checking that variables have been defined
        if (!vm.count("output_dir")) {
            cout << "output_dir was not set.\n";
            cout << desc << "\n";
            return 1;
        }
        {% for name in param_names %}
        if (!vm.count("{{ name}}")) {
            cout << "{{ name }} was not set.\n";
            cout << desc << "\n";
            return 1;
        }
        {% endfor %}

        // Collect variables
        auto output_dir = vm["output_dir"];
        {% for name in param_names %}
        double {{ name }} = vm["{{ name }}"].as<double>();
        {% endfor %}

    }
    catch(exception& e) {
        ExecutableSupport::PrintError(e.GetMessage());
        exit_code = ExecutableSupport::EXIT_ERROR;
    }

    ExecutableSupport::FinalizePetsc();
    return exit_code;

}
