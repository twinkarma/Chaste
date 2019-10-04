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


/**
 * Concrete FieldNoyesReactionSystem class
 */
#ifndef FIELDNOYESREACTIONSYSTEM_HPP_
#define FIELDNOYESREACTIONSYSTEM_HPP_

#include "AbstractOdeSystem.hpp"
#include "OdeSystemInformation.hpp"

class FieldNoyesReactionSystem : public AbstractOdeSystem
{
public:
    FieldNoyesReactionSystem() : AbstractOdeSystem(3)
    {
        mpSystemInfo = OdeSystemInformation<FieldNoyesReactionSystem>::Instance();
    }

    void EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY)
    {
        const double epsilon = 0.05;
        const double p = 6.7;
        const double q = 1e-4;
        const double f = 0.5;

        rDY[0] = (rY[1] - rY[0]*rY[1] + rY[0]*(1 - q*rY[0]))/epsilon;
        rDY[1] = -rY[1] - rY[0]*rY[1] + 2*f*rY[2];
        rDY[2] = (rY[0] - rY[2])/p;
    }
};

template<>
void OdeSystemInformation<FieldNoyesReactionSystem>::Initialise()
{
    this->mVariableNames.push_back("Variable_1");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(1.0);

    this->mVariableNames.push_back("Variable_2");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(1.0);

    this->mVariableNames.push_back("Variable_3");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(1.0);

    this->mInitialised = true;
}

#endif /*FIELDNOYESREACTIONSYSTEM_HPP_*/
