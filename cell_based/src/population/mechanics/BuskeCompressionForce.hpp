/*

Copyright (c) 2005-2013, University of Oxford.
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

#ifndef BUSKECOMPRESSIONFORCE_HPP_
#define BUSKECOMPRESSIONFORCE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractForce.hpp"
#include "TetrahedralMesh.hpp"

/**
 * A force law employed by Buske et al (2011) in their overlapping spheres
 * model of the intestinal crypt (doi:10.1371/journal.pcbi.1001045).
 *
 * Length is scaled by natural length. \todo does this mean natural radius of a cell? If so at what age? (#1764)
 * Time is in hours.
 *
 * This class specifically calculates the compression force which forms part of equation (A6) in the Buske paper.
 */
template<unsigned DIM>
class BuskeCompressionForce : public AbstractForce<DIM>
{
    friend class TestForcesNotForRelease;
private:

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractForce<DIM> >(*this);
        archive & mCompressionEnergyParameter;
    }

    /**
     * Compression energy parameter.
     *
     * Represented by the parameter K in the model by Buske et al (2011) in
     * their off-lattice model of the intestinal crypt
     * (doi:10.1371/journal.pcbi.1001045).
     *
     * Note: K is the bulk modulus of the spheres.
     */
    double mCompressionEnergyParameter;

public:

    /**
     * Constructor.
     */
    BuskeCompressionForce();

    /**
     * Get mCompressionEnergyParameter.
     */
    double GetCompressionEnergyParameter();

    /**
     * Set mCompressionEnergyParameter.
     *
     * @param compressionEnergyParameter the new value of mCompressionEnergyParameter
     */
    void SetCompressionEnergyParameter(double compressionEnergyParameter);

    /**
     * Overridden AddForceContribution() method.
     *
     * @param rCellPopulation a cell population object
     */
    void AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation);

    /**
     * Overridden OutputForceParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputForceParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(BuskeCompressionForce)

#endif /*BUSKECOMPRESSIONFORCE_HPP_*/
