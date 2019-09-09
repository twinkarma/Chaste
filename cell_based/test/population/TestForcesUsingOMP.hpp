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

#ifndef TESTFORCESUSINGOMP_HPP_
#define TESTFORCESUSINGOMP_HPP_

#include <cxxtest/TestSuite.h>

#include "CheckpointArchiveTypes.hpp"
#include "FarhadifarForce.hpp"
#include "WildTypeCellMutationState.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "CellLabel.hpp"
#include "SmartPointers.hpp"
#include "FileComparison.hpp"
#include "SimpleTargetAreaModifier.hpp"
#include "OffLatticeSimulation.hpp"

#include "PetscSetupAndFinalize.hpp"

#include "omp.h"
#include "FarhadifarForceOMP.hpp"

class TestForcesUsingOMP : public AbstractCellBasedTestSuite
{
public:

    void TestFarhadifarForceOMP()
    {
        std::vector<c_vector<double, 2> > applied_forces(6);
        {
            /**
             * Here we test that the forces are applied correctly to individual nodes.
             * We apply the force to something like this:
             *  . ____ . ____ .
             *  |      |      |
             *  |      |      |
             *  . ____ . ____ .
             */
            std::vector<Node<2>*> nodes;
            // the boolean says wether the node is a boundary node or not
            nodes.push_back(new Node<2>(0, true, 0.0, 0.0));
            nodes.push_back(new Node<2>(1, true, 2.0, 0.0));
            nodes.push_back(new Node<2>(2, true, 4.0, 0.0));
            nodes.push_back(new Node<2>(3, true, 4.0, 2.0));
            nodes.push_back(new Node<2>(4, true, 2.0, 2.0));
            nodes.push_back(new Node<2>(5, true, 0.0, 2.0));

            // make two square elements out of these nodes
            std::vector<Node<2>*> nodes_elem_0, nodes_elem_1;
            unsigned node_indices_elem_0[4] = {0, 1, 4, 5};
            unsigned node_indices_elem_1[4] = {1, 2, 3, 4};

            for (unsigned i=0; i<4; i++)
            {
                nodes_elem_0.push_back(nodes[node_indices_elem_0[i]]);
                nodes_elem_1.push_back(nodes[node_indices_elem_1[i]]);
            }

            std::vector<VertexElement<2,2>*> vertex_elements;
            vertex_elements.push_back(new VertexElement<2,2>(0, nodes_elem_0));
            vertex_elements.push_back(new VertexElement<2,2>(1, nodes_elem_1));

            // Make a vertex mesh
            MutableVertexMesh<2,2> vertex_mesh(nodes, vertex_elements);

            TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 2u);
            TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 6u);

            // Get a cell population
            CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
            std::vector<CellPtr> cells;
            cells_generator.GenerateBasic(cells, vertex_mesh.GetNumElements(), std::vector<unsigned>());
            VertexBasedCellPopulation<2> cell_population(vertex_mesh, cells);

            // Set the birth time to -5 such that the target area modifier assigns mature cell target areas
            for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
                    cell_iter != cell_population.End();
                    ++cell_iter)
            {
                cell_iter->SetBirthTime(-5.0);
            }

            MAKE_PTR(SimpleTargetAreaModifier<2>,p_growth_modifier);
            p_growth_modifier->UpdateTargetAreas(cell_population);

            // Now let's make a FarhadifarForce and apply it to the population
            FarhadifarForce<2> force;

            force.AddForceContribution(cell_population);

            for (unsigned int i=0; i<applied_forces.size(); ++i)
            {
                applied_forces[i] = cell_population.rGetMesh().GetNode(i)->rGetAppliedForce();
            }

            c_vector<double, 2> applied_force_0;
            c_vector<double, 2> applied_force_1;
            applied_force_0 = cell_population.rGetMesh().GetNode(0)->rGetAppliedForce();
            applied_force_1 = cell_population.rGetMesh().GetNode(1)->rGetAppliedForce();

            // If this is a Farhadifar force, this will be the force at the vertices
            TS_ASSERT_DELTA(applied_force_0[0], 3.44, 1e-10);
            TS_ASSERT_DELTA(applied_force_0[1], 3.44, 1e-10);
            TS_ASSERT_DELTA(applied_force_1[0], 0.0, 1e-10);
            TS_ASSERT_DELTA(applied_force_1[1], 6.76, 1e-10);
        }
        std::vector<c_vector<double, 2> > applied_forces_OMP(6);
        {
            /**
             * Here we test that the forces are applied correctly to individual nodes.
             * We apply the force to something like this:
             *  . ____ . ____ .
             *  |      |      |
             *  |      |      |
             *  . ____ . ____ .
             */
            std::vector<Node<2>*> nodes;
            // the boolean says wether the node is a boundary node or not
            nodes.push_back(new Node<2>(0, true, 0.0, 0.0));
            nodes.push_back(new Node<2>(1, true, 2.0, 0.0));
            nodes.push_back(new Node<2>(2, true, 4.0, 0.0));
            nodes.push_back(new Node<2>(3, true, 4.0, 2.0));
            nodes.push_back(new Node<2>(4, true, 2.0, 2.0));
            nodes.push_back(new Node<2>(5, true, 0.0, 2.0));

            // make two square elements out of these nodes
            std::vector<Node<2>*> nodes_elem_0, nodes_elem_1;
            unsigned node_indices_elem_0[4] = {0, 1, 4, 5};
            unsigned node_indices_elem_1[4] = {1, 2, 3, 4};

            for (unsigned i=0; i<4; i++)
            {
                nodes_elem_0.push_back(nodes[node_indices_elem_0[i]]);
                nodes_elem_1.push_back(nodes[node_indices_elem_1[i]]);
            }

            std::vector<VertexElement<2,2>*> vertex_elements;
            vertex_elements.push_back(new VertexElement<2,2>(0, nodes_elem_0));
            vertex_elements.push_back(new VertexElement<2,2>(1, nodes_elem_1));

            // Make a vertex mesh
            MutableVertexMesh<2,2> vertex_mesh(nodes, vertex_elements);

            TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 2u);
            TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 6u);

            // Get a cell population
            CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
            std::vector<CellPtr> cells;
            cells_generator.GenerateBasic(cells, vertex_mesh.GetNumElements(), std::vector<unsigned>());
            VertexBasedCellPopulation<2> cell_population(vertex_mesh, cells);

            // Set the birth time to -5 such that the target area modifier assigns mature cell target areas
            for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
                    cell_iter != cell_population.End();
                    ++cell_iter)
            {
                cell_iter->SetBirthTime(-5.0);
            }

            MAKE_PTR(SimpleTargetAreaModifier<2>,p_growth_modifier);
            p_growth_modifier->UpdateTargetAreas(cell_population);

            // Now let's make a FarhadifarForce and apply it to the population
            FarhadifarForceOMP<2> force;

            force.AddForceContribution(cell_population);

            for (unsigned int i=0; i<applied_forces.size(); ++i)
            {
                applied_forces_OMP[i] = cell_population.rGetMesh().GetNode(i)->rGetAppliedForce();
            }
            c_vector<double, 2> applied_force_0;
            c_vector<double, 2> applied_force_1;
            applied_force_0 = cell_population.rGetMesh().GetNode(0)->rGetAppliedForce();
            applied_force_1 = cell_population.rGetMesh().GetNode(1)->rGetAppliedForce();

            // If this is a Farhadifar force, this will be the force at the vertices
            TS_ASSERT_DELTA(applied_force_0[0], 3.44, 1e-10);
            TS_ASSERT_DELTA(applied_force_0[1], 3.44, 1e-10);
            TS_ASSERT_DELTA(applied_force_1[0], 0.0, 1e-10);
            TS_ASSERT_DELTA(applied_force_1[1], 6.76, 1e-10);
        }
        for (unsigned int i=0; i<applied_forces.size(); ++i)
        {
            TS_ASSERT_DELTA(applied_forces[i][0], applied_forces_OMP[i][0], 1e-12);
            TS_ASSERT_DELTA(applied_forces[i][1], applied_forces_OMP[i][1], 1e-12);
        }
    }
};
#endif /*TESTFORCESUSINGOMP_HPP_*/
