//
// Created by twin on 25/01/19.
//

#ifndef EDGE_HPP_
#define EDGE_HPP_

#include <set>
#include <vector>

#include "Node.hpp"
#include "VertexElement.hpp"

/**
 * An edge in a finite element mesh
 */
template<unsigned SPACE_DIM>
class Edge {

private:

    /** Index of this edge within the mesh **/
    unsigned mIndex;

    /** Nodes that form this edge **/
    std::vector<Node<SPACE_DIM>*> mNodes;


    /** Elements that this edge belongs to **/
    std::set<unsigned> mElementIndices;

public:

    Edge(unsigned index)
    {
        this->mIndex = index;
    }

    void SetNodes(Node<SPACE_DIM>* node0, Node<SPACE_DIM>* node1)
    {
        mNodes.clear();
        mNodes.push_back(node0);
        mNodes.push_back(node1);
    }

    unsigned GetNumNodes()
    {
        return mNodes.size();
    }

    void AddElement(unsigned elementIndex)
    {
        mElementIndices.insert(elementIndex);
    }

    void DeleteElement(unsigned elementIndex)
    {
        mElementIndices.erase(elementIndex);
    }

    unsigned GetNumElements()
    {
        return mElementIndices.size();
    }

    bool IsEdgeValid()
    {
        //MUST have 2 existing nodes to form an edge
        if(mNodes.size() != 2)
            return false;


        //Nodes should not be nullptr
        for(auto node: mNodes)
        {
            if(node == nullptr)
                return false;
        }

        //Can't have associated elements if we're less than 2D
        if(SPACE_DIM <= 1 && mElementIndices > 0){
            return false;
        }

        //An ege can only have a maximum of two elements in 2D
        if(SPACE_DIM == 2 && mElementIndices > 2){
            return false;
        }


        return true;
    }



};


#endif //EDGE_HPP_
