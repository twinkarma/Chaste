//
// Created by twin on 30/04/19.
//

#ifndef EDGEREMAPINFO_HPP_
#define EDGEREMAPINFO_HPP_

#include <vector>

/**
 * Storage class contains a mapping to the old local edges index during cell division.
 *
 */
class EdgeRemapInfo {

    /**
     * Contains a mapping to the old local edges index. Negative value means a new edge
     *
     */
    std::vector<long int> mEdgesMapping;

    /**
     * Status
     * 0 Edge has not changed
     * 1 Edge has been split between two elements
     * 2 Completely new edge was created
     */
    std::vector<unsigned char> mEdgeStatus;

public:

    EdgeRemapInfo(std::vector<long int> edgesMapping, std::vector<unsigned char> edgesStatus)
    {
        mEdgesMapping = edgesMapping;
        mEdgeStatus = edgesStatus;
    }

    /**
     * Contains a mapping to the old local edges index. Negative value means a new edge
     *
     */
    std::vector<long int>& GetEdgesMapping()
    {
        return mEdgesMapping;
    }

    /**
     * Status
     * 0 Edge has not changed
     * 1 Edge has been split between two elements
     * 2 Completely new edge was created
     */
    std::vector<unsigned char>& GetEdgesStatus()
    {
        return mEdgeStatus;
    }


};


#endif //EDGEREMAPINFO_HPP_