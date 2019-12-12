#ifndef FUNCTOR_DE_LINK_NODES_H_
#define FUNCTOR_DE_LINK_NODES_H_

#include "SystemStructures.h"


struct functor_de_link_nodes {

	unsigned* globalNeighbors;
	unsigned* currentEdgeCountVec;
	unsigned maxNbrCount;
	double* lengthZero;
	unsigned maxNodeCount;

	__host__ __device__
	functor_de_link_nodes(

		unsigned* _globalNeighbors,
		double* _lengthZero,
		unsigned* _currentEdgeCountVec,
		unsigned _maxNbrCount,
		unsigned _maxNodeCount) :
		globalNeighbors(_globalNeighbors),
		lengthZero(_lengthZero),
		currentEdgeCountVec(_currentEdgeCountVec),
		maxNbrCount(_maxNbrCount),
		maxNodeCount(_maxNodeCount) {}

	__device__
	unsigned operator() (const unsigned& nodeId) {
		unsigned possibleEdgeBegin = maxNbrCount * nodeId;
		unsigned possibleEdgeEnd = maxNbrCount * nodeId + maxNbrCount;//test all possible for now

		unsigned numDelinked=0;

		for (unsigned index = possibleEdgeBegin; index < possibleEdgeEnd; index++) {
		//index is the location of nbrId

			unsigned nbrId = globalNeighbors[index];//loop through all possible neighbors
			if (nbrId < maxNodeCount) {
				bool isNodeNbrOfNodeId = false;//assume not nbr, and check if true.
				unsigned targetPossibleEdgeBegin = maxNbrCount * nbrId;
				unsigned targetPossibleEdgeEnd = maxNbrCount * nbrId + maxNbrCount;
				for (unsigned indexNbr = targetPossibleEdgeBegin; indexNbr < targetPossibleEdgeEnd; indexNbr++ ) {
					if ((globalNeighbors[indexNbr]) == nodeId) {
						isNodeNbrOfNodeId = true;//nbrId is connected to nodeId
						break;
					}
				}
				if (!isNodeNbrOfNodeId) {
					//then nbrId is not connected to nodeId, but nodeId is connected to nbrId.
					//we thus remove nbrId from nodeId.
					//these two mess everything up
					globalNeighbors[index] = maxNodeCount;
					lengthZero[index] = 0.0;
					(currentEdgeCountVec[nodeId])-=1;//subtract from edgecount if nbr removed.
					numDelinked +=1;
				}

			}
		}
		return numDelinked;
	}

};
#endif