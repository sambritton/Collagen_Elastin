#ifndef FUNCTOR_DE_LINK_NODES_H_
#define FUNCTOR_DE_LINK_NODES_H_

#include "system_structures.h"


struct functor_de_link_nodes {

	unsigned* global_neighbors;
	unsigned* currentEdgeCountVec;
	unsigned max_nbr_count;
	double* global_length_zero;
	unsigned max_node_count;

	__host__ __device__
	functor_de_link_nodes(

		unsigned* _global_neighbors,
		double* _global_length_zero,
		unsigned* _currentEdgeCountVec,
		unsigned _max_nbr_count,
		unsigned _max_node_count) :
		global_neighbors(_global_neighbors),
		global_length_zero(_global_length_zero),
		currentEdgeCountVec(_currentEdgeCountVec),
		max_nbr_count(_max_nbr_count),
		max_node_count(_max_node_count) {}

	__device__
	unsigned operator() (const unsigned& nodeId) {
		unsigned possibleEdgeBegin = max_nbr_count * nodeId;
		unsigned possibleEdgeEnd = max_nbr_count * nodeId + max_nbr_count;//test all possible for now

		unsigned numDelinked=0;

		for (unsigned index = possibleEdgeBegin; index < possibleEdgeEnd; index++) {
		//index is the location of nbrId

			unsigned nbrId = global_neighbors[index];//loop through all possible neighbors
			if (nbrId < max_node_count) {
				bool isNodeNbrOfNodeId = false;//assume not nbr, and check if true.
				unsigned targetPossibleEdgeBegin = max_nbr_count * nbrId;
				unsigned targetPossibleEdgeEnd = max_nbr_count * nbrId + max_nbr_count;
				for (unsigned indexNbr = targetPossibleEdgeBegin; indexNbr < targetPossibleEdgeEnd; indexNbr++ ) {
					if ((global_neighbors[indexNbr]) == nodeId) {
						isNodeNbrOfNodeId = true;//nbrId is connected to nodeId
						break;
					}
				}
				if (!isNodeNbrOfNodeId) {
					//then nbrId is not connected to nodeId, but nodeId is connected to nbrId.
					//we thus remove nbrId from nodeId.
					//these two mess everything up
					global_neighbors[index] = max_node_count;
					global_length_zero[index] = 0.0;
					(currentEdgeCountVec[nodeId])-=1;//subtract from edgecount if nbr removed.
					numDelinked +=1;
				}

			}
		}
		return numDelinked;
	}

};
#endif