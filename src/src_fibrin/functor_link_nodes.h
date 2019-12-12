#ifndef FUNCTOR_LINK_NODES_H_
#define FUNCTOR_LINK_NODES_H_

#include "SystemStructures.h"


struct functor_link_nodes {
	double* nodeLocXAddr;
	double* nodeLocYAddr;
	double* nodeLocZAddr;
	unsigned* globalNeighbors;
	unsigned* currentEdgeCountVec;
	double*	lengthZero;

	unsigned* id_value_expanded;
	unsigned* keyBegin;
	unsigned* keyEnd;

	double fiberDiameter;
	unsigned maxNbrCount;
	unsigned maxNodeCount;

	unsigned maxLinksPerIteration;
	unsigned* idTempVecLeft;
	unsigned* idTempVecRight;

	__host__ __device__
		//
		functor_link_nodes(
			double* _nodeLocXAddr,
			double* _nodeLocYAddr,
			double* _nodeLocZAddr,
			unsigned* _globalNeighbors,
			unsigned* _currentEdgeCountVec,
			double* _lengthZero,

			unsigned* _id_value_expanded,
			unsigned* _keyBegin,
			unsigned* _keyEnd,

			double& _fiberDiameter,
			unsigned& _maxNbrCount,
			unsigned& _maxNodeCount,

			unsigned& _maxLinksPerIteration,
			unsigned* _idTempVecLeft,
			unsigned* _idTempVecRight) :

		nodeLocXAddr(_nodeLocXAddr),
		nodeLocYAddr(_nodeLocYAddr),
		nodeLocZAddr(_nodeLocZAddr),
		globalNeighbors(_globalNeighbors),
		currentEdgeCountVec(_currentEdgeCountVec),
		lengthZero(_lengthZero),

		id_value_expanded(_id_value_expanded),
		keyBegin(_keyBegin),
		keyEnd(_keyEnd),

		fiberDiameter(_fiberDiameter),
		maxNbrCount(_maxNbrCount),
		maxNodeCount(_maxNodeCount),

		maxLinksPerIteration(_maxLinksPerIteration),
		idTempVecLeft(_idTempVecLeft),
		idTempVecRight(_idTempVecRight) {}

	__device__
		
	unsigned operator() (const Tuu& u2) {
		unsigned nodeId = thrust::get<0>(u2);//node to attempt link from.
		unsigned bucketId = thrust::get<1>(u2);//bucket containing nodeId
		
		unsigned final_id_left = 0;
		unsigned final_id_right = 0;

		unsigned possibleEdgeBegin = maxNbrCount * nodeId;
		unsigned possibleEdgeEnd = maxNbrCount * nodeId + maxNbrCount;//test all possible for now
		unsigned numPlacedLinks = 0;

		unsigned save_index = maxLinksPerIteration * nodeId;//changes after adding an edge
		unsigned save_index_end = maxLinksPerIteration * nodeId + maxLinksPerIteration;


		if ( nodeId < maxNodeCount )  {
			//tsafety

			//beginning and end of attempted linking end id's in id_value_expanded
			unsigned beginIndex = keyBegin[bucketId];
			unsigned endIndex = keyEnd[bucketId];

			for (unsigned iter = beginIndex; iter < endIndex; iter++ ) {

				unsigned candidateId = id_value_expanded[iter];//test id

				if ( (candidateId < maxNodeCount) && (save_index < save_index_end) ) {
					//then candidateId is not a dpd particle.
					bool candidateIsNew = true;

					//we now attempt to link nodeId and candidateEdge.
					if ((candidateId != nodeId) ){//&& (currentEdgeCountVec[nodeId] < maxNbrCount)) {
						double dist = sqrt(
							((nodeLocXAddr[nodeId] - nodeLocXAddr[candidateId]) * (nodeLocXAddr[nodeId] - nodeLocXAddr[candidateId])) +
							((nodeLocYAddr[nodeId] - nodeLocYAddr[candidateId]) * (nodeLocYAddr[nodeId] - nodeLocYAddr[candidateId])) +
							((nodeLocZAddr[nodeId] - nodeLocZAddr[candidateId]) * (nodeLocZAddr[nodeId] - nodeLocZAddr[candidateId])));

						if ((dist < fiberDiameter)) {
							//then we have a possible link if no previous link was placed.

							//make sure placement is possible for nodeID
							//other thread will take care of candidateID
							for (unsigned k = possibleEdgeBegin; k < possibleEdgeEnd; k++) {
								if (candidateId == (globalNeighbors[k])) {
									candidateIsNew = false;
									break;
								}
							}

							if (candidateIsNew == true) {

								//must regenerate for each new candidate
								//unsigned currentEdgeCount = currentEdgeCountVec[nodeId];

								//then the new edge can be placed
								unsigned location;
								//decide placement of nbr only in row of nodeId,
								//another thread will place for candidateId
								for (unsigned place = possibleEdgeBegin; place < possibleEdgeEnd; place++ ) {
									if (maxNodeCount == globalNeighbors[place]) {
										//maxNodeCount is the default value of the matrix
										//then we have found location where we can place a node scine global neighbors default is maxNodeCount
										location = place;
										break;
									}
								}
								globalNeighbors[ location ] = candidateId;
								lengthZero[ location ] = dist;
								//only if a node is placed do we record it.
	
								final_id_left = min(nodeId, candidateId);
								final_id_right = max(nodeId, candidateId);
								idTempVecLeft[save_index] = final_id_left;
								idTempVecRight[save_index] = final_id_right;
								(currentEdgeCountVec[nodeId])+=1;
								save_index += 1;
								numPlacedLinks += 1;
								
							}
						}
					}
				}
			}
		}
		return (numPlacedLinks);
	}
};


#endif