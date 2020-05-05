#ifndef FUNCTOR_LINK_NODES_H_
#define FUNCTOR_LINK_NODES_H_

#include "system_structures.h"


struct functor_link_nodes {
	double* nodeLocXAddr;
	double* nodeLocYAddr;
	double* nodeLocZAddr;
	bool* node_is_collagen_addr;
	bool* node_is_elastin_addr;

	unsigned* currentEdgeCountVec;
	unsigned* global_neighbors;
	double*	global_length_zero;
	bool* global_isedge_collagen;
	bool* global_isedge_elastin;

	unsigned* id_value_expanded;
	unsigned* key_begin_net_intc;
	unsigned* key_end_net_intc;

	double collagen_diameter;
	double elastin_diameter;
	unsigned max_nbr_count;
	unsigned max_node_count;

	unsigned max_links_per_iteration;
	unsigned* idTempVecLeft;
	unsigned* idTempVecRight;

	__host__ __device__
		//
		functor_link_nodes(
			double* _nodeLocXAddr,
			double* _nodeLocYAddr,
			double* _nodeLocZAddr,
			bool* _node_is_collagen_addr,
			bool* _node_is_elastin_addr,

			unsigned* _currentEdgeCountVec,
			unsigned* _global_neighbors,
			double* _global_length_zero,
			bool* _global_isedge_collagen,
			bool* _global_isedge_elastin,

			unsigned* _id_value_expanded,
			unsigned* _key_begin,
			unsigned* _key_end,

			double& _collagen_diameter,
			double& _elastin_diameter,
			unsigned& _max_nbr_count,
			unsigned& _max_node_count,

			unsigned& _max_links_per_iteration,
			unsigned* _idTempVecLeft,
			unsigned* _idTempVecRight) :

		nodeLocXAddr(_nodeLocXAddr),
		nodeLocYAddr(_nodeLocYAddr),
		nodeLocZAddr(_nodeLocZAddr),
		node_is_collagen_addr(_node_is_collagen_addr),
		node_is_elastin_addr(_node_is_elastin_addr),

		currentEdgeCountVec(_currentEdgeCountVec),
		global_neighbors(_global_neighbors),
		global_length_zero(_global_length_zero),
		global_isedge_collagen(_global_isedge_collagen),
		global_isedge_elastin(_global_isedge_elastin),

		id_value_expanded(_id_value_expanded),
		key_begin_net_intc(_key_begin),
		key_end_net_intc(_key_end),

		collagen_diameter(_collagen_diameter),
		elastin_diameter(_elastin_diameter),
		max_nbr_count(_max_nbr_count),
		max_node_count(_max_node_count),

		max_links_per_iteration(_max_links_per_iteration),
		idTempVecLeft(_idTempVecLeft),
		idTempVecRight(_idTempVecRight) {}

	__device__

	unsigned operator() (const Tuu& u2) {
		unsigned nodeId = thrust::get<0>(u2);//node to attempt link from.
		bool node_is_collagen = node_is_collagen_addr[nodeId];
		bool node_is_elastin = node_is_elastin_addr[nodeId];
		unsigned bucketId = thrust::get<1>(u2);//bucket containing nodeId

		unsigned final_id_left = 0;
		unsigned final_id_right = 0;

		unsigned possibleEdgeBegin = max_nbr_count * nodeId;
		unsigned possibleEdgeEnd = max_nbr_count * nodeId + max_nbr_count;//test all possible for now
		unsigned numPlacedLinks = 0;

		unsigned save_index = max_links_per_iteration * nodeId;//changes after adding an edge
		unsigned save_index_end = max_links_per_iteration * nodeId + max_links_per_iteration;


		if ( nodeId < max_node_count )  {
			//beginning and end of attempted linking end id's in id_value_expanded
			unsigned beginIndex = key_begin_net_intc[bucketId];
			unsigned endIndex = key_end_net_intc[bucketId];

			for (unsigned iter = beginIndex; iter < endIndex; iter++ ) {

				unsigned candidateId = id_value_expanded[iter];//test id

				if ( (candidateId < max_node_count) && (save_index < save_index_end) ) {
					//after candidate is chosen, determine link type
					bool candidate_is_collagen = node_is_collagen_addr[candidateId];
					bool candidate_is_elastin = node_is_elastin_addr[candidateId];
					bool connection_is_collagen = false;
					bool connection_is_elastin = false;
					double fiber_diameter = 0.0;//fiber diameter changes depending

					if (node_is_elastin && candidate_is_elastin){
						connection_is_elastin = true;
						fiber_diameter = elastin_diameter;
					}
					else if (node_is_collagen && candidate_is_collagen){
						connection_is_collagen = true;
						fiber_diameter = collagen_diameter;
					}
					else{
						//then one is collagen and another is elasin
						connection_is_elastin = true;
						fiber_diameter = (elastin_diameter + collagen_diameter) / 2.0;
					}

					//then candidateId is not a dpd particle.
					bool candidateIsNew = true;

					//we now attempt to link nodeId and candidateEdge.
					if ((candidateId != nodeId) ){
						double dist = sqrt(
							((nodeLocXAddr[nodeId] - nodeLocXAddr[candidateId]) * (nodeLocXAddr[nodeId] - nodeLocXAddr[candidateId])) +
							((nodeLocYAddr[nodeId] - nodeLocYAddr[candidateId]) * (nodeLocYAddr[nodeId] - nodeLocYAddr[candidateId])) +
							((nodeLocZAddr[nodeId] - nodeLocZAddr[candidateId]) * (nodeLocZAddr[nodeId] - nodeLocZAddr[candidateId])));


						if ((dist < fiber_diameter)) {
							//then we have a possible link if no previous link was placed.

							//make sure placement is possible for nodeID
							//other thread will take care of candidateID
							for (unsigned k = possibleEdgeBegin; k < possibleEdgeEnd; k++) {
								if (candidateId == (global_neighbors[k])) {
									candidateIsNew = false;
									break;
								}
							}

							if (candidateIsNew == true) {

								//then the new edge can be placed
								unsigned location;
								//decide placement of nbr only in row of nodeId,
								//another thread will place for candidateId
								for (unsigned place = possibleEdgeBegin; place < possibleEdgeEnd; place++ ) {
									if (max_node_count == global_neighbors[place]) {
										//max_node_count is the default value of the matrix
										//then we have found location where we can place a node scine global neighbors default is max_node_count
										location = place;
										break;
									}
								}

								//create updated lengths
								global_neighbors[ location ] = candidateId;
								global_length_zero[ location ] = dist;

								//depending on the type of nodes connecting, create a different link.
								if (connection_is_elastin){
									global_isedge_collagen [ location ] = false;
									global_isedge_elastin [ location ] = true;
								}
								else if (connection_is_collagen){
										global_isedge_collagen [ location ] = true;
										global_isedge_elastin [ location ] = false;
								}
								
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
