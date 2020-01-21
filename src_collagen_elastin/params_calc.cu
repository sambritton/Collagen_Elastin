#include "system_structures.h"
#include "functor_misc.h"

#include "system.h"
#include "Params_Calc.h"
#include "functor_calc_strain_params.h"


void Params_Calc(
    EdgeInfoVecs& edgeInfoVecs,
    NodeInfoVecs& nodeInfoVecs,
    GeneralParams& generalParams,
    PltInfoVecs& pltInfoVecs) {

		//count positive and negative strains for edges that are not added. If an edge is added, a zero is placed on that strain.
		//notice that each thread will count edges twice, so divide by two at the end
	
		thrust::fill(nodeInfoVecs.discretized_edges_strain.begin(), nodeInfoVecs.discretized_edges_strain.end(),0.0);
		thrust::fill(nodeInfoVecs.discretized_edges_alignment.begin(), nodeInfoVecs.discretized_edges_alignment.end(),0.0);	

		//copy current host information to device for strain calculation. 
		thrust::copy(nodeInfoVecs.hostEdgeLeft.begin(),
			nodeInfoVecs.hostEdgeLeft.begin() + generalParams.currentEdgeCount,
			nodeInfoVecs.device_edge_left.begin());

		thrust::copy(nodeInfoVecs.hostEdgeRight.begin(),
			nodeInfoVecs.hostEdgeRight.begin() + generalParams.currentEdgeCount,
			nodeInfoVecs.device_edge_right.begin());

		thrust::transform(
			thrust::make_zip_iterator(
				thrust::make_tuple(
					nodeInfoVecs.device_edge_left.begin(),
					nodeInfoVecs.device_edge_right.begin())),
					 
			thrust::make_zip_iterator(
				thrust::make_tuple(
					nodeInfoVecs.device_edge_left.begin(),
					nodeInfoVecs.device_edge_right.begin())) + generalParams.currentEdgeCount,
					
			//outputs discretized strain etc			
			thrust::make_zip_iterator(
				thrust::make_tuple(
					nodeInfoVecs.discretized_edges_strain.begin(),
					nodeInfoVecs.discretized_edges_alignment.begin())),
					
			functor_calc_strain_params(
				generalParams.originLinkCount,
				generalParams.originEdgeCount,
				generalParams.originNodeCount,
				generalParams.max_node_count,
				generalParams.max_nbr_count,
				thrust::raw_pointer_cast(nodeInfoVecs.node_loc_x.data()),
				thrust::raw_pointer_cast(nodeInfoVecs.node_loc_y.data()),
				thrust::raw_pointer_cast(nodeInfoVecs.node_loc_z.data()),
				thrust::raw_pointer_cast(edgeInfoVecs.num_origin_nbr_per_node_vec.data()),
				thrust::raw_pointer_cast(edgeInfoVecs.current_node_edge_count_vec.data()),
				thrust::raw_pointer_cast(edgeInfoVecs.global_neighbors.data()),
				thrust::raw_pointer_cast(edgeInfoVecs.global_length_zero.data()) ));
		
			thrust::transform(
				thrust::make_zip_iterator(
					thrust::make_tuple(
						nodeInfoVecs.node_force_x.begin(),
						nodeInfoVecs.node_force_y.begin(),
						nodeInfoVecs.node_force_z.begin())),
				thrust::make_zip_iterator(
					thrust::make_tuple(
						nodeInfoVecs.node_force_x.begin(),
						nodeInfoVecs.node_force_y.begin(),
						nodeInfoVecs.node_force_z.begin())) + generalParams.max_node_count,
				nodeInfoVecs.sum_forces_on_node.begin(),//save vector
				
                functor_norm());

			//platelets
			thrust::transform(
				thrust::make_zip_iterator(
					thrust::make_tuple(
						pltInfoVecs.pltForceX.begin(),
						pltInfoVecs.pltForceY.begin(),
						pltInfoVecs.pltForceZ.begin())),
				thrust::make_zip_iterator(
					thrust::make_tuple(
						pltInfoVecs.pltForceX.begin(),
						pltInfoVecs.pltForceY.begin(),
						pltInfoVecs.pltForceZ.begin())) + generalParams.maxPltCount,
				pltInfoVecs.sumForcesOnPlt.begin(),//save vector
				functor_norm());
};