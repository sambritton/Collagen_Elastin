#include "system_structures.h"
#include "functor_misc.h"

#include "system.h"
#include "params_calc.h"
#include "functor_calc_strain_params.h"


void params_calc(
    EdgeInfoVecs& edgeInfoVecs,
    NodeInfoVecs& nodeInfoVecs,
    GeneralParams& generalParams,
    PltInfoVecs& pltInfoVecs) {

		//count positive and negative strains for edges that are not added. If an edge is added, a zero is placed on that strain.
		//notice that each thread will count edges twice, so divide by two at the end
	
		thrust::fill(nodeInfoVecs.discretized_edges_strain.begin(), nodeInfoVecs.discretized_edges_strain.end(),0.0);
		thrust::fill(nodeInfoVecs.discretized_edges_alignment.begin(), nodeInfoVecs.discretized_edges_alignment.end(),0.0);	

		//copy current host information to device for strain calculation. 
		thrust::copy(nodeInfoVecs.host_edge_left.begin(),
			nodeInfoVecs.host_edge_left.begin() + generalParams.current_edge_count,
			nodeInfoVecs.device_edge_left.begin());

		thrust::copy(nodeInfoVecs.host_edge_right.begin(),
			nodeInfoVecs.host_edge_right.begin() + generalParams.current_edge_count,
			nodeInfoVecs.device_edge_right.begin());

		thrust::transform(
			thrust::make_zip_iterator(
				thrust::make_tuple(
					nodeInfoVecs.device_edge_left.begin(),
					nodeInfoVecs.device_edge_right.begin())),
					 
			thrust::make_zip_iterator(
				thrust::make_tuple(
					nodeInfoVecs.device_edge_left.begin(),
					nodeInfoVecs.device_edge_right.begin())) + generalParams.current_edge_count,
					
			//outputs discretized strain etc			
			thrust::make_zip_iterator(
				thrust::make_tuple(
					nodeInfoVecs.discretized_edges_strain.begin(),
					nodeInfoVecs.discretized_edges_alignment.begin())),
					
			functor_calc_strain_params(
				generalParams.origin_node_count,
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

};