#include "functor_advance_pos.h"

#include "system_structures.h"
#include "system.h"
#include "advance_positions.h"


void advance_positions(
	NodeInfoVecs& nodeInfoVecs,
	GeneralParams& generalParams,
	EdgeInfoVecs& edgeInfoVecs,
	RandVecs& randVecs) {


		//At this point, the previous node location is the same as the current node,
		//we can therefore use previous node locations to update nodeLoc.
		unsigned _seed = rand();
    	thrust::counting_iterator<unsigned> index_sequence_begin(_seed);

    	thrust::transform(thrust::device, index_sequence_begin, index_sequence_begin + (generalParams.max_node_count),
        	randVecs.gaussianData.begin(), psrunifgen(-1.0, 1.0));

		thrust::counting_iterator<unsigned> nodeIndexBegin(0);

		thrust::transform(
			thrust::make_zip_iterator(
				thrust::make_tuple(
					nodeIndexBegin,
					nodeInfoVecs.node_loc_x.begin(),
					nodeInfoVecs.node_loc_y.begin(),
					nodeInfoVecs.node_loc_z.begin(),
					nodeInfoVecs.node_upper_selection_pull.begin(),
					nodeInfoVecs.node_lower_selection_pull.begin())),
			thrust::make_zip_iterator(
				thrust::make_tuple(
					nodeIndexBegin,
					nodeInfoVecs.node_loc_x.begin(),
					nodeInfoVecs.node_loc_y.begin(),
					nodeInfoVecs.node_loc_z.begin(),
					nodeInfoVecs.node_upper_selection_pull.begin(),
					nodeInfoVecs.node_lower_selection_pull.begin())) + generalParams.max_node_count,
			//second vector begin
			thrust::make_zip_iterator(
				thrust::make_tuple(
					randVecs.gaussianData.begin(),
					nodeInfoVecs.node_force_x.begin(),
					nodeInfoVecs.node_force_y.begin(),
					nodeInfoVecs.node_force_z.begin())),
			//save result in third vector to test values
			thrust::make_zip_iterator(
				thrust::make_tuple(
					nodeInfoVecs.node_loc_x.begin(),
					nodeInfoVecs.node_loc_y.begin(),
					nodeInfoVecs.node_loc_z.begin(),
					nodeInfoVecs.node_vel.begin())),
			functor_advance_pos(
				generalParams.dt,
				edgeInfoVecs.viscosity_collagen,
				edgeInfoVecs.viscosity_elastin,
				edgeInfoVecs.temperature,
				edgeInfoVecs.kB,
				edgeInfoVecs.node_mass,
				generalParams.max_node_count,
				thrust::raw_pointer_cast(nodeInfoVecs.node_is_collagen.data()),
				thrust::raw_pointer_cast(nodeInfoVecs.node_is_elastin.data()),
				thrust::raw_pointer_cast(nodeInfoVecs.is_node_fixed.data())));

}
