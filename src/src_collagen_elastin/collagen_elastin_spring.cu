#include "system.h"
#include "system_structures.h"
#include "WLC_Force.h"
#include "functor_wlc.h"

/*
the structure of global_length_zero_index is
0  1  2  3
4  5  6  7
8  9  10 11
12 13 14 15 for a 4 node system.
index/4 = row,
index%4 = col. If you apply force to column node always or row node always then
each thread will apply opposing forces to springs.
if you decide to apply force to column instead of rows, you'll need sign change
LengthZero_value is symmetric, so values line up correctly.
*/

void calc_spring_force(
	NodeInfoVecs& nodeInfoVecs,
	EdgeInfoVecs& edgeInfoVecs,
	GeneralParams& generalParams) {


	thrust::counting_iterator<unsigned> startEdgeIter(0);

	//
	thrust::for_each(
		thrust::make_zip_iterator(
			thrust::make_tuple(startEdgeIter,
								nodeInfoVecs.is_node_fixed.begin() )),
		thrust::make_zip_iterator(
			thrust::make_tuple(startEdgeIter,
								nodeInfoVecs.is_node_fixed.begin() )) + generalParams.max_node_count,
		functor_collagen_elastin(
			thrust::raw_pointer_cast(nodeInfoVecs.node_loc_x.data()),
			thrust::raw_pointer_cast(nodeInfoVecs.node_loc_y.data()),
			thrust::raw_pointer_cast(nodeInfoVecs.node_loc_z.data()),
			thrust::raw_pointer_cast(nodeInfoVecs.node_force_x.data()),
			thrust::raw_pointer_cast(nodeInfoVecs.node_force_y.data()),
			thrust::raw_pointer_cast(nodeInfoVecs.node_force_z.data()),

 			edgeInfoVecs.collagen_spring_constant,
			edgeInfoVecs.kB,
			generalParams.persistence_len_monomer,
			edgeInfoVecs.CLM,
			edgeInfoVecs.temperature,
			generalParams.max_nbr_count,
			generalParams.max_node_count,
			generalParams.nummonfiberarea,

			thrust::raw_pointer_cast(edgeInfoVecs.global_length_zero.data()),
			thrust::raw_pointer_cast(edgeInfoVecs.global_neighbors.data()),
			thrust::raw_pointer_cast(edgeInfoVecs.global_isedge_collagen.data()),
			thrust::raw_pointer_cast(edgeInfoVecs.global_isedge_elastin.data()),
			thrust::raw_pointer_cast(edgeInfoVecs.current_node_edge_count_vec.data()),
			thrust::raw_pointer_cast(edgeInfoVecs.num_origin_nbr_per_node_vec.data()) ) );
};
