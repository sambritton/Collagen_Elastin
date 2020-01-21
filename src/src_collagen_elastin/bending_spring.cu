#include "system.h"
#include "system_structures.h"
#include "bending_spring.h"

#include "functor_misc.h"
#include "functor_torsion.h"


void calc_bending_spring_force(
	NodeInfoVecs& nodeInfoVecs,
	BendInfoVecs& bendInfoVecs,
	GeneralParams& generalParams)  {

const double PI = 3.14159265358979323846;
if (bendInfoVecs.total_bend_count>0) {

		thrust::counting_iterator<unsigned> startTorsionIter(0);
		thrust::counting_iterator<unsigned> endTorsionIter(bendInfoVecs.total_bend_count);

		//for_each guarrantees order. This is needed for iter count and saving to torsion force vectors.
		thrust::for_each(
			thrust::make_zip_iterator(
				thrust::make_tuple(
					startTorsionIter,
					bendInfoVecs.leftIndex.begin(),
					bendInfoVecs.centerIndex.begin(),
					bendInfoVecs.rightIndex.begin(),
					bendInfoVecs.angleZero.begin())),
			thrust::make_zip_iterator(
				thrust::make_tuple(
					endTorsionIter,
					bendInfoVecs.leftIndex.end(),
					bendInfoVecs.centerIndex.end(),
					bendInfoVecs.rightIndex.end(),
					bendInfoVecs.angleZero.end())),
			functor_torsion(
				thrust::raw_pointer_cast(nodeInfoVecs.node_loc_x.data()),
				thrust::raw_pointer_cast(nodeInfoVecs.node_loc_y.data()),
				thrust::raw_pointer_cast(nodeInfoVecs.node_loc_z.data()),
				thrust::raw_pointer_cast(bendInfoVecs.forceX.data()),
				thrust::raw_pointer_cast(bendInfoVecs.forceY.data()),
				thrust::raw_pointer_cast(bendInfoVecs.forceZ.data()),

				thrust::raw_pointer_cast(nodeInfoVecs.is_node_fixed.data()),
				thrust::raw_pointer_cast(nodeInfoVecs.node_is_collagen.data()),
				thrust::raw_pointer_cast(nodeInfoVecs.node_is_elastin.data()),

				bendInfoVecs.bend_stiffness_collagen,
				bendInfoVecs.bend_stiffness_elastin,
				generalParams.max_node_count,
				bendInfoVecs.total_bend_count,
				PI));

		//reduce by key to get forces.Notice leftIndex is 1/3rd the length of torsion.forceX
		//this vector will be sorted each iteration, so it needs to be recopied.
		thrust::copy(bendInfoVecs.leftIndex.begin(), bendInfoVecs.leftIndex.end(), bendInfoVecs.tempTorIndices.begin());
		thrust::copy(bendInfoVecs.centerIndex.begin(), bendInfoVecs.centerIndex.end(), bendInfoVecs.tempTorIndices.begin() + bendInfoVecs.total_bend_count);
		thrust::copy(bendInfoVecs.rightIndex.begin(), bendInfoVecs.rightIndex.end(), bendInfoVecs.tempTorIndices.begin() + 2 * bendInfoVecs.total_bend_count);


		//key, then value. Each vector returns sorted
		thrust::sort_by_key(bendInfoVecs.tempTorIndices.begin(), bendInfoVecs.tempTorIndices.begin() + 3 * bendInfoVecs.total_bend_count,
			thrust::make_zip_iterator(
				thrust::make_tuple(
					bendInfoVecs.forceX.begin(),
					bendInfoVecs.forceY.begin(),
					bendInfoVecs.forceZ.begin())), thrust::less<unsigned>());


		thrust::fill(bendInfoVecs.tempForceX.begin(), bendInfoVecs.tempForceX.end(), 0);
		thrust::fill(bendInfoVecs.tempForceY.begin(), bendInfoVecs.tempForceY.end(), 0);
		thrust::fill(bendInfoVecs.tempForceZ.begin(), bendInfoVecs.tempForceZ.end(), 0);
		thrust::fill(bendInfoVecs.reducedIds.begin(), bendInfoVecs.reducedIds.end(), 0);

		unsigned endKey = thrust::get<0>(
			thrust::reduce_by_key(
				bendInfoVecs.tempTorIndices.begin(),
				bendInfoVecs.tempTorIndices.begin() + 3*bendInfoVecs.total_bend_count,
			thrust::make_zip_iterator(
				thrust::make_tuple(
					bendInfoVecs.forceX.begin(),
					bendInfoVecs.forceY.begin(),
					bendInfoVecs.forceZ.begin())),
			bendInfoVecs.reducedIds.begin(),
			thrust::make_zip_iterator(
				thrust::make_tuple(
					bendInfoVecs.tempForceX.begin(),
					bendInfoVecs.tempForceY.begin(),
					bendInfoVecs.tempForceZ.begin())),
			thrust::equal_to<unsigned>(), CVec3Add())) - bendInfoVecs.reducedIds.begin();//binary_pred, binary_op


		thrust::for_each(
			thrust::make_zip_iterator(//1st begin
				thrust::make_tuple(
					bendInfoVecs.reducedIds.begin(),
					bendInfoVecs.tempForceX.begin(),
					bendInfoVecs.tempForceY.begin(),
					bendInfoVecs.tempForceZ.begin())),
			thrust::make_zip_iterator(//1st end
				thrust::make_tuple(
					bendInfoVecs.reducedIds.begin(),
					bendInfoVecs.tempForceX.begin(),
					bendInfoVecs.tempForceY.begin(),
					bendInfoVecs.tempForceZ.begin())) + endKey,
			functor_add_UCVec3_CVec3(
				generalParams.max_node_count,
				thrust::raw_pointer_cast(nodeInfoVecs.node_force_x.data()),
				thrust::raw_pointer_cast(nodeInfoVecs.node_force_y.data()),
				thrust::raw_pointer_cast(nodeInfoVecs.node_force_z.data())));

	}


}
