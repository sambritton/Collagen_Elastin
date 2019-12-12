#include "System.h"

#include "SystemStructures.h"

#include "Torsion_Force.h"

#include "functor_misc.h"

#include "functor_torsion.h"


void Torsion_Force(
	NodeInfoVecs& nodeInfoVecs,
	TorsionInfoVecs& torsionInfoVecs,
	GeneralParams& generalParams)  {
	
const double PI = 3.14159265358979323846;  
if (generalParams.totalTorsionCount>0) { 

		thrust::counting_iterator<unsigned> startTorsionIter(0);
		thrust::counting_iterator<unsigned> endTorsionIter(generalParams.totalTorsionCount);
 
		//for_each guarrantees order. This is needed for iter count and saving to torsion force vectors.
		thrust::for_each(
			thrust::make_zip_iterator(
				thrust::make_tuple(
					startTorsionIter,
					torsionInfoVecs.leftIndex.begin(),
					torsionInfoVecs.centerIndex.begin(),
					torsionInfoVecs.rightIndex.begin(),
					torsionInfoVecs.angleZero.begin())),
			thrust::make_zip_iterator(
				thrust::make_tuple(
					endTorsionIter,
					torsionInfoVecs.leftIndex.end(),
					torsionInfoVecs.centerIndex.end(),
					torsionInfoVecs.rightIndex.end(),
					torsionInfoVecs.angleZero.end())),
			functor_torsion(
				thrust::raw_pointer_cast(nodeInfoVecs.nodeLocX.data()),
				thrust::raw_pointer_cast(nodeInfoVecs.nodeLocY.data()),
				thrust::raw_pointer_cast(nodeInfoVecs.nodeLocZ.data()),
				thrust::raw_pointer_cast(torsionInfoVecs.forceX.data()),
				thrust::raw_pointer_cast(torsionInfoVecs.forceY.data()),
				thrust::raw_pointer_cast(torsionInfoVecs.forceZ.data()),
				thrust::raw_pointer_cast(nodeInfoVecs.isNodeFixed.data()),
				generalParams.torsionStiffness,
				generalParams.maxNodeCount,
				generalParams.totalTorsionCount,
				PI));  

		//reduce by key to get forces.Notice leftIndex is 1/3rd the length of torsion.forceX
		//this vector will be sorted each iteration, so it needs to be recopied.
		thrust::copy(torsionInfoVecs.leftIndex.begin(), torsionInfoVecs.leftIndex.end(), torsionInfoVecs.tempTorIndices.begin());
		thrust::copy(torsionInfoVecs.centerIndex.begin(), torsionInfoVecs.centerIndex.end(), torsionInfoVecs.tempTorIndices.begin() + generalParams.totalTorsionCount);
		thrust::copy(torsionInfoVecs.rightIndex.begin(), torsionInfoVecs.rightIndex.end(), torsionInfoVecs.tempTorIndices.begin() + 2 * generalParams.totalTorsionCount);


		//key, then value. Each vector returns sorted		
		thrust::sort_by_key(torsionInfoVecs.tempTorIndices.begin(), torsionInfoVecs.tempTorIndices.begin() + 3 * generalParams.totalTorsionCount,
			thrust::make_zip_iterator(
				thrust::make_tuple(
					torsionInfoVecs.forceX.begin(),
					torsionInfoVecs.forceY.begin(),
					torsionInfoVecs.forceZ.begin())), thrust::less<unsigned>());


		thrust::fill(torsionInfoVecs.tempForceX.begin(), torsionInfoVecs.tempForceX.end(), 0);
		thrust::fill(torsionInfoVecs.tempForceY.begin(), torsionInfoVecs.tempForceY.end(), 0);
		thrust::fill(torsionInfoVecs.tempForceZ.begin(), torsionInfoVecs.tempForceZ.end(), 0);
		thrust::fill(torsionInfoVecs.reducedIds.begin(), torsionInfoVecs.reducedIds.end(), 0);

		unsigned endKey = thrust::get<0>(
			thrust::reduce_by_key(
				torsionInfoVecs.tempTorIndices.begin(), 
				torsionInfoVecs.tempTorIndices.begin() + 3*generalParams.totalTorsionCount,
			thrust::make_zip_iterator(
				thrust::make_tuple(
					torsionInfoVecs.forceX.begin(),
					torsionInfoVecs.forceY.begin(),
					torsionInfoVecs.forceZ.begin())),
			torsionInfoVecs.reducedIds.begin(),
			thrust::make_zip_iterator(
				thrust::make_tuple(
					torsionInfoVecs.tempForceX.begin(),
					torsionInfoVecs.tempForceY.begin(),
					torsionInfoVecs.tempForceZ.begin())),
			thrust::equal_to<unsigned>(), CVec3Add())) - torsionInfoVecs.reducedIds.begin();//binary_pred, binary_op

		
		thrust::for_each(
			thrust::make_zip_iterator(//1st begin
				thrust::make_tuple(
					torsionInfoVecs.reducedIds.begin(),
					torsionInfoVecs.tempForceX.begin(),
					torsionInfoVecs.tempForceY.begin(),
					torsionInfoVecs.tempForceZ.begin())),
			thrust::make_zip_iterator(//1st end
				thrust::make_tuple(
					torsionInfoVecs.reducedIds.begin(),
					torsionInfoVecs.tempForceX.begin(),
					torsionInfoVecs.tempForceY.begin(),
					torsionInfoVecs.tempForceZ.begin())) + endKey,
			functor_add_UCVec3_CVec3(
				generalParams.maxNodeCount,
				thrust::raw_pointer_cast(nodeInfoVecs.nodeForceX.data()),
				thrust::raw_pointer_cast(nodeInfoVecs.nodeForceY.data()),
				thrust::raw_pointer_cast(nodeInfoVecs.nodeForceZ.data())));

	}

	
}