#include "functor_advance_pos.h"

#include "SystemStructures.h"
#include "System.h"
#include "Advance_Positions_Fibrin.h"


void Advance_Positions_Fibrin(
	NodeInfoVecs& nodeInfoVecs,
	GeneralParams& generalParams,
	RandVecs& randVecs) {


		//At this point, the previous node location is the same as the current node,
		//we can therefore use previous node locations to update nodeLoc.
		unsigned _seed = rand();
    	thrust::counting_iterator<unsigned> index_sequence_begin(_seed);

    	thrust::transform(thrust::device, index_sequence_begin, index_sequence_begin + (generalParams.maxNodeCount),
        	randVecs.gaussianData.begin(), psrunifgen(-1.0, 1.0));

		thrust::counting_iterator<unsigned> nodeIndexBegin(0);

		thrust::transform(
			thrust::make_zip_iterator(
				thrust::make_tuple(
					nodeIndexBegin,
					nodeInfoVecs.nodeLocX.begin(),
					nodeInfoVecs.nodeLocY.begin(),
					nodeInfoVecs.nodeLocZ.begin())),
			thrust::make_zip_iterator(
				thrust::make_tuple(
					nodeIndexBegin,
					nodeInfoVecs.nodeLocX.begin(),
					nodeInfoVecs.nodeLocY.begin(),
					nodeInfoVecs.nodeLocZ.begin())) + generalParams.maxNodeCount,
			//second vector begin
			thrust::make_zip_iterator(
				thrust::make_tuple(
					randVecs.gaussianData.begin(),
					nodeInfoVecs.nodeForceX.begin(),
					nodeInfoVecs.nodeForceY.begin(),
					nodeInfoVecs.nodeForceZ.begin())),
			//save result in third vector to test values
			thrust::make_zip_iterator(
				thrust::make_tuple(
					nodeInfoVecs.nodeLocX.begin(),
					nodeInfoVecs.nodeLocY.begin(),
					nodeInfoVecs.nodeLocZ.begin(),
					nodeInfoVecs.nodeVelocity.begin())),
			functor_advance_pos(generalParams.dtTemp,
				generalParams.viscousDamp_Fibrin,
				generalParams.temperature,
				generalParams.kB,
				generalParams.nodeMass,
				generalParams.maxNodeCount,
				thrust::raw_pointer_cast(nodeInfoVecs.isNodeFixed.data())));

}
