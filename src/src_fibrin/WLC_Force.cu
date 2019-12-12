#include "System.h"
#include "SystemStructures.h"
#include "WLC_Force.h" 
#include "functor_wlc.h"

/*
the structure of lengthZero_index is 
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

void WLC_Force(
	NodeInfoVecs& nodeInfoVecs,
	WLCInfoVecs& wlcInfoVecs,  
	GeneralParams& generalParams) {
 
 
	thrust::counting_iterator<unsigned> startEdgeIter(0);
			  
	//
	thrust::for_each( 
		thrust::make_zip_iterator( 
			thrust::make_tuple(startEdgeIter,
								nodeInfoVecs.isNodeFixed.begin() )),
		thrust::make_zip_iterator(
			thrust::make_tuple(startEdgeIter,
								nodeInfoVecs.isNodeFixed.begin() )) + generalParams.maxNodeCount,
		functor_wlc(
			thrust::raw_pointer_cast(nodeInfoVecs.nodeLocX.data()),
			thrust::raw_pointer_cast(nodeInfoVecs.nodeLocY.data()),
			thrust::raw_pointer_cast(nodeInfoVecs.nodeLocZ.data()),
			thrust::raw_pointer_cast(nodeInfoVecs.nodeForceX.data()),
			thrust::raw_pointer_cast(nodeInfoVecs.nodeForceY.data()),
			thrust::raw_pointer_cast(nodeInfoVecs.nodeForceZ.data()),
 
			generalParams.kB,
			generalParams.persistenceLengthMon,
			generalParams.CLM,
			generalParams.temperature,
			generalParams.maxNeighborCount,
			generalParams.maxNodeCount,
			generalParams.nummonfiberarea,

			thrust::raw_pointer_cast(wlcInfoVecs.lengthZero.data()),
			thrust::raw_pointer_cast(wlcInfoVecs.globalNeighbors.data()),
			thrust::raw_pointer_cast(wlcInfoVecs.currentNodeEdgeCountVector.data()),
			thrust::raw_pointer_cast(wlcInfoVecs.numOriginalNeighborsNodeVector.data()) ) );
};

