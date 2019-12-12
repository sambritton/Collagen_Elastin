#include "SystemStructures.h"
#include "functor_misc.h"

#include "System.h"
#include "Params_Calc.h"
#include "functor_calc_strain_params.h"


void Params_Calc(
    WLCInfoVecs& wlcInfoVecs,
    NodeInfoVecs& nodeInfoVecs,
    GeneralParams& generalParams,
    PltInfoVecs& pltInfoVecs) {

		//count positive and negative strains for edges that are not added. If an edge is added, a zero is placed on that strain.
		//notice that each thread will count edges twice, so divide by two at the end
	
		thrust::fill(nodeInfoVecs.discretizedEdgeStrain.begin(), nodeInfoVecs.discretizedEdgeStrain.end(),0.0);
		thrust::fill(nodeInfoVecs.discretizedEdgeAlignment.begin(), nodeInfoVecs.discretizedEdgeAlignment.end(),0.0);	

		//copy current host information to device for strain calculation. 
		thrust::copy(nodeInfoVecs.hostEdgeLeft.begin(),
			nodeInfoVecs.hostEdgeLeft.begin() + generalParams.currentEdgeCount,
			nodeInfoVecs.deviceEdgeLeft.begin());

		thrust::copy(nodeInfoVecs.hostEdgeRight.begin(),
			nodeInfoVecs.hostEdgeRight.begin() + generalParams.currentEdgeCount,
			nodeInfoVecs.deviceEdgeRight.begin());

		thrust::transform(
			thrust::make_zip_iterator(
				thrust::make_tuple(
					nodeInfoVecs.deviceEdgeLeft.begin(),
					nodeInfoVecs.deviceEdgeRight.begin())),
					 
			thrust::make_zip_iterator(
				thrust::make_tuple(
					nodeInfoVecs.deviceEdgeLeft.begin(),
					nodeInfoVecs.deviceEdgeRight.begin())) + generalParams.currentEdgeCount,
					
			//outputs discretized strain etc			
			thrust::make_zip_iterator(
				thrust::make_tuple(
					nodeInfoVecs.discretizedEdgeStrain.begin(),
					nodeInfoVecs.discretizedEdgeAlignment.begin())),
					
			functor_calc_strain_params(
				generalParams.originLinkCount,
				generalParams.originEdgeCount,
				generalParams.originNodeCount,
				generalParams.maxNodeCount,
				generalParams.maxNeighborCount,
				thrust::raw_pointer_cast(nodeInfoVecs.nodeLocX.data()),
				thrust::raw_pointer_cast(nodeInfoVecs.nodeLocY.data()),
				thrust::raw_pointer_cast(nodeInfoVecs.nodeLocZ.data()),
				thrust::raw_pointer_cast(wlcInfoVecs.numOriginalNeighborsNodeVector.data()),
				thrust::raw_pointer_cast(wlcInfoVecs.currentNodeEdgeCountVector.data()),
				thrust::raw_pointer_cast(wlcInfoVecs.globalNeighbors.data()),
				thrust::raw_pointer_cast(wlcInfoVecs.lengthZero.data()) ));
		
			thrust::transform(
				thrust::make_zip_iterator(
					thrust::make_tuple(
						nodeInfoVecs.nodeForceX.begin(),
						nodeInfoVecs.nodeForceY.begin(),
						nodeInfoVecs.nodeForceZ.begin())),
				thrust::make_zip_iterator(
					thrust::make_tuple(
						nodeInfoVecs.nodeForceX.begin(),
						nodeInfoVecs.nodeForceY.begin(),
						nodeInfoVecs.nodeForceZ.begin())) + generalParams.maxNodeCount,
				nodeInfoVecs.sumForcesOnNode.begin(),//save vector
				
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