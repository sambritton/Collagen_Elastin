#ifndef FUNCTOR_STRAIN_H_
#define FUNCTOR_STRAIN_H_

#include "system_structures.h"

struct functor_strain {
	unsigned axis;
	double maxNetworkLength;

	__host__ __device__
	functor_strain(
		unsigned& _axis,
		double& _maxNetworkLength):
		axis(_axis),
		maxNetworkLength(_maxNetworkLength) {}

	__device__
	double operator() (const Tbbdd& b2d2) {
		
		bool is_collagen_node = thrust::get<0>(b2d2);
		bool isStrainNode = thrust::get<1>(b2d2);
		double position = 0.0;
		if (axis == 0){
			//zposition
			position = thrust::get<2>(b2d2);
		}else if (axis == 1){
			//xposition
			position = thrust::get<3>(b2d2);
		}

		if ((isStrainNode) && (is_collagen_node)) {
			return position;	
		}
		return 0.0;
	}
};

#endif
