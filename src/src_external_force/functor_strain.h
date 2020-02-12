#ifndef FUNCTOR_STRAIN_H_
#define FUNCTOR_STRAIN_H_

#include "system_structures.h"

struct functor_strain {
	unsigned maxNodeCount;
	double maxNetworkLength;

	__host__ __device__
	functor_strain(
		unsigned& _maxNodeCount,
		double& _maxNetworkLength):
		maxNodeCount(_maxNodeCount),
		maxNetworkLength(_maxNetworkLength) {}

	__device__
	double operator() (const Tbbd& b2d1) {
		
		bool is_collagen_node = thrust::get<0>(b2d1);
		bool isStrainNode = thrust::get<1>(b2d1);

		double zpos=0.0;
		if ((isStrainNode) && (is_collagen_node)) {
			zpos = thrust::get<2>(b2d1);	
		}
		return zpos;
	}
};

#endif
