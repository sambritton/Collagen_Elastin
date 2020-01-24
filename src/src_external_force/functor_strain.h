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
	double operator() (const Tubd& u1b1d1) {
		unsigned id = thrust::get<0>(u1b1d1);
		bool isStrainNode = thrust::get<1>(u1b1d1);

		if (isStrainNode) {
			double zpos = thrust::get<2>(u1b1d1);
			return zpos;
		}
		else {
			return (0.0);
		}
	}

};

#endif
