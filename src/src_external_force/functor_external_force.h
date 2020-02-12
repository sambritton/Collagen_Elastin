#ifndef FUNCTOR_EXTERNAL_FORCE_H_
#define FUNCTOR_EXTERNAL_FORCE_H_

#include "system_structures.h"

struct functor_external_force {
	bool* isNodeFixedAddr;
	double* forceXAddr;
	double* forceYAddr;
	double* forceZAddr;
	double magForce;
	double originalNetworkLength;
	double strain_proportion_end_sim;
	double averageLowerStrain;
	double averageUpperStrain;

	__host__ __device__
		functor_external_force(
			bool* _isNodeFixedAddr,
			double*	_forceXAddr,
			double*	_forceYAddr,
			double*	_forceZAddr,
			double& _magForce,
			double& _originalNetworkLength,
			double& _strain_proportion_end_sim,
			double& _averageLowerStrain,
			double& _averageUpperStrain) :
		isNodeFixedAddr(_isNodeFixedAddr),
		forceXAddr(_forceXAddr),
		forceYAddr(_forceYAddr),
		forceZAddr(_forceZAddr),
		
		magForce(_magForce),
		originalNetworkLength(_originalNetworkLength),
		strain_proportion_end_sim(_strain_proportion_end_sim),
		averageLowerStrain(_averageLowerStrain),
		averageUpperStrain(_averageUpperStrain) {}

	__device__
	void operator()(const Tudbbbb& u1d1b4) {

		unsigned id = thrust::get<0>(u1d1b4);
		double locZ = thrust::get<1>(u1d1b4);
		bool isFixed = thrust::get<2>(u1d1b4);
		bool is_node_collagen = thrust::get<3>(u1d1b4);
		bool isUpperStrainNode = thrust::get<4>(u1d1b4);
		bool isLowerStrainNode = thrust::get<5>(u1d1b4);

		//pull top
		if ((!isFixed) && (isUpperStrainNode)) {
			//if not fixed, we can apply force unless it is too far away, then we fix it.
			if (locZ > originalNetworkLength * (strain_proportion_end_sim) ) {
				isNodeFixedAddr[id] = true;
			}

			double upperDiff = abs(locZ - averageUpperStrain);

			//only apply force if within 1 of the average.
			if ( (upperDiff < 0.5) ){
				double dirX = 0.0;//tForceX / normForce;
				double dirY = 0.0;//tForceY / normForce;
				double dirZ = 1.0;
				forceXAddr[id] = dirX * (magForce);
				forceYAddr[id] = dirY * (magForce);
				forceZAddr[id] = dirZ * (magForce);
			}
			if ((is_node_collagen) && (upperDiff < 2.0)) {
				double dirX = 0.0;//tForceX / normForce;
				double dirY = 0.0;//tForceY / normForce;
				double dirZ = 1.0;
				forceXAddr[id] = dirX * (magForce);
				forceYAddr[id] = dirY * (magForce);
				forceZAddr[id] = dirZ * (magForce);
			}
		}


		//pull bottom
		else if ((!isFixed) && (isLowerStrainNode) ) {
			//safety fix
			if (locZ < -(originalNetworkLength * (strain_proportion_end_sim))) {
				isNodeFixedAddr[id] = true;
			}


			double lowerDiff = abs(locZ - averageLowerStrain);

			if (lowerDiff < 0.5) {
				double dirX = 0.0;//tForceX / normForce;
				double dirY = 0.0;//tForceY / normForce;
				double dirZ = -1.0;
				forceXAddr[id] = dirX * (magForce);
				forceYAddr[id] = dirY * (magForce);
				forceZAddr[id] = dirZ * (magForce);
			}
			if ((is_node_collagen) && (lowerDiff < 2.0)) {
				double dirX = 0.0;//tForceX / normForce;
				double dirY = 0.0;//tForceY / normForce;
				double dirZ = -1.0;
				forceXAddr[id] = dirX * (magForce);
				forceYAddr[id] = dirY * (magForce);
				forceZAddr[id] = dirZ * (magForce);
			}

		}
	}
};

#endif
