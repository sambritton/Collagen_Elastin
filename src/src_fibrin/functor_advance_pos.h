#ifndef FUNCTOR_ADVANCE_POS_H_
#define FUNCTOR_ADVANCE_POS_H_

#include "SystemStructures.h"

//used to advance position and velocity from known force
//Velocity Verlet Like algorithm used : https://arxiv.org/pdf/1212.1244.pdf
struct functor_advance_pos : public thrust::binary_function<UCVec3, CVec4, CVec4> {
	double dt;
	double viscosity;
	double temperature;
	double kB;
	double mass;
	bool* isNodeFixedAddr;
	unsigned maxNodeCount;

	__host__ __device__
		//
		functor_advance_pos(
			double& _dt,
			double& _viscosity,
			double& _temperature,
			double& _kB,
			double& _mass,
			unsigned& _maxNodeCount,
			bool* _isNodeFixedAddr) :
		dt(_dt),
		viscosity(_viscosity),
		temperature(_temperature),
		kB(_kB),
		mass(_mass),
		maxNodeCount(_maxNodeCount),
		isNodeFixedAddr(_isNodeFixedAddr) {}

	__device__
		CVec4 operator()(const UCVec3 &p3, const CVec4 &g1f3) {

		unsigned id = thrust::get<0>(p3);
		bool isFixed =false;//true if fixed, false if movable.
		if (id < maxNodeCount) {
			isFixed = isNodeFixedAddr[id];
		}

		if (!isFixed) {
			double locX = thrust::get<1>(p3);
			double locY = thrust::get<2>(p3);
			double locZ = thrust::get<3>(p3);

			//random data
			double gaussianData = thrust::get<0>(g1f3);

			//normally you would have x_(n+1) - x(n) / dt = F / eta + F_b / eta, and F_b = sqrt(2*kb*T*eta/dt) after multiplication, additive noise
			//becomes sqrt(2*kb*t*dt/eta) * N(0,1)
			double noise = sqrt(2.0 * kB* temperature * dt  / viscosity) * gaussianData;

			double accX = (thrust::get<1>(g1f3));
			double accY = (thrust::get<2>(g1f3));
			double accZ = (thrust::get<3>(g1f3));


			//update positions
			double xLocRes = locX + (dt/viscosity) * (accX) + noise;
			double yLocRes = locY + (dt/viscosity) * (accY) + noise;

			double zLocRes = locZ + (dt/viscosity) * (accZ) + noise;

			double velocity = sqrt((xLocRes - locX) * (xLocRes - locX) + (yLocRes - locY) * (yLocRes - locY) + (zLocRes - locZ) * (zLocRes - locZ));



			return thrust::make_tuple(xLocRes, yLocRes, zLocRes, velocity);
		}
		else {
			//do not move position.
			return thrust::make_tuple(thrust::get<1>(p3),
			thrust::get<2>(p3), thrust::get<3>(p3),0.0);
		}
	}

};

#endif