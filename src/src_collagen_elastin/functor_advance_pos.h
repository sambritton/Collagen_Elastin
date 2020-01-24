#ifndef FUNCTOR_ADVANCE_POS_H_
#define FUNCTOR_ADVANCE_POS_H_

#include "system_structures.h"

//used to advance position and velocity from known force
//Velocity Verlet Like algorithm used : https://arxiv.org/pdf/1212.1244.pdf
struct functor_advance_pos : public thrust::binary_function<UCVec3, CVec4, CVec4> {
	double dt;
	double viscosity_collagen;
	double viscosity_elastin;
	double temperature;
	double kB;
	double mass;
	unsigned max_node_count;

	bool* isNodeCollagenAddr;
	bool* isNodeElastinAddr;
	bool* isNodeFixedAddr;

	__host__ __device__
		//
		functor_advance_pos(
			double& _dt,
			double& _viscosity_collagen,
			double& _viscosity_elastin,
			double& _temperature,
			double& _kB,
			double& _mass,
			unsigned& _max_node_count,
			bool* _isNodeCollagenAddr,
			bool* _isNodeElastinAddr,
			bool* _isNodeFixedAddr) :
		dt(_dt),
		viscosity_collagen(_viscosity_collagen),
		viscosity_elastin(_viscosity_elastin),
		temperature(_temperature),
		kB(_kB),
		mass(_mass),
		max_node_count(_max_node_count),
		isNodeCollagenAddr(_isNodeCollagenAddr),
		isNodeElastinAddr(_isNodeElastinAddr),
		isNodeFixedAddr(_isNodeFixedAddr) {}

	__device__
		CVec4 operator()(const UCVec3 &id1p3, const CVec4 &g1f3) {

		unsigned id = thrust::get<0>(id1p3);
		bool isFixed = false;//true if fixed, false if movable.
		bool isCollagen = false;
		bool isElastin = false;
		double viscosity = viscosity_collagen;
		if (id < max_node_count) {
			isFixed = isNodeFixedAddr[id];
			//isCollagen = isNodeCollagenAddr[id];
			isElastin = isNodeElastinAddr[id];
		}
		if (isElastin){
			viscosity = viscosity_elastin;
		}

		if (!isFixed) {
			double locX = thrust::get<1>(id1p3);
			double locY = thrust::get<2>(id1p3);
			double locZ = thrust::get<3>(id1p3);

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
			return thrust::make_tuple(thrust::get<1>(id1p3),
			thrust::get<2>(id1p3), thrust::get<3>(id1p3),0.0);
		}
	}

};

#endif
