#ifndef FUNCTOR_BENDING_H_
#define FUNCTOR_BENDING_H_

#include "system_structures.h"
struct functor_bending {
	double * locXAddr;
	double * locYAddr;
	double * locZAddr;
	double * forceXAddr;
	double * forceYAddr;
	double * forceZAddr;

	bool* isFixedVector;
	bool* node_is_collagen;
	bool* node_is_elastin;

	double bend_stiffness_collagen;
	double bend_stiffness_elastin;
	unsigned max_node_count;
	unsigned total_bending_spring_count;

	__host__ __device__
		functor_bending(
			double* _locXAddr,
			double* _locYAddr,
			double* _locZAddr,
			double* _forceXAddr,
			double* _forceYAddr,
			double* _forceZAddr,

			bool* _isFixedVector,
			bool* _node_is_collagen,
			bool* _node_is_elastin,

			double& _bend_stiffness_collagen,
			double& _bend_stiffness_elastin,
			unsigned& _max_node_count,
			unsigned& _total_bending_spring_count) :

		locXAddr(_locXAddr),
		locYAddr(_locYAddr),
		locZAddr(_locZAddr),
		forceXAddr(_forceXAddr),
		forceYAddr(_forceYAddr),
		forceZAddr(_forceZAddr),

		isFixedVector(_isFixedVector),
		node_is_collagen(_node_is_collagen),
		node_is_elastin(_node_is_elastin),

		bend_stiffness_collagen(_bend_stiffness_collagen),
		bend_stiffness_elastin(_bend_stiffness_elastin),
		max_node_count(_max_node_count),
		total_bending_spring_count(_total_bending_spring_count) {}

	__device__
		//maybe hand in id and vector of neighbors???
		void operator() (const Tuuuud& u4d1) {
		unsigned counter = thrust::get<0>(u4d1);
		unsigned indexLeft = thrust::get<1>(u4d1);
		unsigned indexCenter = thrust::get<2>(u4d1);
		unsigned indexRight = thrust::get<3>(u4d1);
		double angleZero = thrust::get<4>(u4d1);

		const double PI = 3.14159265358979323846;

		bool is_left_collagen = node_is_collagen[indexLeft];
		bool is_center_collagen = node_is_collagen[indexCenter];
		bool is_right_collagen = node_is_collagen[indexRight];

		bool is_left_elastin = node_is_elastin[indexLeft];
		bool is_center_elastin = node_is_elastin[indexCenter];
		bool is_right_elastin = node_is_elastin[indexRight];

		bool is_spring_collagen = false;
		bool is_spring_elastin = false;

		unsigned collagen_count = is_left_collagen + is_center_collagen + is_right_collagen;
		unsigned elastin_count = is_left_elastin + is_center_elastin + is_right_elastin;

		//depending on the count of nodes, we need to use a different bending constant.
		double bending_stiffness = 0.0;
		if (collagen_count > elastin_count){
			bending_stiffness = bend_stiffness_collagen;
		}
		if (elastin_count >= collagen_count){
			bending_stiffness = bend_stiffness_elastin;
		}

		unsigned indexLeftApply = counter;
		unsigned indexCenterApply = counter + total_bending_spring_count;
		unsigned indexRightApply = counter + 2*total_bending_spring_count;


		double distLCX = locXAddr[indexLeft] - locXAddr[indexCenter];
		double distLCY = locYAddr[indexLeft] - locYAddr[indexCenter];
		double distLCZ = locZAddr[indexLeft] - locZAddr[indexCenter];
		double distRCX = locXAddr[indexRight] - locXAddr[indexCenter];
		double distRCY = locYAddr[indexRight] - locYAddr[indexCenter];
		double distRCZ = locZAddr[indexRight] - locZAddr[indexCenter];

		//lengths between left & center, right & center
		double lenLC = sqrt(distLCX*distLCX + distLCY*distLCY + distLCZ*distLCZ);
		double lenRC = sqrt(distRCX*distRCX + distRCY*distRCY + distRCZ*distRCZ);


		//normalized dot product
		double cosTheta = (distLCX*distRCX + distLCY*distRCY + distLCZ*distRCZ) / (lenLC * lenRC);

		//rounding errors
		if (cosTheta < -1.0) {
			cosTheta = -1.0;
		}
		else if (cosTheta > 1.0) {
			cosTheta = 1.0;
		}

		//should be zero to pi.
		double currentAngle = acos(cosTheta);


		double magForce = -1.0 * bending_stiffness * (angleZero - currentAngle);

		//calculate force vectors
		//careful here. If the vectors are parallel, axb = 0;
		CVec3 axb = cross_product(thrust::make_tuple(distLCX, distLCY, distLCZ),
			thrust::make_tuple(distRCX, distRCY, distRCZ));

		CVec3 c1 = cross_product(axb, thrust::make_tuple(distLCX, distLCY, distLCZ));
		CVec3 c2 = cross_product(thrust::make_tuple(distRCX, distRCY, distRCZ), axb);

		double divisorA = sqrt(dot_product(c1, c1));
		double divisorB = sqrt(dot_product(c2, c2));

		CVec3 a = thrust::make_tuple(magForce * thrust::get<0>(c1) / divisorA,
				magForce * thrust::get<1>(c1) / divisorA,
				magForce * thrust::get<2>(c1) / divisorA);

			//calculate force on right node
		CVec3 b = thrust::make_tuple(magForce*thrust::get<0>(c2) / divisorB,
				magForce*thrust::get<1>(c2) / divisorB,
				magForce*thrust::get<2>(c2) / divisorB);

		if (isfinite(thrust::get<0>(a))
			&& isfinite(thrust::get<1>(a))
			&& isfinite(thrust::get<2>(a))
			&& isfinite(thrust::get<0>(b))
			&& isfinite(thrust::get<1>(b))
			&& isfinite(thrust::get<2>(b))) {

				//notice we write force before applying it. We thus rewrite and do not add to the current value.
			if (!isFixedVector[indexCenter]){
				forceXAddr[indexCenterApply] = -thrust::get<0>(a) - thrust::get<0>(b);
				forceYAddr[indexCenterApply] = -thrust::get<1>(a) - thrust::get<1>(b);
				forceZAddr[indexCenterApply] = -thrust::get<2>(a) - thrust::get<2>(b);
			}
			else {
				forceXAddr[indexCenterApply] = 0;
				forceYAddr[indexCenterApply] = 0;
				forceZAddr[indexCenterApply] = 0;
			}

			if (!isFixedVector[indexLeft]) {
				forceXAddr[indexLeftApply] = thrust::get<0>(a);
				forceYAddr[indexLeftApply] = thrust::get<1>(a);
				forceZAddr[indexLeftApply] = thrust::get<2>(a);
			}
			else {
				forceXAddr[indexLeftApply] = 0;
				forceYAddr[indexLeftApply] = 0;
				forceZAddr[indexLeftApply] = 0;
			}

			if (!isFixedVector[indexRight]){
				forceXAddr[indexRightApply] = thrust::get<0>(b);
				forceYAddr[indexRightApply] = thrust::get<1>(b);
				forceZAddr[indexRightApply] = thrust::get<2>(b);
			}
			else {
				forceXAddr[indexRightApply] = 0;
				forceYAddr[indexRightApply] = 0;
				forceZAddr[indexRightApply] = 0;
			}
		}
		else  {
			forceXAddr[indexRightApply] = 0;
			forceYAddr[indexRightApply] = 0;
			forceZAddr[indexRightApply] = 0;
			forceXAddr[indexCenterApply] = 0;
			forceYAddr[indexCenterApply] = 0;
			forceZAddr[indexCenterApply] = 0;
			forceXAddr[indexLeftApply] = 0;
			forceYAddr[indexLeftApply] = 0;
			forceZAddr[indexLeftApply] = 0;
		}

	}
};

#endif
