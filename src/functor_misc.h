#ifndef FUNCTOR_MISC_H_
#define FUNCTOR_MISC_H_

#include "SystemStructures.h"

__device__
inline CVec3 cross_product(CVec3 v1, CVec3 v2) {
	return thrust::make_tuple(thrust::get<1>(v1)*thrust::get<2>(v2) - thrust::get<2>(v1)*thrust::get<1>(v2),
		-(thrust::get<0>(v1)*thrust::get<2>(v2) - thrust::get<2>(v1)*thrust::get<0>(v2)),
		thrust::get<0>(v1)*thrust::get<1>(v2) - thrust::get<1>(v1)*thrust::get<0>(v2));
};

__device__
inline double dot_product(CVec3 v1, CVec3 v2) {
	return (thrust::get<0>(v1)*thrust::get<0>(v2) +
		thrust::get<1>(v1)*thrust::get<1>(v2) +
		thrust::get<2>(v1)*thrust::get<2>(v2));
}


struct functor_add_UCVec3_CVec3 {//same as torsion
	unsigned maxNodeCount;
	double* forceXAddr;
	double* forceYAddr;
	double* forceZAddr;

	__host__ __device__
	//
		functor_add_UCVec3_CVec3(
				unsigned& _maxNodeCount,
				double* _forceXAddr,
				double* _forceYAddr,
				double* _forceZAddr) :
			maxNodeCount(maxNodeCount),
			forceXAddr(_forceXAddr),
			forceYAddr(_forceYAddr),
			forceZAddr(_forceZAddr) {}

	__device__
	void operator() (const Tuddd& u1d3) {
			unsigned idToAssign = thrust::get<0>(u1d3);
			
			if (idToAssign < maxNodeCount) {
				if (!isnan(thrust::get<1>(u1d3)) && !isnan(thrust::get<2>(u1d3)) && !isnan(thrust::get<3>(u1d3))) {

					forceXAddr[idToAssign] += thrust::get<1>(u1d3);
					forceYAddr[idToAssign] += thrust::get<2>(u1d3);
					forceZAddr[idToAssign] += thrust::get<3>(u1d3);
				}
			}
	}

};


struct functor_add_UCVec3_CVec3_pltVol {
	unsigned maxNodeCount;
	double* forceXAddr;
	double* forceYAddr;
	double* forceZAddr;
	bool* isNodeInPltVolAddr;

	__host__ __device__
	//
		functor_add_UCVec3_CVec3_pltVol(
				unsigned& _maxNodeCount,
				double* _forceXAddr,
				double* _forceYAddr,
				double* _forceZAddr,
				bool* _isNodeInPltVolAddr) :
			maxNodeCount(maxNodeCount),
			forceXAddr(_forceXAddr),
			forceYAddr(_forceYAddr),
			forceZAddr(_forceZAddr),
			isNodeInPltVolAddr(_isNodeInPltVolAddr) {}

	__device__
	void operator() (const Tuddd& u1d3) {
			unsigned idToAssign = thrust::get<0>(u1d3);
			double f1 = thrust::get<1>(u1d3); 
			double f2 = thrust::get<2>(u1d3); 
			double f3 = thrust::get<3>(u1d3); 

			if (idToAssign < maxNodeCount) {
				if ((f1 != 0.0) || (f2 != 0.0) || (f3 !=0.0)) {

					forceXAddr[idToAssign] += thrust::get<1>(u1d3);
					forceYAddr[idToAssign] += thrust::get<2>(u1d3);
					forceZAddr[idToAssign] += thrust::get<3>(u1d3); 
					isNodeInPltVolAddr[idToAssign] = true;
				}
				if( (f1 == 0.0) && (f2 == 0.0) && (f3 == 0.0) ){
					isNodeInPltVolAddr[idToAssign] = false;
				}
			}


	}

};

struct functor_norm {

	__host__ __device__
		double operator() (const CVec3& vec) {
		//divide force by fiber cross section to get stress
		double result = sqrt(thrust::get<0>(vec) * thrust::get<0>(vec) +
			thrust::get<1>(vec) * thrust::get<1>(vec) +
			thrust::get<2>(vec) * thrust::get<2>(vec));

		return result;


	}
};


struct functor_prob_detach {//same as torsion
	double dt;
	double P;
	unsigned maxIdCountFlag;

	__host__ __device__
		//
		functor_prob_detach(
			double& _dt,
			double& _P,
			unsigned& _maxIdCountFlag):
		dt(_dt),
		P(_P),
		maxIdCountFlag(_maxIdCountFlag) {}

	__device__
		unsigned operator() (const Tud& u1d1) {
		unsigned current_node_plt_attached = thrust::get<0>(u1d1);
		double prob = thrust::get<1>(u1d1);

		if (prob < dt * P) {
			return maxIdCountFlag;
		}
		else {
			return current_node_plt_attached;
		}

	}

};

#endif
