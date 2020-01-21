#ifndef FUNCTOR_CALC_STRAIN_PARAMS_H_
#define FUNCTOR_CALC_STRAIN_PARAMS_H_

#include "system_structures.h"
struct functor_calc_strain_params : public thrust::unary_function<Tuu, Tdd> {
	unsigned originLinkCount;
	unsigned originEdgeCount;
	unsigned originNodeCount;
	unsigned max_node_count;
	unsigned max_nbr_count;

	double* locXAddr;
	double* locYAddr;
	double* locZAddr;

	unsigned* originNbrVec;
	unsigned* currentNbrVec;
	unsigned* global_neighbors;
	double* lenZero;

	__host__ __device__
	functor_calc_strain_params(
		unsigned& _originLinkCount,
		unsigned& _originEdgeCount,
		unsigned& _originNodeCount,
		unsigned& _max_node_count,
		unsigned& _max_nbr_count,	

		double* _locXAddr, 
		double* _locYAddr, 
		double* _locZAddr,

		unsigned* _originNbrVec,
		unsigned* _currentNbrVec,
		unsigned* _global_neighbors,
		double* _lenZero) :
		originLinkCount(_originLinkCount),
		originEdgeCount(_originEdgeCount),
		originNodeCount(_originNodeCount),
		max_node_count(_max_node_count),
		max_nbr_count(_max_nbr_count),

		locXAddr(_locXAddr), 
		locYAddr(_locYAddr), 
		locZAddr(_locZAddr),
		
		originNbrVec(_originNbrVec),
		currentNbrVec(_currentNbrVec),
		global_neighbors(_global_neighbors),
		lenZero(_lenZero){}
 
	__device__
	Tdd operator() (const Tuu& u2) {

		unsigned idL = thrust::get<0>(u2);//use IDL as row
		unsigned idR = thrust::get<1>(u2);//use IDR as col

		//identify matrix location
		unsigned idMat;
		for (unsigned i = idL; i < max_nbr_count; i++) {
			unsigned idCol = global_neighbors[i];//represents nbr id
			if (idR == idCol) {
				idMat = idCol;
				break;
			}
		}		

		double global_length_zero = lenZero[idMat];

		//unsigned originNbrCount = originNbrVec[idL];//number of original neighbors on idRow
		//unsigned currentNbrCount = currentNbrVec[idL];

		double edgeStrain = 0.0;
		double edgeAlignment = 0.0;
		
		if ((global_length_zero > 0.0) ) {
			double posXA_XB = locXAddr[idL] - locXAddr[idR];
			double posYA_YB = locYAddr[idL] - locYAddr[idR];
			double posZA_ZB = locZAddr[idL] - locZAddr[idR];
			double currentLength = sqrt(
				(posXA_XB) * (posXA_XB) +
				(posYA_YB) * (posYA_YB) +
				(posZA_ZB) * (posZA_ZB));
			edgeStrain = ((currentLength - global_length_zero) / global_length_zero);
			edgeAlignment = abs(1.0 * posZA_ZB / currentLength);
		}


		
		
		return thrust::make_tuple(edgeStrain, edgeAlignment);
	}
};

#endif