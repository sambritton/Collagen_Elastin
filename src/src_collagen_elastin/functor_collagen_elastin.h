#ifndef FUNCTOR_COLLAGEN_ELASTIN_H_
#define FUNCTOR_COLLAGEN_ELASTIN_H_

#include "system_structures.h"

struct functor_collagen_elastin {
	double* locXAddr;
	double* locYAddr;
	double* locZAddr;
	double* forceXAddr;
	double* forceYAddr;
	double* forceZAddr;

	double collagen_spring_constant;
	double Kb; //convert to nN and microns
	double PLengthMon;
	double CLM;
	double Temp;
	unsigned max_nbr_count;
	unsigned max_node_count;
	double num_mon_elastin_area;

	double* lenZero;
	unsigned* edgeCountVec;
	unsigned* global_neighbors;
	bool* global_is_edge_collagen;
	bool* global_is_edge_elastin;
	unsigned* numOriginalNeighborsVec;

	__host__ __device__

		functor_collagen_elastin(
			double* _locXAddr,
			double* _locYAddr,
			double* _locZAddr,
			double* _forceXAddr,
			double* _forceYAddr,
			double* _forceZAddr,

			double& _collagen_spring_constant,
			double& _Kb,
			double& _PLengthMon,
			double& _CLM,
			double& _Temp,
			unsigned& _max_nbr_count,
			unsigned& _max_node_count,
			double& _num_mon_elastin_area,

			double* _lenZero,
			unsigned* _edgeCountVec,
			unsigned* _global_neighbors,
			bool* _global_is_edge_collagen,
			bool* _global_is_edge_elastin,
			unsigned* _numOriginalNeighborsVec) :

		locXAddr(_locXAddr),
		locYAddr(_locYAddr),
		locZAddr(_locZAddr),
		forceXAddr(_forceXAddr),
		forceYAddr(_forceYAddr),
		forceZAddr(_forceZAddr),

		collagen_spring_constant(_collagen_spring_constant),
		Kb(_Kb),
		PLengthMon(_PLengthMon),
		CLM(_CLM),
		Temp(_Temp),
		max_nbr_count(_max_nbr_count),
		max_node_count(_max_node_count),
		num_mon_elastin_area(_num_mon_elastin_area),

		lenZero(_lenZero),
		edgeCountVec(_edgeCountVec),
		global_neighbors(_global_neighbors),
		global_is_edge_collagen(_global_is_edge_collagen),
		global_is_edge_elastin(_global_is_edge_elastin),
		numOriginalNeighborsVec(_numOriginalNeighborsVec) {}

	__device__
	void operator()(const Tub& u1b1) {
		//idA represents row.
		unsigned idA = thrust::get<0>(u1b1);
		//unsigned numOriginalNeighbors = numOriginalNeighborsVec[idA];

		bool isFixed = thrust::get<1>(u1b1);
		double sumForceX = 0;
		double sumForceY = 0;
		double sumForceZ = 0;

		if (!isFixed) {
			//only apply force if not fixed.

			unsigned beginIndex = idA * max_nbr_count;
			unsigned endIndex = beginIndex + max_nbr_count;


			for (unsigned i = beginIndex; i < endIndex; i++) {//currentSpringCount is the length of index and value vectors
				unsigned idB = global_neighbors[i];//look through possible neighbors. May contain ULONG_MAX
				if (idB < max_node_count){
					if (idA == 0 || idA == 1){
						double no_var=1;
					}
					double length_zero = lenZero[i];
					bool is_edge_collagen = global_is_edge_collagen[i];
					bool is_edge_elastin = global_is_edge_elastin[i];
					if (length_zero > 0) {

						double posXA_XB = locXAddr[idB] - locXAddr[idA];
						double posYA_YB = locYAddr[idB] - locYAddr[idA];
						double posZA_ZB = locZAddr[idB] - locZAddr[idA];

						double length_current = sqrt(
							(posXA_XB) * (posXA_XB)+
							(posYA_YB) * (posYA_YB)+
							(posZA_ZB) * (posZA_ZB));

						double magForceX = 0.0;
						double magForceY = 0.0;
						double magForceZ = 0.0;

						if (is_edge_elastin){
							//apply wlc force
							double strain = ((length_current - length_zero) / length_zero);

							double dL_norm = strain / ( CLM);//CLM is unitless since it was already normalized.
							double magForce = (num_mon_elastin_area*(Kb*Temp) / PLengthMon) * ( 0.25 * pow(1.0 - dL_norm, -2.0) - 0.25 + dL_norm);

							magForceX = (posXA_XB / length_current) * magForce;
							magForceY = (posYA_YB / length_current) * magForce;
							magForceZ = (posZA_ZB / length_current) * magForce;
						}
						else if (is_edge_collagen){
							//apply linear spring force

							double strain = ((length_current - length_zero) / length_zero);

							double dL_norm = strain / ( CLM);//CLM is unitless since it was already normalized.
							double plength = PLengthMon/5.0;
							double magForce = (num_mon_elastin_area*(Kb*Temp) / plength) * ( 0.25 * pow(1.0 - dL_norm, -2.0) - 0.25 + dL_norm);

							magForceX = (posXA_XB / length_current) * magForce;
							magForceY = (posYA_YB / length_current) * magForce;
							magForceZ = (posZA_ZB / length_current) * magForce;
							//double magnitude = collagen_spring_constant * (length_current - length_zero);
							//magForceX = magnitude * (posXA_XB/length_current);
							//magForceY = magnitude * (posYA_YB/length_current);
							//magForceZ = magnitude * (posZA_ZB/length_current);
						}

						sumForceX += magForceX;
						sumForceY += magForceY;
						sumForceZ += magForceZ;

					}
				}
			}


			if (isfinite(sumForceX))
				forceXAddr[idA] += sumForceX;

			if (isfinite(sumForceY))
				forceYAddr[idA] += sumForceY;

			if (isfinite(sumForceY))
				forceZAddr[idA] += sumForceZ;
		}

	}
};
#endif
