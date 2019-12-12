#ifndef FUNCTOR_WLC_H_
#define FUNCTOR_WLC_H_

#include "SystemStructures.h"
	
struct functor_wlc {
	double* locXAddr;
	double* locYAddr;
	double* locZAddr;
	double* forceXAddr;
	double* forceYAddr;
	double* forceZAddr;

	double Kb; //convert to nN and microns
	double PLengthMon;
	double CLM;
	double Temp;
	unsigned maxNeighborCount;
	unsigned maxNodeCount;
	unsigned numMonFiberArea;

	double* lenZero;
	unsigned* edgeCountVec;
	unsigned* globalNeighbors;
	unsigned* numOriginalNeighborsVec;

	__host__ __device__

		functor_wlc(
			double* _locXAddr, 
			double* _locYAddr, 
			double* _locZAddr,
			double* _forceXAddr, 
			double* _forceYAddr, 
			double* _forceZAddr, 

			double& _Kb, 
			double& _PLengthMon, 
			double& _CLM, 
			double& _Temp,
			unsigned& _maxNeighborCount,
			unsigned& _maxNodeCount,
			unsigned& _numMonFiberArea,

			double* _lenZero,
			unsigned* _globalNeighbors,
			unsigned* _edgeCountVec,
			unsigned* _numOriginalNeighborsVec) :

		locXAddr(_locXAddr),
		locYAddr(_locYAddr),
		locZAddr(_locZAddr),
		forceXAddr(_forceXAddr),
		forceYAddr(_forceYAddr),
		forceZAddr(_forceZAddr),

		Kb(_Kb), 
		PLengthMon(_PLengthMon), 
		CLM(_CLM), 
		Temp(_Temp),
		maxNeighborCount(_maxNeighborCount),
		maxNodeCount(_maxNodeCount),
		numMonFiberArea(_numMonFiberArea),

		lenZero(_lenZero),
		globalNeighbors(_globalNeighbors),
		edgeCountVec(_edgeCountVec),
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
			
			unsigned beginIndex = idA * maxNeighborCount;
			unsigned endIndex = beginIndex + maxNeighborCount;
			

			for (unsigned i = beginIndex; i < endIndex; i++) {//currentSpringCount is the length of index and value vectors
				unsigned idB = globalNeighbors[i];//look through possible neighbors. May contain ULONG_MAX
				if (idB < maxNodeCount){
				
					double lengthZero = lenZero[i];
					if (lengthZero > 0) {
						
						double posXA_XB = locXAddr[idB] - locXAddr[idA];
						double posYA_YB = locYAddr[idB] - locYAddr[idA];
						double posZA_ZB = locZAddr[idB] - locZAddr[idA];
		
						double currentLength = sqrt(
							(posXA_XB) * (posXA_XB)+
							(posYA_YB) * (posYA_YB)+
							(posZA_ZB) * (posZA_ZB));
	
						double strain = ((currentLength - lengthZero) / lengthZero);
						if (strain>2){
							strain=2;
						}

						double dL_norm = strain / ( CLM);//CLM is unitless since it was already normalized. 
						double magForce = (numMonFiberArea*(Kb*Temp) / PLengthMon) * ( 0.25 * pow(1.0 - dL_norm, -2.0) - 0.25 + dL_norm);
				
						double magForceX = (posXA_XB / currentLength) * magForce;
						double magForceY = (posYA_YB / currentLength) * magForce;
						double magForceZ = (posZA_ZB / currentLength) * magForce;
						
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
