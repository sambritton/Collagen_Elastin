#ifndef FUNCTOR_EXTERNAL_PULL_H_
#define FUNCTOR_EXTERNAL_PULL_H_

#include "system_structures.h"

struct functor_external_pull {
	bool* isNodeFixedAddr;
	double* node_loc_x;
	double* node_loc_z;
	unsigned axis;
	double pull_ammount;
	double originalNetworkLength;
	double strain_proportion_end_sim;
	double averageLowerStrain;
	double averageUpperStrain;

	__host__ __device__
		functor_external_pull(
			bool* _isNodeFixedAddr,
			double*	_node_loc_x,
			double*	_node_loc_z,

			double& _pull_ammount,
			unsigned& _axis,
			double& _originalNetworkLength,
			double& _strain_proportion_end_sim,
			double& _averageLowerStrain,
			double& _averageUpperStrain) :
		isNodeFixedAddr(_isNodeFixedAddr),
		node_loc_x(_node_loc_x),
		node_loc_z(_node_loc_z),
		
		pull_ammount(_pull_ammount),
		axis(_axis),
		originalNetworkLength(_originalNetworkLength),
		strain_proportion_end_sim(_strain_proportion_end_sim),
		averageLowerStrain(_averageLowerStrain),
		averageUpperStrain(_averageUpperStrain) {}

	__device__
	void operator()(const Tuddbbbb& u1d2b4) {

		unsigned id = thrust::get<0>(u1d2b4);
		double locZ = thrust::get<1>(u1d2b4);
		double locX = thrust::get<2>(u1d2b4);
		bool isFixed = thrust::get<3>(u1d2b4);
		bool is_node_collagen = thrust::get<4>(u1d2b4);
		bool isUpperStrainNode = thrust::get<5>(u1d2b4);
		bool isLowerStrainNode = thrust::get<6>(u1d2b4);

		double location=0.0;
		double dirX = 0.0;//tForceX / normForce;
		double dirY = 0.0;//tForceY / normForce;
		double dirZ = 0.0;

		//pull top
		if ( isUpperStrainNode ) {
            if (axis==0){
                double new_z_pos = locZ + pull_ammount;
                node_loc_z[id] = new_z_pos;
            }            
            else if (axis!=0){
                double new_x_pos = locX + pull_ammount;
                node_loc_x[id] = new_x_pos;
            }
		}

		//pull bottom
		else if (isLowerStrainNode ) {
            if (axis==0){
                double new_z_pos = locZ - pull_ammount;
                node_loc_z[id] = new_z_pos;
            }            
            else if (axis!=0){
                double new_x_pos = locX - pull_ammount;
                node_loc_x[id] = new_x_pos;
            }
		}
	}
};

#endif
