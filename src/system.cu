#include <thrust/system_error.h>
#include <thrust/binary_search.h>
#include <thrust/reduce.h>
#include <algorithm>
#include <thrust/replace.h>
#include <thrust/unique.h>
#include <thrust/gather.h>
#include <ostream>
#include <stdio.h>
#include <curand.h>
#include <curand_kernel.h>
#include <thrust/sort.h>
#include <thrust/transform_reduce.h>
#include <math.h>

#include "storage.h"
#include "system_builder.h"
#include "collagen_elastin_spring.h"
#include "bending_spring.h"
#include "advance_positions.h"
#include "bucket_scheme.h"
#include "link_nodes.h"
#include "external_force.h"
#include "system.h"
#include "functor_misc.h"

using namespace thrust::placeholders;

void System::set_bucket_scheme(){
	init_dim_general(nodeInfoVecs, domainParams, auxVecs, generalParams);
	init_net_inct_bucket(nodeInfoVecs, domainParams, auxVecs, generalParams);
	build_net_inct_bucket(nodeInfoVecs, domainParams, auxVecs, generalParams);
	extend_net_inct_bucket(nodeInfoVecs, domainParams, auxVecs, generalParams);
}

void System::solve_forces() {

	thrust::fill(nodeInfoVecs.node_force_x.begin(),nodeInfoVecs.node_force_x.end(),0);
	thrust::fill(nodeInfoVecs.node_force_y.begin(),nodeInfoVecs.node_force_y.end(),0);
	thrust::fill(nodeInfoVecs.node_force_z.begin(),nodeInfoVecs.node_force_z.end(),0);
	
	if (generalParams.linking == true) {
		link_nodes(nodeInfoVecs, edgeInfoVecs, auxVecs, generalParams);
	}
		
	calc_bending_spring_force(nodeInfoVecs, bendInfoVecs, generalParams);
	  
	calc_spring_force(nodeInfoVecs, edgeInfoVecs, generalParams);
	
	
    //std::cout<<"pre norm: " << std::flush;
	extensionParams.totalAppliedForce = thrust::transform_reduce(
		thrust::make_zip_iterator(
			thrust::make_tuple(
				nodeInfoVecs.node_force_x.begin(),
				nodeInfoVecs.node_force_y.begin(),
				nodeInfoVecs.node_force_z.begin())),
		thrust::make_zip_iterator(
			thrust::make_tuple(
				nodeInfoVecs.node_force_x.begin(),
				nodeInfoVecs.node_force_y.begin(),
				nodeInfoVecs.node_force_z.begin())) + generalParams.max_node_count,
			functor_norm(), 0.0, thrust::plus<double>() );
			thrust::transform(
				thrust::make_zip_iterator(
					thrust::make_tuple(
						nodeInfoVecs.node_force_x.begin(),
						nodeInfoVecs.node_force_y.begin(),
						nodeInfoVecs.node_force_z.begin())),
				thrust::make_zip_iterator(
					thrust::make_tuple(
						nodeInfoVecs.node_force_x.begin(),
						nodeInfoVecs.node_force_y.begin(),
						nodeInfoVecs.node_force_z.begin())) + generalParams.max_node_count,
				nodeInfoVecs.sum_forces_on_node.begin(),//save vector
				functor_norm());
	
	extensionParams.applied_force_upper = thrust::transform_reduce(
												thrust::make_zip_iterator(
													thrust::make_tuple(
														nodeInfoVecs.node_upper_selection_pull.begin(),
														nodeInfoVecs.sum_forces_on_node.begin())),
												thrust::make_zip_iterator(
													thrust::make_tuple(
														nodeInfoVecs.node_upper_selection_pull.begin(),
														nodeInfoVecs.sum_forces_on_node.begin())) + generalParams.max_node_count,
												functor_sum_pulled_forces(), 0.0, thrust::plus<double>());
	
	extensionParams.applied_force_lower = thrust::transform_reduce(
												thrust::make_zip_iterator(
													thrust::make_tuple(
														nodeInfoVecs.node_lower_selection_pull.begin(),
														nodeInfoVecs.sum_forces_on_node.begin())),
													thrust::make_zip_iterator(
														thrust::make_tuple(
															nodeInfoVecs.node_lower_selection_pull.begin(),
															nodeInfoVecs.sum_forces_on_node.begin())) + generalParams.max_node_count,
												functor_sum_pulled_forces(), 0.0, thrust::plus<double>());
};


void System::solve_system() {

	double lastTime = 0.0;
	bool runIters = true;
	std::cout << " setting initial bucket scheme "  << std::endl;
	std::cout<<"magnitudeForce at beginning of simulation: "<< generalParams.magnitudeForce<<std::endl;

	//set initial epsilon
	generalParams.epsilon = (generalParams.epsilon_factor) *
		sqrt(6.0 * edgeInfoVecs.kB * edgeInfoVecs.temperature * generalParams.dt / edgeInfoVecs.viscosity_elastin);
	std::cout<< "Initial Epsilon: " << generalParams.epsilon << std::endl;
	//initialize images
	storage->print_VTK_file();
	storage->save_params();
	
	set_bucket_scheme();
	external_force(
		nodeInfoVecs,
		generalParams,
		extensionParams,
		domainParams);//set initial step and strain parameters
	std::cout<<"starting system" << std::flush;
	while (runIters == true) {

		generalParams.iterationCounter++;
		generalParams.currentTime += generalParams.dt;
		//if (generalParams.iterationCounter % 50 == 0){
			//std::cout << "current iter: " <<generalParams.iterationCounter<<  std::endl;
		set_bucket_scheme();
		//}
		advance_positions(
			nodeInfoVecs,
			generalParams,
			edgeInfoVecs,
      		randVecs);
		
		solve_forces(); //resets and solves forces for next time step
		double maxVel = *(thrust::max_element(nodeInfoVecs.node_vel.begin(), nodeInfoVecs.node_vel.end()));

		thrust::device_vector<double>::iterator iter = thrust::max_element(nodeInfoVecs.node_vel.begin(), nodeInfoVecs.node_vel.end());
		unsigned position = iter - nodeInfoVecs.node_vel.begin();
		double max_val = *iter;
		
		thrust::device_vector<double>::iterator iterx = thrust::max_element(nodeInfoVecs.node_loc_x.begin(), nodeInfoVecs.node_loc_x.end());
		double max_x = *iterx;
		std::cout<<"max x: " << max_x << std::endl;

		if ((generalParams.iterationCounter % 1000) == 0) {
			double currentStrain = (extensionParams.averageUpperStrain - extensionParams.averageLowerStrain) /
			(extensionParams.originAverageUpperStrain - extensionParams.originAverageLowerStrain ) - 1.0;
			if (currentStrain>4.0){
				runIters=false;
			}		
			std::cout<<" current strain: " << currentStrain << std::endl;
			std::cout<<" max velocity: " << maxVel << std::endl;
			std::cout<<" epsilon: " << generalParams.epsilon << std::endl;
			std::cout << " extensionParams.averageUpperStrain: " << extensionParams.averageUpperStrain << std::endl;
			std::cout << " extensionParams.averageLowerStrain: " << extensionParams.averageLowerStrain << std::endl;
			
			std::cout << " extensionParams.applied_force_lower: " << extensionParams.applied_force_lower << std::endl;
			std::cout << " extensionParams.applied_force_upper: " << extensionParams.applied_force_upper << std::endl;
		
		}
		//store sum of all forces on each node. Used in stress calculations

		if ((generalParams.iterationCounter % 5000) == 0) {								
			storage->print_VTK_file();
			storage->save_params();
		}

		if ((maxVel < generalParams.epsilon) && (generalParams.iterationCounter % 50 == 0)) {
			//perform pulling
			external_force(
				nodeInfoVecs,
				generalParams,
				extensionParams,
				domainParams);

			std::cout<<"Maximum vel: "<< maxVel <<std::endl;
			generalParams.magnitudeForce += generalParams.df;
			
		}
		///////////////////////////////////////////////////////////////////////////////
		//EQUILIBRIUM END
		//////////////////////////////////////////////////////////////////////

	}

};

System::System()  {};

void System::assign_storage(std::shared_ptr<Storage> _storage) {
	storage = _storage;
}

void System::initialize_system(HostNodeInfoVecs& hostNodeInfoVecs) {

	std::cout<< "total Edge Count: "<< generalParams.origin_edge_count << std::endl;
	std::cout << "max num nodes: " << generalParams.max_node_count << std::endl;

	nodeInfoVecs.origin_edge_left = hostNodeInfoVecs.host_spring_edge_left;
	nodeInfoVecs.origin_edge_right = hostNodeInfoVecs.host_spring_edge_right;

	set_node_vecs(//calls initDimensionBucketScheme
		hostNodeInfoVecs);

	set_bend_vecs(hostNodeInfoVecs);

	set_edge_vecs(hostNodeInfoVecs);

	set_extras();
};


void System::set_node_vecs(
	HostNodeInfoVecs& hostNodeInfoVecs) {

	randVecs.gaussianData.resize(generalParams.max_node_count);

	nodeInfoVecs.id_edges_made_temp.resize(generalParams.max_node_count * generalParams.max_links_per_iteration);//corresponds to upperAdj vector size plus a single value to hold number of added nodes
	thrust::fill(nodeInfoVecs.id_edges_made_temp.begin(), nodeInfoVecs.id_edges_made_temp.end(), 0);

	nodeInfoVecs.sum_forces_on_node.resize(generalParams.max_node_count);

	nodeInfoVecs.node_upper_selection_pull.resize(generalParams.max_node_count);
	nodeInfoVecs.node_lower_selection_pull.resize(generalParams.max_node_count);

	nodeInfoVecs.node_vel.resize(generalParams.max_node_count);

  	nodeInfoVecs.node_is_collagen.resize(generalParams.max_node_count);
  	nodeInfoVecs.node_is_elastin.resize(generalParams.max_node_count);

	nodeInfoVecs.node_loc_x.resize(generalParams.max_node_count);
	nodeInfoVecs.node_loc_y.resize(generalParams.max_node_count);
	nodeInfoVecs.node_loc_z.resize(generalParams.max_node_count);
	nodeInfoVecs.node_vel_x.resize(generalParams.max_node_count);
	nodeInfoVecs.node_vel_y.resize(generalParams.max_node_count);
	nodeInfoVecs.node_vel_z.resize(generalParams.max_node_count);


	nodeInfoVecs.node_force_x.resize(generalParams.max_node_count);
	nodeInfoVecs.node_force_y.resize(generalParams.max_node_count);
	nodeInfoVecs.node_force_z.resize(generalParams.max_node_count);

	nodeInfoVecs.discretized_edges_strain.resize(generalParams.max_node_count * generalParams.max_nbr_count);
	nodeInfoVecs.discretized_edges_alignment.resize(generalParams.max_node_count * generalParams.max_nbr_count);

	//sized larger for input later
	
	nodeInfoVecs.device_edge_left.resize(generalParams.max_node_count * generalParams.max_nbr_count);
	nodeInfoVecs.device_edge_right.resize(generalParams.max_node_count * generalParams.max_nbr_count);

	nodeInfoVecs.host_edge_left.resize(generalParams.max_node_count * generalParams.max_nbr_count);
	nodeInfoVecs.host_edge_right.resize(generalParams.max_node_count * generalParams.max_nbr_count);


	thrust::fill(nodeInfoVecs.discretized_edges_strain.begin(), nodeInfoVecs.discretized_edges_strain.end(),0.0);
	thrust::fill(nodeInfoVecs.host_edge_right.begin(), nodeInfoVecs.host_edge_right.end(), 0);	//fill force and velocity with zeros for computation.
	thrust::fill(nodeInfoVecs.host_edge_left.begin(), nodeInfoVecs.host_edge_left.end(), 0);	//fill force and velocity with zeros for computation.
	thrust::fill(nodeInfoVecs.id_edges_made_temp.begin(), nodeInfoVecs.id_edges_made_temp.end(), 0);

	thrust::fill(nodeInfoVecs.sum_forces_on_node.begin(), nodeInfoVecs.sum_forces_on_node.end(), 0);

	thrust::fill(nodeInfoVecs.node_upper_selection_pull.begin(),
		nodeInfoVecs.node_upper_selection_pull.end(),false);

	thrust::fill(nodeInfoVecs.node_lower_selection_pull.begin(),
		nodeInfoVecs.node_lower_selection_pull.end(),false);

	thrust::copy(hostNodeInfoVecs.host_node_is_collagen.begin(), hostNodeInfoVecs.host_node_is_collagen.end(), nodeInfoVecs.node_is_collagen.begin());
	thrust::copy(hostNodeInfoVecs.host_node_is_elastin.begin(), hostNodeInfoVecs.host_node_is_elastin.end(), nodeInfoVecs.node_is_elastin.begin());
	thrust::copy(hostNodeInfoVecs.host_pos_x.begin(), hostNodeInfoVecs.host_pos_x.end(), nodeInfoVecs.node_loc_x.begin());
	thrust::copy(hostNodeInfoVecs.host_pos_y.begin(), hostNodeInfoVecs.host_pos_y.end(), nodeInfoVecs.node_loc_y.begin());
	thrust::copy(hostNodeInfoVecs.host_pos_z.begin(), hostNodeInfoVecs.host_pos_z.end(), nodeInfoVecs.node_loc_z.begin());

	nodeInfoVecs.links_made_individual_thread.resize(generalParams.max_node_count);

	nodeInfoVecs.id_temp_linked_left.resize(generalParams.max_node_count * generalParams.max_links_per_iteration);
	nodeInfoVecs.id_temp_linked_right.resize(generalParams.max_node_count * generalParams.max_links_per_iteration);
	//copy fixed positions
	nodeInfoVecs.host_id_left.resize(generalParams.max_node_count * generalParams.max_links_per_iteration);
	nodeInfoVecs.host_id_right.resize(generalParams.max_node_count * generalParams.max_links_per_iteration);

	nodeInfoVecs.is_node_fixed.resize(generalParams.max_node_count);
	thrust::fill(nodeInfoVecs.is_node_fixed.begin(), nodeInfoVecs.is_node_fixed.end(), false);

	//now that all the nodes are loaded in, choose the top to apply strain, and fix the bottom

	determine_bounds();

	//at this point all nodes are filled, so we can generate domainParams before seeding dpd particles.
	init_dim_general(
		nodeInfoVecs,
		domainParams,
		auxVecs,
		generalParams);

	//set original parameters for domain. others will be reset as simulation takes place.
	domainParams.origin_min_x = domainParams.min_x;
	domainParams.origin_max_x = domainParams.max_x;
	domainParams.origin_min_y = domainParams.min_y;
	domainParams.origin_max_y = domainParams.max_y;
	domainParams.origin_min_z = domainParams.min_z;
	domainParams.origin_max_z = domainParams.max_z;
	std::cout<< "node count : " <<nodeInfoVecs.node_loc_y.size()<< std::endl;


	auxVecs.id_bucket_net_intc.resize(generalParams.max_node_count);
	auxVecs.id_value_net_intc.resize(generalParams.max_node_count);
	auxVecs.id_value_expanded_net_intc.resize(27 * (generalParams.max_node_count));
	auxVecs.id_bucket_expanded_net_intc.resize(27 *( generalParams.max_node_count));

};

void System::determine_bounds() {
	//determin z positions of nodes to be pulled and fixed.

	thrust::device_vector<double> pos_temp;
	pos_temp.resize(generalParams.max_node_count);
	if (extensionParams.axis == 0) {
		thrust::copy(nodeInfoVecs.node_loc_z.begin(), nodeInfoVecs.node_loc_z.end(), pos_temp.begin());
	}else{
		thrust::copy(nodeInfoVecs.node_loc_x.begin(), nodeInfoVecs.node_loc_x.end(), pos_temp.begin());
	}


	thrust::device_vector<double>::iterator iter_max = thrust::max_element(pos_temp.begin(), pos_temp.end());
	thrust::device_vector<double>::iterator iter_min = thrust::min_element(pos_temp.begin(), pos_temp.end());
	//unsigned position_max = iter - pos_temp.begin();
	double max_val = *iter_max;
	double min_val = *iter_min;
	
	
	double length = max_val - min_val;
	std::cout<<"start end pos_temp: "<< max_val << " "<< min_val <<std::endl;

	//upperLevelAlt pulls 10% default. Set in main.cpp using input
	if (generalParams.pull_percent >= 1.0 || generalParams.pull_percent < 0.0) {
		generalParams.pull_percent=1.0;
		std::cout<<"ERROR PULL PERCENT MUST BE LESS THAN ONE AND LARGER THAN 0.0"<<std::endl;
	}
	double pull_width = generalParams.pull_percent * length;
	double upperLevelAlt = max_val - pull_width;


	double lowerLevel = min_val + pull_width;

	std::cout<<"minimal level final choice for strain choice: " << lowerLevel <<std::endl;

	std::cout<<"maximal level final choice for strain choice: " << upperLevelAlt <<std::endl;

	//apply strain only to original nodes and not added edge subdivision nodes. Set top and bottom

	thrust::replace_if(nodeInfoVecs.node_upper_selection_pull.begin(), nodeInfoVecs.node_upper_selection_pull.begin() + generalParams.origin_node_count,
						pos_temp.begin(),
						IsGreaterThanLevel( upperLevelAlt ), true);

	thrust::replace_if(nodeInfoVecs.node_lower_selection_pull.begin(), nodeInfoVecs.node_lower_selection_pull.begin() + generalParams.origin_node_count,
						pos_temp.begin(),
						IsLessThanLevel( lowerLevel ), true);

	generalParams.numUpperStrainNodes_collagen = thrust::transform_reduce(		
													thrust::make_zip_iterator(
														thrust::make_tuple(
															nodeInfoVecs.node_upper_selection_pull.begin(),
															nodeInfoVecs.node_is_collagen.begin())),
															
													thrust::make_zip_iterator(
														thrust::make_tuple(
															nodeInfoVecs.node_upper_selection_pull.begin(),
															nodeInfoVecs.node_is_collagen.begin())) + generalParams.max_node_count,
													IsEqualToOne_and_node_type(),					
													0, thrust::plus<unsigned>());

	generalParams.numUpperStrainNodes_elastin = thrust::transform_reduce(		
													thrust::make_zip_iterator(
														thrust::make_tuple(
															nodeInfoVecs.node_upper_selection_pull.begin(),
															nodeInfoVecs.node_is_elastin.begin())),
															
													thrust::make_zip_iterator(
														thrust::make_tuple(
															nodeInfoVecs.node_upper_selection_pull.begin(),
															nodeInfoVecs.node_is_elastin.begin())) + generalParams.max_node_count,
													IsEqualToOne_and_node_type(),					
													0, thrust::plus<unsigned>());								


	generalParams.numLowerStrainNodes_collagen = thrust::transform_reduce(		
													thrust::make_zip_iterator(
														thrust::make_tuple(
															nodeInfoVecs.node_lower_selection_pull.begin(),
															nodeInfoVecs.node_is_collagen.begin())),
															
													thrust::make_zip_iterator(
														thrust::make_tuple(
															nodeInfoVecs.node_lower_selection_pull.begin(),
															nodeInfoVecs.node_is_collagen.begin())) + generalParams.max_node_count,
													IsEqualToOne_and_node_type(),					
													0, thrust::plus<unsigned>());

	generalParams.numLowerStrainNodes_elastin = thrust::transform_reduce(		
													thrust::make_zip_iterator(
														thrust::make_tuple(
															nodeInfoVecs.node_lower_selection_pull.begin(),
															nodeInfoVecs.node_is_elastin.begin())),
															
													thrust::make_zip_iterator(
														thrust::make_tuple(
															nodeInfoVecs.node_lower_selection_pull.begin(),
															nodeInfoVecs.node_is_elastin.begin())) + generalParams.max_node_count,
													IsEqualToOne_and_node_type(),					
													0, thrust::plus<unsigned>());		

	std::cout<<"first try numLowerStrainNodes_elastin: " << generalParams.numLowerStrainNodes_elastin <<std::endl;
	std::cout<<"first try numUpperStrainNodes_elastin: " << generalParams.numUpperStrainNodes_elastin <<std::endl;
	std::cout<<"first try numLowerStrainNodes_collagen: " << generalParams.numLowerStrainNodes_collagen <<std::endl;
	std::cout<<"first try numUpperStrainNodes_collagen: " << generalParams.numUpperStrainNodes_collagen <<std::endl;
													
	//We need to pull the same number of collagen and elastin nodes from the top and bottom. 
	unsigned num_more_collagen = 0;
	if (generalParams.numLowerStrainNodes_collagen > generalParams.numUpperStrainNodes_collagen){
		num_more_collagen = generalParams.numLowerStrainNodes_collagen - generalParams.numUpperStrainNodes_collagen;
	}else { num_more_collagen = generalParams.numUpperStrainNodes_collagen - generalParams.numLowerStrainNodes_collagen; }

	for (unsigned i = 0; i < num_more_collagen; i++){
		for (unsigned j = 0; j < generalParams.max_node_count; j++){
			bool is_collagen = nodeInfoVecs.node_is_collagen[j];
			if (is_collagen){
				if (generalParams.numLowerStrainNodes_collagen > generalParams.numUpperStrainNodes_collagen){		
					bool is_pulled = nodeInfoVecs.node_lower_selection_pull[j];
					if (is_pulled){
						nodeInfoVecs.node_lower_selection_pull[j]=false;
						break;
					}
				}else{
					bool is_pulled = nodeInfoVecs.node_upper_selection_pull[j];
					if (is_pulled){
						nodeInfoVecs.node_upper_selection_pull[j]=false;
						break;
					}
				}
			}
		}
	}
	unsigned num_more_elastin = 0;
	if (generalParams.numLowerStrainNodes_elastin > generalParams.numUpperStrainNodes_elastin){
		num_more_elastin = generalParams.numLowerStrainNodes_elastin - generalParams.numUpperStrainNodes_elastin;
	}else { num_more_elastin = generalParams.numUpperStrainNodes_elastin - generalParams.numLowerStrainNodes_elastin; }
	

	for (unsigned i = 0; i < num_more_elastin; i++){
		for (unsigned j = 0; j < generalParams.max_node_count; j++){
			bool is_elastin = nodeInfoVecs.node_is_elastin[j];
			if (is_elastin){
				if (generalParams.numLowerStrainNodes_elastin > generalParams.numUpperStrainNodes_elastin){		
					bool is_pulled = nodeInfoVecs.node_lower_selection_pull[j];
					if (is_pulled){
						nodeInfoVecs.node_lower_selection_pull[j]=false;
						break;
					}
				}else{
					bool is_pulled = nodeInfoVecs.node_upper_selection_pull[j];
					if (is_pulled){
						nodeInfoVecs.node_upper_selection_pull[j]=false;
						break;
					}
				}
			}
		}
	}
	
	//Test Fix pulled selection and move them 
	/*for (unsigned i = 0; i < nodeInfoVecs.node_lower_selection_pull.size(); i++){
		bool is_upper_pulled = nodeInfoVecs.node_upper_selection_pull[i];
		bool is_lower_pulled = nodeInfoVecs.node_lower_selection_pull[i];
		if (is_upper_pulled || is_lower_pulled){
			nodeInfoVecs.is_node_fixed[i] = true;
			std::cout<<"fixing: " << i << std::endl;
		}
	}*/
	generalParams.numUpperStrainNodes_collagen = thrust::transform_reduce(		
		thrust::make_zip_iterator(
			thrust::make_tuple(
				nodeInfoVecs.node_upper_selection_pull.begin(),
				nodeInfoVecs.node_is_collagen.begin())),
				
		thrust::make_zip_iterator(
			thrust::make_tuple(
				nodeInfoVecs.node_upper_selection_pull.begin(),
				nodeInfoVecs.node_is_collagen.begin())) + generalParams.max_node_count,
		IsEqualToOne_and_node_type(),					
		0, thrust::plus<unsigned>());

	generalParams.numUpperStrainNodes_elastin = thrust::transform_reduce(		
			thrust::make_zip_iterator(
				thrust::make_tuple(
					nodeInfoVecs.node_upper_selection_pull.begin(),
					nodeInfoVecs.node_is_elastin.begin())),
					
			thrust::make_zip_iterator(
				thrust::make_tuple(
					nodeInfoVecs.node_upper_selection_pull.begin(),
					nodeInfoVecs.node_is_elastin.begin())) + generalParams.max_node_count,
			IsEqualToOne_and_node_type(),					
			0, thrust::plus<unsigned>());								


	generalParams.numLowerStrainNodes_collagen = thrust::transform_reduce(		
			thrust::make_zip_iterator(
				thrust::make_tuple(
					nodeInfoVecs.node_lower_selection_pull.begin(),
					nodeInfoVecs.node_is_collagen.begin())),
					
			thrust::make_zip_iterator(
				thrust::make_tuple(
					nodeInfoVecs.node_lower_selection_pull.begin(),
					nodeInfoVecs.node_is_collagen.begin())) + generalParams.max_node_count,
			IsEqualToOne_and_node_type(),					
			0, thrust::plus<unsigned>());

	generalParams.numLowerStrainNodes_elastin = thrust::transform_reduce(		
			thrust::make_zip_iterator(
				thrust::make_tuple(
					nodeInfoVecs.node_lower_selection_pull.begin(),
					nodeInfoVecs.node_is_elastin.begin())),
					
			thrust::make_zip_iterator(
				thrust::make_tuple(
					nodeInfoVecs.node_lower_selection_pull.begin(),
					nodeInfoVecs.node_is_elastin.begin())) + generalParams.max_node_count,
			IsEqualToOne_and_node_type(),					
			0, thrust::plus<unsigned>());		
	
	generalParams.numUpperStrainNodes=generalParams.numUpperStrainNodes_elastin + generalParams.numUpperStrainNodes_collagen;
	generalParams.numLowerStrainNodes=generalParams.numLowerStrainNodes_elastin + generalParams.numLowerStrainNodes_collagen;

	std::cout<<"numLowerStrainNodes_elastin: " << generalParams.numLowerStrainNodes_elastin <<std::endl;
	std::cout<<"numUpperStrainNodes_elastin: " << generalParams.numUpperStrainNodes_elastin <<std::endl;
	std::cout<<"numLowerStrainNodes_collagen: " << generalParams.numLowerStrainNodes_collagen <<std::endl;
	std::cout<<"numUpperStrainNodes_collagen: " << generalParams.numUpperStrainNodes_collagen <<std::endl;
	
	unsigned numFixed = thrust::count_if(nodeInfoVecs.is_node_fixed.begin(),nodeInfoVecs.is_node_fixed.end(), IsEqualToOne() );
	std::cout<<"number of nodes fixed: " << numFixed <<std::endl;
	pos_temp.resize(0);

}

void System::set_bend_vecs(
	HostNodeInfoVecs& hostNodeInfoVecs) {

	bendInfoVecs.leftIndex.resize(bendInfoVecs.bend_factor * bendInfoVecs.total_bend_count);
	bendInfoVecs.centerIndex.resize(bendInfoVecs.bend_factor * bendInfoVecs.total_bend_count);
	bendInfoVecs.rightIndex.resize(bendInfoVecs.bend_factor * bendInfoVecs.total_bend_count);
	bendInfoVecs.angleZero.resize(bendInfoVecs.bend_factor * bendInfoVecs.total_bend_count);

	thrust::fill(bendInfoVecs.leftIndex.begin(),bendInfoVecs.leftIndex.end(),ULONG_MAX);
	thrust::fill(bendInfoVecs.centerIndex.begin(),bendInfoVecs.centerIndex.end(),ULONG_MAX);
	thrust::fill(bendInfoVecs.rightIndex.begin(),bendInfoVecs.rightIndex.end(),ULONG_MAX);

	//after default value is set, set the real id's
	thrust::copy(hostNodeInfoVecs.host_torsion_index_left.begin(), hostNodeInfoVecs.host_torsion_index_left.end(), bendInfoVecs.leftIndex.begin());
	thrust::copy(hostNodeInfoVecs.host_torsion_index_center.begin(), hostNodeInfoVecs.host_torsion_index_center.end(), bendInfoVecs.centerIndex.begin());
	thrust::copy(hostNodeInfoVecs.host_torsion_index_right.begin(), hostNodeInfoVecs.host_torsion_index_right.end(), bendInfoVecs.rightIndex.begin());

	thrust::transform(
		thrust::make_zip_iterator(
			thrust::make_tuple(
				bendInfoVecs.leftIndex.begin(),
				bendInfoVecs.centerIndex.begin(),
				bendInfoVecs.rightIndex.begin())),
		thrust::make_zip_iterator(
			thrust::make_tuple(
				bendInfoVecs.leftIndex.begin(),
				bendInfoVecs.centerIndex.begin(),
				bendInfoVecs.rightIndex.begin())) + bendInfoVecs.total_bend_count,
			bendInfoVecs.angleZero.begin(),//save vector
		functor_initial_angle(
			thrust::raw_pointer_cast(nodeInfoVecs.node_loc_x.data()),
			thrust::raw_pointer_cast(nodeInfoVecs.node_loc_y.data()),
			thrust::raw_pointer_cast(nodeInfoVecs.node_loc_z.data())));

	//		std::cout<<" in NSD device values"<<std::endl;
	/*for (unsigned i = 0; i<bendInfoVecs.total_bend_count; i++) {
		unsigned n0 = bendInfoVecs.leftIndex[i];
		unsigned n1 = bendInfoVecs.centerIndex[i];
		unsigned n2 = bendInfoVecs.rightIndex[i];
		std::cout<< "angle : "<< n0<< " " << n1<< " " << n2<< " " << bendInfoVecs.angleZero[i]<<std::endl;
	}*/

	//3x bigger since each spring affects 3 nodes.
	bendInfoVecs.forceX.resize(bendInfoVecs.bend_factor * bendInfoVecs.total_bend_count);
	bendInfoVecs.forceY.resize(bendInfoVecs.bend_factor * bendInfoVecs.total_bend_count);
	bendInfoVecs.forceZ.resize(bendInfoVecs.bend_factor * bendInfoVecs.total_bend_count);
	bendInfoVecs.tempForceX.resize(bendInfoVecs.bend_factor * bendInfoVecs.total_bend_count);
	bendInfoVecs.tempForceY.resize(bendInfoVecs.bend_factor * bendInfoVecs.total_bend_count);
	bendInfoVecs.tempForceZ.resize(bendInfoVecs.bend_factor * bendInfoVecs.total_bend_count);


	thrust::fill(bendInfoVecs.forceX.begin(), bendInfoVecs.forceX.end(), 0.0);
	thrust::fill(bendInfoVecs.forceY.begin(), bendInfoVecs.forceY.end(), 0.0);
	thrust::fill(bendInfoVecs.forceZ.begin(), bendInfoVecs.forceZ.end(), 0.0);

	bendInfoVecs.tempTorIndices.resize(bendInfoVecs.bend_factor * bendInfoVecs.total_bend_count);
	bendInfoVecs.reducedIds.resize(bendInfoVecs.bend_factor * bendInfoVecs.total_bend_count);


};

void System::set_edge_vecs(
	HostNodeInfoVecs& hostNodeInfoVecs ) {

	edgeInfoVecs.global_neighbors.resize(generalParams.max_node_count * generalParams.max_nbr_count);
  	edgeInfoVecs.global_isedge_collagen.resize(generalParams.max_node_count * generalParams.max_nbr_count);
  	edgeInfoVecs.global_isedge_elastin.resize(generalParams.max_node_count * generalParams.max_nbr_count);

	edgeInfoVecs.current_node_edge_count_vec.resize(generalParams.max_node_count);

	edgeInfoVecs.global_length_zero.resize(generalParams.max_node_count * generalParams.max_nbr_count);
	edgeInfoVecs.num_origin_nbr_per_node_vec.resize(generalParams.max_node_count);


  	thrust::fill(edgeInfoVecs.global_neighbors.begin(), edgeInfoVecs.global_neighbors.end(), generalParams.max_node_count);
  	thrust::fill(edgeInfoVecs.global_isedge_collagen.begin(), edgeInfoVecs.global_isedge_collagen.end(), false);
  	thrust::fill(edgeInfoVecs.global_isedge_elastin.begin(), edgeInfoVecs.global_isedge_elastin.end(), false);

  	thrust::fill(edgeInfoVecs.current_node_edge_count_vec.begin(), edgeInfoVecs.current_node_edge_count_vec.end(),0);
	thrust::fill(edgeInfoVecs.global_length_zero.begin(), edgeInfoVecs.global_length_zero.end(), 0.0);



	nodeInfoVecs.host_edge_left = hostNodeInfoVecs.host_spring_edge_left;
	nodeInfoVecs.host_edge_right = hostNodeInfoVecs.host_spring_edge_right;
	//scan through hostAdj and put in device.
	for (unsigned id = 0; id < hostNodeInfoVecs.host_spring_length_zero.size(); id++) {
		generalParams.totalNumberOfEdges++;
		unsigned idL = hostNodeInfoVecs.host_spring_edge_left[id];
		unsigned idR = hostNodeInfoVecs.host_spring_edge_right[id];

     	bool is_idL_collagen = hostNodeInfoVecs.host_node_is_collagen[idL];
		bool is_idR_collagen = hostNodeInfoVecs.host_node_is_collagen[idR];
		bool is_edge_collagen = false;
		bool is_edge_elastin = false;
		if (is_idL_collagen && is_idR_collagen) {
			is_edge_collagen=true;
		}
		else { is_edge_elastin = true;}
		//std::cout<< "linking " << idL << " to " <<idR << std::endl;

		 double edgeLen = hostNodeInfoVecs.host_spring_length_zero[id];
				//we use the global_length_zero vector to identify edges as well.

    	//node id is row, column node is connected to row node.
		//add edge for left node
		unsigned edgeNumL = edgeInfoVecs.current_node_edge_count_vec[idL]; //number of edges on (nodeId = row)	is that entry in cECV
		unsigned indexL = idL*generalParams.max_nbr_count + edgeNumL;
		edgeInfoVecs.global_length_zero[indexL] = edgeLen;
		edgeInfoVecs.global_neighbors[indexL] = idR;
		edgeInfoVecs.global_isedge_collagen[indexL] = is_edge_collagen;
		edgeInfoVecs.global_isedge_elastin[indexL] = is_edge_elastin;

		(edgeInfoVecs.current_node_edge_count_vec[idL])++; //right connects to left

		//add edge for right node
		unsigned edgeNumR = edgeInfoVecs.current_node_edge_count_vec[idR]; //number of edges on (nodeId = row)	is that entry in cECV
		unsigned indexR = idR*generalParams.max_nbr_count + edgeNumR;
		edgeInfoVecs.global_length_zero[indexR] = edgeLen;
		edgeInfoVecs.global_neighbors[indexR] = idL;
		edgeInfoVecs.global_isedge_collagen[indexR] = is_edge_collagen;
		edgeInfoVecs.global_isedge_elastin[indexR] = is_edge_elastin;

		(edgeInfoVecs.current_node_edge_count_vec[idR])++; //left connects to right
		generalParams.current_edge_count += 1;
	}
	//at this point current_node_edge_count_vec holds the number of edges, copy this to
	thrust::copy(edgeInfoVecs.current_node_edge_count_vec.begin(), edgeInfoVecs.current_node_edge_count_vec.end(), edgeInfoVecs.num_origin_nbr_per_node_vec.begin());
};

void System::set_extras() {
	if (extensionParams.axis==0){
		extensionParams.originalNetworkLength = domainParams.max_z; //compression along x extensionParams.axis
	}else{
		extensionParams.originalNetworkLength = domainParams.max_x; }
};
