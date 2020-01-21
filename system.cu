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

#include "collagen_elastin_spring.h"
#include "bending_spring.h"
#include "advance_position.h"
#include "bucket_scheme.h"
#include "system.h"
#include "functors_misc.h"

using namespace thrust::placeholders;

system::solve_forces() {

	thrust::fill(nodeInfoVecs.node_force_x.begin(),nodeInfoVecs.node_force_x.end(),0);
	thrust::fill(nodeInfoVecs.node_force_y.begin(),nodeInfoVecs.node_force_y.end(),0);
	thrust::fill(nodeInfoVecs.node_force_z.begin(),nodeInfoVecs.node_force_z.end(),0);

	init_dim_general(nodeInfoVecs, domainParams, auxVecs, generalParams);

	init_net_inct_bucket(nodeInfoVecs, domainParams, auxVecs, generalParams);

	extend_net_inct_bucket(nodeInfoVecs, domainParams, auxVecs, generalParams);

	double addedLinks = generalParams.currentEdgeCount - generalParams.originEdgeCount;

	if (generalParams.linking == true) {
			link_nodes(nodeInfoVecs, edgeInfoVecs, auxVecs, generalParams);
	}

	//apply external force.
	increment_external_force(nodeInfoVecs,
		generalParams,
		extensionParams,
		domainParams);

	//only counts external force on network nodes since force has been reset.
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

  calc_bending_spring_force(nodeInfoVecs, bendInfoVecs, generalParams);
	calc_spring_force(nodeInfoVecs, edgeInfoVecs, generalParams);
};


void system::solve_system() {

	double lastTime = 0.0;
	storage->update_storage();//initial position storage
	bool runIters = true;

	//set initial epsilon
	//generalParams.epsilon = (1.0) *
	//	sqrt(6.0*edgeInfoVecs.kB * edgeInfoVecs.temperature * generalParams.dtTemp / edgeInfoVecs.viscosity_collagen);

	while (runIters == true) {

		generalParams.iterationCounter++;
		generalParams.currentTime += generalParams.dtTemp;
		//std::cout << "current time: " << std::endl;

		advance_positions(
			nodeInfoVecs,
		 	generalParams,
      randVecs);


		if ((generalParams.iterationCounter % 20000) == 0) {
			storage->print_VTK_file();
		}

		solve_forces(); //resets and solves forces for next time step

		if ((generalParams.iterationCounter % 20000) == 0) {
			double currentStrain = (extensionParams.averageUpperStrain - extensionParams.averageLowerStrain) /
			(extensionParams.originAverageUpperStrain - extensionParams.originAverageLowerStrain ) - 1.0;
			if (currentStrain>4.0){
				runIters=false;
			}
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

			GetStrainParameters(nodeInfoVecs,
				edgeInfoVecs,
				generalParams,
				domainParams);
			storage->updateTotalStrain();
		}


		double maxVel = *(thrust::max_element(nodeInfoVecs.node_vel.begin(), nodeInfoVecs.node_vel.end()));
		//std::cout<<"maxvelocity: "<< maxVel<< std::endl;

		//difference in time
 		if (abs(generalParams.currentTime - lastTime) > (generalParams.lagTime)) {
			 //move epsilon. It will be reset

			generalParams.epsilon += 0.01;
			lastTime = generalParams.currentTime;

			std::cout<<"updating epsilon: "<< generalParams.epsilon<<std::endl;

			double addedEdges = generalParams.currentEdgeCount - generalParams.originEdgeCount;
			std::cout<<"added edges: "<< addedEdges <<std::endl;
			std::cout<<"Minz: "<< domainParams.min_z<<std::endl;
			std::cout<<"Maxz: "<< domainParams.max_z<<std::endl;
			std::cout<<"Miny: "<< domainParams.min_y<<std::endl;
			std::cout<<"Maxy: "<< domainParams.max_y<<std::endl;
			std::cout<<"Minx: "<< domainParams.min_x<<std::endl;
			std::cout<<"Maxx: "<< domainParams.max_x<<std::endl;
		}

		if (maxVel < generalParams.epsilon) {
			//store sum of all forces on each node. Used in stress calculations
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

			storage->updateStorage();

			generalParams.totalNumberOfEdges += nodeInfoVecs.id_edges_made_host.size();
			nodeInfoVecs.id_edges_made_host.resize(0);


			generalParams.epsilon = (1.0) *
				sqrt(6.0 * edgeInfoVecs.kB * edgeInfoVecs.temperature * generalParams.dtTemp / edgeInfoVecs.viscosity_collagen);

			std::cout<<"Maximum vel: "<< maxVel <<std::endl;
			std::cout<<"updating epsilon back to original: "<< generalParams.epsilon<<std::endl;
			generalParams.magnitudeForce += generalParams.df;
			std::cout<<"magnitudeForce: "<< generalParams.magnitudeForce<<std::endl;

		}
		///////////////////////////////////////////////////////////////////////////////
		//EQUILIBRIUM END
		//////////////////////////////////////////////////////////////////////

	}

};

system::system()  {};

void system::assign_storage(std::shared_ptr<Storage> _storage) {
	storage = _storage;
}

void system::initialize_system(HostNodeInfoVecs& hostNodeInfoVecs) {

	std::cout<< "total Edge Count: "<< generalParams.originEdgeCount << std::endl;
	std::cout << "max num nodes: " << generalParams.max_node_count << std::endl;

	nodeInfoVecs.origin_edge_left = hostNodeInfoVecs.host_spring_edge_left;
	nodeInfoVecs.origin_edge_right = hostNodeInfoVecs.host_spring_edge_right;

	set_node_vecs(//calls initDimensionBucketScheme
		hostNodeInfoVecs);

	set_bend_vecs(hostNodeInfoVecs);

	set_edge_vecs(hostNodeInfoVecs);

	set_extras();
};


void system::set_node_vecs(
	thrust::host_vector<bool>& hostNodeInfoVecs.hostNodeInfoVecs.host_is_node_fixed,
	thrust::host_vector<double>& hostNodeInfoVecs.host_pos_x,
	thrust::host_vector<double>& hostNodeInfoVecs.host_pos_y,
	thrust::host_vector<double>& hostNodeInfoVecs.host_pos_z,
	thrust::host_vector<unsigned>& hostSpringDivisionCount) {


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


	thrust::fill(nodeInfoVecs.discretized_edges_strain.begin(), nodeInfoVecs.discretized_edges_strain.end(),0.0);
	thrust::fill(nodeInfoVecs.device_edge_right.begin(), nodeInfoVecs.device_edge_right.end(), 0);	//fill force and velocity with zeros for computation.
	thrust::fill(nodeInfoVecs.device_edge_left.begin(), nodeInfoVecs.device_edge_left.end(), 0);	//fill force and velocity with zeros for computation.
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


	//copy fixed positions
	nodeInfoVecs.is_node_fixed.resize(generalParams.max_node_count);
	thrust::fill(nodeInfoVecs.is_node_fixed.begin(), nodeInfoVecs.is_node_fixed.end(), false);

	//now that all the nodes are loaded in, choose the top to apply strain, and fix the bottom

	determine_bounds();

	//at this point all nodes are filled, so we can generate domainParams before seeding dpd particles.
	init_dim_general(
		nodeInfoVecs,
		domainParams,
		auxVecs,
		generalParams,
		extensionParams);

	//set original parameters for domain. others will be reset as simulation takes place.
	domainParams.origin_min_x = domainParams.min_x;
	domainParams.origin_max_x = domainParams.max_x;
	domainParams.origin_min_y = domainParams.min_y;
	domainParams.origin_max_y = domainParams.max_y;
	domainParams.origin_min_z = domainParams.min_z;
	domainParams.origin_max_z = domainParams.max_z;
	std::cout<< "node count : " <<nodeInfoVecs.node_loc_y.size()<< std::endl;


	auxVecs.id_bucket_net_intc.resize(generalParams.max_node_count + dpdParticleVariables.particleCount);
	auxVecs.id_value_net_intc.resize(generalParams.max_node_count + dpdParticleVariables.particleCount);
	auxVecs.id_value_expanded_net_intc.resize(27 * (generalParams.max_node_count + dpdParticleVariables.particleCount));
	auxVecs.id_bucket_expanded_net_intc.resize(27 *( generalParams.max_node_count + dpdParticleVariables.particleCount));

};

void system::determine_bounds() {
	//determin z positions of nodes to be pulled and fixed.

	thrust::device_vector<double> zPosTemp;
	zPosTemp.resize(generalParams.max_node_count);
	thrust::copy(nodeInfoVecs.node_loc_z.begin(), nodeInfoVecs.node_loc_z.end(), zPosTemp.begin());

	//not used
	//pull at least 10% of nodes.
	unsigned tempNodeAmmount = static_cast<unsigned>( 0.25 * generalParams.max_node_count ); //pull 10% of top nodes

	//sort in increasing order
	thrust::sort(zPosTemp.begin(), zPosTemp.end(), thrust::less<double>());
	double length = zPosTemp[ zPosTemp.size()-1 ];
	std::cout<<"start end ZposTemp: "<< zPosTemp[0] << " "<< zPosTemp[zPosTemp.size()-1]<<std::endl;

	//upperLevelAlt pulls 10% default. Set in main.cpp using input
	if (generalParams.pullPercent >= 1.0) {
		std::cout<<"ERROR PULL PERCENT MUST BE LESS THAN ONE"<<std::endl;;
	}
	double upperLevelAlt = (1.0 - generalParams.pullPercent) * length;


	double lowerLevel = abs (upperLevelAlt - (zPosTemp[zPosTemp.size()-1]));

	std::cout<<"minimal level final choice for strain choice: " << lowerLevel <<std::endl;

	std::cout<<"maximal level final choice for strain choice: " << upperLevelAlt <<std::endl;

	//apply strain only to original nodes and not added edge subdivision nodes. Set top and bottom

	thrust::replace_if(nodeInfoVecs.node_upper_selection_pull.begin(), nodeInfoVecs.node_upper_selection_pull.begin() + generalParams.originNodeCount,
						nodeInfoVecs.node_loc_z.begin(),
						IsGreaterThanLevel( upperLevelAlt ), true);

	thrust::replace_if(nodeInfoVecs.node_lower_selection_pull.begin(), nodeInfoVecs.node_lower_selection_pull.begin() + generalParams.originNodeCount,
						nodeInfoVecs.node_loc_z.begin(),
						IsLessThanLevel( lowerLevel ), true);

	generalParams.numUpperStrainNodes = thrust::count_if(nodeInfoVecs.node_upper_selection_pull.begin(),nodeInfoVecs.node_upper_selection_pull.end(), IsEqualToOne( ) );
	generalParams.numLowerStrainNodes = thrust::count_if(nodeInfoVecs.node_lower_selection_pull.begin(),nodeInfoVecs.node_lower_selection_pull.end(), IsEqualToOne( ) );

	std::cout<<"number of nodes pulled for strain: " << generalParams.numLowerStrainNodes + generalParams.numUpperStrainNodes <<std::endl;

	unsigned numFixed = thrust::count_if(nodeInfoVecs.is_node_fixed.begin(),nodeInfoVecs.is_node_fixed.end(), IsEqualToOne() );
	std::cout<<"number of nodes fixed: " << numFixed <<std::endl;
	zPosTemp.resize(0);

}

void system::set_bend_vecs(
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
	for (unsigned i = 0; i<bendInfoVecs.total_bend_count; i++) {
		unsigned n0 = bendInfoVecs.leftIndex[i];
		unsigned n1 = bendInfoVecs.centerIndex[i];
		unsigned n2 = bendInfoVecs.rightIndex[i];
		std::cout<< "angle : "<< n0<< " " << n1<< " " << n2<< " " << bendInfoVecs.angleZero[i]<<std::endl;
	}

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

void system::set_edge_vecs(
	HostNodeInfoVecs& hostNodeInfoVecs ) {

	edgeInfoVecs.global_neighbors.resize(generalParams.max_node_count * generalParams.max_nbr_count);
  edgeInfoVecs.global_isedge_collagen.resize(generalParams.max_node_count * generalParams.max_nbr_count);
  edgeInfoVecs.global_isedge_elastin.resize(generalParams.max_node_count * generalParams.max_nbr_count);

	edgeInfoVecs.current_node_edge_count_vec.resize(generalParams.max_node_count);

	edgeInfoVecs.global_length_zero.resize(generalParams.max_node_count * generalParams.max_nbr_count);
	edgeInfoVecs.num_origin_nbr_per_node_vec.resize(generalParams.max_node_count);


  thrust::fill(edgeInfoVecs.global_neighbors.begin(), edgeInfoVecs.global_neighbors.end(), ULONG_MAX);
  thrust::fill(edgeInfoVecs.global_isedge_collagen.begin(), edgeInfoVecs.global_isedge_collagen.end(), false);
  thrust::fill(edgeInfoVecs.global_isedge_elastin.begin(), edgeInfoVecs.global_isedge_elastin.end(), false);

  thrust::fill(edgeInfoVecs.current_node_edge_count_vec.begin(), edgeInfoVecs.current_node_edge_count_vec.end(),0);
	thrust::fill(edgeInfoVecs.global_length_zero.begin(), edgeInfoVecs.global_length_zero.end(), 0.0);



	nodeInfoVecs.device_edge_left = host_sub_edge_left;
	nodeInfoVecs.device_edge_right = host_sub_edge_right;
	//scan through hostAdj and put in device.
	for (unsigned id = 0; id < host_sub_len_zero.size(); id++) {
		generalParams.totalNumberOfEdges++;
		 unsigned idL = host_sub_edge_left[id];
		 unsigned idR = host_sub_edge_right[id];
     bool is_edge_collagen = host_sub_edge_is_collagen[id];
     bool is_edge_elastin = host_sub_edge_is_elastin[id];
		//std::cout<< "linking " << idL << " to " <<idR << std::endl;

		 double edgeLen = host_sub_len_zero[id];
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
		generalParams.currentEdgeCount += 1;
	}
	//at this point current_node_edge_count_vec holds the number of edges, copy this to
	thrust::copy(edgeInfoVecs.current_node_edge_count_vec.begin(), edgeInfoVecs.current_node_edge_count_vec.end(), edgeInfoVecs.num_origin_nbr_per_node_vec.begin());
};

void system::set_extras() {
	extensionParams.originalNetworkLength = domainParams.max_z; //compression along x extensionParams.axis
	extensionParams.originalNetworkWidth = domainParams.max_x;  //strain along z extensionParams.axis.
};