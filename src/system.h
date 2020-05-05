

#ifndef SYSTEM_H_
#define SYSTEM_H_

#pragma once
//this file is included by NSI.h

#include <thrust/system_error.h>
#include <memory>
#include <math.h>
#include <thrust/extrema.h>
#include <thrust/sort.h>
#include <thrust/copy.h>
#include <thrust/for_each.h>
#include <thrust/transform.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/tuple.h>
#include <thrust/device_vector.h>
#include <thrust/device_ptr.h>
#include <thrust/remove.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/copy.h>
#include <thrust/execution_policy.h>
#include <thrust/pair.h>
#include <stdint.h>

struct RandVecs {

  thrust::device_vector<double> gaussianData;

};


//Data Structure for node location. velocity and force
struct NodeInfoVecs {
	thrust::device_vector<bool> node_upper_selection_pull;//holds id's of top 10% for strain test
	thrust::device_vector<bool> node_lower_selection_pull;//holds id's of top 10% for strain test

	thrust::device_vector<unsigned> spring_division_count;
	//holds sum of forces for node on given time_step
	thrust::host_vector<unsigned> id_edges_made_host;
	thrust::host_vector<unsigned> host_id_left;//left id linked
	thrust::host_vector<unsigned> host_id_right;//right id linked
	

	thrust::host_vector<unsigned> host_edge_left;
	thrust::host_vector<unsigned> host_edge_right;

	thrust::device_vector<unsigned> device_edge_left;
	thrust::device_vector<unsigned> device_edge_right;

	thrust::device_vector<unsigned> origin_edge_left;
	thrust::device_vector<unsigned> origin_edge_right;

	thrust::device_vector<unsigned> id_edges_made_temp;
	thrust::device_vector<unsigned> links_made_individual_thread;
	thrust::device_vector<unsigned> id_temp_linked_left;
	thrust::device_vector<unsigned> id_temp_linked_right;

	thrust::device_vector<double> sum_forces_on_node;

	thrust::device_vector<double> discretized_edges_strain; //counts strain of edge
	
	thrust::device_vector<double> discretized_edges_alignment;

//true if fixed, false if movable. Default false.
	thrust::device_vector<bool> is_node_fixed;
  thrust::device_vector<bool> node_is_collagen;
  thrust::device_vector<bool> node_is_elastin;

	// X,Y,Z, location, velocity and force of all nodes
	thrust::device_vector<double> node_loc_x;
	thrust::device_vector<double> node_loc_y;
	thrust::device_vector<double> node_loc_z;
	thrust::device_vector<double> node_vel_x;
	thrust::device_vector<double> node_vel_y;
	thrust::device_vector<double> node_vel_z;

	thrust::device_vector<double> node_vel;

//holds forces to advance position and velocity
	thrust::device_vector<double> node_force_x;
	thrust::device_vector<double> node_force_y;
	thrust::device_vector<double> node_force_z;


};

//struct used for linking of nodes in network
struct AuxVecs {
  thrust::device_vector<unsigned> id_bucket_net_intc; //bucket id
  thrust::device_vector<unsigned> id_value_net_intc; //node id

  // bucket id expanded denotes what  the bucket IDs are for the neighbors of a certain point
  // id value expanded means each point ( represented by its global rank) will have multiple copies including self
  thrust::device_vector<unsigned> id_bucket_expanded_net_intc;
  thrust::device_vector<unsigned> id_value_expanded_net_intc;

  //entry key_begin_net_intc[bucketKey] returns start of node indices to search for interaction
  //entry key_end_net_intc[bucketKey] returns end of node indices to search for interaction
  thrust::device_vector<unsigned> key_begin_net_intc;
  thrust::device_vector<unsigned> key_end_net_intc;

	unsigned end_index_bucket_keys_net_intc;
};

struct DomainParams {
	double min_x;
	double max_x;
	double min_y;
	double max_y;
	double min_z;
	double max_z;
	double origin_min_x;
	double origin_max_x;
	double origin_min_y;
	double origin_max_y;
	double origin_min_z;
	double origin_max_z;
	double grid_spacing_net_intc = 0.5;
	unsigned bucket_count_x;
	unsigned bucket_count_y;
	unsigned bucket_count_z;
	unsigned total_bucket_count_net_intc=0;
};

//Data for edge node id's
struct EdgeInfoVecs {
	double collagen_spring_constant = 10.0; //linear spring for collagen

  	//wlc spring for elastin
	double kB, CLM, temperature;
  	double viscosity_collagen, viscosity_elastin;
	double node_mass = 1;
	double persistence_len_monomer;
  	double num_mon_elastin_area;
  	double collagen_diameter = 0.5;
	double elastin_diameter = 0.1;

	thrust::device_vector<unsigned> num_origin_nbr_per_node_vec;//holds how many original edges a node was connected to.

	//note: flattened matrix formatting
	thrust::device_vector<unsigned> global_neighbors;
	thrust::device_vector<bool> global_isedge_collagen;
	thrust::device_vector<bool> global_isedge_elastin;

	thrust::device_vector<unsigned> current_node_edge_count_vec;
	//thrust::device_vector<unsigned> currentBindCountPerOriginalEdgeVector;

	thrust::device_vector<double> global_length_zero;

};

struct BendInfoVecs {
	unsigned currentSpringCount = 0;
	//thrust::device_vector<unsigned> angleZero_index;//saved as row
	//thrust::device_vector<double> angleZero_value;
	thrust::device_vector<unsigned> leftIndex;
	thrust::device_vector<unsigned> centerIndex;
	thrust::device_vector<unsigned> rightIndex;
	thrust::device_vector<double> angleZero;
	thrust::device_vector<double> forceX;//triple leftIndexLength
	thrust::device_vector<double> forceY;
	thrust::device_vector<double> forceZ;
	thrust::device_vector<double> tempForceX;//triple leftIndexLength
	thrust::device_vector<double> tempForceY;
	thrust::device_vector<double> tempForceZ;

	thrust::device_vector<unsigned> tempTorIndices;
	thrust::device_vector<unsigned> reducedIds;

	unsigned bend_factor = 3;
	unsigned total_bend_count;//total bending springs
	double bend_stiffness_collagen;
	double bend_stiffness_elastin;
};

struct ExtensionParams {
	//reset and used to calculate current equilibrium state.
	double totalAppliedForce=0;
	double applied_force_upper=0;
	double applied_force_lower=0;

	unsigned currentNumberPulled = 0;
	unsigned nextNumberPulled = 0;
	double target_strain = 0.1;
	double current_strain = 1;

	unsigned axis = 0; //default axis for boundary condition is zero.
	double currentNetworkLength;
	double originalNetworkLength = 0;
	double originalNetworkWidth = 0;
	double currentStrainNetworkWidthMininum = 0.0;

	double strain_proportion_end_sim = 4.0;//if applying strain, how far to run L_new/L_original

	double averageLowerStrain = 1.0;
	double averageUpperStrain = 1.0;
	double originAverageLowerStrain = 0.0;
	double originAverageUpperStrain = 1.0;

	double percentOriginalEdgesUnderStrain1 = 0.0;
	double percentOriginalEdgesExtended = 0.0;
	double percentOriginalEdgesCompressed = 0.0;
	double percentAddedEdgesUnderStrain1 = 0.0;
	double percentAddedEdgesExtended = 0.0;
	double percentAddedEdgesCompressed = 0.0;
	double averageStrainAddedEdges = 0.0;
	double averageStrainOriginalEdges = 0.0;

};

struct GeneralParams{
	//general computation
	bool run_sim = true; //default true to begin sim. Sim ends when runSim == false
	bool strain_sim = false;

	double pull_ammount;

	double pull_percent; // set in main.cpp

	unsigned max_nbr_count = 100;
	unsigned max_node_count;//after discretize
	unsigned origin_node_count;//pre discretize

	unsigned origin_edge_count = 0; //total links set at beginning. Constants
	unsigned current_edge_count = 0;//total non constant if links are made

	unsigned sub_node_count = 0;//maximal subnode division for longest edge


	//parameters for advancing timestep and determining equilibrium
	double df, dt, epsilon, maxForce, epsilon_factor;
	double magnitudeForce = 0.0;
	double currentTime = 0.0;

	double currentLength=5;

	//total equilibrium iters and linking determiner
	unsigned iterationCounter = 0;
	bool linking = true;

	unsigned numUpperStrainNodes_elastin;
	unsigned numUpperStrainNodes_collagen;
	unsigned numUpperStrainNodes = 0;
	
	unsigned numLowerStrainNodes_elastin;
	unsigned numLowerStrainNodes_collagen;
	unsigned numLowerStrainNodes = 0;

	unsigned totalNumberOfEdges = 0;//updated after linking
	unsigned max_links_per_iteration = 5;
};

class HostNodeInfoVecs;
class Storage;

class System {
public:
	DomainParams domainParams;
	NodeInfoVecs nodeInfoVecs;
	AuxVecs auxVecs;
	EdgeInfoVecs edgeInfoVecs;
	BendInfoVecs bendInfoVecs;
	ExtensionParams extensionParams;
	GeneralParams generalParams;
	RandVecs randVecs;

	std::shared_ptr<Storage> storage;

public:

	System();

  //use from cpu side to begin solving system.
	void solve_system();

	void solve_forces(); //use in system.cu

	void set_bucket_scheme();

	void assign_storage(std::shared_ptr<Storage> _storage);

	void initialize_system(HostNodeInfoVecs& hostNodeInfoVecs);

	void determine_bounds();//used for strain.

	void set_node_vecs(HostNodeInfoVecs& hostNodeInfoVecs);

	void set_bend_vecs(HostNodeInfoVecs& hostNodeInfoVecs);

	void set_edge_vecs(HostNodeInfoVecs& hostNodeInfoVecs);

	void set_extras();
};


#endif /*NODESYSTEMIMPLDEVICE_H_*/
