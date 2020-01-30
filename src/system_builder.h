#ifndef SYSTEMBUILDER_H_
#define SYSTEMBUILDER_H_

#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include "glm.hpp"
#include <list>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>

struct HostNodeInfoVecs {
	thrust::host_vector<unsigned> host_node_id;
	thrust::host_vector<unsigned> host_spring_edge_left;
	thrust::host_vector<unsigned> host_spring_edge_right;
	thrust::host_vector<double> host_spring_length_zero;
	
	thrust::host_vector<unsigned> host_sub_spring_edge_left;
	thrust::host_vector<unsigned> host_sub_spring_edge_right;
	thrust::host_vector<double> host_sub_spring_length_zero;

	thrust::host_vector<unsigned> host_torsion_index_left;
	thrust::host_vector<unsigned> host_torsion_index_center;
	thrust::host_vector<unsigned> host_torsion_index_right;
	thrust::host_vector<double> host_torsion_angle_zero;

	thrust::host_vector<bool> host_node_is_collagen;
	thrust::host_vector<bool> host_node_is_elastin;
	thrust::host_vector<double> host_pos_x;
	thrust::host_vector<double> host_pos_y;
	thrust::host_vector<double> host_pos_z;

	thrust::host_vector<bool> host_is_node_fixed;
};

class System;

class System_Builder {
public:
	System_Builder(double _epsilon, double _dt, double _df, double _targetStrain);
	~System_Builder();


	//the reason we need to have each buildNode contain previous and next node id's is
	//that we will use this structure to set preferred edge angle to incorporate bending.
	struct BuildNode {

		std::vector<unsigned> next;
		std::vector<unsigned> prev;
		unsigned id;

		BuildNode(unsigned id_) : id(id_) {}	 //the only thing set is the main id
	};



public:
	HostNodeInfoVecs hostNodeInfoVecs;
	
	double default_pull_percent = 0.0;
	double epsilon, dt, df, default_target_strain;

	double default_bend_stiffness_collagen=0.5;
	double default_bend_stiffness_elastin=0.1;
	double default_collagen_spring_constant=10;

	double default_persistence_len_monomer = 1.0;
	double default_temperature = 300.0; // 300' kelvin ~ 27' celsius
	double default_CLM = 1.0;
	double default_kB = 1.3806488e-8;//converted from nN microns
	double default_num_mon_elastin_area=1100;
	double default_viscosity_collagen=1;
	double default_viscosity_elastin=1;

	double defaultUnitsPerExtraNode = 1.0;//??
	double default_collagen_diameter = 0.1; //used for fiber-fiber linking threshold distance.
	double default_elastin_diameter = 0.1; //used for fiber-fiber linking threshold distance.

	bool use_extra_nodes = false;
	bool use_constant_number_extra_nodes = false;

	unsigned default_extra_nodes_per_edge = 0;

	bool default_strain_sim = false;
	bool default_linking = true;
	unsigned axis = 0; //x-axis default. Used for BoundaryForce set in input file. and used in ESOD.h


	std::vector<std::shared_ptr<BuildNode>> buildNodes;

	std::vector<glm::dvec3> nodePositions;

	unsigned add_elastin_node(glm::dvec3);
	unsigned add_collagen_node(glm::dvec3);


	void put_spring(unsigned, unsigned);
	void put_sub_spring(unsigned, unsigned);
	//void put_collagen_spring(unsigned, unsigned);
	//void put_elastin_spring(unsigned, unsigned);
	void put_bending_spring(unsigned, unsigned, unsigned);
	void add_sub_nodes(void);
	std::list<glm::dvec3> fill_space(glm::dvec3, glm::dvec3, unsigned);
	void generate_build_node_triplets(void);
	void fix_node(unsigned);
	//void setSystemForParallelComputation();
	std::shared_ptr<System> create();


};

#endif
