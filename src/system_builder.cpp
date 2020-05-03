/*
 * System_Builder.cpp
 *
 *  Created on: 25 авг. 2014 г.
 *      Author: yan & SRB
 */

#include <set>
#include <list>
#include <vector>
#include <memory>
#include "system_builder.h"
#include "system.h"
# define M_PI 3.14159265358979323846  /* pi */



System_Builder::System_Builder(
	double _epsilon,
	double _dt,
	double _df,
	double _targetStrain):
	epsilon(_epsilon),
	dt(_dt),
	df(_df),
	default_target_strain(_targetStrain) {}

System_Builder::~System_Builder() {
}

//currently unused
/*
void System_Builder::put_elastin_spring(unsigned n1, unsigned n2) {
	double global_length_zero = glm::length(nodePositions[n1] - nodePositions[n2]);

	//store for parallel solving.
	hostNodeInfoVecs.host_spring_length_zero.push_back(global_length_zero);
	hostNodeInfoVecs.host_spring_edge_left.push_back(n1);
	hostNodeInfoVecs.host_spring_edge_right.push_back(n2);
	//hostNodeInfoVecs.host_sub_edge_is_collagen.push_back(false)
	//hostNodeInfoVecs.host_sub_edge_is_elastin.push_back(true)
}

void System_Builder::put_collagen_spring(unsigned n1, unsigned n2) {


	double global_length_zero = glm::length(nodePositions[n1] - nodePositions[n2]);

	//store for parallel solving.
	hostNodeInfoVecs.host_spring_length_zero.push_back(global_length_zero);
	hostNodeInfoVecs.host_spring_edge_left.push_back(n1);
	hostNodeInfoVecs.host_spring_edge_right.push_back(n2);
	//hostNodeInfoVecs.host_sub_edge_is_collagen.push_back(true)
	//hostNodeInfoVecs.host_sub_edge_is_elastin.push_back(false)
};*/

void System_Builder::put_spring(unsigned n1, unsigned n2) {

	//std::cout<<"adding spring between: " << n1 << " " << n2 << std::endl;
	double global_length_zero = glm::length(nodePositions[n1] - nodePositions[n2]);

	hostNodeInfoVecs.host_spring_edge_left.push_back(n1);
	hostNodeInfoVecs.host_spring_edge_right.push_back(n2);
	
	hostNodeInfoVecs.host_spring_length_zero.push_back(global_length_zero);
};


void System_Builder::put_sub_spring(unsigned n1, unsigned n2) {

	
	//std::cout<<"adding sub spring between: " << n1 << " " << n2 << std::endl;
	double global_length_zero = glm::length(nodePositions[n1] - nodePositions[n2]);

	hostNodeInfoVecs.host_sub_spring_edge_left.push_back(n1);
	hostNodeInfoVecs.host_sub_spring_edge_right.push_back(n2);
	hostNodeInfoVecs.host_sub_spring_length_zero.push_back(global_length_zero);
};


void System_Builder::put_bending_spring(unsigned n1, unsigned n2, unsigned n3) {

	glm::dvec3 p1 = nodePositions[n1] - nodePositions[n2];
	glm::dvec3 p2 = nodePositions[n3] - nodePositions[n2];

	double cos_theta = glm::dot(p1, p2) / (glm::length(p1) * glm::length(p2));
	double preferredAngle;
	if (cos_theta > 1) {
		cos_theta = 1.0;
	}
	else if (cos_theta < -1) {
		cos_theta = -1.0;
	}

	preferredAngle = std::acos(cos_theta);


	hostNodeInfoVecs.host_torsion_index_left.push_back(n1);
	hostNodeInfoVecs.host_torsion_index_center.push_back(n2);
	hostNodeInfoVecs.host_torsion_index_right.push_back(n3);
	hostNodeInfoVecs.host_torsion_angle_zero.push_back(preferredAngle);

}


//all node pointers made here. We take this opportunity to
unsigned System_Builder::add_collagen_node(glm::dvec3 pos) {

	unsigned newId = buildNodes.size();
	//std::cout<<"adding new collagen node: " << newId<< std::endl;
	//notice the node and buildnode have the same corresponding id's.
	std::shared_ptr<BuildNode> ptr1(new BuildNode(newId));
	buildNodes.push_back(ptr1);

	nodePositions.push_back(pos);

	hostNodeInfoVecs.host_node_id.push_back(newId);

	hostNodeInfoVecs.host_pos_x.push_back(pos.x);
	hostNodeInfoVecs.host_pos_y.push_back(pos.y);
	hostNodeInfoVecs.host_pos_z.push_back(pos.z);

	hostNodeInfoVecs.host_node_is_collagen.push_back(true);
	hostNodeInfoVecs.host_node_is_elastin.push_back(false);

	hostNodeInfoVecs.host_is_node_fixed.push_back(false);
	return newId;

} 

unsigned System_Builder::add_elastin_node(glm::dvec3 pos) {

	unsigned newId = buildNodes.size();
	
	//std::cout<<"adding new elastin node: " << newId<< std::endl;
	//notice the node and buildnode have the same corresponding id's.
	std::shared_ptr<BuildNode> ptr1(new BuildNode(newId));
	buildNodes.push_back(ptr1);

	nodePositions.push_back(pos);

	hostNodeInfoVecs.host_node_id.push_back(newId);

	hostNodeInfoVecs.host_pos_x.push_back(pos.x);
	hostNodeInfoVecs.host_pos_y.push_back(pos.y);
	hostNodeInfoVecs.host_pos_z.push_back(pos.z);

	hostNodeInfoVecs.host_node_is_collagen.push_back(false);
	hostNodeInfoVecs.host_node_is_elastin.push_back(true);

	hostNodeInfoVecs.host_is_node_fixed.push_back(false);
	return newId;

}



//list holds positions of subnodes
std::list<glm::dvec3> System_Builder::fill_space(glm::dvec3 from, glm::dvec3 to, unsigned subNodes) {

	std::list<glm::dvec3> list;
	glm::dvec3 segment = (to - from) / (subNodes + 1.0);
	for (unsigned i = 1; i <= subNodes; ++i)
		list.push_back(from + segment * (double) i);
	return list;
}



//this method is made to fix an error in how the buildNodes are generated.
//Each node's neighbors are scanned, and depending on how many there are, a different number of springs
//is added.
//notice that we sort the neighbors so torsion springs will be placed in a specific ordering. This is used when preferred angles are stored in a matrix.
void System_Builder::generate_build_node_triplets() {
	
	//std::cout<<"node positions" << std::endl;
	for (unsigned i = 0; i < nodePositions.size(); i++) {
		//std::cout<<nodePositions[i]<<std::endl;
		auto ptrBN = buildNodes[i];
		unsigned center = ptrBN->id;

		//find neighbors of center
		std::vector<unsigned> neighbors;
		for (unsigned i = 0; i < hostNodeInfoVecs.host_spring_edge_left.size(); ++i) {
			unsigned idLeft = hostNodeInfoVecs.host_spring_edge_left[i];
			unsigned idRight = hostNodeInfoVecs.host_spring_edge_right[i];
			if (idLeft == center) {
				neighbors.push_back(idRight);
			}
			if (idRight == center) {
				neighbors.push_back(idLeft);
			}
		}
		//std::cout<<"nbrs size"<<neighbors.size() <<std::endl;
		/*for (unsigned temp = 0; temp < neighbors.size(); temp++){
			std::cout<<neighbors[temp] << std::endl;
		}*/
		//now that we have the neighbors, we'll add pairs to prev and next
		//we'll sort them before adding.

		std::sort(neighbors.begin(), neighbors.end());

		//for the buildNode related to center
		//we only continue if there are more than one neighbor.
		//if two neightbors, we need one torsion spring
		if (neighbors.size() == 2) {
			ptrBN->prev.push_back(neighbors[0]);
			ptrBN->next.push_back(neighbors[1]);
			continue;
		}

		//if n>2 neighbors we need n torsion springs
		if (neighbors.size() > 2) {
			for (unsigned j = 0; j < neighbors.size(); ++j) {

				ptrBN->prev.push_back(neighbors[j]);
				unsigned index = (j + 1) % neighbors.size();
				ptrBN->next.push_back(neighbors[index]);
			}
		}
		
		//std::cout<<buildNodes[i]->id<<std::endl;
		//std::cout<<buildNodes[i]->id<<std::endl;
		//std::cout<<buildNodes[i]->id<<std::endl;

	}
}

void System_Builder::fix_node(unsigned id) {
	hostNodeInfoVecs.host_is_node_fixed[id] = true;
}

void System_Builder::add_sub_nodes() {
	std::cout << "setting subnodes. Total edges: "<< hostNodeInfoVecs.host_spring_edge_left.size() << std::endl;
	//if use Extra nodes is true, we complete this section
	if (use_extra_nodes) {
		//std::cout << "use extra nodes: " << useExtraNodes << std::endl;
		for (unsigned i = 0; i < hostNodeInfoVecs.host_spring_edge_left.size(); ++i) {

			unsigned idLeft = hostNodeInfoVecs.host_spring_edge_left[i];
			unsigned idRight = hostNodeInfoVecs.host_spring_edge_right[i];
			bool id_left_is_collagen = hostNodeInfoVecs.host_node_is_collagen[idLeft];
			bool id_left_is_elastin = hostNodeInfoVecs.host_node_is_elastin[idLeft];
			bool id_right_is_collagen = hostNodeInfoVecs.host_node_is_collagen[idRight];
			bool id_right_is_elastin = hostNodeInfoVecs.host_node_is_elastin[idRight];
			bool sub_nodes_are_collagen=false;
			bool sub_nodes_are_elastin=false;
			if (id_left_is_collagen && id_right_is_collagen){
				sub_nodes_are_collagen=true;
			}
			else {
				//in this case, one or both of the sub nodes on the edge is elastin.
				sub_nodes_are_elastin = true;
			}
			unsigned sub_node_count = 0; //default subnode amount

			if (use_constant_number_extra_nodes) {
				sub_node_count = default_extra_nodes_per_edge;
			}
			else {
				double length = glm::length(nodePositions[idLeft] - nodePositions[idRight]);
				sub_node_count = (unsigned)glm::round(length / defaultUnitsPerExtraNode);

				//if we are not using a constant number of nodes, make the variable
				//hold the largest amount of subnodes per edge.
				/*if (sub_node_count > default_extra_nodes_per_edge) {
					default_extra_nodes_per_edge = sub_node_count;
				}*/
			}

			if (sub_node_count == 0)
				std::cout << "subnodes cannot be placed" << "\n";

			// fill space between two nodes by adding a number of extra nodes
			// you can choose these or they will be rounded.
			auto points = fill_space(nodePositions[idLeft], nodePositions[idRight], sub_node_count);
			std::vector<unsigned> subNodeIds;

			//notice that each subnode is added to the vector of no des
			for (glm::dvec3& point : points) {

				//this addnode method assures that if subnodes are used, then the id and pos are still in hostposX,Y,Z.
				//this must be done before using put_spring method
				//number of edges made by adding n nodes is n+1.
				if (sub_nodes_are_collagen){
					subNodeIds.push_back(add_collagen_node(point));
				}
				else if (sub_nodes_are_elastin) {
					subNodeIds.push_back(add_elastin_node(point));
				}
			}

			if (sub_node_count == 0)
			{
				std::cout << "subnodesCount = 0 after fill_space" << "\n";
				continue;
			}

			//Now that extra nodes are placed, place extra springs along each divided edge.
			//link head
			put_sub_spring(idLeft, subNodeIds.front());
			// link all nodes to make a chain
			for (unsigned i = 0; i < subNodeIds.size() - 1; ++i) {
				put_sub_spring(subNodeIds[i], subNodeIds[i + 1]);
			}
			// link main nodes with and tail of list of sub-nodes

			put_sub_spring(subNodeIds.back(), idRight);


		}

	}

}

//idea: need to use subnode id's for linking during simulation.


//adds all constraints to the nodesystem model so that it can use the constraints.
std::shared_ptr<System> System_Builder::create() {
	//before adding subnodes generate original numbers of nodes and links
	unsigned origin_edge_count = hostNodeInfoVecs.host_spring_edge_left.size();
	unsigned origin_node_count = hostNodeInfoVecs.host_pos_x.size();

	//first add all subnodes
	add_sub_nodes();
	/*for (unsigned edge = 0; edge < hostNodeInfoVecs.host_spring_edge_left.size(); edge++ ){
		std::cout<< "preleft " <<hostNodeInfoVecs.host_spring_edge_left[edge] << std::endl;
		std::cout<< "preright " <<hostNodeInfoVecs.host_spring_edge_right[edge] << std::endl;
	}
	for (unsigned edge = 0; edge < hostNodeInfoVecs.host_sub_spring_edge_left.size(); edge++ ){
		std::cout<< "sub left " <<hostNodeInfoVecs.host_sub_spring_edge_left[edge] << std::endl;
		std::cout<< "sub right " <<hostNodeInfoVecs.host_sub_spring_edge_right[edge] << std::endl;
	}*/
	if (hostNodeInfoVecs.host_sub_spring_edge_left.size() > hostNodeInfoVecs.host_spring_edge_left.size()){
		hostNodeInfoVecs.host_spring_edge_left = hostNodeInfoVecs.host_sub_spring_edge_left;
		hostNodeInfoVecs.host_spring_edge_right = hostNodeInfoVecs.host_sub_spring_edge_right;
		hostNodeInfoVecs.host_spring_length_zero = hostNodeInfoVecs.host_sub_spring_length_zero;
	}
	/*for (unsigned edge = 0; edge < hostNodeInfoVecs.host_spring_edge_left.size(); edge++ ){
		std::cout<< "left " <<hostNodeInfoVecs.host_spring_edge_left[edge] << std::endl;
		std::cout<< "right " <<hostNodeInfoVecs.host_spring_edge_right[edge] << std::endl;
	}*/

	std::cout << "total edge count: " << hostNodeInfoVecs.host_spring_length_zero.size() << std::endl;

	unsigned numNodes = hostNodeInfoVecs.host_pos_x.size();
	unsigned numEdges = hostNodeInfoVecs.host_spring_edge_left.size();
	std::cout << "node count before adding: " << origin_node_count << std::endl;
	std::cout << "node count after adding: " << numNodes << std::endl;


	//preferred edge angle is always taken into account.
	std::cout << "Torsion springs mounting" << std::endl;
	generate_build_node_triplets(); //sets all nodes in buildNodes neighbors for linking in next 10 lines
	for (auto& center : buildNodes) {
		//if (center->id >= originNodeCount) {
			//only use bending springs along fibers

			for (unsigned i = 0; i < (center->next).size(); ++i) {
				unsigned left = center->prev[i];
				unsigned right = center->next[i];
				
				//std::cout << "bending spring: " << left << " " << center->id<< " " << right << std::endl;
				put_bending_spring(left, center->id, right);//set theta value here
			}
	}


	//now all the edges and variables are set.
	//so set the system and return a pointer to main.
	std::shared_ptr<System> system_ptr = std::make_shared<System>();

	system_ptr->generalParams.max_node_count = hostNodeInfoVecs.host_pos_x.size();
	system_ptr->bendInfoVecs.total_bend_count = hostNodeInfoVecs.host_torsion_angle_zero.size();
	system_ptr->bendInfoVecs.bend_stiffness_collagen = default_bend_stiffness_collagen;
	system_ptr->bendInfoVecs.bend_stiffness_elastin = default_bend_stiffness_elastin;

	system_ptr->generalParams.origin_node_count = origin_node_count;
	system_ptr->generalParams.origin_edge_count = origin_edge_count;

	system_ptr->generalParams.sub_node_count = default_extra_nodes_per_edge;
	system_ptr->generalParams.epsilon_factor = epsilon;
	system_ptr->generalParams.dt = dt;
	system_ptr->generalParams.df = df;
	system_ptr->generalParams.pull_ammount = pull_ammount;
	

	system_ptr->edgeInfoVecs.collagen_spring_constant = default_collagen_spring_constant;
	system_ptr->edgeInfoVecs.kB = default_kB;
	system_ptr->edgeInfoVecs.CLM = default_CLM;
	system_ptr->edgeInfoVecs.viscosity_collagen = default_viscosity_collagen;
	system_ptr->edgeInfoVecs.viscosity_elastin = default_viscosity_elastin;
	system_ptr->edgeInfoVecs.collagen_diameter = default_collagen_diameter;
	system_ptr->edgeInfoVecs.elastin_diameter = default_elastin_diameter;
	system_ptr->edgeInfoVecs.temperature = default_temperature;
	system_ptr->edgeInfoVecs.persistence_len_monomer = default_persistence_len_monomer;
	system_ptr->edgeInfoVecs.num_mon_elastin_area=default_num_mon_elastin_area;

	system_ptr->generalParams.linking = default_linking;
	system_ptr->generalParams.strain_sim = default_strain_sim;

	system_ptr->extensionParams.target_strain = default_target_strain;
	system_ptr->extensionParams.axis = axis;

	system_ptr->generalParams.pull_percent = default_pull_percent;
	std::cout<<"pull percent in builder.cpp: " << system_ptr->generalParams.pull_percent<< std::endl;

	system_ptr->initialize_system(hostNodeInfoVecs);

	return system_ptr;

}
