#include <sys/stat.h>

#include <iomanip> // setprecision
#include <sstream> // stringstream

#include "system.h"
#include "system_builder.h"

#include "system_structures.h"
#include "storage.h"
#include <numeric>

Storage::Storage(std::weak_ptr<System> a_system,
	std::weak_ptr<System_Builder> b_system , const std::string& a_fileName) {
	
	system = a_system;
	builder = b_system;
	bn = a_fileName; //this will be used later to open files
	//std::ofstream statesOutput(a_fileName + ".sta");
	//std::ofstream statesOutputStrain(a_fileName + "_Strain.sta");

	std::shared_ptr<System> sys = system.lock();

	if ( sys ){
		std::stringstream stream_min;
		std::stringstream stream_max;

		unsigned domain_size = ceil((sys->domainParams.max_x + 
			sys->domainParams.max_y + 
			sys->domainParams.max_z) / 3.0);

		unsigned max_nodes = sys->generalParams.max_node_count;
		unsigned max_z = sys->domainParams.max_z;
		unsigned max_x = sys->domainParams.max_x;
		unsigned axis = sys->extensionParams.axis;
		int pull_ammount = int(100*sys->generalParams.pull_ammount);
		//std::stringstream tmp;
		//tmp << std::setprecision(3) << std::fixed << pull_ammount;
		//pull_ammount = stod(tmp.str());
		
		int epsilon = int(1000*sys->generalParams.epsilon_factor);
		std::stringstream tmp1;
		tmp1 << std::setprecision(4) << std::fixed << epsilon;
		epsilon = stod(tmp1.str());

		std::string str_nodes = "_max_nodes_";
		std::string str_z = "_max_z_";
		std::string str_x = "_max_x_";
		std::string str_axis = "_axis_";
		std::string str_pull = "_pullwidth_";
		std::string str_eps = "_eps_";

		std::string str_a = "Animation_";
		std::string str_p = "Params_";
		
		str_animation = str_a +
			str_nodes + std::to_string(max_nodes)+
			str_z + std::to_string(max_z)+
			str_x + std::to_string(max_x)+
			str_axis + std::to_string(axis)+
			str_pull + std::to_string(pull_ammount)+
			str_eps + std::to_string(epsilon);

		const int dir_err_anim = mkdir(str_animation.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		if (-1 == dir_err_anim)
		{printf("Error creating directory animation test!n");}
		else {printf("making folder!n"); printf(str_animation.c_str());}

		str_params = str_p + 			
			str_nodes + std::to_string(max_nodes)+
			str_z + std::to_string(max_z)+
			str_x + std::to_string(max_x)+
			str_axis + std::to_string(axis)+
			str_pull + std::to_string(pull_ammount)+
			str_eps + std::to_string(epsilon);

		const int dir_err_params = mkdir(str_params.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		if (-1 == dir_err_params)
		{printf("Error creating directory params!n");}
		else {printf("making folder!n"); printf(str_params.c_str());}
	}
};

void Storage::save_params(void) {
	std::shared_ptr<System> sys = system.lock();
	if (sys) {

		double currentStrain = (sys->extensionParams.averageUpperStrain - sys->extensionParams.averageLowerStrain) /
			 (sys->extensionParams.originAverageUpperStrain - sys->extensionParams.originAverageLowerStrain ) - 1.0;
		//first create a new file using the current network strain

		std::string format = ".sta";
		
		std::string strain =  std::to_string(sys->generalParams.currentTime);
		std::string initial = str_params+"/Param_";
		std::ofstream ofs;
		std::string Filename = initial + strain + format;
		ofs.open(Filename.c_str());


		unsigned max_nbr_count = sys->generalParams.max_nbr_count;
		unsigned max_node_count = sys->generalParams.max_node_count;
		unsigned originalNodeCount = sys->generalParams.origin_node_count;
		unsigned originalEdgeCount = sys->generalParams.origin_edge_count;
		unsigned edgeCountDiscretize = sys->generalParams.current_edge_count;
		//Now first place strain
		ofs << std::setprecision(5) <<std::fixed<< "time " << sys->generalParams.currentTime<<std::endl;
		ofs << std::setprecision(5) <<std::fixed<< "network_strain " << currentStrain<<std::endl;
		ofs << std::setprecision(5) <<std::fixed<< "min_x " << sys->domainParams.min_x<<std::endl;
		ofs << std::setprecision(5) <<std::fixed<< "max_x " << sys->domainParams.max_x<<std::endl;
		ofs << std::setprecision(5) <<std::fixed<< "min_y " << sys->domainParams.min_y<<std::endl;
		ofs << std::setprecision(5) <<std::fixed<< "max_y " << sys->domainParams.max_y<<std::endl;
		ofs << std::setprecision(5) <<std::fixed<< "min_z " << sys->domainParams.min_x<<std::endl;
		ofs << std::setprecision(5) <<std::fixed<< "max_z " << sys->domainParams.max_x<<std::endl;

		ofs << std::setprecision(5) <<std::fixed<< "force_upper " << sys->extensionParams.applied_force_upper<<std::endl;
		ofs << std::setprecision(5) <<std::fixed<< "force_lower " << sys->extensionParams.applied_force_lower<<std::endl;


		//ofs << std::setprecision(5) <<std::fixed<< "total_applied_force " << sys->extensionParams.totalAppliedForce<<std::endl;
		ofs << std::setprecision(5) <<std::fixed<< "original_node_count " << originalNodeCount <<std::endl;
		ofs << std::setprecision(5) <<std::fixed<< "node_count_discretize " << max_node_count <<std::endl;
		ofs << std::setprecision(5) <<std::fixed<< "original_edge_count " << originalEdgeCount <<std::endl;
		ofs << std::setprecision(5) <<std::fixed<< "edge_count_discretize " << edgeCountDiscretize <<std::endl;

		//place nodes
		for (unsigned i = 0; i < sys->nodeInfoVecs.node_loc_x.size(); i++) {
			double x = sys->nodeInfoVecs.node_loc_x[i];
			double y = sys->nodeInfoVecs.node_loc_y[i];
			double z = sys->nodeInfoVecs.node_loc_z[i];
			ofs << std::setprecision(5) <<std::fixed<< "node " << x << " " << y << " " << z <<std::endl;

		}
		//place force node is experiencing
		for (unsigned i = 0; i < sys->nodeInfoVecs.node_loc_x.size(); i++) {
			ofs << std::setprecision(5) <<std::fixed<< "force_on_node " << sys->nodeInfoVecs.sum_forces_on_node[i]<<std::endl;

		}

		//place original edges
		for (unsigned edge = 0; edge < sys->generalParams.origin_edge_count; edge++) {
			unsigned idL = sys->nodeInfoVecs.host_edge_left[edge];
			unsigned idR = sys->nodeInfoVecs.host_edge_right[edge];
			ofs <<"original_edge_discretized " <<idL <<" "<< idR <<std::endl;

		}

		//place added edges
		for (unsigned edge = sys->generalParams.origin_edge_count; edge < sys->generalParams.current_edge_count; edge++) {
			unsigned idL = sys->nodeInfoVecs.host_edge_left[edge];
			unsigned idR = sys->nodeInfoVecs.host_edge_right[edge];
			ofs <<"added_edge " <<idL <<" "<< idR <<std::endl;

		}

		//original edge strain
		for (unsigned i = 0; i < sys->generalParams.origin_edge_count; i++ ){
			double val = sys->nodeInfoVecs.discretized_edges_strain[i];

			ofs << std::setprecision(5)<< std::fixed<<"original_edge_strain " << val <<std::endl;
		}

		//original edge alignment
		for (unsigned i = 0; i < sys->generalParams.origin_edge_count; i++ ){
			double val = sys->nodeInfoVecs.discretized_edges_alignment[i];
			ofs << std::setprecision(5)<< std::fixed<<"original_edge_alignment " << val <<std::endl;
		}

		//added edge strain
		for (unsigned i = sys->generalParams.origin_edge_count; i < sys->generalParams.current_edge_count; i++ ){
			double val = sys->nodeInfoVecs.discretized_edges_strain[i];
			ofs << std::setprecision(5)<< std::fixed<<"added_edge_strain " << val <<std::endl;
		}

		//added links per node.
		for (unsigned i = 0; i < sys->generalParams.max_node_count; i++ ){
			unsigned val = sys->edgeInfoVecs.current_node_edge_count_vec[i] -
				sys->edgeInfoVecs.num_origin_nbr_per_node_vec[i];
			ofs << std::setprecision(5)<< std::fixed<<"bind_sites_per_node " << val <<std::endl;
		}



	}
}


void Storage::print_VTK_file() {

	std::shared_ptr<System> sys = system.lock();
	if (sys) {

		unsigned max_node_count = sys->generalParams.max_node_count;
		unsigned max_nbr_count = sys->generalParams.max_nbr_count;
		unsigned num_collagen_edges = 0;
		unsigned num_elastin_edges=0;
		unsigned numEdges = sys->nodeInfoVecs.host_edge_left.size();
		for (unsigned edge = 0; edge < numEdges; edge++) {	
			unsigned idA = sys->nodeInfoVecs.host_edge_left[edge];
			unsigned idB = sys->nodeInfoVecs.host_edge_right[edge];
			bool is_A_collagen = sys->nodeInfoVecs.node_is_collagen[idA];
			bool is_B_collagen = sys->nodeInfoVecs.node_is_collagen[idB];
			if (is_A_collagen && is_B_collagen){ num_collagen_edges+=1;}
			else{num_elastin_edges+=1;}
		}
		
		iteration+=1;
		unsigned digits = ceil(log10(iteration + 1));
		std::string format = ".vtk";
		std::string Number;
		std::string initial = str_animation + "/Collagen_Network_";
		std::ofstream ofs;
		if (digits == 1 || digits == 0) {
			Number = "0000" + std::to_string(iteration);
		}
		else if (digits == 2) {
			Number = "000" + std::to_string(iteration);
		}
		else if (digits == 3) {
			Number = "00" + std::to_string(iteration);
		}
		else if (digits == 4) {
			Number = "0" + std::to_string(iteration);
		}

		std::string Filename = initial + Number + format;

		ofs.open(Filename.c_str());

		ofs << "# vtk DataFile Version 3.0" << std::endl;
		ofs << "Point representing Sub_cellular elem model" << std::endl;
		ofs << "ASCII" << std::endl << std::endl;
		ofs << "DATASET UNSTRUCTURED_GRID" << std::endl;


		ofs << "POINTS " << max_node_count << " float" << std::endl;
		for (unsigned i = 0; i< max_node_count; i++) {
			double xPos = sys->nodeInfoVecs.node_loc_x[i];
			double yPos = sys->nodeInfoVecs.node_loc_y[i];
			double zPos = sys->nodeInfoVecs.node_loc_z[i];

			ofs << std::setprecision(5) <<std::fixed<< xPos << " " << yPos << " " << zPos << " " << '\n'<< std::fixed;
		}
		//now plot particles
		unsigned numCells = num_collagen_edges;
		unsigned numNumsInCells = 3 * num_collagen_edges;

		ofs << "CELLS " << numCells << " " << numNumsInCells << std::endl;

		for (unsigned edge = 0; edge < numEdges; edge++) {
			
			unsigned idA = sys->nodeInfoVecs.host_edge_left[edge];
			unsigned idB = sys->nodeInfoVecs.host_edge_right[edge];
			bool is_A_collagen = sys->nodeInfoVecs.node_is_collagen[idA];
			bool is_B_collagen = sys->nodeInfoVecs.node_is_collagen[idB];
			if (is_A_collagen && is_B_collagen){
				ofs<< 2 << " " << idA << " " << idB << std::endl;
			}
		}

		ofs << "CELL_TYPES " << numCells << std::endl;
		for (unsigned i = 0; i<num_collagen_edges; i++) {
			ofs << 3 << std::endl;
		}

		
		ofs << "CELL_DATA " << numCells << std::endl;
		ofs << "SCALARS Fiber_Strain double " << std::endl;
		ofs << "LOOKUP_TABLE default "  << std::endl;
		for (unsigned edge = 0; edge < numEdges; edge++) {
			unsigned idA = sys->nodeInfoVecs.host_edge_left[edge];
			unsigned idB = sys->nodeInfoVecs.host_edge_right[edge];
			bool is_A_collagen = sys->nodeInfoVecs.node_is_collagen[idA];
			bool is_B_collagen = sys->nodeInfoVecs.node_is_collagen[idB];
			if (is_A_collagen && is_B_collagen){
				unsigned begin = idA * max_nbr_count;
				unsigned end = begin + max_nbr_count;
				double L0;
				for (unsigned i = begin; i < end; i++) {
					unsigned idTemp = sys->edgeInfoVecs.global_neighbors[i];
					if (idTemp == idB){
						L0 = sys->edgeInfoVecs.global_length_zero[i];
					}
				}
				double xL = sys->nodeInfoVecs.node_loc_x[idA];
				double yL = sys->nodeInfoVecs.node_loc_y[idA];
				double zL = sys->nodeInfoVecs.node_loc_z[idA];
				double xR = sys->nodeInfoVecs.node_loc_x[idB];
				double yR = sys->nodeInfoVecs.node_loc_y[idB];
				double zR = sys->nodeInfoVecs.node_loc_z[idB];

				double L1 = std::sqrt( (xL - xR)*(xL - xR)+(yL - yR)*(yL - yR)+(zL - zR)*(zL - zR));
				double strain = (L1 - L0) / L0;
				ofs << std::fixed << strain   << std::endl;
			}

		}

		ofs.close();


		//Now print elastin
		
		initial = str_animation + "/Elastin_Network_";
		Filename = initial + Number + format;

		ofs.open(Filename.c_str());

		ofs << "# vtk DataFile Version 3.0" << std::endl;
		ofs << "Point representing Sub_cellular elem model" << std::endl;
		ofs << "ASCII" << std::endl << std::endl;
		ofs << "DATASET UNSTRUCTURED_GRID" << std::endl;


		ofs << "POINTS " << max_node_count << " float" << std::endl;
		for (unsigned i = 0; i< max_node_count; i++) {
			double xPos = sys->nodeInfoVecs.node_loc_x[i];
			double yPos = sys->nodeInfoVecs.node_loc_y[i];
			double zPos = sys->nodeInfoVecs.node_loc_z[i];

			ofs << std::setprecision(5) <<std::fixed<< xPos << " " << yPos << " " << zPos << " " << '\n'<< std::fixed;
		}
		//now plot particles
		numCells = num_elastin_edges;
		numNumsInCells = 3 * num_elastin_edges;

		ofs << "CELLS " << numCells << " " << numNumsInCells << std::endl;

		for (unsigned edge = 0; edge < numEdges; edge++) {
			
			unsigned idA = sys->nodeInfoVecs.host_edge_left[edge];
			unsigned idB = sys->nodeInfoVecs.host_edge_right[edge];
			bool is_A_collagen = sys->nodeInfoVecs.node_is_collagen[idA];
			bool is_B_collagen = sys->nodeInfoVecs.node_is_collagen[idB];
			if ((!is_A_collagen) || (!is_B_collagen)){
				ofs<< 2 << " " << idA << " " << idB << std::endl;
			}
		}

		ofs << "CELL_TYPES " << numCells << std::endl;
		for (unsigned i = 0; i<num_elastin_edges; i++) {
			ofs << 3 << std::endl;
		}

		
		ofs << "CELL_DATA " << numCells << std::endl;
		ofs << "SCALARS Fiber_Strain double " << std::endl;
		ofs << "LOOKUP_TABLE default "  << std::endl;
		for (unsigned edge = 0; edge < numEdges; edge++) {
			unsigned idA = sys->nodeInfoVecs.host_edge_left[edge];
			unsigned idB = sys->nodeInfoVecs.host_edge_right[edge];
			bool is_A_collagen = sys->nodeInfoVecs.node_is_collagen[idA];
			bool is_B_collagen = sys->nodeInfoVecs.node_is_collagen[idB];
			if ((!is_A_collagen) || (!is_B_collagen) ){
				unsigned begin = idA * max_nbr_count;
				unsigned end = begin + max_nbr_count;
				double L0;
				for (unsigned i = begin; i < end; i++) {
					unsigned idTemp = sys->edgeInfoVecs.global_neighbors[i];
					if (idTemp == idB){
						L0 = sys->edgeInfoVecs.global_length_zero[i];
					}
				}
				double xL = sys->nodeInfoVecs.node_loc_x[idA];
				double yL = sys->nodeInfoVecs.node_loc_y[idA];
				double zL = sys->nodeInfoVecs.node_loc_z[idA];
				double xR = sys->nodeInfoVecs.node_loc_x[idB];
				double yR = sys->nodeInfoVecs.node_loc_y[idB];
				double zR = sys->nodeInfoVecs.node_loc_z[idB];

				double L1 = std::sqrt( (xL - xR)*(xL - xR)+(yL - yR)*(yL - yR)+(zL - zR)*(zL - zR));
				double strain = (L1 - L0) / L0;
				ofs << std::fixed << strain   << std::endl;
			}

		}

		ofs.close();

	}
};
