

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
	std::ofstream statesOutput(a_fileName + ".sta");
	std::ofstream statesOutputStrain(a_fileName + "_Strain.sta");

	std::shared_ptr<System> sysA = system.lock();
	std::shared_ptr<System_Builder> sysB = builder.lock();

	if ((sysA) && (sysB) ){
		unsigned max_node_count = sysA->generalParams.max_node_count;
		unsigned max_nbr_count = sysA->generalParams.max_nbr_count;

		statesOutput << "node_count " << max_node_count << '\n';
		statesOutput << "origin_node_count " << sysA->generalParams.origin_node_count << '\n';
		statesOutput << "origin_link_count " << sysA->generalParams.origin_edge_count << '\n';
		statesOutput << "sub_node_count " << sysA->generalParams.sub_node_count << std::endl;//system->getSubNodesSize() << '\n';
		statesOutput << "link_count " << sysA->generalParams.current_edge_count << '\n';

		for (unsigned edge = 0; edge < sysB->hostNodeInfoVecs.host_spring_edge_left.size(); edge++) {
			unsigned idLeft = sysB->hostNodeInfoVecs.host_spring_edge_left[edge];
			unsigned idRight = sysB->hostNodeInfoVecs.host_spring_edge_right[edge];
			statesOutput << '\n' << idLeft << ' ' << idRight;
		}

	}


	statesOutput.close();
}

void Storage::updateStrain() {

/*	std::shared_ptr<system> sys = system.lock();
	if (sys) {

	statesOutputStrain.open(bn + "_Strain.sta", std::ofstream::out | std::ofstream::app);
		statesOutputStrain << "\ntime " << sys->generalParams.currentTime;
		statesOutputStrain << "\nforce " << sys->extensionParams.totalAppliedForce;

		statesOutputStrain << "\nupper_XPos " << sys->domainParams.max_x;
		statesOutputStrain << "\nlower_XPos " << sys->domainParams.min_x;

		statesOutputStrain << "\nupper_YPos " << sys->domainParams.max_y;
		statesOutputStrain << "\nlower_YPos " << sys->domainParams.min_y;

		statesOutputStrain << "\nupper_ZPosAve " << sys->extensionParams.averageUpperStrain;
		statesOutputStrain << "\nlower_ZPosAve " << sys->extensionParams.averageLowerStrain;
		statesOutputStrain << "\noriginal_extended_percent " << sys->edgeInfoVecs.percentOriginalEdgesExtended;
		statesOutputStrain << "\noriginal_compressed_percent " << sys->edgeInfoVecs.percentOriginalEdgesCompressed;
		statesOutputStrain << "\noriginal_average_strain " << sys->edgeInfoVecs.averageStrainOriginalEdges;


		for (unsigned i = 0; i < sys->edgeInfoVecs.strainBucketOriginalNeg.size(); i++ ) {
			statesOutputStrain << " \noriginal_strain_neg " << sys->edgeInfoVecs.strainBucketOriginalNeg[i] / (2.0 * sys->generalParams.originEdgeCount);

		}
		for (unsigned i = 0; i < sys->edgeInfoVecs.strainBucketOriginalPos.size(); i++ ) {
			statesOutputStrain << " \noriginal_strain_pos " << sys->edgeInfoVecs.strainBucketOriginalPos[i] /  (2.0 * sys->generalParams.originEdgeCount);
		}




		statesOutputStrain << "\nadded_extended_percent " << sys->edgeInfoVecs.percentAddedEdgesExtended;
		statesOutputStrain << "\nadded_compressed_percent " << sys->edgeInfoVecs.percentAddedEdgesCompressed;
		statesOutputStrain << "\nadded_average_strain " << sys->edgeInfoVecs.averageStrainAddedEdges;

		double sumOfNumsAdded = std::accumulate(sys->edgeInfoVecs.strainBucketAddedNeg.begin(),
			sys->edgeInfoVecs.strainBucketAddedNeg.end(),0.0);
		sumOfNumsAdded += std::accumulate(sys->edgeInfoVecs.strainBucketAddedPos.begin(),
			sys->edgeInfoVecs.strainBucketAddedPos.end(),0.0);
		for (unsigned i = 0; i < sys->edgeInfoVecs.strainBucketAddedNeg.size(); i++ ) {
			statesOutputStrain << " \nadded_strain_neg " << sys->edgeInfoVecs.strainBucketAddedNeg[i]/sumOfNumsAdded;

		}
		for (unsigned i = 0; i < sys->edgeInfoVecs.strainBucketAddedPos.size(); i++ ) {
			statesOutputStrain << " \nadded_strain_pos " << sys->edgeInfoVecs.strainBucketAddedPos[i]/sumOfNumsAdded;

		}


		for (unsigned i = 0; i < sys->edgeInfoVecs.alignmentAverage.size(); i++ ) {
			double numEdgesInBin = sys->edgeInfoVecs.numberOfEdgesAlignment[i];
			double val = 0.0;

			if (numEdgesInBin != 0.0) {
				val  = sys->edgeInfoVecs.alignmentAverage[i]/numEdgesInBin;
			}
			statesOutputStrain << " \nslice_alignment " << val;

		}


	}
	statesOutputStrain.flush();
	statesOutputStrain.close();*/

};

void Storage::updateTotalStrain(void) {
	std::shared_ptr<System> sys = system.lock();
	if (sys) {

		double currentStrain = (sys->extensionParams.averageUpperStrain - sys->extensionParams.averageLowerStrain) /
			 (sys->extensionParams.originAverageUpperStrain - sys->extensionParams.originAverageLowerStrain ) - 1.0;
		//first create a new file using the current network strain

		std::string format = ".sta";
		std::string strain =  std::to_string(currentStrain);
		std::string initial = "StrainTest/Strain_";
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


void Storage::print_VTK_File() {

	std::shared_ptr<System> sys = system.lock();
	if (sys) {
		iteration+=1;
		unsigned digits = ceil(log10(iteration + 1));
		std::string format = ".vtk";
		std::string Number;
		std::string initial = "AnimationTest/Network_";
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


		unsigned max_node_count = sys->generalParams.max_node_count;
		unsigned max_nbr_count = sys->generalParams.max_nbr_count;
		unsigned numEdges = sys->generalParams.current_edge_count;

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
		unsigned numCells = numEdges;
		unsigned numNumsInCells = 3 * numEdges;


		ofs << "CELLS " << numCells << " " << numNumsInCells << std::endl;

		for (unsigned idA = 0; idA < max_node_count; idA++ ){

			unsigned beginIndex = idA * max_nbr_count;
			unsigned endIndex = beginIndex + max_nbr_count;
			for (unsigned i = beginIndex; i < endIndex; i++) {//currentSpringCount is the length of index and value vectors
				unsigned idB = sys->edgeInfoVecs.global_neighbors[i];//look through possible neighbors. May contain ULONG_MAX

				if ((idA < idB) && (idB < max_node_count) ) {
					ofs<< 2 << " " << idA << " " << idB << std::endl;
				}
			}
		}

		ofs << "CELL_TYPES " << numCells << std::endl;
		for (unsigned i = 0; i<numEdges; i++) {
			ofs << 3 << std::endl;
		}

		//
		ofs << "CELL_DATA " << numCells << std::endl;
		ofs << "SCALARS Fiber_Strain double " << std::endl;
		ofs << "LOOKUP_TABLE default "  << std::endl;
		for (unsigned edge = 0; edge < numEdges; edge++) {
			unsigned idA = sys->nodeInfoVecs.host_edge_left[edge];
			unsigned idB = sys->nodeInfoVecs.host_edge_right[edge];

			unsigned begin = idA * sys->generalParams.max_nbr_count;
			unsigned end = begin + sys->generalParams.max_nbr_count;
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

		ofs.close();

	}
};
