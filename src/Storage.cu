

#include "NodeSystemDevice.h"
#include "NodeSystemBuilder.h"

#include "NodeSystemDeviceFunctors.h"
#include "ForceDiagramStorage.h"
#include <numeric>

ForceDiagramStorage::ForceDiagramStorage(std::weak_ptr<NodeSystemDevice> a_system,
	std::weak_ptr<NodeSystemBuilder> b_system , const std::string& a_fileName) {
	//std::cout << "FDM constructor" << std::endl;
	system = a_system;
	builder = b_system;
	bn = a_fileName; //this will be used later to open files
	std::ofstream statesOutput(a_fileName + ".sta");
	std::ofstream statesOutputStrain(a_fileName + "_Strain.sta");

	std::shared_ptr<NodeSystemDevice> sysA = system.lock();
	std::shared_ptr<NodeSystemBuilder> sysB = builder.lock();
	
	if ((sysA) && (sysB) ){
		unsigned maxNodeCount = sysA->generalParams.maxNodeCount;
		unsigned maxNeighborCount = sysA->generalParams.maxNeighborCount;
		
		statesOutput << "node_count " << maxNodeCount << '\n';
		statesOutput << "origin_node_count " << sysA->generalParams.originNodeCount << '\n';
		statesOutput << "origin_link_count " << sysA->generalParams.originLinkCount << '\n';
		statesOutput << "sub_node_count " << sysA->generalParams.subNodeCount << std::endl;//system->getSubNodesSize() << '\n';
		statesOutput << "link_count " << sysA->generalParams.originEdgeCount << '\n';
		
		for (unsigned edge = 0; edge < sysB->hostWLCEdgeLeft.size(); edge++) {
			unsigned idLeft = sysB->hostWLCEdgeLeft[edge];
			unsigned idRight = sysB->hostWLCEdgeRight[edge];
			statesOutput << '\n' << idLeft << ' ' << idRight;
		}

	}


	statesOutput.close();
}

void ForceDiagramStorage::updateStrain() {
	
/*	std::shared_ptr<NodeSystemDevice> sys = system.lock();
	if (sys) {
		
	statesOutputStrain.open(bn + "_Strain.sta", std::ofstream::out | std::ofstream::app);
		statesOutputStrain << "\ntime " << sys->generalParams.currentTime;
		statesOutputStrain << "\nforce " << sys->compressionParams.totalAppliedForce;
		
		statesOutputStrain << "\nupper_XPos " << sys->domainParams.maxX;
		statesOutputStrain << "\nlower_XPos " << sys->domainParams.minX;
		
		statesOutputStrain << "\nupper_YPos " << sys->domainParams.maxY;
		statesOutputStrain << "\nlower_YPos " << sys->domainParams.minY;
		
		statesOutputStrain << "\nupper_ZPosAve " << sys->compressionParams.averageUpperStrain;
		statesOutputStrain << "\nlower_ZPosAve " << sys->compressionParams.averageLowerStrain;
		statesOutputStrain << "\noriginal_extended_percent " << sys->wlcInfoVecs.percentOriginalEdgesExtended;
		statesOutputStrain << "\noriginal_compressed_percent " << sys->wlcInfoVecs.percentOriginalEdgesCompressed;
		statesOutputStrain << "\noriginal_average_strain " << sys->wlcInfoVecs.averageStrainOriginalEdges;


		for (unsigned i = 0; i < sys->wlcInfoVecs.strainBucketOriginalNeg.size(); i++ ) {
			statesOutputStrain << " \noriginal_strain_neg " << sys->wlcInfoVecs.strainBucketOriginalNeg[i] / (2.0 * sys->generalParams.originEdgeCount);
	
		}		
		for (unsigned i = 0; i < sys->wlcInfoVecs.strainBucketOriginalPos.size(); i++ ) {
			statesOutputStrain << " \noriginal_strain_pos " << sys->wlcInfoVecs.strainBucketOriginalPos[i] /  (2.0 * sys->generalParams.originEdgeCount);
		}
	

 

		statesOutputStrain << "\nadded_extended_percent " << sys->wlcInfoVecs.percentAddedEdgesExtended;
		statesOutputStrain << "\nadded_compressed_percent " << sys->wlcInfoVecs.percentAddedEdgesCompressed;
		statesOutputStrain << "\nadded_average_strain " << sys->wlcInfoVecs.averageStrainAddedEdges;
		
		double sumOfNumsAdded = std::accumulate(sys->wlcInfoVecs.strainBucketAddedNeg.begin(),
			sys->wlcInfoVecs.strainBucketAddedNeg.end(),0.0);
		sumOfNumsAdded += std::accumulate(sys->wlcInfoVecs.strainBucketAddedPos.begin(),
			sys->wlcInfoVecs.strainBucketAddedPos.end(),0.0);
		for (unsigned i = 0; i < sys->wlcInfoVecs.strainBucketAddedNeg.size(); i++ ) {
			statesOutputStrain << " \nadded_strain_neg " << sys->wlcInfoVecs.strainBucketAddedNeg[i]/sumOfNumsAdded;
	
		}		
		for (unsigned i = 0; i < sys->wlcInfoVecs.strainBucketAddedPos.size(); i++ ) {
			statesOutputStrain << " \nadded_strain_pos " << sys->wlcInfoVecs.strainBucketAddedPos[i]/sumOfNumsAdded;
	
		}


		for (unsigned i = 0; i < sys->wlcInfoVecs.alignmentAverage.size(); i++ ) {
			double numEdgesInBin = sys->wlcInfoVecs.numberOfEdgesAlignment[i];
			double val = 0.0;
				
			if (numEdgesInBin != 0.0) {
				val  = sys->wlcInfoVecs.alignmentAverage[i]/numEdgesInBin;
			}
			statesOutputStrain << " \nslice_alignment " << val;
	
		}		
		
		
	} 
	statesOutputStrain.flush();
	statesOutputStrain.close();*/
	
};

void ForceDiagramStorage::updateTotalStrain(void) {
	std::shared_ptr<NodeSystemDevice> sys = system.lock();
	if (sys) {

		double currentStrain = (sys->compressionParams.averageUpperStrain - sys->compressionParams.averageLowerStrain) /
			 (sys->compressionParams.originAverageUpperStrain - sys->compressionParams.originAverageLowerStrain ) - 1.0;
		//first create a new file using the current network strain
		
		std::string format = ".sta";
		std::string strain =  std::to_string(currentStrain);
		std::string initial = "StrainTest/Strain_";
		std::ofstream ofs;
		std::string Filename = initial + strain + format;
		ofs.open(Filename.c_str());



		unsigned maxNeighborCount = sys->generalParams.maxNeighborCount;
		unsigned maxNodeCount = sys->generalParams.maxNodeCount;
		unsigned originalNodeCount = sys->generalParams.originNodeCount;
		unsigned originalEdgeCount = sys->generalParams.originLinkCount;
		unsigned edgeCountDiscretize = sys->generalParams.originEdgeCount;
		//Now first place strain
		ofs << std::setprecision(5) <<std::fixed<< "time " << sys->generalParams.currentTime<<std::endl;
		ofs << std::setprecision(5) <<std::fixed<< "network_strain " << currentStrain<<std::endl;
		ofs << std::setprecision(5) <<std::fixed<< "minX " << sys->domainParams.minX<<std::endl;
		ofs << std::setprecision(5) <<std::fixed<< "maxX " << sys->domainParams.maxX<<std::endl;
		ofs << std::setprecision(5) <<std::fixed<< "minY " << sys->domainParams.minY<<std::endl;
		ofs << std::setprecision(5) <<std::fixed<< "maxY " << sys->domainParams.maxY<<std::endl;
		ofs << std::setprecision(5) <<std::fixed<< "minZ " << sys->domainParams.minX<<std::endl;
		ofs << std::setprecision(5) <<std::fixed<< "maxZ " << sys->domainParams.maxX<<std::endl;
		
		
		ofs << std::setprecision(5) <<std::fixed<< "total_applied_force " << sys->compressionParams.totalAppliedForce<<std::endl;
		ofs << std::setprecision(5) <<std::fixed<< "original_node_count " << originalNodeCount <<std::endl;
		ofs << std::setprecision(5) <<std::fixed<< "node_count_discretize " << maxNodeCount <<std::endl;
		ofs << std::setprecision(5) <<std::fixed<< "original_edge_count " << originalEdgeCount <<std::endl;
		ofs << std::setprecision(5) <<std::fixed<< "edge_count_discretize " << edgeCountDiscretize <<std::endl;
		
		//place nodes
		for (unsigned i = 0; i < sys->nodeInfoVecs.nodeLocX.size(); i++) {
			double x = sys->nodeInfoVecs.nodeLocX[i];
			double y = sys->nodeInfoVecs.nodeLocY[i];
			double z = sys->nodeInfoVecs.nodeLocZ[i];
			ofs << std::setprecision(5) <<std::fixed<< "node " << x << " " << y << " " << z <<std::endl;
		
		}
		//place force node is experiencing
		for (unsigned i = 0; i < sys->nodeInfoVecs.nodeLocX.size(); i++) {
			ofs << std::setprecision(5) <<std::fixed<< "force_on_node " << sys->nodeInfoVecs.sumForcesOnNode[i]<<std::endl;
		
		}

		//place original edges
		for (unsigned edge = 0; edge < sys->generalParams.originEdgeCount; edge++) {
			unsigned idL = sys->nodeInfoVecs.deviceEdgeLeft[edge];
			unsigned idR = sys->nodeInfoVecs.deviceEdgeRight[edge];
			ofs <<"original_edge_discretized " <<idL <<" "<< idR <<std::endl;
			
		}
				 
		//place added edges
		for (unsigned edge = sys->generalParams.originEdgeCount; edge < sys->generalParams.currentEdgeCount; edge++) {
			unsigned idL = sys->nodeInfoVecs.deviceEdgeLeft[edge];
			unsigned idR = sys->nodeInfoVecs.deviceEdgeRight[edge];
			ofs <<"added_edge " <<idL <<" "<< idR <<std::endl;
			
		}

		//original edge strain
		for (unsigned i = 0; i < sys->generalParams.originEdgeCount; i++ ){
			double val = sys->nodeInfoVecs.discretizedEdgeStrain[i];

			ofs << std::setprecision(5)<< std::fixed<<"original_edge_strain " << val <<std::endl;
		}
				
		//original edge alignment
		for (unsigned i = 0; i < sys->generalParams.originEdgeCount; i++ ){
			double val = sys->nodeInfoVecs.discretizedEdgeAlignment[i];
			ofs << std::setprecision(5)<< std::fixed<<"original_edge_alignment " << val <<std::endl;
		}

		//added edge strain
		for (unsigned i = sys->generalParams.originEdgeCount; i < sys->generalParams.currentEdgeCount; i++ ){
			double val = sys->nodeInfoVecs.discretizedEdgeStrain[i];
			ofs << std::setprecision(5)<< std::fixed<<"added_edge_strain " << val <<std::endl;
		}
		
		//added links per node.
		for (unsigned i = 0; i < sys->generalParams.maxNodeCount; i++ ){
			unsigned val = sys->wlcInfoVecs.currentNodeEdgeCountVector[i] - 
				sys->wlcInfoVecs.numOriginalNeighborsNodeVector[i];
			ofs << std::setprecision(5)<< std::fixed<<"bind_sites_per_node " << val <<std::endl;
		}



	}
}


void ForceDiagramStorage::print_VTK_File() {
	
	std::shared_ptr<NodeSystemDevice> sys = system.lock();
	if (sys) {	
		iteration+=1;
		unsigned digits = ceil(log10(iteration + 1));
		std::string format = ".vtk";
		std::string Number;
		std::string initial = "AnimationTest/FibrinNetwork_";
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
		
	
		unsigned maxNodeCount = sys->generalParams.maxNodeCount;
		unsigned maxNeighborCount = (sys->generalParams).maxNeighborCount;
		
		unsigned numEdges = sys->generalParams.currentEdgeCount;//sys->nodeInfoVecs.hostEdgeRight.size();
		
		ofs << "# vtk DataFile Version 3.0" << std::endl;
		ofs << "Point representing Sub_cellular elem model" << std::endl;
		ofs << "ASCII" << std::endl << std::endl;
		ofs << "DATASET UNSTRUCTURED_GRID" << std::endl;
		
		 
		ofs << "POINTS " << maxNodeCount << " float" << std::endl;
		for (unsigned i = 0; i< maxNodeCount; i++) {
			double xPos = sys->nodeInfoVecs.nodeLocX[i];
			double yPos = sys->nodeInfoVecs.nodeLocY[i];
			double zPos = sys->nodeInfoVecs.nodeLocZ[i];
			
			ofs << std::setprecision(5) <<std::fixed<< xPos << " " << yPos << " " << zPos << " " << '\n'<< std::fixed;
		}
		//now plot particles

		
		unsigned numCells = numEdges;
		unsigned numNumsInCells = 3 * numEdges;
		
		
		ofs << "CELLS " << numCells << " " << numNumsInCells << std::endl;
		/*for (unsigned i = 0; i< numEdges; i++) {
			unsigned idL = sys->nodeInfoVecs.hostEdgeLeft[i];
			
			unsigned idR = sys->nodeInfoVecs.hostEdgeRight[i];
			ofs<< 2 << " " << idL << " " << idR << std::endl;
		}*/

		
		
		for (unsigned idA = 0; idA < maxNodeCount; idA++ ){
			
			unsigned beginIndex = idA * maxNeighborCount;
			unsigned endIndex = beginIndex + maxNeighborCount;
			for (unsigned i = beginIndex; i < endIndex; i++) {//currentSpringCount is the length of index and value vectors
				unsigned idB = sys->wlcInfoVecs.globalNeighbors[i];//look through possible neighbors. May contain ULONG_MAX
	
				if ((idA < idB) && (idB < maxNodeCount) ) {
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
		ofs << "SCALARS magnitude double " << std::endl;
		ofs << "LOOKUP_TABLE default "  << std::endl;
		
		for (unsigned idA = 0; idA < maxNodeCount; idA++ ){
			
			unsigned beginIndex = idA * maxNeighborCount;
			unsigned endIndex = beginIndex + maxNeighborCount;
			for (unsigned i = beginIndex; i < endIndex; i++) {//currentSpringCount is the length of index and value vectors
				unsigned idB = sys->wlcInfoVecs.globalNeighbors[i];//look through possible neighbors. May contain ULONG_MAX
	
				if ((idA < idB) && (idB < maxNodeCount) ) {
					unsigned index = idA * maxNeighborCount + idB;
					double L0 = sys->wlcInfoVecs.lengthZero[i];
					double xL = sys->nodeInfoVecs.nodeLocX[idA];
					double yL = sys->nodeInfoVecs.nodeLocY[idA];
					double zL = sys->nodeInfoVecs.nodeLocZ[idA];
					double xR = sys->nodeInfoVecs.nodeLocX[idB];
					double yR = sys->nodeInfoVecs.nodeLocY[idB];
					double zR = sys->nodeInfoVecs.nodeLocZ[idB];
					
				
					
					double L1 = std::sqrt( (xL - xR)*(xL - xR)+(yL - yR)*(yL - yR)+(zL - zR)*(zL - zR));
					double strain = (L1 - L0) / L0;
					ofs << std::fixed << strain   << std::endl;
				}
			}
		}	

		ofs.close();
	
	}
}

void ForceDiagramStorage::updateStorage() {

	//currentAddedEdges = (system->getDynamicLinks().size() - previousAddedEdges);
	
	std::shared_ptr<NodeSystemDevice> sys = system.lock();
	if (sys) {
		statesOutput.open(bn + ".sta", std::ofstream::out | std::ofstream::app);
		//output.open(bn + ".grm", std::ofstream::out | std::ofstream::app);
		statesOutput << "\nextended percent " << sys->wlcInfoVecs.percentOriginalEdgesExtended;
		statesOutput << "\nforce " << sys->compressionParams.totalAppliedForce;
		statesOutput << "\ntime " << sys->generalParams.currentTime;
		statesOutput << "\nadded edges " << ((sys->nodeInfoVecs.idEdgesMadeHost).size());

		unsigned maxNodeCount = sys->generalParams.maxNodeCount;

		//print new added edges	for current time step recording
		
		for (unsigned i = 0; i < (sys->nodeInfoVecs.idEdgesMadeHost.size()); i++) {
			unsigned idUpper = sys->nodeInfoVecs.idEdgesMadeHost[i];
			if (idUpper != 0) {
				unsigned first = idUpper - maxNodeCount*(idUpper / maxNodeCount); //represents column
				unsigned second = (idUpper / maxNodeCount); //represents row

					statesOutput << '\n' << first << ' ' << second;
			}

		}
		



		for (unsigned i = 0; i < maxNodeCount; ++i) {


			double xPos = sys->nodeInfoVecs.nodeLocX[i];
			double yPos = sys->nodeInfoVecs.nodeLocY[i];
			double zPos = sys->nodeInfoVecs.nodeLocZ[i];
			double xForce = sys->nodeInfoVecs.nodeVelX[i];
			double yForce = sys->nodeInfoVecs.nodeVelY[i];
			double zForce = sys->nodeInfoVecs.nodeVelZ[i];
			double sumOfForces = sys->nodeInfoVecs.sumForcesOnNode[i];
			statesOutput << '\n' << i;

			//auto pos = node->getPosition();
			//auto vel = node->getVelocity();
 
			//for (int k = 0; k < 3; ++k)
			statesOutput << ' ' << xPos << ' ' << yPos << ' ' << zPos;

			//for (int k = 0; k < 3; ++k)
			statesOutput << ' ' << xForce << ' ' << yForce << ' ' << zForce << ' ' << sumOfForces;

		}
	}

	output << magnitudeForce << ' ' << std::endl;
	statesOutput.flush();
	output.flush();
	statesOutput.close();
	output.close();

	std::cout << "*** one step completed ***\n\n";
					 
}




