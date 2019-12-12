/*
* NodeSystemImplDevice.cu
*
* Created on 8/1/2017
* 		Author: SRB
*/


 
//#include <thrust/version.h>

//#include <cuda_runtime.h> 
//#include <cuda.h>
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
#include "ForceDiagramStorage.h"
#include "IncrementExternalForceOnDevice.h"
#include "LinkNodesOnDevice.h"
#include "CalculateEquilibrium.h"
#include "WLCSolveOnDevice.h"
#include "TorsionSolveOnDevice.h"
#include "AdvancePositionOnDevice.h"
#include "BucketSchemeOnDevice.h"
#include "DPDParticle.h"
#include "NodeSystemDevice.h"
#include "NodeSystemDeviceFunctors.h"




using namespace thrust::placeholders;

										  
double NodeSystemDevice::solveForcesOnDevice() {
	//std::cout<<"force: " << nodeInfoVecs.nodeForceX[10] <<" " << nodeInfoVecs.nodeForceY[10] <<" "<<  nodeInfoVecs.nodeForceZ[10]<<std::endl;
	//std::cout<<"loc: " << nodeInfoVecs.nodeLocX[10] <<" " << nodeInfoVecs.nodeLocY[10] <<" "<<  nodeInfoVecs.nodeLocZ[10]<<std::endl;
	//////////////////////////////////////////////////////////////////////////////////////////
	//RESET FORCE TO ZERO AT BEGINNING/////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////
	thrust::fill(nodeInfoVecs.nodeForceX.begin(),nodeInfoVecs.nodeForceX.end(),0);
	thrust::fill(nodeInfoVecs.nodeForceY.begin(),nodeInfoVecs.nodeForceY.end(),0);
	thrust::fill(nodeInfoVecs.nodeForceZ.begin(),nodeInfoVecs.nodeForceZ.end(),0);
	
	//std::cout<<"post fill zeros" << std::endl;
	//std::cout<<"force: " << nodeInfoVecs.nodeForceX[10] <<" " << nodeInfoVecs.nodeForceY[10] <<" "<<  nodeInfoVecs.nodeForceZ[10]<<std::endl;
	///////////////////////////////////////////////////////////////////////////////////////
	//RESET FORCE TO ZERO AT BEGINNING////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////

	try{initDimensionBucketScheme(
			nodeInfoVecs,
			domainParams, 
			auxVecs, 
			generalParams,
			dpdParticleVariables,
			compressionParams);} //reset dimensions before bucketting domain. Possibly replace with larger domain?}
	catch(thrust::system_error &e){std::cerr << "Error initializing buckets: " << e.what() << std::endl; exit(-1);}
	
	try{buildBucketScheme(nodeInfoVecs, domainParams, 
			auxVecs, generalParams, dpdParticleVariables);}
	catch(thrust::system_error &e){std::cerr << "Error building buckets: " << e.what() << std::endl; exit(-1);}
	
	try{extendBucketScheme(nodeInfoVecs, domainParams, auxVecs);}
	catch(thrust::system_error &e){std::cerr << "Error extending buckets: " << e.what() << std::endl; exit(-1);}
	cudaThreadSynchronize();

	/////////////////////////////////////////////////////////////////////////////
	//safety feature begin //////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////
/*	double tempSizes[] = {std::abs(domainParams.maxX),std::abs(domainParams.minX),
					std::abs(domainParams.maxY), std::abs(domainParams.minY),
					std::abs(domainParams.maxY), std::abs(domainParams.minY)};

  	// using default comparison:
  	double tempMax = *std::max_element(tempSizes,tempSizes+6);

	if ( (tempMax > 10 * compressionParams.originalNetworkLength) ||
		(tempMax > 10 * compressionParams.originalNetworkWidth) ) {
		generalParams.runSim = false;
		std::cout<<"safety feature invoked" << std::endl;
	}
	
	//stop sim if length is larger than strain proportion times original length
	if (abs(compressionParams.averageUpperStrain - compressionParams.averageLowerStrain) > (compressionParams.strainProportion * compressionParams.originalNetworkLength)) {
		generalParams.runSim = false;                                                                             
		std::cout<<"stopping sim, maximumsize reached"<<std::endl;
		 
	}*/
	////////////////////////////////////////////////////////////////////////////
	//safety feature end////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////
	
	/////////////////////////////////////////////////////////////////////
	//LINKING BEGIN//////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////	
	//test turning linking off. 
	double addedLinks = generalParams.currentEdgeCount - generalParams.originEdgeCount;
	

	if (generalParams.linking == true) {
		//default is linking turned on, i.e. Linking = false
		try {
			LinkNodesOnDevice(
					nodeInfoVecs,
					wlcInfoVecs,
					auxVecs,
					torsionInfoVecs,
					generalParams);
			}
		catch(thrust::system_error &e) { 

			std::cerr << "Error linking: " << e.what() << std::endl;		
		}
	}
	cudaThreadSynchronize();

	//apply external force. 
	try {				 									
	IncrementExternalForceOnDevice(nodeInfoVecs, 
		generalParams,
		compressionParams,
		domainParams);
		cudaThreadSynchronize();
	 }
	
	catch(thrust::system_error &e){std::cerr << "Error applying Force : " << e.what() << std::endl; exit(-1);}

	//only counts external force on network nodes since force has been reset. 
	try { compressionParams.totalAppliedForce = (thrust::transform_reduce(
			thrust::make_zip_iterator(
				thrust::make_tuple(
					nodeInfoVecs.nodeForceX.begin(),
					nodeInfoVecs.nodeForceY.begin(),
					nodeInfoVecs.nodeForceZ.begin())),
			thrust::make_zip_iterator(				
				thrust::make_tuple(
					nodeInfoVecs.nodeForceX.begin(),
					nodeInfoVecs.nodeForceY.begin(),
					nodeInfoVecs.nodeForceZ.begin())) + generalParams.maxNodeCount,
				NormFunctor(), 0.0, thrust::plus<double>() ) ); } 
	catch(thrust::system_error &e){
		std::cerr << "Error reduce in total force applied: " << e.what() << std::endl; exit(-1);}
					
				//std::cout<<"totalApplied Force: "<< compressionParams.totalAppliedForce<<std::endl;
	try { TorsionSolveOnDevice(nodeInfoVecs, torsionInfoVecs, generalParams); 
		cudaThreadSynchronize(); }
	catch(thrust::system_error &e){std::cerr << "Error Torsion: " << e.what() << std::endl; exit(-1);}
	

	try { WLCSolveOnDevice(nodeInfoVecs, wlcInfoVecs, generalParams); 
		cudaThreadSynchronize(); }
	catch(thrust::system_error &e){std::cerr << "Error WLC : " << e.what() << std::endl; exit(-1);}
	
	return 0.0;//compressionParams.currentNetworkLength;
};


void NodeSystemDevice::solveSystemDevice() {
	//(generalParams.magnitudeForce <= generalParams.maxForce)
	//(compressionParams.targetStrain <= generalParams.maxForce)
	
	double lastTime = 0.0;
	storage->updateStorage();//initial position storage
	bool runIters = true;

	//set initial epsilon
	generalParams.epsilon = (1.0) * 
		sqrt(6.0*generalParams.kB * generalParams.temperature * generalParams.dtTemp / generalParams.viscousDamp);
	std::cout<<"new gen eps begin: "<< generalParams.epsilon<<std::endl;

	while (runIters == true) {
		
		generalParams.iterationCounter++;
		generalParams.currentTime += generalParams.dtTemp;
		//std::cout << "current time: " << std::endl;

		double unused = AdvancePositionOnDevice(
			nodeInfoVecs,
		 	generalParams,
			dpdParticleVariables);
		

		if ((generalParams.iterationCounter % 20000) == 0) {
			storage->print_VTK_File();
			//unsigned maxLinked = *( thrust::max_element(wlcInfoVecs.currentNodeEdgeCountVector.begin(),wlcInfoVecs.currentNodeEdgeCountVector.end()) );
			//std::cout<<"maximum neighbors currently: "<< maxLinked<<std::endl;
			//std::cout<<"maximum neighbors: "<< generalParams.maxNeighborCount <<std::endl;
			
		}

		generalParams.currentLength = solveForcesOnDevice(); //resets and solves forces for next time step
		cudaThreadSynchronize();
		if ((generalParams.iterationCounter % 20000) == 0) {
			double currentStrain = (compressionParams.averageUpperStrain - compressionParams.averageLowerStrain) /
			(compressionParams.originAverageUpperStrain - compressionParams.originAverageLowerStrain ) - 1.0;
			if (currentStrain>4.0){
				runIters=false;
			} 
			thrust::transform( 
				thrust::make_zip_iterator(
					thrust::make_tuple( 
						nodeInfoVecs.nodeForceX.begin(),
						nodeInfoVecs.nodeForceY.begin(),
						nodeInfoVecs.nodeForceZ.begin())),
				thrust::make_zip_iterator(
					thrust::make_tuple(
						nodeInfoVecs.nodeForceX.begin(),
						nodeInfoVecs.nodeForceY.begin(),
						nodeInfoVecs.nodeForceZ.begin())) + generalParams.maxNodeCount,
				nodeInfoVecs.sumForcesOnNode.begin(),//save vector
				NormFunctor()); 
			
			GetStrainParameters(nodeInfoVecs,
				wlcInfoVecs,  
				generalParams,
				domainParams);
			storage->updateTotalStrain();
		}


		double maxVel = *(thrust::max_element(nodeInfoVecs.nodeVelocity.begin(), nodeInfoVecs.nodeVelocity.end()));
		//std::cout<<"maxvelocity: "<< maxVel<< std::endl;

		//difference in time 
 		if (abs(generalParams.currentTime - lastTime) > (generalParams.lagTime)) {
			 //move epsilon. It will be reset 
 
			generalParams.epsilon += 0.01;
			lastTime = generalParams.currentTime; 
			
			std::cout<<"updating epsilon: "<< generalParams.epsilon<<std::endl;

			double addedEdges = generalParams.currentEdgeCount - generalParams.originEdgeCount;
			std::cout<<"added edges: "<< addedEdges <<std::endl;
			std::cout<<"Minz: "<< domainParams.minZ<<std::endl;
			std::cout<<"Maxz: "<< domainParams.maxZ<<std::endl;
			std::cout<<"Miny: "<< domainParams.minY<<std::endl;
			std::cout<<"Maxy: "<< domainParams.maxY<<std::endl;
			std::cout<<"Minx: "<< domainParams.minX<<std::endl;
			std::cout<<"Maxx: "<< domainParams.maxX<<std::endl;
		} 
	
		if (maxVel < generalParams.epsilon) {
		//only reached if dtMax == dtTemp
			//lastTime = generalParams.currentTime;
			
		
			//store sum of all forces on each node. Used in stress calculations
			thrust::transform( 
				thrust::make_zip_iterator(
					thrust::make_tuple( 
						nodeInfoVecs.nodeForceX.begin(),
						nodeInfoVecs.nodeForceY.begin(),
						nodeInfoVecs.nodeForceZ.begin())),
				thrust::make_zip_iterator(
					thrust::make_tuple(
						nodeInfoVecs.nodeForceX.begin(),
						nodeInfoVecs.nodeForceY.begin(),
						nodeInfoVecs.nodeForceZ.begin())) + generalParams.maxNodeCount,
				nodeInfoVecs.sumForcesOnNode.begin(),//save vector
				NormFunctor()); 
				
			
			
			storage->updateStorage();
			
			generalParams.totalNumberOfEdges += nodeInfoVecs.idEdgesMadeHost.size();
			nodeInfoVecs.idEdgesMadeHost.resize(0);
	
			
			generalParams.epsilon = (1.0) * 
				sqrt(6.0 * generalParams.kB * generalParams.temperature * generalParams.dtTemp / generalParams.viscousDamp);

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



NodeSystemDevice::NodeSystemDevice()  {};

void NodeSystemDevice::assignForceDiagramStorage(std::shared_ptr<ForceDiagramStorage> _storage) {
	storage = _storage;
} 

//__host__ __device__
void NodeSystemDevice::initializeSystem(
	thrust::host_vector<bool>& hostIsNodeFixed,
	thrust::host_vector<double>& hostPosX,
	thrust::host_vector<double>& hostPosY,
	thrust::host_vector<double>& hostPosZ,
	thrust::host_vector<unsigned>& hostWLCEdgeLeft,
	thrust::host_vector<unsigned>& hostWLCEdgeRight,
	thrust::host_vector<double>& hostWLCLenZero,

	thrust::host_vector<unsigned>& hostWLCSubEdgeLeft,
	thrust::host_vector<unsigned>& hostWLCSubEdgeRight,
	thrust::host_vector<double>& hostWLCSubLenZero,
	thrust::host_vector<unsigned>& hostSpringDivisionCount,
	thrust::host_vector<unsigned>& hostTorsionIndexLeft,
	thrust::host_vector<unsigned>& hostTorsionIndexCenter,
	thrust::host_vector<unsigned>& hostTorsionIndexRight,
	thrust::host_vector<double>& hostTorsionAngleZero) {
	
	std::cout<< "total Edge Count: "<< generalParams.originEdgeCount << std::endl;
	std::cout << "max num nodes: " << generalParams.maxNodeCount << std::endl;

	nodeInfoVecs.hostOriginalEdgeLeft = hostWLCEdgeLeft;
	nodeInfoVecs.hostOriginalEdgeRight = hostWLCEdgeRight;

	setNodeVecs(//calls initDimensionBucketScheme
		hostIsNodeFixed,
		hostPosX,
		hostPosY,
		hostPosZ,
		hostSpringDivisionCount);
		
	

	setTorsionVecs(
		hostTorsionIndexLeft,
		hostTorsionIndexCenter,
		hostTorsionIndexRight,
		hostTorsionAngleZero);

	setWLCVecs(	hostWLCSubEdgeLeft,
				hostWLCSubEdgeRight,
				hostWLCSubLenZero );

	setExtras();
};


void NodeSystemDevice::setNodeVecs(
	thrust::host_vector<bool>& hostIsNodeFixed,
	thrust::host_vector<double>& hostPosX,
	thrust::host_vector<double>& hostPosY,
	thrust::host_vector<double>& hostPosZ,
	thrust::host_vector<unsigned>& hostSpringDivisionCount) {


	nodeInfoVecs.idEdgesMadeTemp.resize(generalParams.maxNodeCount * generalParams.maxLinksPerIteration);//corresponds to upperAdj vector size plus a single value to hold number of added nodes
	thrust::fill(nodeInfoVecs.idEdgesMadeTemp.begin(), nodeInfoVecs.idEdgesMadeTemp.end(), 0);

	nodeInfoVecs.sumForcesOnNode.resize(generalParams.maxNodeCount);

	nodeInfoVecs.nodeUpperChoiceForStrain.resize(generalParams.maxNodeCount);
	nodeInfoVecs.nodeLowerChoiceForStrain.resize(generalParams.maxNodeCount);

	nodeInfoVecs.springDivisionCount.resize(generalParams.maxNodeCount);

	
	nodeInfoVecs.prevNodeLocX.resize(generalParams.maxNodeCount);
	nodeInfoVecs.prevNodeLocY.resize(generalParams.maxNodeCount);
	nodeInfoVecs.prevNodeLocZ.resize(generalParams.maxNodeCount);
	nodeInfoVecs.prevNodeVelX.resize(generalParams.maxNodeCount);
	nodeInfoVecs.prevNodeVelY.resize(generalParams.maxNodeCount);
	nodeInfoVecs.prevNodeVelZ.resize(generalParams.maxNodeCount);
	nodeInfoVecs.prevNodeForceX.resize(generalParams.maxNodeCount);
	nodeInfoVecs.prevNodeForceY.resize(generalParams.maxNodeCount);
	nodeInfoVecs.prevNodeForceZ.resize(generalParams.maxNodeCount);

	nodeInfoVecs.nodeVelocity.resize(generalParams.maxNodeCount);
	
	nodeInfoVecs.nodeLocX.resize(generalParams.maxNodeCount);
	nodeInfoVecs.nodeLocY.resize(generalParams.maxNodeCount);
	nodeInfoVecs.nodeLocZ.resize(generalParams.maxNodeCount);
	nodeInfoVecs.nodeVelX.resize(generalParams.maxNodeCount);
	nodeInfoVecs.nodeVelY.resize(generalParams.maxNodeCount);
	nodeInfoVecs.nodeVelZ.resize(generalParams.maxNodeCount);
	
	
	nodeInfoVecs.nodeForceX.resize(generalParams.maxNodeCount); 
	nodeInfoVecs.nodeForceY.resize(generalParams.maxNodeCount);
	nodeInfoVecs.nodeForceZ.resize(generalParams.maxNodeCount);

	nodeInfoVecs.discretizedEdgeStrain.resize(generalParams.maxNodeCount * generalParams.maxNeighborCount);
	nodeInfoVecs.discretizedEdgeAlignment.resize(generalParams.maxNodeCount * generalParams.maxNeighborCount);
	
	//sized larger for input later
	nodeInfoVecs.deviceEdgeLeft.resize(generalParams.maxNodeCount * generalParams.maxNeighborCount);
	nodeInfoVecs.deviceEdgeRight.resize(generalParams.maxNodeCount * generalParams.maxNeighborCount);


	thrust::fill(nodeInfoVecs.discretizedEdgeStrain.begin(), nodeInfoVecs.discretizedEdgeStrain.end(),0.0);
	thrust::fill(nodeInfoVecs.deviceEdgeRight.begin(), nodeInfoVecs.deviceEdgeRight.end(), 0);	//fill force and velocity with zeros for computation.
	thrust::fill(nodeInfoVecs.deviceEdgeLeft.begin(), nodeInfoVecs.deviceEdgeLeft.end(), 0);	//fill force and velocity with zeros for computation.
	thrust::fill(nodeInfoVecs.idEdgesMadeTemp.begin(), nodeInfoVecs.idEdgesMadeTemp.end(), 0);
	
	thrust::fill(nodeInfoVecs.sumForcesOnNode.begin(), nodeInfoVecs.sumForcesOnNode.end(), 0);

	thrust::fill(nodeInfoVecs.nodeUpperChoiceForStrain.begin(), 
		nodeInfoVecs.nodeUpperChoiceForStrain.end(),false);
		
	thrust::fill(nodeInfoVecs.nodeLowerChoiceForStrain.begin(), 
		nodeInfoVecs.nodeLowerChoiceForStrain.end(),false);

	thrust::copy(hostSpringDivisionCount.begin(),hostSpringDivisionCount.end(), nodeInfoVecs.springDivisionCount.begin());
 



	thrust::copy(hostPosX.begin(), hostPosX.end(), nodeInfoVecs.prevNodeLocX.begin());
	thrust::copy(hostPosY.begin(), hostPosY.end(), nodeInfoVecs.prevNodeLocY.begin());
	thrust::copy(hostPosZ.begin(), hostPosZ.end(), nodeInfoVecs.prevNodeLocZ.begin());
	thrust::copy(hostPosX.begin(), hostPosX.end(), nodeInfoVecs.nodeLocX.begin());
	thrust::copy(hostPosY.begin(), hostPosY.end(), nodeInfoVecs.nodeLocY.begin());
	thrust::copy(hostPosZ.begin(), hostPosZ.end(), nodeInfoVecs.nodeLocZ.begin());
	
	
	//copy fixed positions
	nodeInfoVecs.isNodeFixed.resize(generalParams.maxNodeCount);
	thrust::fill(nodeInfoVecs.isNodeFixed.begin(), nodeInfoVecs.isNodeFixed.end(), false);

	//now that all the nodes are loaded in, choose the top to apply strain, and fix the bottom
	
	determineBounds();
	
	//at this point all nodes are filled, so we can generate domainParams before seeding dpd particles. 
	initDimensionBucketScheme(
		nodeInfoVecs, 
		domainParams, 
		auxVecs, 
		generalParams, 
		dpdParticleVariables,
		compressionParams); 
	
	//set original parameters for domain. others will be reset as simulation takes place. 
	domainParams.originMinX = domainParams.minX;
	domainParams.originMaxX = domainParams.maxX;
	domainParams.originMinY = domainParams.minY; 
	domainParams.originMaxY = domainParams.maxY;
	domainParams.originMinZ = domainParams.minZ;
	domainParams.originMaxZ = domainParams.maxZ;
	std::cout<< "node count : " <<nodeInfoVecs.nodeLocY.size()<< std::endl;


	auxVecs.bucketKeys.resize(generalParams.maxNodeCount + dpdParticleVariables.particleCount);
	auxVecs.bucketValues.resize(generalParams.maxNodeCount + dpdParticleVariables.particleCount);
	auxVecs.bucketValuesIncludingNeighbor.resize(27 * (generalParams.maxNodeCount + dpdParticleVariables.particleCount));
	auxVecs.bucketKeysExpanded.resize(27 *( generalParams.maxNodeCount + dpdParticleVariables.particleCount));

};

void NodeSystemDevice::determineBounds() {
	//determin z positions of nodes to be pulled and fixed. 
	
	thrust::device_vector<double> zPosTemp;
	zPosTemp.resize(generalParams.maxNodeCount);
	thrust::copy(nodeInfoVecs.nodeLocZ.begin(), nodeInfoVecs.nodeLocZ.end(), zPosTemp.begin());

	//not used
	//pull at least 10% of nodes. 
	unsigned tempNodeAmmount = static_cast<unsigned>( 0.25 * generalParams.maxNodeCount ); //pull 10% of top nodes
	
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

	thrust::replace_if(nodeInfoVecs.nodeUpperChoiceForStrain.begin(), nodeInfoVecs.nodeUpperChoiceForStrain.begin() + generalParams.originNodeCount, 
						nodeInfoVecs.nodeLocZ.begin(),  
						IsGreaterThanLevel( upperLevelAlt ), true);
						
	thrust::replace_if(nodeInfoVecs.nodeLowerChoiceForStrain.begin(), nodeInfoVecs.nodeLowerChoiceForStrain.begin() + generalParams.originNodeCount, 
						nodeInfoVecs.nodeLocZ.begin(),  
						IsLessThanLevel( lowerLevel ), true);
		
	generalParams.numUpperStrainNodes = thrust::count_if(nodeInfoVecs.nodeUpperChoiceForStrain.begin(),nodeInfoVecs.nodeUpperChoiceForStrain.end(), IsEqualToOne( ) );
	generalParams.numLowerStrainNodes = thrust::count_if(nodeInfoVecs.nodeLowerChoiceForStrain.begin(),nodeInfoVecs.nodeLowerChoiceForStrain.end(), IsEqualToOne( ) );
	
	std::cout<<"number of nodes pulled for strain: " << generalParams.numLowerStrainNodes + generalParams.numUpperStrainNodes <<std::endl;

	unsigned numFixed = thrust::count_if(nodeInfoVecs.isNodeFixed.begin(),nodeInfoVecs.isNodeFixed.end(), IsEqualToOne() );
	std::cout<<"number of nodes fixed: " << numFixed <<std::endl;
	zPosTemp.resize(0);

}

void NodeSystemDevice::setTorsionVecs(
	thrust::host_vector<unsigned>& hostTorsionIndexLeft,
	thrust::host_vector<unsigned>& hostTorsionIndexCenter,
	thrust::host_vector<unsigned>& hostTorsionIndexRight,
	thrust::host_vector<double>& hostTorsionAngleZero) {

	unsigned torsion_factor = 500;

	torsionInfoVecs.leftIndex.resize(torsion_factor * generalParams.totalTorsionCount);
	torsionInfoVecs.centerIndex.resize(torsion_factor * generalParams.totalTorsionCount);
	torsionInfoVecs.rightIndex.resize(torsion_factor * generalParams.totalTorsionCount);
	torsionInfoVecs.angleZero.resize(torsion_factor * generalParams.totalTorsionCount);

	thrust::fill(torsionInfoVecs.leftIndex.begin(),torsionInfoVecs.leftIndex.end(),ULONG_MAX);
	thrust::fill(torsionInfoVecs.centerIndex.begin(),torsionInfoVecs.centerIndex.end(),ULONG_MAX);
	thrust::fill(torsionInfoVecs.rightIndex.begin(),torsionInfoVecs.rightIndex.end(),ULONG_MAX);

	//after default value is set, set the real id's
	thrust::copy(hostTorsionIndexLeft.begin(), hostTorsionIndexLeft.end(), torsionInfoVecs.leftIndex.begin());
	thrust::copy(hostTorsionIndexCenter.begin(), hostTorsionIndexCenter.end(), torsionInfoVecs.centerIndex.begin());
	thrust::copy(hostTorsionIndexRight.begin(), hostTorsionIndexRight.end(), torsionInfoVecs.rightIndex.begin());
	
	thrust::transform( 
		thrust::make_zip_iterator(
			thrust::make_tuple( 
				torsionInfoVecs.leftIndex.begin(),
				torsionInfoVecs.centerIndex.begin(),
				torsionInfoVecs.rightIndex.begin())),
		thrust::make_zip_iterator(
			thrust::make_tuple(
				torsionInfoVecs.leftIndex.begin(),
				torsionInfoVecs.centerIndex.begin(),
				torsionInfoVecs.rightIndex.begin())) + generalParams.totalTorsionCount,
			torsionInfoVecs.angleZero.begin(),//save vector
		TorsionAngleFunctor(
			thrust::raw_pointer_cast(nodeInfoVecs.nodeLocX.data()),
			thrust::raw_pointer_cast(nodeInfoVecs.nodeLocY.data()),
			thrust::raw_pointer_cast(nodeInfoVecs.nodeLocZ.data())));
	
	//		std::cout<<" in NSD device values"<<std::endl;
	for (unsigned i = 0; i<generalParams.totalTorsionCount; i++) {
		unsigned n0 = torsionInfoVecs.leftIndex[i];
		unsigned n1 = torsionInfoVecs.centerIndex[i];
		unsigned n2 = torsionInfoVecs.rightIndex[i];
		std::cout<< "angle : "<< n0<< " " << n1<< " " << n2<< " " << torsionInfoVecs.angleZero[i]<<std::endl; 
	}  

	//3x bigger since each spring affects 3 nodes. 
	torsionInfoVecs.forceX.resize(torsion_factor * 3 * generalParams.totalTorsionCount);
	torsionInfoVecs.forceY.resize(torsion_factor * 3 * generalParams.totalTorsionCount);
	torsionInfoVecs.forceZ.resize(torsion_factor * 3 * generalParams.totalTorsionCount);
	torsionInfoVecs.tempForceX.resize(torsion_factor * 3 * generalParams.totalTorsionCount);
	torsionInfoVecs.tempForceY.resize(torsion_factor * 3 * generalParams.totalTorsionCount);
	torsionInfoVecs.tempForceZ.resize(torsion_factor * 3 * generalParams.totalTorsionCount);


	thrust::fill(torsionInfoVecs.forceX.begin(), torsionInfoVecs.forceX.end(), 0.0);
	thrust::fill(torsionInfoVecs.forceY.begin(), torsionInfoVecs.forceY.end(), 0.0);
	thrust::fill(torsionInfoVecs.forceZ.begin(), torsionInfoVecs.forceZ.end(), 0.0);

	torsionInfoVecs.tempTorIndices.resize(torsion_factor * 3 * generalParams.totalTorsionCount);
	torsionInfoVecs.reducedIds.resize(torsion_factor * 3 * generalParams.totalTorsionCount);

 
};

void NodeSystemDevice::setWLCVecs(
	thrust::host_vector<unsigned>& hostWLCSubEdgeLeft,
	thrust::host_vector<unsigned>& hostWLCSubEdgeRight,
	thrust::host_vector<double>& hostWLCSubLenZero ) {

	wlcInfoVecs.globalNeighbors.resize(generalParams.maxNodeCount * generalParams.maxNeighborCount);
	wlcInfoVecs.currentNodeEdgeCountVector.resize(generalParams.maxNodeCount);

	wlcInfoVecs.lengthZero.resize(generalParams.maxNodeCount * generalParams.maxNeighborCount);
	wlcInfoVecs.numOriginalNeighborsNodeVector.resize(generalParams.maxNodeCount);


	thrust::fill(wlcInfoVecs.globalNeighbors.begin(), wlcInfoVecs.globalNeighbors.end(), ULONG_MAX);
	thrust::fill(wlcInfoVecs.currentNodeEdgeCountVector.begin(), wlcInfoVecs.currentNodeEdgeCountVector.end(),0);
	thrust::fill(wlcInfoVecs.lengthZero.begin(), wlcInfoVecs.lengthZero.end(), 0.0);

			   

	nodeInfoVecs.deviceEdgeLeft = hostWLCSubEdgeLeft;
	nodeInfoVecs.deviceEdgeRight = hostWLCSubEdgeRight;
	//scan through hostAdj and put in device.
	for (unsigned id = 0; id < hostWLCSubLenZero.size(); id++) {
		generalParams.totalNumberOfEdges++;
		 unsigned idL = hostWLCSubEdgeLeft[id];
		 unsigned idR = hostWLCSubEdgeRight[id]; 
		 
		//std::cout<< "linking " << idL << " to " <<idR << std::endl;
		
		 double edgeLen = hostWLCSubLenZero[id];	
				//we use the lengthZero vector to identify edges as well.
				//node id is row, column node is connected to row node.
				
		//add edge for left node 		
		unsigned edgeNumL = wlcInfoVecs.currentNodeEdgeCountVector[idL]; //number of edges on (nodeId = row)	is that entry in cECV
		unsigned indexL = idL*generalParams.maxNeighborCount + edgeNumL;
		wlcInfoVecs.lengthZero[indexL] = edgeLen;
		wlcInfoVecs.globalNeighbors[indexL] = idR;
		(wlcInfoVecs.currentNodeEdgeCountVector[idL])++; //right connects to left
  
		//add edge for right node
		unsigned edgeNumR = wlcInfoVecs.currentNodeEdgeCountVector[idR]; //number of edges on (nodeId = row)	is that entry in cECV
		unsigned indexR = idR*generalParams.maxNeighborCount + edgeNumR;
		wlcInfoVecs.lengthZero[indexR] = edgeLen;
		wlcInfoVecs.globalNeighbors[indexR] = idL;
		(wlcInfoVecs.currentNodeEdgeCountVector[idR])++; //left connects to right
		generalParams.currentEdgeCount += 1; 
	} 
	//at this point currentNodeEdgeCountVector holds the number of edges, copy this to 
	thrust::copy(wlcInfoVecs.currentNodeEdgeCountVector.begin(), wlcInfoVecs.currentNodeEdgeCountVector.end(), wlcInfoVecs.numOriginalNeighborsNodeVector.begin());
};

void NodeSystemDevice::setExtras() {
	compressionParams.originalNetworkLength = domainParams.maxZ; //compression along x compressionParams.axis
	compressionParams.originalNetworkWidth = domainParams.maxX;  //strain along z compressionParams.axis.
};


