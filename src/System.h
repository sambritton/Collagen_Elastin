

#ifndef NODESYSTEMIMPLDEVICE_H_
#define NODESYSTEMIMPLDEVICE_H_

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


 
//Data Structure for node location. velocity and force
struct NodeInfoVecs {
	thrust::device_vector<bool> nodeUpperChoiceForStrain;//holds id's of top 10% for strain test
	thrust::device_vector<bool> nodeLowerChoiceForStrain;//holds id's of top 10% for strain test
	 
	thrust::device_vector<unsigned> springDivisionCount;
	//holds sum of forces for node on given time_step
	thrust::host_vector<unsigned> idEdgesMadeHost;
	
	thrust::device_vector<unsigned> deviceEdgeLeft;
	thrust::device_vector<unsigned> deviceEdgeRight;
	
	thrust::device_vector<unsigned> hostOriginalEdgeLeft;
	thrust::device_vector<unsigned> hostOriginalEdgeRight;
	
	thrust::device_vector<unsigned> idEdgesMadeTemp;
 
	thrust::device_vector<double> sumForcesOnNode;
	 
	thrust::device_vector<double> discretizedEdgeStrain; //counts strain of edge
	//thrust::device_vector<double> addedEdgeStrain;
	thrust::device_vector<double> discretizedEdgeAlignment;

//true if fixed, false if movable. Default false.
	thrust::device_vector<bool> isNodeFixed;


	// X,Y,Z, location, velocity and force of all nodes
	thrust::device_vector<double> prevNodeLocX;
	thrust::device_vector<double> prevNodeLocY;
	thrust::device_vector<double> prevNodeLocZ;	
	thrust::device_vector<double> prevNodeVelX;
	thrust::device_vector<double> prevNodeVelY;
	thrust::device_vector<double> prevNodeVelZ;	
	thrust::device_vector<double> prevNodeForceX;
	thrust::device_vector<double> prevNodeForceY;
	thrust::device_vector<double> prevNodeForceZ;
 
//now holding node and dpd particles. DPD has indices >=maxNodeCount. Total vector length is maxNodeCount + DPDParticleCount
	thrust::device_vector<double> nodeLocX;
	thrust::device_vector<double> nodeLocY;
	thrust::device_vector<double> nodeLocZ;
	thrust::device_vector<double> nodeVelX;
	thrust::device_vector<double> nodeVelY;
	thrust::device_vector<double> nodeVelZ;

	
	thrust::device_vector<double> nodeVelocity;
	
//holds forces to advance position and velocity
	thrust::device_vector<double> nodeForceX;
	thrust::device_vector<double> nodeForceY;
	thrust::device_vector<double> nodeForceZ;
	

};


struct DPDParticleVariables {
	unsigned particleCount = 0;//stores count
	double alpha = .75;//strength of conservative force 
	double R_c = 0.5;//cutoff length for interaction
	double particleMass = 1;
	double gridSpacing = 0.75; //same as R_c/2?
	

};

//struct used for linking of nodes in network 
struct AuxVecs {
	// bucket key means which bucket ID does a certain point fit into
	thrust::device_vector<unsigned> bucketKeys;	//bucket id
	// bucket value means global rank of a certain point
	thrust::device_vector<unsigned> bucketValues;//node id
	// bucket key expanded means what are the bucket IDs are the neighbors of a certain point
	thrust::device_vector<unsigned> bucketKeysExpanded;
	// bucket value expanded means each point ( represented by its global rank) will have multiple copies
	thrust::device_vector<unsigned> bucketValuesIncludingNeighbor;
	
	// begin position of a keys in bucketKeysExpanded and bucketValuesIncludingNeighbor
	//entry keyBegin[bucketKey] returns start of indices to link
	thrust::device_vector<unsigned> keyBegin;
	// end position of a keys in bucketKeysExpanded and bucketValuesIncludingNeighbor
	thrust::device_vector<unsigned> keyEnd;
	
	unsigned endIndexBucketKeys; 
};



struct DomainParams {
	double minX;
	double maxX;
	double minY;
	double maxY;
	double minZ;
	double maxZ;
	double originMinX;
	double originMaxX;
	double originMinY;
	double originMaxY;
	double originMinZ;
	double originMaxZ;
	double dpdMinX;
	double dpdMaxX;
	double dpdMinY;
	double dpdMaxY;
	double dpdMinZ;
	double dpdMaxZ;
	double gridSpacing = 0.5;
	unsigned XBucketCount;
	unsigned YBucketCount;
	unsigned ZBucketCount;
	unsigned totalBucketCount;
};

 
//Data for edge node id's
struct WLCInfoVecs {
	double percentOriginalEdgesUnderStrain1 = 0.0;
	double percentOriginalEdgesExtended = 0.0;
	double percentOriginalEdgesCompressed = 0.0;
	double percentAddedEdgesUnderStrain1 = 0.0;
	double percentAddedEdgesExtended = 0.0;
	double percentAddedEdgesCompressed = 0.0;
	double averageStrainAddedEdges = 0.0;
	double averageStrainOriginalEdges = 0.0;
 
    thrust::device_vector<unsigned> numOriginalNeighborsNodeVector;//holds how many original edges a node was connected to. 
	thrust::device_vector<unsigned> globalNeighbors;
	thrust::device_vector<unsigned> currentNodeEdgeCountVector;
	thrust::device_vector<unsigned> currentBindCountPerOriginalEdgeVector;

	thrust::device_vector<double> lengthZero;

};

struct TorsionInfoVecs {
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

};

struct CompressionParams {

	bool compressionTurnedStrain = false;//when true, df is incremented and shear strain is applied
	double totalAppliedForce = 0.0;
	
	double averageStrainPositions;
	//strain parameters
	bool rheometerSim = false;
	double compressionPercent = 1.0;//determines percent of compression before strain is applied. Default at 1 means it does nothing unless user set

	//parameters for springs set by constructor

	
	//reset and used to calculate current equilibrium state.
	unsigned currentNumberPulled = 0;
	unsigned nextNumberPulled = 0;
	double targetStrain = 0.1;
	double currentStrain = 1;
	
	unsigned axis = 0; //default axis for boundary condition is zero. 
	double currentNetworkLength;
	double originalNetworkLength = 0;
	double originalNetworkWidth = 0;
	double currentStrainNetworkWidthMininum = 0.0;
	
	double strainProportion = 4.0;//if applying strain, how far to run L_new/L_original
	
	double averageLowerStrain = 0.0;
	double averageUpperStrain = 0.0;
	double originAverageLowerStrain = 0.0;
	double originAverageUpperStrain = 0.0;
	
};

struct GeneralParams{
	//general computation
	bool runSim = true; //default true to begin sim. Sim ends when runSim == false
	bool useDPDParticles = false;
	bool strainSim = false;

	double lagTime = 1.0;//set in main.cpp usin
	double pullPercent = 10.0; // set in main.cpp
 
	unsigned maxNeighborCount = 50;
	unsigned maxNodeCount;//after discretize
	unsigned originNodeCount;//pre discretize

	unsigned originLinkCount;//constant unsubdivided count of edges
	unsigned originEdgeCount = 0; //total links set at beginning. Constants
	unsigned currentEdgeCount = 0;//total non constant if links are made

	unsigned totalTorsionCount;//total bending springs
	unsigned subNodeCount = 0;//maximal subnode division for longest edge
	
	double kB, CLM, temperature, torsionStiffness, viscousDamp;
	double nodeMass = 1;
	double persistenceLengthMon;

	double fiberDiameter = 0.1;
	
	//parameters for advancing timestep and determining equilibrium
	double df, dtMax, dtTemp, epsilon, maxForce;
	double magnitudeForce = 0.0;
	double currentTime = 0.0; 

	double currentLength=5;

	//total equilibrium iters and linking determiner
	unsigned iterationCounter = 0;
	bool linking = true;
	bool boundaryCondition = false;

	unsigned numUpperStrainNodes = 0;
	unsigned numLowerStrainNodes = 0;
	
	double proportionOfStrainToCompressed;//one if all strain, 
	unsigned totalNumberOfEdges = 0;//updated after linking

	unsigned maxLinksPerIteration = 5;
};


class ForceDiagramStorage;

class NodeSystemDevice {
public:
	DomainParams domainParams;
	NodeInfoVecs nodeInfoVecs;
	AuxVecs auxVecs;
	WLCInfoVecs wlcInfoVecs;
	TorsionInfoVecs torsionInfoVecs;
	CompressionParams compressionParams;
	GeneralParams generalParams;
	DPDParticleVariables dpdParticleVariables;
	
	std::shared_ptr<ForceDiagramStorage> storage;
	
public:

	NodeSystemDevice();

	double solveForcesOnDevice();


	void assignForceDiagramStorage(std::shared_ptr<ForceDiagramStorage> _storage);

	void initializeSystem(
		thrust::host_vector<bool>& _hostIsNodeFixed,
		thrust::host_vector<double>& _hostPosX,
		thrust::host_vector<double>& _hostPosY,
		thrust::host_vector<double>& _hostPosZ,
		thrust::host_vector<unsigned>& hostWLCEdgeLeft,
		thrust::host_vector<unsigned>& hostWLCEdgeRight,
		thrust::host_vector<double>& hostWLCLenZero,
		thrust::host_vector<unsigned>& hostWLCEdgeSubLeft,
		thrust::host_vector<unsigned>& hostWLCEdgeSubRight,
		thrust::host_vector<double>& hostWLCLenSubZero,
		thrust::host_vector<unsigned>& hostSpringDivisionCount,
		thrust::host_vector<unsigned>& _hostTorsionIndexLeft,
		thrust::host_vector<unsigned>& _hostTorsionIndexCenter,
		thrust::host_vector<unsigned>& _hostTorsionIndexRight,
		thrust::host_vector<double>& _hostTorsionAngleZero);


	//use from cpu side to begin solving system.
	void solveSystemDevice();
	

	void determineBounds();//used for strain.

	void setNodeVecs(
		thrust::host_vector<bool>& hostIsNodeFixed,
		thrust::host_vector<double>& hostPosX, 
		thrust::host_vector<double>& hostPosY, 
		thrust::host_vector<double>& hostPosZ,
		thrust::host_vector<unsigned>& hostSpringDivisionCount);
		
	void setTorsionVecs(	
		thrust::host_vector<unsigned>& hostTorsionIndexLeft,	
		thrust::host_vector<unsigned>& hostTorsionIndexCenter,
		thrust::host_vector<unsigned>& hostTorsionIndexRight,
		thrust::host_vector<double>& hostTorsionAngleZero);
	
	void setWLCVecs(	
		thrust::host_vector<unsigned>& hostWLCSubEdgeLeft,
		thrust::host_vector<unsigned>& hostWLCSubEdgeRight,
		thrust::host_vector<double>& hostWLCSubLenZero );
		
	void setExtras();
};


#endif /*NODESYSTEMIMPLDEVICE_H_*/