/*
* ForceDiagramSorage.h
*
*  created on 6/1/2017 SRB
*
*/
//During graph deformation, this file stores position and velocity of nodes at a given time step


#ifndef FORCEDIAGRAMSTORAGE_H_
#define FORCEDIAGRAMSTORAGE_H_

#include <fstream>	     
#include <memory>
#include <iomanip>

class NodeSystemDevice;
class NodeSystemBuilder;

class ForceDiagramStorage {
	
	std::weak_ptr<NodeSystemDevice> system;
	std::weak_ptr<NodeSystemBuilder> builder;
	//std::shared_ptr<ExternalForce> grip;
	std::ofstream output;
	std::ofstream statesOutput;
	
	std::ofstream statesOutputStrain;
	std::string bn;

	unsigned stepCounter = 0;
	unsigned stepInterval = 10;

	double forceStep = 0.0;	 //increments in which force will be increased. 
	double magnitudeForce = 0.0;  //how much force we are currently at.
	int currentAddedEdges = 0;
	int previousAddedEdges = 0;
	unsigned iteration = 0;

public: 
	ForceDiagramStorage(std::weak_ptr<NodeSystemDevice> a_system,
		std::weak_ptr<NodeSystemBuilder> b_system, const std::string& a_filename);

	void updateStrain(void);
	void updateTotalStrain(void);

	void updateStorage(void);
	void print_VTK_File(void);
};

#endif /*FORCEDIAGRAMSTORAGE_H_*/