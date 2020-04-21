#ifndef STORAGE_H_
#define STORAGE_H_

#include <fstream>
#include <memory>
#include <iomanip>

class System;
class System_Builder;

class Storage {

	std::weak_ptr<System> system;
	std::weak_ptr<System_Builder> builder;
	//std::shared_ptr<ExternalForce> grip;
	//std::ofstream output;
	//std::ofstream statesOutput;

	std::ofstream statesOutputStrain;
	std::string bn;
	std::string str_animation;
	std::string str_params; 

	unsigned stepCounter = 0;
	unsigned stepInterval = 10;

	double forceStep = 0.0;	 //increments in which force will be increased.
	double magnitudeForce = 0.0;  //how much force we are currently at.
	int currentAddedEdges = 0;
	int previousAddedEdges = 0;
	unsigned iteration = 0;

public:
	Storage(std::weak_ptr<System> a_system,
		std::weak_ptr<System_Builder> b_system, const std::string& a_filename);

	void save_params(void);

	void print_VTK_file(void);
};

#endif
