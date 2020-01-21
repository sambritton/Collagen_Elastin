#ifndef STORAGE_H_
#define STORAGE_H_

#include <fstream>
#include <memory>
#include <iomanip>

class system;
class system_builder;

class Storage {

	std::weak_ptr<system> system;
	std::weak_ptr<system_builder> builder;
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
	Storage(std::weak_ptr<system> a_system,
		std::weak_ptr<system_builder> b_system, const std::string& a_filename);

	void updateStrain(void);
	void updateTotalStrain(void);

	void updateStorage(void);
	void print_VTK_File(void);
};

#endif
