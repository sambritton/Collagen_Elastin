

#include <iomanip>
#include <string>
#include <memory>
#include <fstream>							 
#include <ctime>
#include <stdio.h>
#include <inttypes.h>
#include <cstddef>

#include "pugixml/include/pugixml.hpp"

#include "NodeSystemDevice.h"

#include "NodeSystemBuilder.h"
#include "ForceDiagramStorage.h"




enum UsageError {
	UE_NOT_ENOUGH_ARGUMENTS = 1,
	UE_UNKNOWN_ARGUMENT,
	UE_INCOMPLETE_SYSTEM,
};

glm::dvec3 fDir;
bool InitExt = true;
bool LogState = false;
bool firstStep = true;
bool newChoiceOfDataField = false; //to input different values for P,C,Bs etc use data fields. 
bool noLinking = false;
//input -data_field=int will give index of vector for corresponding values saved in 
//newChoiceosDataField
int choiceOfDataField;
double Extension;
double initialSize;
unsigned updateEachStep = 10;
double largestForce = 0.0;//used for while loop in solve method. Set in Force input file 
double targetStrain = 1.0;//used to set target for how much network should be compressed using rheometerSim. Set between 0.05 and 1.  

double equilibriumLagTime = 1.0;
double pullPercent = 10.0;
//Set rheometerSim true by command --rheometerSim=(some input fraction for targetStrain).
//we'll push 'down' along the axis defined in builder->axis
//we do not use the distributed force per node, instead we input the total force to be 
//imparted on the network

bool rheomSim = false;



void printUsage() {
	std::cout << "Usage1: <program-name> [params] <scheme-file>" << std::endl;
}

//
std::shared_ptr<NodeSystemDevice> createNodeSystem(const char* schemeFile, std::shared_ptr<NodeSystemBuilder> builder)	{
	pugi::xml_document doc;
	pugi::xml_parse_result parseResult = doc.load_file(schemeFile);

	if (!parseResult) {
		std::cout << "parse error in createNodeSystem: " << parseResult.description() << std::endl;
		return nullptr;
	}
	pugi::xml_node root = doc.child("data");
	pugi::xml_node nodes = root.child("nodes"); 
	pugi::xml_node links = root.child("links");
	pugi::xml_node props = root.child("settings");

	//first, we'll input settings

	if (!(root && nodes && links)) {
		std::cout << "couldn't find nessesary data\n";
		//return false;
	}

	if (auto p = props.child("resistance"))
		builder->defaultResistance = (p.text().as_double());

	if (auto p = props.child("spring-stiffness"))
		builder->defaultSpringStiffness = (p.text().as_double());

	if (auto p = props.child("torsion-stiffness"))
		builder->defaultTorsionSpringStiffness = (p.text().as_double());

	if (auto p = props.child("persistance-length"))
		builder->defaultPersistanceLength = (p.text().as_double());

	if (auto p = props.child("absolute-temperature"))
		builder->defaultTemperature = (p.text().as_double());

	if (auto p = props.child("contour-length-multiplier"))
		builder->defaultContourLengthMultiplier = (p.text().as_double());

	if (auto p = props.child("units-per-extra-node"))
		builder->defaultUnitsPerExtraNode = (p.text().as_double());
	
	if (auto p = props.child("dpd-particles"))
		builder->dpdParticles = (p.text().as_bool());

	if (auto p = props.child("extra-nodes-per-edge"))
		builder->defaultExtraNodesPerEdge = (p.text().as_uint());

	if (auto p = props.child("use-extra-nodes"))
		builder->useExtraNodes = (p.text().as_bool());

	if (auto p = props.child("constant-extra-nodes"))
		builder->useConstantNumberOfExtraNodes = (p.text().as_bool());

	if (auto p = props.child("link-diameter"))
		builder->defaultLinkDiameter = (p.text().as_double());

	if (auto p = props.child("worm-like"))
		builder->wormLikeEnabled = (p.text().as_bool());

	if (auto p = props.child("use-linking"))
		builder->linking = (p.text().as_bool());
	
	if (auto p = props.child("strain-test")) {
		builder->compressionPercent = (p.text().as_double());
		builder->strainSim = true;	
	}
	
	if (auto p = props.child("boundary-condition"))
		builder->boundaryConditionEnabled = (p.text().as_bool());

	if (auto p = props.child("axis"))
		builder->axis = (p.text().as_uint());

	if (rheomSim == true) {
		builder->defaultRheometerSim = true;//sets RheometerSim=true in both NSB and NSI
	}

		//extra data fields for use	units are in nN and micrometers.
	std::vector<double> dataEModulus(5);
	std::vector<double> dataPLength(5);
	std::vector<double> dataBStiffness(5);
	
	//I = 1/4*pi*r^4
	double inertia = 0.25 * 3.1415926535897* std::pow(((builder->defaultLinkDiameter)/2), 4);
	//MPa, .1, 1, 10, 25, 50
	dataEModulus[0] = 100;	 // units: nN/micron^2
	dataEModulus[1] = 1000;
	dataEModulus[2] = 10000;
	dataEModulus[3] = 25000;
	dataEModulus[4] = 50000;

    //Bs = EI;
	dataBStiffness[0] = dataEModulus[0] * inertia; // units: nN*micron^2
	dataBStiffness[1] = dataEModulus[1] * inertia;
	dataBStiffness[2] = dataEModulus[2] * inertia;
	dataBStiffness[3] = dataEModulus[3] * inertia;
	dataBStiffness[4] = dataEModulus[4] * inertia;
	
	//P = EI/KBT
	dataPLength[0] = dataBStiffness[0] / ((builder->defaultTemperature)*(builder->defaultBoltzmannConstant));	
	dataPLength[1] = dataBStiffness[1] / ((builder->defaultTemperature)*(builder->defaultBoltzmannConstant));
	dataPLength[2] = dataBStiffness[2] / ((builder->defaultTemperature)*(builder->defaultBoltzmannConstant));
	dataPLength[3] = dataBStiffness[3] / ((builder->defaultTemperature)*(builder->defaultBoltzmannConstant));
	dataPLength[4] = dataBStiffness[4] / ((builder->defaultTemperature)*(builder->defaultBoltzmannConstant));


	//auto builder = std::make_shared<NodeSystemBuilder>();
	std::cout << "builder ptr address: " << builder << std::endl;
	std::vector<unsigned> originNodes;
	
	double mass;
	double x, y, z, targx, targy, targz, rad; //variables to be used reading in data.
	double defaultMass = nodes.attribute("default-mass").as_double(-1.0);
	builder->defaultMass = defaultMass;
	//debuggin printl
	//std::cout<<"default mass: " << defaultMass << std::endl;
	for (auto node = nodes.child("node"); node; node = node.next_sibling("node")) {
		mass = node.attribute("mass").as_double(defaultMass);
		if (mass < 0.0) {
			std::cout << "parse error: node mass is undefined\n";
			return 0;
		}
		const char* text = node.text().as_string();

		//std::cout<<"mass: " << mass << std::endl;
		if (3 != sscanf(text, "%lf %lf %lf", &x, &y, &z)) {
			std::cout << "parse node error\n";
			return 0;
		}
		int unused = builder->addNode(mass, glm::dvec3(x, y, z), glm::dvec3(0,0,0),1);
		
		//after the node is created, assign an external force
		//default force is zero and default target is dbl_max;
		//these will be updated below if they are set in parameter input
		//std::cout << "origin nodes size: " << originNodes.size() << std::endl;
	}

	unsigned from, to;	
	for (auto link = links.child("link"); link; link = link.next_sibling("link")) {
		if (2 != sscanf(link.text().as_string(""), "%u %u" , &from, &to)) {
			std::cout << "parse link error\n";
			return 0;
		}
		//std::cout << "putting spring between: " << from << ' ' <<to<<  std::endl;
		builder->putSpring(from, to); //adds edges into saved vectors
		
	}

	std::cout << "post springs" << std::endl;


	//forces is used to apply mutliple forces in different directions.
	//each node will have an associated externalforce object.
	//else we set the force to zero.

	pugi::xml_node forces = root.child("forces");
	std::cout << "preforces" << std::endl;
	if (forces) {

		for (auto node : forces.children()) {
			if (std::string("grip") == node.name()) {

				//auto grip = std::make_shared<ExternalForce>();
				//if we apply force
				if (3 == sscanf(node.attribute("value").as_string(), "%lf %lf %lf", &x, &y, &z)) { 
 					//always set force
					for (auto n : node.children("node")) {
						std::cout << "id of node in constraint: " << n.text().as_uint() << std::endl;
						//totalforce is set in builder
						unsigned idTemp = n.text().as_uint();
						//now update hostTargetForces.
						builder->hostTargetForceX[idTemp] = x;
						builder->hostTargetForceY[idTemp] = y;
						builder->hostTargetForceZ[idTemp] = z;
						//builder->updateExtConstraintForce(n.text().as_uint(), glm::dvec3(x, y, z));

						//if we have target and radius, we set them. 
						if (3 == sscanf(node.attribute("target").as_string(), "%lf %lf %lf", &targx, &targy, &targz)) {
							for (auto n : node.children("node")) {
								//builder->updateExtConstraintTarget(n.text().as_uint(), glm::dvec3(targx, targy, targz));
							}

							if (1 == sscanf(node.attribute("radius").as_string(), "%lf", &rad)) {
								for (auto n : node.children("node")) {
									//builder->updateExtConstraintRadius(n.text().as_uint(), double(rad));
								}
							}
						}

						//std::cout << " reading in data: " << targx << " " << targy << " " << targz << " " << rad << std::endl;

					}
					//find the input force with largest norm. 
					double largestForceTemp = glm::sqrt(glm::dot(glm::dvec3(x, y, z), glm::dvec3(x, y, z)));
					if (largestForceTemp > largestForce) {
						largestForce = largestForceTemp;
						
					}

				}				 
				else {
					std::cout << "Input Data Parse Error" << std::endl;
					
				}
			}
		}
	} 
	std::cout << "largest force: " << largestForce << std::endl;
	builder->maxForce = largestForce;



	

	//if we want to use different data fields, we'll rewrite them 
	//This must be done before create is called since create uses the values.
	std::cout << "new choice for data in main.cpp true or false: " << newChoiceOfDataField << std::endl;
	if (newChoiceOfDataField == true) {

		builder->defaultTorsionSpringStiffness = dataBStiffness[choiceOfDataField];
		//builder->defaultPersistanceLength = dataPLength[choiceOfDataField];

	}

	//replace false entries with true if node is fixed.
	pugi::xml_node fixedRoot = root.child("fixed");
	if (fixedRoot) {
		for (auto node = fixedRoot.child("node"); node; node = node.next_sibling("node"))
			builder->fixNode(node.text().as_uint());
	}


	std::cout << "post fixed" << std::endl;

	//last, set and add non resistance and non external constraints.
	auto model = builder->create();  

std::cout << "model built" << "\n";

return model;
}


std::string generateOutputFileName(std::string inputFileName)
{
	time_t now;																   
	const int MAX_DATE = 64;
	char theDate[MAX_DATE];
	theDate[0] = '\0';

	now = time(nullptr);

	if (now != -1) {
		strftime(theDate, MAX_DATE, "_%Y.%m.%d_%H-%M-%S", gmtime(&now));
		return inputFileName + theDate;
	}
	return "";
}



int diagram_main(int argc, char** argv) {

	time_t t0,t1;
	t0 = time(0);
	std::cout << "inside diagram_main " << asctime(localtime(&t0)) << std::endl;
	if (argc < 3) {
		std::cout << "not enough parameters" << std::endl;
		printUsage();
		return UE_NOT_ENOUGH_ARGUMENTS;
	}

	std::cout << "diagram_main started" << std::endl;

	double forceStep = 0.0;
	double epsilon = 0.0;
	double timeStep = 0.01;
	

	bool forceStepEncountered = false;
	bool epsilonEncountered = false;
	bool dtEncountered = false;

	for (int i = 0; i < argc-1; ++i) {


		std::string arg = argv[i];
		unsigned pos = arg.find('=');

		if (pos == std::string::npos)
		{

			std::cout << "UNKNOWN ARG " << arg << std::endl;
			printUsage();
			return UE_UNKNOWN_ARGUMENT;
		}
		std::string key = arg.substr(0, pos);
		std::string val = arg.substr(pos + 1);

		if (key == "-df") {
			forceStep = std::atof(val.c_str());
			forceStepEncountered = true;
			continue;
		}
		if (key == "-eps") {
			epsilon = std::atof(val.c_str());
			epsilonEncountered = true;
			continue;
		}
		if (key == "-dt") {
			timeStep = std::atof(val.c_str());
			dtEncountered = true;
			continue;
		}									  


		if (key == "--lagTime") {
			equilibriumLagTime = std::atof(val.c_str());
			continue;
		}


		if (key == "--pullPercent") {
			pullPercent = std::atof(val.c_str());
			continue;
		}
		
		if (key == "--dataField") {
			choiceOfDataField = std::atof(val.c_str());
			newChoiceOfDataField = true;
			continue;
		}
		
		printUsage();
		return UE_UNKNOWN_ARGUMENT;
	}

	if (!(forceStepEncountered && epsilonEncountered && dtEncountered)) {
		printUsage();
		std::cout << "not enough parameters: df, dt and eps requred" << std::endl;
		return UE_NOT_ENOUGH_ARGUMENTS;
	};
	



	auto builder = std::make_shared<NodeSystemBuilder>(epsilon, timeStep, forceStep, targetStrain);
	builder->lagTime = equilibriumLagTime;
	builder->pullPercent = pullPercent;
	//sets all parameters and edges on device side
	auto system = createNodeSystem(argv[argc-1], builder);

	std::cout<<"pull percent: "<< system->generalParams.pullPercent <<std::endl;
	std::cout<<"lagTime : "<< system->generalParams.lagTime <<std::endl;
	


	//fDir = glm::normalize(system->getGrip()->getForce());
	auto outputFileName = generateOutputFileName(argv[argc-1]);
	
	//once the system is set, we'll store the initial values via the ptr system.
	//ForceDiagramStocrrage storage( system, outputFileName);
	auto storage = std::make_shared<ForceDiagramStorage>(system, builder, outputFileName);
	//pass 
	std::cout << "assigning fdiagram in main" << std::endl;
	system->assignForceDiagramStorage(storage);
	std::cout << "post fdiagram in main" << std::endl;

	std::cout << "solving system in main" << std::endl;
	system->solveSystemDevice();


	t1 = time(0);  //current time at the end of solving the system.
	int total,hours,min,sec;
	total = difftime(t1,t0);
	hours = total / 3600;
	min = (total % 3600) / 60;
	sec = (total % 3600) % 60;
	// текущее время расчета
	std::cout << "Total time hh: " << hours << " mm:" << min << " ss:" << sec <<"\n";
	return 0;
}


int main(int argc, char** argv)
{
	std::cout << argc;
	std::cout << std::endl;

	if (argc < 2) {
		printUsage();
		return UE_NOT_ENOUGH_ARGUMENTS;
	}

	// choose calculation mode
	std::string mode = argv[1];

	std::cout << mode << std::endl;
	if (mode == "diagram") {
		return diagram_main(argc - 2, argv + 2);
	}
	return 0;
}
