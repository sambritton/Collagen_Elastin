

#include <iomanip>
#include <string>
#include <memory>
#include <fstream>
#include <ctime>
#include <stdio.h>
#include <inttypes.h>
#include <cstddef>

#include "pugixml.hpp"

#include "system.h"

#include "system_builder.h"
#include "storage.h"




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

int choiceOfDataField;
double Extension;
double initialSize;
unsigned updateEachStep = 10;
double largestForce = 0.0;//used for while loop in solve method. Set in Force input file
double targetStrain = 1.0;//used to set target for how much network should be compressed using rheometerSim. Set between 0.05 and 1.

double equilibriumLagTime = 1.0;
double pullPercent = 0.1;

void printUsage() {
	std::cout << "Usage1: [params] <scheme-file>" << std::endl;
}

//
std::shared_ptr<System> create_system(const char* schemeFile, std::shared_ptr<System_Builder> builder)	{
	pugi::xml_document doc;
	pugi::xml_parse_result parseResult = doc.load_file(schemeFile);

	if (!parseResult) {
		std::cout << "parse error in create_system: " << parseResult.description() << std::endl;
		return nullptr;
	}
	pugi::xml_node root = doc.child("data");
	pugi::xml_node nodes = root.child("nodes");
	pugi::xml_node edges = root.child("edges");
	pugi::xml_node props = root.child("settings");

	//first, we'll input settings

	if (!(root && nodes && edges)) {
		std::cout << "couldn't find nessesary data\n";
	}

	if (auto p = props.child("viscosity-elastin"))
		builder->default_viscosity_elastin = (p.text().as_double());

	if (auto p = props.child("viscosity-collagen"))
		builder->default_viscosity_collagen = (p.text().as_double());

	if (auto p = props.child("collagen-spring-constant"))
		builder->default_collagen_spring_constant = (p.text().as_double());

	if (auto p = props.child("bend-stiffness-collagen"))
		builder->default_bend_stiffness_collagen = (p.text().as_double());
	if (auto p = props.child("bend-stiffness-elastin"))
		builder->default_bend_stiffness_elastin = (p.text().as_double());

	if (auto p = props.child("persistance-len-monomer"))
		builder->default_persistence_len_monomer = (p.text().as_double());

	if (auto p = props.child("temperature"))
		builder->default_temperature = (p.text().as_double());

	if (auto p = props.child("contour-length-multiplier"))
		builder->default_CLM = (p.text().as_double());

	if (auto p = props.child("units-per-extra-node"))
		builder->defaultUnitsPerExtraNode = (p.text().as_double());

	if (auto p = props.child("extra-nodes-per-edge"))
		builder->default_extra_nodes_per_edge = (p.text().as_uint());
		std::cout<<"extra nodes: "<< builder->default_extra_nodes_per_edge << std::endl;

	if (auto p = props.child("use-extra-nodes"))
		builder->use_extra_nodes = (p.text().as_bool());

	if (auto p = props.child("constant-number-extra-nodes"))
		builder->use_constant_number_extra_nodes = (p.text().as_bool());

	if (auto p = props.child("collagen-diameter"))
		builder->default_collagen_diameter = (p.text().as_double());
	if (auto p = props.child("elastin-diameter"))
		builder->default_elastin_diameter = (p.text().as_double());

	if (auto p = props.child("use-linking"))
		builder->default_linking = (p.text().as_bool());

	if (auto p = props.child("strain-test")) {
		builder->default_strain_sim = true;
	}

	if (auto p = props.child("axis")){
		builder->axis = (p.text().as_uint());
		std::cout<<"pulling axis: " << builder->axis << std::endl;
	}

	std::vector<unsigned> originNodes;

	double x, y, z, rad; //variables to be used reading in data.

	//first add collagen points
	for (auto node = nodes.child("node_collagen"); node; node = node.next_sibling("node_collagen")) {

		const char* text = node.text().as_string();

		if (3 != sscanf(text, "%lf %lf %lf", &x, &y, &z)) {
			std::cout << "parse node error\n";
			return 0;
		}
		builder->add_collagen_node(glm::dvec3(x, y, z));
		//std::cout<<"setting collagen node" <<std::endl;
	}
	//now add elastin points
	for (auto node = nodes.child("node_elastin"); node; node = node.next_sibling("node_elastin")) {

		const char* text = node.text().as_string();

		if (3 != sscanf(text, "%lf %lf %lf", &x, &y, &z)) {
			std::cout << "parse node error\n";
			return 0;
		}
		builder->add_elastin_node(glm::dvec3(x, y, z));
		
		//std::cout<<"setting elastin node" <<std::endl;
	}

	unsigned from, to;
	for (auto edge = edges.child("edge_collagen"); edge; edge = edge.next_sibling("edge_collagen")) {
		if (2 != sscanf(edge.text().as_string(""), "%u %u" , &from, &to)) {
			std::cout << "parse edge error\n";
			return 0;
		}
		//std::cout << "putting spring between: " << from << ' ' <<to<<  std::endl;
		builder->put_spring(from, to); //adds edges into saved vectors
		
		//std::cout<<"setting collagen spring" <<std::endl;

	}
	for (auto edge = edges.child("edge_elastin"); edge; edge = edge.next_sibling("edge_elastin")) {
		if (2 != sscanf(edge.text().as_string(""), "%u %u" , &from, &to)) {
			std::cout << "parse edge error\n";
			return 0;
		}
		//std::cout << "putting spring between: " << from << ' ' <<to<<  std::endl;
		builder->put_spring(from, to); //adds edges into saved vectors
		
		//std::cout<<"setting elastin spring" <<std::endl;

	}

	std::cout << "post springs" << std::endl;

	//replace false entries with true if node is fixed.
	pugi::xml_node fixedRoot = root.child("fixed");
	if (fixedRoot) {
		for (auto node = fixedRoot.child("node"); node; node = node.next_sibling("node"))
			builder->fix_node(node.text().as_uint());
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






int main(int argc, char** argv)
{
	std::cout << argc << std::endl;
	std::cout<< argv[0] << std::endl;
	std::cout<< argv[1] << std::endl;
	if (argc < 1) {
		printUsage();
		return UE_NOT_ENOUGH_ARGUMENTS;
	}


	time_t t0,t1;
	t0 = time(0);

	double forceStep = 0.0;
	double epsilon = 1.0;
	double timeStep = 0.001;
	double pull_ammount=0.01;


	bool forceStepEncountered = false;
	bool epsilonEncountered = false;
	bool dtEncountered = false;

	for (int i = 0; i < argc; ++i) {
		std::string arg = argv[i];
		unsigned pos = arg.find('=');

		std::string key = arg.substr(0, pos);
		std::string val = arg.substr(pos + 1);

		if (key == "-df") {
			forceStep = std::atof(val.c_str());
			forceStepEncountered = true;
			std::cout<<"force: "<< forceStep << std::endl;
			continue;
		}
		if (key == "-dpull") {
			pull_ammount = std::atof(val.c_str());
			std::cout<<"pull: "<< pull_ammount << std::endl;
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
	}


	std::cout<<"pre builder"<<std::endl;

	auto builder = std::make_shared<System_Builder>(epsilon, timeStep, forceStep, targetStrain);
	builder->default_pull_percent = pullPercent;
	builder->pull_ammount = pull_ammount;
	//sets all parameters and edges on device side
	auto system = create_system(argv[argc-1], builder);

	std::cout<<"pull percent: "<< system->generalParams.pull_percent <<std::endl;
	
	auto outputFileName = generateOutputFileName(argv[argc-1]);

	//once the system is set, we'll store the initial values via the ptr system.
	auto storage = std::make_shared<Storage>(system, builder, outputFileName);
	//pass
	std::cout << "assigning fdiagram in main" << std::endl;
	system->assign_storage(storage);
	std::cout << "post fdiagram in main" << std::endl;

	std::cout << "solving system in main" << std::endl;
	system->solve_system();

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
