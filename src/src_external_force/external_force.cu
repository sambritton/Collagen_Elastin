#include "system_structures.h"
#include "system.h"
#include "external_force.h"
#include "functor_strain.h"
#include "functor_external_force.h"


void external_force(
	NodeInfoVecs& nodeInfoVecs,
	GeneralParams& generalParams,
	ExtensionParams& extensionParams,
	DomainParams& domainParams){
/*
	for (unsigned i = 0; i < nodeInfoVecs.node_is_collagen.size(); i++) {
		std::cout<< "is collagen: " << nodeInfoVecs.node_is_collagen[i] << std::endl;
		std::cout<< "is pulled: " << nodeInfoVecs.node_upper_selection_pull[i] << std::endl;
		std::cout<< "z: " << nodeInfoVecs.node_loc_z[i] << std::endl;
	}*/
	try{
		//try only counting collagen
		extensionParams.averageUpperStrain = (thrust::transform_reduce(
			thrust::make_zip_iterator(
				thrust::make_tuple(
					nodeInfoVecs.node_is_collagen.begin(),
					nodeInfoVecs.node_upper_selection_pull.begin(),
					nodeInfoVecs.node_loc_z.begin())),
			thrust::make_zip_iterator(
				thrust::make_tuple(
					nodeInfoVecs.node_is_collagen.begin(),
					nodeInfoVecs.node_upper_selection_pull.begin(),
					nodeInfoVecs.node_loc_z.begin())) + generalParams.max_node_count,
			functor_strain(generalParams.max_node_count, extensionParams.originalNetworkLength),
				0.0,
			thrust::plus<double>())) / generalParams.numUpperStrainNodes_collagen;
			

			extensionParams.averageLowerStrain = (thrust::transform_reduce(
				thrust::make_zip_iterator(
					thrust::make_tuple(
						nodeInfoVecs.node_is_collagen.begin(),
						nodeInfoVecs.node_lower_selection_pull.begin(),
						nodeInfoVecs.node_loc_z.begin())),
				thrust::make_zip_iterator(
					thrust::make_tuple(
						nodeInfoVecs.node_is_collagen.begin(),
						nodeInfoVecs.node_lower_selection_pull.begin(),
						nodeInfoVecs.node_loc_z.begin())) + generalParams.max_node_count,
				functor_strain(generalParams.max_node_count, extensionParams.originalNetworkLength),
					0.0,
				thrust::plus<double>())) / generalParams.numLowerStrainNodes_collagen;

	if (generalParams.iterationCounter == 1) {
		extensionParams.originAverageUpperStrain = extensionParams.averageUpperStrain;
		extensionParams.originAverageLowerStrain = extensionParams.averageLowerStrain;
	}

	//Apply External Force. 
	//Currently, we apply forces to all nodes withing a range of the average hight of those chosen for the upper/lower section. 
	//We always apply force to the collagen though. maybe remove that feature? currently collagen is allowed a window of 2micron, elastin 0.5
	thrust::counting_iterator<unsigned> indexBeginA(0);

	thrust::for_each(
		thrust::make_zip_iterator(
			thrust::make_tuple(
				indexBeginA,
				nodeInfoVecs.node_loc_z.begin(),
				nodeInfoVecs.is_node_fixed.begin(),
				nodeInfoVecs.node_is_collagen.begin(),
				nodeInfoVecs.node_upper_selection_pull.begin(),
				nodeInfoVecs.node_lower_selection_pull.begin())),
		thrust::make_zip_iterator(
			thrust::make_tuple(
				indexBeginA,
				nodeInfoVecs.node_loc_z.begin(),
				nodeInfoVecs.is_node_fixed.begin(),
				nodeInfoVecs.node_is_collagen.begin(),
				nodeInfoVecs.node_upper_selection_pull.begin(),
				nodeInfoVecs.node_lower_selection_pull.begin())) + generalParams.max_node_count,
		functor_external_force(
			thrust::raw_pointer_cast(nodeInfoVecs.is_node_fixed.data()),
			thrust::raw_pointer_cast(nodeInfoVecs.node_force_x.data()),
			thrust::raw_pointer_cast(nodeInfoVecs.node_force_y.data()),
			thrust::raw_pointer_cast(nodeInfoVecs.node_force_z.data()),

			generalParams.magnitudeForce,
			extensionParams.originalNetworkLength,
			extensionParams.strain_proportion_end_sim,
			extensionParams.averageLowerStrain,
			extensionParams.averageUpperStrain));

		}
		catch(int e){std::cout<<"test"<< e <<std::flush;}
  };
