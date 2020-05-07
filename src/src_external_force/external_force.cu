#include "system_structures.h"
#include "system.h"
#include "external_force.h"
#include "functor_strain.h"
#include "functor_external_force.h"
#include "functor_external_pull.h"


void external_force(
	NodeInfoVecs& nodeInfoVecs,
	GeneralParams& generalParams,
	ExtensionParams& extensionParams,
	DomainParams& domainParams){

	if ((generalParams.numUpperStrainNodes_collagen > 0) && (generalParams.numLowerStrainNodes_collagen > 0)){
		//try only counting collagen
		extensionParams.averageUpperStrain = (thrust::transform_reduce(
			thrust::make_zip_iterator(
				thrust::make_tuple(
					nodeInfoVecs.node_is_collagen.begin(),
					nodeInfoVecs.node_upper_selection_pull.begin(),
					nodeInfoVecs.node_loc_z.begin(),
					nodeInfoVecs.node_loc_x.begin())),
			thrust::make_zip_iterator(
				thrust::make_tuple(
					nodeInfoVecs.node_is_collagen.begin(),
					nodeInfoVecs.node_upper_selection_pull.begin(),
					nodeInfoVecs.node_loc_z.begin(),
					nodeInfoVecs.node_loc_x.begin())) + generalParams.max_node_count,
			functor_strain(extensionParams.axis, extensionParams.originalNetworkLength),
				0.0,
			thrust::plus<double>())) / generalParams.numUpperStrainNodes_collagen;
			

		extensionParams.averageLowerStrain = (thrust::transform_reduce(
			thrust::make_zip_iterator(
				thrust::make_tuple(
					nodeInfoVecs.node_is_collagen.begin(),
					nodeInfoVecs.node_lower_selection_pull.begin(),
					nodeInfoVecs.node_loc_z.begin(),
					nodeInfoVecs.node_loc_x.begin())),
			thrust::make_zip_iterator(
				thrust::make_tuple(
					nodeInfoVecs.node_is_collagen.begin(),
					nodeInfoVecs.node_lower_selection_pull.begin(),
					nodeInfoVecs.node_loc_z.begin(),
					nodeInfoVecs.node_loc_x.begin())) + generalParams.max_node_count,
			functor_strain(extensionParams.axis, extensionParams.originalNetworkLength),
				0.0,
			thrust::plus<double>())) / generalParams.numLowerStrainNodes_collagen;

			thrust::counting_iterator<unsigned> indexBeginA(0);

			thrust::for_each(
				thrust::make_zip_iterator(
					thrust::make_tuple(
						indexBeginA,
						nodeInfoVecs.node_loc_z.begin(),
						nodeInfoVecs.node_loc_x.begin(),
						nodeInfoVecs.is_node_fixed.begin(),
						nodeInfoVecs.node_is_collagen.begin(),
						nodeInfoVecs.node_upper_selection_pull.begin(),
						nodeInfoVecs.node_lower_selection_pull.begin())),
				thrust::make_zip_iterator(
					thrust::make_tuple(
						indexBeginA,
						nodeInfoVecs.node_loc_z.begin(),
						nodeInfoVecs.node_loc_x.begin(),
						nodeInfoVecs.is_node_fixed.begin(),
						nodeInfoVecs.node_is_collagen.begin(),
						nodeInfoVecs.node_upper_selection_pull.begin(),
						nodeInfoVecs.node_lower_selection_pull.begin())) + generalParams.max_node_count,
				functor_external_pull(
					thrust::raw_pointer_cast(nodeInfoVecs.is_node_fixed.data()),
					thrust::raw_pointer_cast(nodeInfoVecs.node_loc_x.data()),
					thrust::raw_pointer_cast(nodeInfoVecs.node_loc_z.data()),
					generalParams.pull_ammount,
					extensionParams.axis,
					extensionParams.originalNetworkLength,
					extensionParams.strain_proportion_end_sim,
					extensionParams.averageLowerStrain,
					extensionParams.averageUpperStrain));
			
	}
	
	//also pull elastin for now. 
	if ((generalParams.numUpperStrainNodes_elastin >0) && ( generalParams.numLowerStrainNodes_elastin > 0)){
		
		extensionParams.averageUpperStrain = (thrust::transform_reduce(
			thrust::make_zip_iterator(
				thrust::make_tuple(
					nodeInfoVecs.node_is_elastin.begin(),
					nodeInfoVecs.node_upper_selection_pull.begin(),
					nodeInfoVecs.node_loc_z.begin(),
					nodeInfoVecs.node_loc_x.begin())),
			thrust::make_zip_iterator(
				thrust::make_tuple(
					nodeInfoVecs.node_is_elastin.begin(),
					nodeInfoVecs.node_upper_selection_pull.begin(),
					nodeInfoVecs.node_loc_z.begin(),
					nodeInfoVecs.node_loc_x.begin())) + generalParams.max_node_count,
			functor_strain(extensionParams.axis, extensionParams.originalNetworkLength),
				0.0,
			thrust::plus<double>())) / generalParams.numUpperStrainNodes_elastin;
			
		extensionParams.averageLowerStrain = (thrust::transform_reduce(
			thrust::make_zip_iterator(
				thrust::make_tuple(
					nodeInfoVecs.node_is_elastin.begin(),
					nodeInfoVecs.node_lower_selection_pull.begin(),
					nodeInfoVecs.node_loc_z.begin(),
					nodeInfoVecs.node_loc_x.begin())),
			thrust::make_zip_iterator(
				thrust::make_tuple(
					nodeInfoVecs.node_is_elastin.begin(),
					nodeInfoVecs.node_lower_selection_pull.begin(),
					nodeInfoVecs.node_loc_z.begin(),
					nodeInfoVecs.node_loc_x.begin())) + generalParams.max_node_count,
			functor_strain(extensionParams.axis, extensionParams.originalNetworkLength),
				0.0,
			thrust::plus<double>())) / generalParams.numLowerStrainNodes_elastin;
			thrust::counting_iterator<unsigned> indexBeginA(0);

			thrust::for_each(
				thrust::make_zip_iterator(
					thrust::make_tuple(
						indexBeginA,
						nodeInfoVecs.node_loc_z.begin(),
						nodeInfoVecs.node_loc_x.begin(),
						nodeInfoVecs.is_node_fixed.begin(),
						nodeInfoVecs.node_is_elastin.begin(),
						nodeInfoVecs.node_upper_selection_pull.begin(),
						nodeInfoVecs.node_lower_selection_pull.begin())),
				thrust::make_zip_iterator(
					thrust::make_tuple(
						indexBeginA,
						nodeInfoVecs.node_loc_z.begin(),
						nodeInfoVecs.node_loc_x.begin(),
						nodeInfoVecs.is_node_fixed.begin(),
						nodeInfoVecs.node_is_elastin.begin(),
						nodeInfoVecs.node_upper_selection_pull.begin(),
						nodeInfoVecs.node_lower_selection_pull.begin())) + generalParams.max_node_count,
				functor_external_pull(
					thrust::raw_pointer_cast(nodeInfoVecs.is_node_fixed.data()),
					thrust::raw_pointer_cast(nodeInfoVecs.node_loc_x.data()),
					thrust::raw_pointer_cast(nodeInfoVecs.node_loc_z.data()),
					generalParams.pull_ammount,
					extensionParams.axis,
					extensionParams.originalNetworkLength,
					extensionParams.strain_proportion_end_sim,
					extensionParams.averageLowerStrain,
					extensionParams.averageUpperStrain));
	}
	if (generalParams.iterationCounter < 2) {
		extensionParams.originAverageUpperStrain = extensionParams.averageUpperStrain;
		extensionParams.originAverageLowerStrain = extensionParams.averageLowerStrain;
	}

	//Apply External Force. 
	//Currently, we apply forces to all nodes withing a range of the average hight of those chosen for the upper/lower section. 
	//We always apply force to the collagen though. maybe remove that feature? currently collagen is allowed a window of 2micron, elastin 0.5


  };
