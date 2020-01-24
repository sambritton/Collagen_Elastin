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

    thrust::counting_iterator<unsigned> index_begin_upper(0);
		thrust::counting_iterator<unsigned> index_begin_lower(0);

		extensionParams.averageUpperStrain = (thrust::transform_reduce(
			thrust::make_zip_iterator(
				thrust::make_tuple(
					index_begin_upper,
					nodeInfoVecs.node_upper_selection_pull.begin(),
					nodeInfoVecs.node_loc_z.begin())),
			thrust::make_zip_iterator(
				thrust::make_tuple(
					index_begin_upper,
					nodeInfoVecs.node_upper_selection_pull.begin(),
					nodeInfoVecs.node_loc_z.begin())) + generalParams.max_node_count,
			functor_strain(generalParams.max_node_count, extensionParams.originalNetworkLength),
				0.0,
			thrust::plus<double>())) / generalParams.numUpperStrainNodes;

			extensionParams.averageLowerStrain = (thrust::transform_reduce(
				thrust::make_zip_iterator(
					thrust::make_tuple(
						index_begin_lower,
						nodeInfoVecs.node_lower_selection_pull.begin(),
						nodeInfoVecs.node_loc_z.begin())),
				thrust::make_zip_iterator(
					thrust::make_tuple(
						index_begin_lower,
						nodeInfoVecs.node_lower_selection_pull.begin(),
						nodeInfoVecs.node_loc_z.begin())) + generalParams.max_node_count,
				functor_strain(generalParams.max_node_count, extensionParams.originalNetworkLength),
					0.0,
				thrust::plus<double>())) / generalParams.numLowerStrainNodes;

	if (generalParams.iterationCounter == 1) {
		extensionParams.originAverageUpperStrain = extensionParams.averageUpperStrain;
		extensionParams.originAverageLowerStrain = extensionParams.averageLowerStrain;
	}


 
	//Apply External Force
	thrust::counting_iterator<unsigned> indexBeginA(0);

	thrust::for_each(
		thrust::make_zip_iterator(
			thrust::make_tuple(
				indexBeginA,
				nodeInfoVecs.node_loc_z.begin(),
				nodeInfoVecs.is_node_fixed.begin(),
				nodeInfoVecs.node_upper_selection_pull.begin(),
				nodeInfoVecs.node_lower_selection_pull.begin())),
		thrust::make_zip_iterator(
			thrust::make_tuple(
				indexBeginA,
				nodeInfoVecs.node_loc_z.begin(),
				nodeInfoVecs.is_node_fixed.begin(),
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


  };
