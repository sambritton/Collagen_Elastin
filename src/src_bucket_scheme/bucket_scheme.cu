//This file sets the grid for network self interaction

#include "system_structures.h"
#include "bucket_scheme.h"
#include "system.h"

#include "functor_neighbor.h"
#include "functor_bucket_indexer.h"
#include "function_extend.h"

//take domain and discretize into square buckets of size gridspace
void init_dim_general(
	NodeInfoVecs& nodeInfoVecs,
	DomainParams& domainParams,
	AuxVecs& auxVecs,
	GeneralParams& generalParams) {

	double minXTemp = (*(thrust::min_element(nodeInfoVecs.node_loc_x.begin(), nodeInfoVecs.node_loc_x.end())));
	double maxXTemp = (*(thrust::max_element(nodeInfoVecs.node_loc_x.begin(), nodeInfoVecs.node_loc_x.end())));
	double minYTemp = (*(thrust::min_element(nodeInfoVecs.node_loc_y.begin(), nodeInfoVecs.node_loc_y.end())));
	double maxYTemp = (*(thrust::max_element(nodeInfoVecs.node_loc_y.begin(), nodeInfoVecs.node_loc_y.end())));
	double minZTemp = (*(thrust::min_element(nodeInfoVecs.node_loc_z.begin(), nodeInfoVecs.node_loc_z.end())));
	double maxZTemp = (*(thrust::max_element(nodeInfoVecs.node_loc_z.begin(), nodeInfoVecs.node_loc_z.end())));


	double space = 0.0;
	domainParams.min_x = minXTemp - space;
	domainParams.max_x = maxXTemp + space;
	domainParams.min_y = minYTemp - space;
	domainParams.max_y = maxYTemp + space;
	domainParams.min_z = minZTemp - space;
	domainParams.max_z = maxZTemp + space;
};

void init_net_inct_bucket(
	NodeInfoVecs& nodeInfoVecs,
	DomainParams& domainParams,
	AuxVecs& auxVecs,
	GeneralParams& generalParams) {

	unsigned padding = 1;
	if (generalParams.iterationCounter == 0) {
		padding = 1;
	}
	else {
		padding = 1;
	}

	//on the first iteration, we allocate more, we don't plan on using it.
	//always set bucket count. Update total if different.
	domainParams.bucket_count_x = padding * ceil((domainParams.max_x - domainParams.min_x) / domainParams.grid_spacing_net_intc) + 1;
	domainParams.bucket_count_y = padding * ceil((domainParams.max_y - domainParams.min_y) / domainParams.grid_spacing_net_intc) + 1;
	domainParams.bucket_count_z = padding * ceil((domainParams.max_z - domainParams.min_z) / domainParams.grid_spacing_net_intc) + 1;

	if ( (domainParams.bucket_count_x * domainParams.bucket_count_x * domainParams.bucket_count_z) > domainParams.total_bucket_count_net_intc	) {
		std::cout<<"resetting grid for network interact" << std::endl;
		std::cout<<"x-bucket: "<< domainParams.bucket_count_x<<std::endl;
		std::cout<<"y-bucket: "<< domainParams.bucket_count_y<<std::endl;
		std::cout<<"z-bucket: "<< domainParams.bucket_count_z<<std::endl;

		//double amount of buckets in case of resizing networks
		domainParams.total_bucket_count_net_intc = 
			domainParams.bucket_count_x * 
			domainParams.bucket_count_y * 
			domainParams.bucket_count_z;
		std::cout<<"grid: "<< domainParams.grid_spacing_net_intc << std::endl;
		std::cout<<"total bucket count: "<< domainParams.total_bucket_count_net_intc<<std::endl;

		std::cout<<"min_x: " << domainParams.min_x << std::endl;
		std::cout<<"max_x: " << domainParams.max_x << std::endl;
		std::cout<<"min_y: " << domainParams.min_y << std::endl;
		std::cout<<"max_y: " << domainParams.max_y << std::endl;
		std::cout<<"min_z: " << domainParams.min_z << std::endl;
		std::cout<<"max_z: " << domainParams.max_z << std::endl;

		auxVecs.key_begin_net_intc.resize(domainParams.total_bucket_count_net_intc);
		auxVecs.key_end_net_intc.resize(domainParams.total_bucket_count_net_intc);

	}

	thrust::fill(auxVecs.key_begin_net_intc.begin(),auxVecs.key_begin_net_intc.end(),0);
	thrust::fill(auxVecs.key_end_net_intc.begin(),auxVecs.key_end_net_intc.end(),0);

};

//convert buckets into neighboring scheme
void extend_net_inct_bucket(
	NodeInfoVecs& nodeInfoVecs,
	DomainParams& domainParams,
	AuxVecs& auxVecs,
	GeneralParams& generalParams) {

	//memory is already allocated.
	unsigned endIndexExpanded = (auxVecs.end_index_bucket_keys_net_intc) * 27;


	//test for removing copies.
	unsigned valuesCount = auxVecs.id_value_net_intc.size();
	thrust::fill(auxVecs.id_bucket_expanded_net_intc.begin(),auxVecs.id_bucket_expanded_net_intc.end(),0);
	thrust::fill(auxVecs.id_value_expanded_net_intc.begin(),auxVecs.id_value_expanded_net_intc.end(),0);


	/*
	* beginning of constant iterator
	*/
	thrust::constant_iterator<unsigned> first(27);
	/*
	* end of constant iterator.
	* the plus sign only indicate movement of position, not value.
	* e.g. movement is 5 and first iterator is initialized as 9
	* result array is [9,9,9,9,9];
	*/
	thrust::constant_iterator<unsigned> last = first + (auxVecs.end_index_bucket_keys_net_intc); // this is NOT numerical addition!

	expand(first, last,
		thrust::make_zip_iterator(
			thrust::make_tuple(
				auxVecs.id_bucket_net_intc.begin(),
				auxVecs.id_value_net_intc.begin())),
		thrust::make_zip_iterator(
			thrust::make_tuple(
				auxVecs.id_bucket_expanded_net_intc.begin(),
				auxVecs.id_value_expanded_net_intc.begin())));

	thrust::counting_iterator<unsigned> countingBegin(0);

	thrust::transform(
		thrust::make_zip_iterator(
			thrust::make_tuple(
				auxVecs.id_bucket_expanded_net_intc.begin(),
				countingBegin)),
		thrust::make_zip_iterator(
			thrust::make_tuple(
				auxVecs.id_bucket_expanded_net_intc.begin(),
				countingBegin)) + endIndexExpanded,

		auxVecs.id_bucket_expanded_net_intc.begin(),
		functor_neighbor(
			domainParams.bucket_count_x,
			domainParams.bucket_count_y,
			domainParams.bucket_count_z));

	thrust::stable_sort_by_key(auxVecs.id_bucket_expanded_net_intc.begin(),
		auxVecs.id_bucket_expanded_net_intc.end(),
		auxVecs.id_value_expanded_net_intc.begin());


	thrust::counting_iterator<unsigned> search_begin(0);
	thrust::counting_iterator<unsigned> search_begin1(0);

	thrust::lower_bound(auxVecs.id_bucket_expanded_net_intc.begin(),
		auxVecs.id_bucket_expanded_net_intc.end(), search_begin,
		search_begin + domainParams.total_bucket_count_net_intc,
		auxVecs.key_begin_net_intc.begin());

	thrust::upper_bound(auxVecs.id_bucket_expanded_net_intc.begin(),
		auxVecs.id_bucket_expanded_net_intc.end(), search_begin1,
		search_begin1 + domainParams.total_bucket_count_net_intc,
		auxVecs.key_end_net_intc.begin());

		
	/*
	unsigned choice = 0;
	unsigned bucket = auxVecs.id_bucket_net_intc[choice];
	std::cout<<"bucketplt 0: "<< bucket<<std::endl;
	std::cout<<"plt pos: "<<nodeInfoVecs.node_loc_x[0]<<" "<<nodeInfoVecs.node_loc_y[0]<<" "<<nodeInfoVecs.node_loc_z[0]<<std::endl;
	std::cout<<"key len: "<< auxVecs.key_begin_net_intc.size() << std::endl;
	unsigned begin = auxVecs.key_begin_net_intc[bucket];
	unsigned end = auxVecs.key_end_net_intc[bucket];
	
	std::cout<<"from bucket scheme:"<<std::endl;
	for (unsigned i = begin; i < end; i++) {
		
		unsigned nbr = auxVecs.id_value_expanded_net_intc[i];
		unsigned buck = auxVecs.id_bucket_net_intc[nbr];
		double x_dist = nodeInfoVecs.node_loc_x[choice] - nodeInfoVecs.node_loc_x[nbr];
		double y_dist = nodeInfoVecs.node_loc_y[choice] - nodeInfoVecs.node_loc_y[nbr];
		double z_dist = nodeInfoVecs.node_loc_z[choice] - nodeInfoVecs.node_loc_z[nbr];
		double dist = std::sqrt(std::pow(x_dist,2.0)+std::pow(y_dist,2.0)+std::pow(z_dist,2.0));
		if ((dist < 19.0) && (dist>0)){
			std::cout<<"dist: "<< dist<< " between: "<< choice << " and nbr: "<< nbr<<std::endl; 
			std::cout<<"nbr: "<< nbr<< " is in bucket: "<< buck <<std::endl;
		}
	}*/

}


//takes nodes and places in buckets.
void build_net_inct_bucket(
	NodeInfoVecs& nodeInfoVecs,
	DomainParams& domainParams,
	AuxVecs& auxVecs,
	GeneralParams& generalParams) {


	thrust::counting_iterator<unsigned> indexBucketBegin(0);
	// takes counting iterator and coordinates
	// return tuple of keys and values
	// transform the points to their bucket indices

	thrust::transform(
		thrust::make_zip_iterator(
			thrust::make_tuple(
				nodeInfoVecs.node_loc_x.begin(),
				nodeInfoVecs.node_loc_y.begin(),
				nodeInfoVecs.node_loc_z.begin(),
				indexBucketBegin)),
		thrust::make_zip_iterator(
			thrust::make_tuple(
				nodeInfoVecs.node_loc_x.begin(),
				nodeInfoVecs.node_loc_y.begin(),
				nodeInfoVecs.node_loc_z.begin(),
				indexBucketBegin)) + generalParams.max_node_count,
		thrust::make_zip_iterator(
			thrust::make_tuple(
				auxVecs.id_bucket_net_intc.begin(),
				auxVecs.id_value_net_intc.begin())),
		functor_bucket_indexer(
			domainParams.min_x, domainParams.max_x, domainParams.min_y,
			domainParams.max_y, domainParams.min_z, domainParams.max_z,
			domainParams.bucket_count_x,
			domainParams.bucket_count_y,
			domainParams.bucket_count_z,
			domainParams.grid_spacing_net_intc));

//test sorting by node instaed of bucket index
thrust::sort_by_key(auxVecs.id_value_net_intc.begin(),
		auxVecs.id_value_net_intc.begin() + generalParams.max_node_count,
		auxVecs.id_bucket_net_intc.begin());

auxVecs.end_index_bucket_keys_net_intc = generalParams.max_node_count;
};
