//This file sets the grid for network self interaction

#include "SystemStructures.h"
#include "Bucket_Net.h"
#include "System.h"

#include "functor_neighbor.h"
#include "functor_bucket_indexer.h"
#include "function_extend.h"



//take domain and discretize into square buckets of size gridspace
void init_dim_general(
	NodeInfoVecs& nodeInfoVecs,
	PltInfoVecs& pltInfoVecs,
	DomainParams& domainParams,
	AuxVecs& auxVecs,
	GeneralParams& generalParams) {

	double minXTemp = (*(thrust::min_element(nodeInfoVecs.nodeLocX.begin(), nodeInfoVecs.nodeLocX.end())));
	double maxXTemp = (*(thrust::max_element(nodeInfoVecs.nodeLocX.begin(), nodeInfoVecs.nodeLocX.end())));
	double minYTemp = (*(thrust::min_element(nodeInfoVecs.nodeLocY.begin(), nodeInfoVecs.nodeLocY.end())));
	double maxYTemp = (*(thrust::max_element(nodeInfoVecs.nodeLocY.begin(), nodeInfoVecs.nodeLocY.end())));
	double minZTemp = (*(thrust::min_element(nodeInfoVecs.nodeLocZ.begin(), nodeInfoVecs.nodeLocZ.end())));
	double maxZTemp = (*(thrust::max_element(nodeInfoVecs.nodeLocZ.begin(), nodeInfoVecs.nodeLocZ.end())));

	//platelets
	if (generalParams.maxPltCount != 0) {
		domainParams.pltminX = (*(thrust::min_element(pltInfoVecs.pltLocX.begin(), pltInfoVecs.pltLocX.end())));
		domainParams.pltmaxX = (*(thrust::max_element(pltInfoVecs.pltLocX.begin(), pltInfoVecs.pltLocX.end())));
		domainParams.pltminY = (*(thrust::min_element(pltInfoVecs.pltLocY.begin(), pltInfoVecs.pltLocY.end())));
		domainParams.pltmaxY = (*(thrust::max_element(pltInfoVecs.pltLocY.begin(), pltInfoVecs.pltLocY.end())));
		domainParams.pltminZ = (*(thrust::min_element(pltInfoVecs.pltLocZ.begin(), pltInfoVecs.pltLocZ.end())));
		domainParams.pltmaxZ = (*(thrust::max_element(pltInfoVecs.pltLocZ.begin(), pltInfoVecs.pltLocZ.end())));
	}
	else {
		domainParams.pltminX = minXTemp;
		domainParams.pltmaxX = maxXTemp;
		domainParams.pltminY = minYTemp;
		domainParams.pltmaxY = maxYTemp;
		domainParams.pltminZ = minZTemp;
		domainParams.pltmaxZ = maxZTemp;
	}

	double space = 0.0;
	domainParams.minX = min(minXTemp, domainParams.pltminX) - space;
	domainParams.maxX = max(maxXTemp, domainParams.pltmaxX) + space;
	domainParams.minY = min(minYTemp, domainParams.pltminY) - space;
	domainParams.maxY = max(maxYTemp, domainParams.pltmaxY) + space;
	domainParams.minZ = min(minZTemp, domainParams.pltminZ) - space;
	domainParams.maxZ = max(maxZTemp, domainParams.pltmaxZ) + space;
};

void init_net_inct_bucket(
	NodeInfoVecs& nodeInfoVecs,
	PltInfoVecs& pltInfoVecs,
	DomainParams& domainParams,
	AuxVecs& auxVecs,
	GeneralParams& generalParams) {

	unsigned padding = 1;
	if (generalParams.iterationCounter == 0) {
		padding = 2;
	}
	else {
		padding = 1;
	}
	
	//on the first iteration, we allocate more, we don't plan on using it. 
	//always set bucket count. Update total if different. 
	domainParams.XBucketCount_net_intc = padding * ceil((domainParams.maxX - domainParams.minX) / domainParams.gridSpacing_net_intc) + 1;
	domainParams.YBucketCount_net_intc = padding * ceil((domainParams.maxY - domainParams.minY) / domainParams.gridSpacing_net_intc) + 1;
	domainParams.ZBucketCount_net_intc = padding * ceil((domainParams.maxZ - domainParams.minZ) / domainParams.gridSpacing_net_intc) + 1;

	if ( (domainParams.XBucketCount_net_intc * domainParams.YBucketCount_net_intc * domainParams.ZBucketCount_net_intc) > domainParams.totalBucketCount_net_intc	) {
		std::cout<<"resetting grid for network interact" << std::endl;
		std::cout<<"x-bucket: "<< domainParams.XBucketCount_net_intc<<std::endl;
		std::cout<<"y-bucket: "<< domainParams.YBucketCount_net_intc<<std::endl;
		std::cout<<"z-bucket: "<< domainParams.ZBucketCount_net_intc<<std::endl;

		//double amount of buckets in case of resizing networks
		domainParams.totalBucketCount_net_intc = domainParams.XBucketCount_net_intc * domainParams.YBucketCount_net_intc * domainParams.ZBucketCount_net_intc;
		std::cout<<"grid: "<< domainParams.gridSpacing_net_intc << std::endl;
		std::cout<<"total bucket count: "<< domainParams.totalBucketCount_net_intc<<std::endl;

		std::cout<<"minX: " << domainParams.minX << std::endl;
		std::cout<<"maxX: " << domainParams.maxX << std::endl;
		std::cout<<"minY: " << domainParams.minY << std::endl;
		std::cout<<"maxY: " << domainParams.maxY << std::endl;
		std::cout<<"minZ: " << domainParams.minZ << std::endl;
		std::cout<<"maxZ: " << domainParams.maxZ << std::endl;

		auxVecs.keyBegin_net_intc.resize(domainParams.totalBucketCount_net_intc);
		auxVecs.keyEnd_net_intc.resize(domainParams.totalBucketCount_net_intc);
 
	}

	thrust::fill(auxVecs.keyBegin_net_intc.begin(),auxVecs.keyBegin_net_intc.end(),0);
	thrust::fill(auxVecs.keyEnd_net_intc.begin(),auxVecs.keyEnd_net_intc.end(),0);

};

//convert buckets into neighboring scheme
void extend_net_inct_bucket(
	NodeInfoVecs& nodeInfoVecs,
	PltInfoVecs& pltInfoVecs,
	DomainParams& domainParams,
	AuxVecs& auxVecs,
	GeneralParams& generalParams) {

	//memory is already allocated.
	unsigned endIndexExpanded = (auxVecs.endIndexBucketKeys_net_intc) * 27;
	

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
	thrust::constant_iterator<unsigned> last = first + (auxVecs.endIndexBucketKeys_net_intc); // this is NOT numerical addition!

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
			domainParams.XBucketCount_net_intc,
			domainParams.YBucketCount_net_intc,
			domainParams.ZBucketCount_net_intc)); 

	thrust::stable_sort_by_key(auxVecs.id_bucket_expanded_net_intc.begin(),
		auxVecs.id_bucket_expanded_net_intc.end(),
		auxVecs.id_value_expanded_net_intc.begin());


	thrust::counting_iterator<unsigned> search_begin(0);

	thrust::lower_bound(auxVecs.id_bucket_expanded_net_intc.begin(),
		auxVecs.id_bucket_expanded_net_intc.end(), search_begin,
		search_begin + domainParams.totalBucketCount_net_intc,
		auxVecs.keyBegin_net_intc.begin());

	thrust::upper_bound(auxVecs.id_bucket_expanded_net_intc.begin(),
		auxVecs.id_bucket_expanded_net_intc.end(),search_begin,
		search_begin + domainParams.totalBucketCount_net_intc,
		auxVecs.keyEnd_net_intc.begin());


	/*
	unsigned choice = 0;

	unsigned bucket = auxVecs.idPlt_bucket[choice];
	std::cout<<"bucketplt 0: "<< bucket<<std::endl;
	std::cout<<"plt pos: "<<pltInfoVecs.pltLocX[0]<<" "<<pltInfoVecs.pltLocY[0]<<" "<<pltInfoVecs.pltLocZ[0]<<std::endl;
	std::cout<<"key len: "<< auxVecs.keyBegin.size() << std::endl;
	unsigned begin = auxVecs.keyBegin[bucket];
	unsigned end = auxVecs.keyEnd[bucket];
	
	std::cout<<"from bucket scheme:"<<std::endl;
	for (unsigned i = begin; i < end; i++) {
		
		unsigned nbr = auxVecs.id_value_expanded[i];
		unsigned buck = auxVecs.id_bucket[nbr];
		double x_dist = pltInfoVecs.pltLocX[choice] - nodeInfoVecs.nodeLocX[nbr];
		double y_dist = pltInfoVecs.pltLocY[choice] - nodeInfoVecs.nodeLocY[nbr];
		double z_dist = pltInfoVecs.pltLocZ[choice] - nodeInfoVecs.nodeLocZ[nbr];
		double dist = std::sqrt(std::pow(x_dist,2.0)+std::pow(y_dist,2.0)+std::pow(z_dist,2.0));
		if (dist < 1.0){
			std::cout<<"dist: "<< dist<< " between: "<< choice << " and nbr: "<< nbr<<std::endl; 
			std::cout<<"nbr: "<< nbr<< " is in bucket: "<< buck <<std::endl;
		}
	}*/

	/*
	std::cout<<"from all plt:"<<std::endl;
	for (unsigned i = 0; i < generalParams.maxNodeCount; i++) {
		unsigned nbr = i;//auxVecs.id_value_expanded[i];
		unsigned buck = auxVecs.id_bucket[nbr];
		double x_dist = pltInfoVecs.pltLocX[choice] - nodeInfoVecs.nodeLocX[nbr];
		double y_dist = pltInfoVecs.pltLocY[choice] - nodeInfoVecs.nodeLocY[nbr];
		double z_dist = pltInfoVecs.pltLocZ[choice] - nodeInfoVecs.nodeLocZ[nbr];
		double dist = std::sqrt(std::pow(x_dist,2.0)+std::pow(y_dist,2.0)+std::pow(z_dist,2.0));
		if (dist < 1.0){
			std::cout<<"dist: "<< dist<< " between: "<< choice << " and nbr: "<< nbr<<std::endl; 
			std::cout<<"nbr: "<< nbr<< " is in bucket: "<< buck <<std::endl;
		} 
	}*/


}


//takes nodes and places in buckets.
void build_net_inct_bucket(
	NodeInfoVecs& nodeInfoVecs,
	PltInfoVecs& pltInfoVecs,
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
				nodeInfoVecs.nodeLocX.begin(),
				nodeInfoVecs.nodeLocY.begin(),
				nodeInfoVecs.nodeLocZ.begin(),
				indexBucketBegin)),
		thrust::make_zip_iterator(
			thrust::make_tuple(
				nodeInfoVecs.nodeLocX.begin(),
				nodeInfoVecs.nodeLocY.begin(),
				nodeInfoVecs.nodeLocZ.begin(),
				indexBucketBegin)) + generalParams.maxNodeCount,
		thrust::make_zip_iterator(
			thrust::make_tuple(
				auxVecs.id_bucket_net_intc.begin(),
				auxVecs.id_value_net_intc.begin())),
		functor_bucket_indexer(
			domainParams.minX, domainParams.maxX, domainParams.minY,
			domainParams.maxY, domainParams.minZ, domainParams.maxZ,
			domainParams.XBucketCount_net_intc,
			domainParams.YBucketCount_net_intc,
			domainParams.ZBucketCount_net_intc,
			domainParams.gridSpacing_net_intc));

//test sorting by node instaed of bucket index
thrust::sort_by_key(auxVecs.id_value_net_intc.begin(),
		auxVecs.id_value_net_intc.begin() + generalParams.maxNodeCount,
		auxVecs.id_bucket_net_intc.begin());

auxVecs.endIndexBucketKeys_net_intc = generalParams.maxNodeCount;

 
};
