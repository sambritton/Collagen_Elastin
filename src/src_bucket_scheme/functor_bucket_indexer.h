#ifndef FUNCTOR_BUCKET_INDEXER_H_
#define FUNCTOR_BUCKET_INDEXER_H_

#include "system_structures.h"

struct functor_bucket_indexer {
	double min_x;
	double max_x;
	double min_y;
	double max_y;
	double min_z;
	double max_z;

	unsigned bucket_count_x;
	unsigned bucket_count_y;
	unsigned bucket_count_z;
	double unitLen;

	__host__ __device__

	functor_bucket_indexer(
		double _minX,
		double _maxX,
		double _minY,
		double _maxY,
		double _minZ,
		double _maxZ,
		unsigned _XBucketCount,
		unsigned _YBucketCount,
		unsigned _ZBucketCount,
		double _unitLen):

		min_x(_minX),
		max_x(_maxX),
		min_y(_minY),
		max_y(_maxY),
		min_z(_minZ),
		max_z(_maxZ),
		bucket_count_x(_XBucketCount),
		bucket_count_y(_YBucketCount),
		bucket_count_z(_ZBucketCount),
		unitLen(_unitLen) {}

	__device__
	Tuu operator()(const Tdddu& v) {

			unsigned id = thrust::get<3>(v);
			if (id == 0){
				//comment
			}
			unsigned x = static_cast<unsigned>((thrust::get<0>(v) - min_x) / unitLen);
			unsigned y = static_cast<unsigned>((thrust::get<1>(v) - min_y) / unitLen);
			unsigned z = static_cast<unsigned>((thrust::get<2>(v) - min_z) / unitLen);

			unsigned total_bucket_count = bucket_count_x * bucket_count_y * bucket_count_z;
			// return the bucket's linear index and node's global index
			//return thrust::make_tuple(z * XSize * YSize + y * XSize + x, thrust::get<4>(v));
			unsigned bucket = z * bucket_count_x * bucket_count_y + y * bucket_count_x + x;
			//try to make it so bucket does not return unsigned32Max
			if (bucket > total_bucket_count) {
				bucket = 0;
			}
			return thrust::make_tuple(bucket, id);

	}
};

#endif
