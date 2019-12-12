#ifndef FUNCTOR_BUCKET_INDEXER_H_
#define FUNCTOR_BUCKET_INDEXER_H_

#include "SystemStructures.h"

struct functor_bucket_indexer {
	double minX;
	double maxX;
	double minY;
	double maxY;
	double minZ;
	double maxZ;

	unsigned XBucketCount;
	unsigned YBucketCount;
	unsigned ZBucketCount;
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
		double _unitLen) :

		minX(_minX),
		maxX(_maxX),
		minY(_minY),
		maxY(_maxY),
		minZ(_minZ),
		maxZ(_maxZ),
		XBucketCount(_XBucketCount),
		YBucketCount(_YBucketCount),
		ZBucketCount(_ZBucketCount),
		unitLen(_unitLen) {}

	__device__ 
	Tuu operator()(const Tdddu& v) {

			unsigned id = thrust::get<3>(v);
			if (id == 0){
				//comment
			}
			unsigned x = static_cast<unsigned>((thrust::get<0>(v) - minX) / unitLen);
			unsigned y = static_cast<unsigned>((thrust::get<1>(v) - minY) / unitLen);
			unsigned z = static_cast<unsigned>((thrust::get<2>(v) - minZ) / unitLen);


			// return the bucket's linear index and node's global index
			//return thrust::make_tuple(z * XSize * YSize + y * XSize + x, thrust::get<4>(v));
			unsigned bucket = z * XBucketCount * YBucketCount + y * XBucketCount + x;
			//try to make it so bucket does not return unsigned32Max
			if (bucket == ULONG_MAX) {
				bucket = 0;
			}
			return thrust::make_tuple(bucket, id);

	}
};

#endif