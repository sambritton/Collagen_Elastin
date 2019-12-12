#ifndef FUNCTOR_NEIGHBOR_H_
#define FUNCTOR_NEIGHBOR_H_

#include "SystemStructures.h"

/*
* Functor to compute neighbor buckets(center bucket included) of a node.
* @param input1 bucket index of node
* @param input2 pick from the sequence, which is also global rank of the node
*
* @return output1 bucket indices of node ( all neighbors and the center bucket) of node
* @return output2 global rank of node

example with 3x3x3 grid. with node at position (1,1,1). Notice there are 27 possibilities.
0: (1,1,1)
1: (0,0,1)
2: (1,0,1)
The words left & right denote x change, top and bottom denote y change and upper & lower denote z change
 */

struct functor_neighbor : public thrust::unary_function<Tuu, Tuu> {
	unsigned numOfBucketsInXDim;
	unsigned numOfBucketsInYDim;
	unsigned numOfBucketsInZDim;

	__host__ __device__ functor_neighbor(
		unsigned _numOfBucketsInXDim,
		unsigned _numOfBucketsInYDim,
		unsigned _numOfBucketsInZDim ) :
		numOfBucketsInXDim(_numOfBucketsInXDim),
		numOfBucketsInYDim(_numOfBucketsInYDim),
		numOfBucketsInZDim(_numOfBucketsInZDim) {}

	__device__ unsigned operator()(const Tuu &v) {
		unsigned relativeRank = thrust::get<1>(v) % 27;	//27 = 3^3. Takes global node id for calculation

		//takes global bucket id for calculation
		if ((thrust::get<1>(v) == (27*479+1))) {
			//coment

		}
		unsigned area = numOfBucketsInXDim * numOfBucketsInYDim;
		//unsigned volume = area * numOfBucketsInZDim;

		unsigned xPos = thrust::get<0>(v) % numOfBucketsInXDim;	//col
		unsigned xPosLeft = xPos - 1;
		unsigned xPosRight = xPos + 1;
		if (xPos == 0) {
			//wraparound unsigned
			xPosLeft = numOfBucketsInXDim-1;
		}
		if (xPosRight >= numOfBucketsInXDim) {
			xPosRight = 0;
		}


		unsigned zPos = thrust::get<0>(v) / area; //z divide by area
		unsigned zPosUp = zPos + 1;
		unsigned zPosLow = zPos - 1;
		if (zPos == 0 ) {
			//wraparound unsigned
			zPosLow = numOfBucketsInZDim-1;
		}
		if (zPosUp >= numOfBucketsInZDim) {
			zPosUp = 0;
		}

		unsigned yPos = (thrust::get<0>(v) - zPos * area) / numOfBucketsInXDim;	//row
		unsigned yPosTop = yPos + 1;
		unsigned yPosBottom = yPos - 1;

		if (yPos == 0) {
			//wraparound unsigend
			yPosBottom = numOfBucketsInYDim-1;
		}
		if (yPosTop >= numOfBucketsInYDim) {
			yPosTop = 0;
		}

		switch (relativeRank) {
		//middle cases
		case 0:
			return thrust::get<0>(v);
			//break;
		case 1:{
				unsigned topLeft = xPosLeft + yPosTop * numOfBucketsInXDim + zPos * area;
				return (topLeft);
				//break;
		}
		case 2:{
				unsigned top = xPos + yPosTop * numOfBucketsInXDim + zPos * area;
			return (top);
			//break;
		}
		case 3:{
				unsigned topRight = xPosRight + yPosTop * numOfBucketsInXDim + zPos * area;
			return topRight;
			//break;
		}
		case 4:{
				unsigned right = xPosRight + yPos * numOfBucketsInXDim + zPos * area;
			return right;
			//break;
		}
		case 5:{
				unsigned bottomRight = xPosRight + yPosBottom * numOfBucketsInXDim + zPos * area;
			return bottomRight;
			//break;
		}
		case 6:{
				unsigned bottom = xPos + yPosBottom * numOfBucketsInXDim + zPos * area;
			return bottom;
			//break;
		}
		case 7:{
				unsigned bottomLeft = xPosLeft + yPosBottom * numOfBucketsInXDim + zPos * area;
			return bottomLeft;
			//break;
		}
		case 8:{
				unsigned left = xPosLeft + yPos * numOfBucketsInXDim + zPos * area;
			return left;
			//break;
		}
		//lower Z cases
		case 9:{
				unsigned lowerCenter = xPos + yPos * numOfBucketsInXDim +  zPosLow * area;
			return lowerCenter;
			//break;
		}
		case 10:{
				unsigned lowerTopLeft = xPosLeft + yPosTop * numOfBucketsInXDim + zPosLow* area;
			return lowerTopLeft;
			//break;
		}
		case 11:{
				unsigned lowerTop = xPos + yPosTop * numOfBucketsInXDim + zPosLow * area;
			return (lowerTop);
			//break;
		}
		case 12:{
				unsigned lowerTopRight = xPosRight + yPosTop * numOfBucketsInXDim  + zPosLow * area;
			return lowerTopRight;
			//break;
		}
		case 13:{
				unsigned lowerRight = xPosRight + yPos * numOfBucketsInXDim + zPosLow * area;
			return (lowerRight);
			//break;
		}
		case 14:{
				unsigned lowerBottomRight = xPosRight + yPosBottom * numOfBucketsInXDim + zPosLow * area;
			return (lowerBottomRight);
			//break;
		}
		case 15:{
				unsigned lowerBottom = xPos + yPosBottom * numOfBucketsInXDim + zPosLow * area;
			return lowerBottom;
			//break;
		}
		case 16:{
				unsigned lowerBottomLeft = xPosLeft + yPosBottom * numOfBucketsInXDim  + zPosLow * area;
			return lowerBottomLeft;
			//break;
		}
		case 17:{
				unsigned lowerLeft = xPosLeft + yPos * numOfBucketsInXDim + zPosLow * area;
			return lowerLeft;
			//break;
		}
		//upper Z cases
		case 18:{
				unsigned upperCenter = xPos + yPos * numOfBucketsInXDim +  zPosUp * area;
			return (upperCenter);
			//break;
		}
		case 19:{
				unsigned upperTopLeft = xPosLeft + yPosTop * numOfBucketsInXDim + zPosUp * area;
			return (upperTopLeft);
			//break;
		}
		case 20:{
				unsigned upperTop = xPos + yPosTop * numOfBucketsInXDim + zPosUp * area;
			return (upperTop);
			//break;
		}
		case 21:{
				unsigned upperTopRight = xPosRight + yPosTop * numOfBucketsInXDim  + zPosUp * area;
			return (upperTopRight);
			//break;
		}
		case 22:{
				unsigned upperRight = xPosRight + yPos * numOfBucketsInXDim + zPosUp * area;
			return (upperRight);
			//break;
		}
		case 23:{
				unsigned upperBottomRight = xPosRight + yPosBottom * numOfBucketsInXDim + zPosUp * area;
			return (upperBottomRight);
			//break;
		}
		case 24:{
				unsigned upperBottom = xPos + yPosBottom * numOfBucketsInXDim + zPosUp * area;
			return (upperBottom);
			//break;
		}
		case 25:{
				unsigned upperBottomLeft = xPosLeft + yPosBottom * numOfBucketsInXDim  + zPosUp * area;
			return (upperBottomLeft);
			//break;
		}
		case 26:{
				unsigned upperLeft = xPosLeft + yPos * numOfBucketsInXDim + zPosUp * area;
			return (upperLeft);
			//break;
		}
		default:{
			unsigned default_Id = ULONG_MAX;
			return (default_Id);
		}

		}
	}
};

#endif