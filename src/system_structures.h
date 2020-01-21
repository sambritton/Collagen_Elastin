#ifndef SYSTEM_STRUCTURES_H_
#define SYSTEM_STRUCTURES_H_

#include <memory>
#include <cmath>

#include <thrust/random.h>
#include <thrust/count.h>
#include <thrust/extrema.h>
#include <thrust/sort.h>
#include <thrust/for_each.h>
#include <thrust/transform.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/tuple.h>
#include <thrust/device_vector.h>
#include <thrust/device_ptr.h>
#include <thrust/remove.h>
#include <thrust/copy.h>
#include <thrust/execution_policy.h>
#include <thrust/pair.h>
#include <thrust/unique.h>
#include <thrust/remove.h>
#include <thrust/binary_search.h>
#include <thrust/reduce.h>
#include <thrust/replace.h>
#include <thrust/gather.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/random/uniform_real_distribution.h>
#include <stdint.h>
#include <thrust/sequence.h>


struct NodeInfoVecs;
struct GeneralParams;
struct DomainParams;
struct EdgeInfoVecs;
struct BendInfoVecs;
struct AuxVecs;


typedef thrust::tuple<unsigned, bool, double> Tubd;
typedef thrust::tuple<unsigned, bool> Tub;
typedef thrust::tuple<unsigned, double> Tud;
typedef thrust::tuple<bool, double> Tbd;


typedef thrust::tuple<unsigned, unsigned, double> Tuud;

typedef thrust::tuple<unsigned, unsigned, unsigned, unsigned, double> Tuuuud;
typedef thrust::tuple<unsigned, unsigned, unsigned,unsigned> Tuuuu;
typedef thrust::tuple<unsigned, unsigned, unsigned> Tuuu;
typedef thrust::tuple<unsigned, unsigned, bool> Tuub;
typedef thrust::tuple<unsigned, unsigned> Tuu;

typedef thrust::tuple<unsigned, double, double, double> Tuddd;
typedef thrust::tuple<double, double, double, unsigned> Tdddu;
typedef thrust::tuple<double, double> Tdd;

typedef thrust::tuple<bool, double, double, double, double, double, double> BoolCVec6;
typedef thrust::tuple<unsigned, double, double, double, double, double, double,double, double, double> UCVec9;

typedef thrust::tuple<unsigned, unsigned, double, double, double, double, double, double, double> U2CVec7;
typedef thrust::tuple<unsigned, unsigned, double, double, double, double, double, double> U2CVec6;
typedef thrust::tuple<unsigned, double, double, double, double, double, double> UCVec6;
typedef thrust::tuple<unsigned, unsigned, double, double, double> U2CVec3;
typedef thrust::tuple<unsigned, bool, double, double, double> UBCVec3;
typedef thrust::tuple<unsigned, double, double, double> UCVec3;

typedef thrust::tuple<double, double, double, double, double, double, double, double, double> CVec9;
typedef thrust::tuple<double, double, double, double, double, double, double, double> CVec8;
typedef thrust::tuple<double, double, double, double, double, double, double> CVec7;
typedef thrust::tuple<double, double, double, double, double, double> CVec6;
typedef thrust::tuple<double, double, double, double, double> CVec5;
typedef thrust::tuple<double, double, double, double> CVec4;
typedef thrust::tuple<double, double, double> CVec3;
typedef thrust::tuple<double, double> CVec2;

struct CVec3Add : public thrust::binary_function<CVec3, CVec3, CVec3> {
	__host__ __device__
		CVec3 operator()(const CVec3 &vec1, const CVec3 &vec2) {
		return thrust::make_tuple(
            thrust::get<0>(vec1) + thrust::get<0>(vec2),
			thrust::get<1>(vec1) + thrust::get<1>(vec2),
			thrust::get<2>(vec1) + thrust::get<2>(vec2));
	}
};

struct CVec3NormBinary {
	__host__ __device__

		double operator() (const CVec3& vec1, const CVec3& vec2) {
		//divide force by fiber cross section to get stress
		return fabs(
			((thrust::get<0>(vec1) - thrust::get<0>(vec2))) +
			((thrust::get<1>(vec1) - thrust::get<1>(vec2))) +
			((thrust::get<2>(vec1) - thrust::get<2>(vec2))));
	}
};

struct CVec3NormUnary {
	__host__ __device__
		double operator() (const CVec3& vec) {
		//divide force by fiber cross section to get stress
		return (
			sqrt(thrust::get<0>(vec) * thrust::get<0>(vec) +
			thrust::get<1>(vec) * thrust::get<1>(vec) +
			thrust::get<2>(vec) * thrust::get<2>(vec)));
	}
};

struct CVec3InnerProduct {
    __host__ __device__
    double operator() (CVec3& v1 ) {
        return (thrust::get<0>(v1) * thrust::get<0>(v1) +
                thrust::get<1>(v1) * thrust::get<1>(v1) +
                thrust::get<2>(v1) * thrust::get<2>(v1) );
    }
};


///////////////////////////////////////////////
//random number generators
struct psrnormgen {

    double a, b;

    __host__ __device__
	psrnormgen(
		double _a,
		double _b) :
		a(_a),
		b(_b) {}

    __host__ __device__
	double operator()(const unsigned n) const
    {
        thrust::default_random_engine rng(n);
        thrust::normal_distribution<float> dist(a, b);
        rng.discard(n);
        return dist(rng);
    }

};
struct psrunifgen {

    double a, b;

    __host__ __device__
	psrunifgen(
		double _a,
		double _b) :
		a(_a),
		b(_b) {}

    __host__ __device__
	double operator()(const unsigned n) const
    {
        thrust::default_random_engine rng(n);
        thrust::uniform_real_distribution<float> dist(a, b);
        rng.discard(n);
        return dist(rng);
    }

};


struct functor_initial_angle {
	double* locXAddr;
	double* locYAddr;
	double* locZAddr;
	__host__ __device__
		functor_initial_angle(
			double* _locXAddr,
			double* _locYAddr,
			double* _locZAddr) :
			locXAddr(_locXAddr),
			locYAddr(_locYAddr),
			locZAddr(_locZAddr) {}

	__device__
	double operator() (const Tuuu &u3) {
		unsigned indexLeft = thrust::get<0>(u3);
		unsigned indexCenter = thrust::get<1>(u3);
		unsigned indexRight = thrust::get<2>(u3);

		double distLCX = locXAddr[indexLeft] - locXAddr[indexCenter];
		double distLCY = locYAddr[indexLeft] - locYAddr[indexCenter];
		double distLCZ = locZAddr[indexLeft] - locZAddr[indexCenter];
		double distRCX = locXAddr[indexRight] - locXAddr[indexCenter];
		double distRCY = locYAddr[indexRight] - locYAddr[indexCenter];
		double distRCZ = locZAddr[indexRight] - locZAddr[indexCenter];
//		CVec3 distLC =

		//lengths between left & center, right & center
		double lenLC = sqrt(distLCX*distLCX + distLCY*distLCY + distLCZ*distLCZ);
		double lenRC = sqrt(distRCX*distRCX + distRCY*distRCY + distRCZ*distRCZ);


		//normalized dot product
		double cosTheta = (distLCX*distRCX + distLCY*distRCY + distLCZ*distRCZ) / (lenLC * lenRC);

		//rounding errors
		if (cosTheta < -1.0) {
			cosTheta = -1.0;
		}
		else if (cosTheta > 1.0) {
			cosTheta = 1.0;
		}

		//double currentAngle = acos(cosTheta);

		return acos(cosTheta);
	}
};



//used to calculate strain
struct AveStrainFunctor {
  __host__ __device__
  double operator() (const Tbd &b1d1) {
	  bool isStrainNode = thrust::get<0>(b1d1);
	  if (isStrainNode)
		  return thrust::get<1>(b1d1);
	  else
    	return 0.0;
  }
};




struct IsEqualToOne {
	 __host__ __device__

	bool operator()(const unsigned& x) {
		if ( x == 1 ) {
			return true;
		}
		else {
			return false;
		}
	}
};

//return true if not equal to input value.
struct isNotEqualToZero {

	__host__ __device__
	bool operator() (const unsigned& x) {
		return (x != 0);
	}
};

//return true if equal to input value.
struct isEqualZero {
	__host__ __device__
	bool operator() (const unsigned& x) {
		return (x == 0);
	}
};

struct tupleEqual {

  __host__ __device__
    bool operator()(Tuu x, Tuu y)
    {
      return ( (x.get<0>() == y.get<0>()) && (x.get<1>() == y.get<1>()) );
    }
};

struct is_greater_than {
	unsigned limit;

	 __host__ __device__
	is_greater_than(unsigned& _limit) : limit(_limit) {}
  __device__
  bool operator()(const unsigned& x) {
	if ( x > limit ) {
    	return true;
	}
	else {
		return false;
	}
  }
};

struct is_less_than {
	unsigned limit;

	 __host__ __device__
	is_less_than(unsigned& _limit) : limit(_limit) {}
  __device__
  bool operator()(const unsigned& x) {
	if ( x < limit ) {
    	return true;
	}
	else {
		return false;
	}
  }
};

#endif /* SYSTEMSTRUCTURES_H_*/
