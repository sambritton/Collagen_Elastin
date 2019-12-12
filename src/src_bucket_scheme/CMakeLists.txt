include_directories(${CUDA_INCLUDE_DIRS})
include_directories(${THRUST_INCLUDE_DIR})
include_directories(${SRC}) #include for higher source files include

MESSAGE("SRC in SRC_BUCKET: " ${SRC})
add_library(bucket_schemeLib STATIC
	Bucket_Net.cu
	Bucket_Net.h
	Bucket_Plt.cu
	Bucket_Plt.h
	function_extend.h
	functor_bucket_indexer.h
	functor_neighbor.h)


target_compile_features(bucket_schemeLib PUBLIC cxx_std_11)