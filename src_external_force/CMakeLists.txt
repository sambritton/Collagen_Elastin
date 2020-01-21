include_directories(${CUDA_INCLUDE_DIRS})
include_directories(${THRUST_INCLUDE_DIR})
include_directories(${SRC}) #include for higher source files include

MESSAGE("SRC in SRC_BUCKET: " ${SRC})
add_library(external_forceLib STATIC
	external_force.cu
	external_force.h
	functor_strain.h
	functor_external_force.h)


target_compile_features(external_forceLib PUBLIC cxx_std_11)
