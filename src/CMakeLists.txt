cmake_minimum_required(VERSION 2.9 FATAL_ERROR)

set(src_bucket_scheme "${CMAKE_CURRENT_SOURCE_DIR}/src_bucket_scheme")
set(src_fibrin "${CMAKE_CURRENT_SOURCE_DIR}/src_fibrin")
set(src_plt "${CMAKE_CURRENT_SOURCE_DIR}/src_plt")
set(SRC "${CMAKE_CURRENT_SOURCE_DIR}")
MESSAGE("SRC: " ${SRC})

include_directories(src_bucket_scheme)
include_directories(src_fibrin)
include_directories(src_plt)


project(Fibrin_Platelet LANGUAGES CXX CUDA)

find_package(CUDA REQUIRED)

set(CMAKE_CUDA_FLAGS " -g -G ${CMAKE_CUDA_FLAGS}")
	
MESSAGE("CMAKE_CXX_COMPILER_VERSION: " ${CMAKE_CXX_COMPILER_VERSION})
MESSAGE("CMAKE_CUDA_COMPILER: " ${CMAKE_CUDA_COMPILER})
MESSAGE("CMAKE_CUDA_FLAGS: " ${CMAKE_CUDA_FLAGS})

set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS}; -g -G -Xcompiler -DTHRUST_DEVICE_SYSTEM 
  -gencode=arch=compute_20,code=sm_20
  -gencode=arch=compute_30,code=sm_30
  -gencode=arch=compute_35,code=sm_35
  -gencode=arch=compute_50,code=sm_50
  -gencode=arch=compute_60,code=sm_60)
MESSAGE("NVCC flags are : " ${CUDA_NVCC_FLAGS})

# THIS SECTION GETS INCLUDES FOR THRUST
find_path( THRUST_INCLUDE_DIR
  HINTS
	${CUDA_INCLUDE_DIRS}
	NAMES thrust/version.h
	DOC "Thrust headers"
)

if( THRUST_INCLUDE_DIR )

  list( REMOVE_DUPLICATES THRUST_INCLUDE_DIR )
endif( THRUST_INCLUDE_DIR )
MESSAGE("thrust version : " ${THRUST_INCLUDE_DIR}) 
include_directories(${THRUST_INCLUDE_DIR}) #specifies path that compiler looks in while searching for header files


#make PUGIXML LIB
set (PUGI_LOCATION ${CMAKE_CURRENT_LIST_DIR}/src_ext_lib/pugixml/src)
set (PUGI_INCLUDE ${CMAKE_CURRENT_LIST_DIR}/src_ext_lib/pugixml/include)
MESSAGE("PUGI_LOCATION: " ${PUGI_LOCATION})
MESSAGE("CMAKE_CURRENT_LIST_DIR: " ${CMAKE_CURRENT_LIST_DIR})
include_directories(${PUGI_LOCATION})
include_directories(${PUGI_INCLUDE})

add_subdirectory( ${CMAKE_CURRENT_LIST_DIR}/src_ext_lib/pugixml)  
add_library (pugixmlLib STATIC
	${PUGI_LOCATION}/pugiconfig.hpp  ${PUGI_LOCATION}/pugixml.cpp ${PUGI_LOCATION}/pugixml.hpp )

MESSAGE("pugixmlLib: " ${pugixmlLib})
MESSAGE("CMAKE_SOURCE_DIR: " ${CMAKE_SOURCE_DIR})

#END PUGI

#glm
#Get glm includes
set (GLM_LOCATION ${CMAKE_CURRENT_LIST_DIR}/src_ext_lib/glm)
MESSAGE("GLM_LOCATION: " ${GLM_LOCATION})
include_directories(${GLM_LOCATION})
#end GLM


#To include libraries in another folder (src_dg) add them
add_subdirectory(${src_bucket_scheme})
add_subdirectory(${src_fibrin})
add_subdirectory(${src_plt})

add_library(storageLib STATIC
	Storage.cu
	Storage.h
	System.h
	System_Builder.h
	SystemStructures.h)

	
target_compile_features(storageLib PUBLIC cxx_std_11)

add_library(systemLib STATIC
	System.cu
	System.h
	System_Builder.cpp
	System_Builder.h
	functor_misc.h
	SystemStructures.h)

	
target_compile_features(systemLib PUBLIC cxx_std_11)

add_dependencies(systemLib bucket_schemeLib)
add_dependencies(systemLib fibrinLib)
add_dependencies(systemLib pltLib)
add_dependencies(systemLib storageLib)


target_link_libraries(systemLib bucket_schemeLib)
target_link_libraries(systemLib fibrinLib)
target_link_libraries(systemLib pltLib)
target_link_libraries(systemLib storageLib)

set(MAIN_SOURCES main.cpp)
set(TARGET bend-model)
add_executable(${TARGET} ${MAIN_SOURCES})

MESSAGE("TARGET: " ${TARGET})

add_dependencies(${TARGET} pugixmlLib)


target_link_libraries(${TARGET} systemLib )
target_link_libraries(${TARGET} pugixmlLib )

