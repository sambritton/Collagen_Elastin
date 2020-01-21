
include_directories(${CUDA_INCLUDE_DIRS})
include_directories(${THRUST_INCLUDE_DIR})
include_directories(${SRC}) #include for higher source files include

MESSAGE("SRC in SRC_FIBRIN: " ${SRC})
MESSAGE("SRC in CUDA_INCLUDE_DIRS: " ${CUDA_INCLUDE_DIRS})

add_library(col_elaLib STATIC
	advance_positions.cu
	advance_positions.h
	functor_advance_pos.h
	link_nodes.cu
	link_nodes.h
	functor_link_nodes.h
	functor_de_link_nodes.h
	collagen_elastin_spring.cu
	collagen_elastin_spring.h
	functor_collagen_elastin.h
	bending_spring.cu
	bending_spring.h
	functor_bending.h
	params_calc.cu
	params_calc.h
	functor_calc_strain_params.h)


target_compile_features(col_elaLib PUBLIC cxx_std_11)