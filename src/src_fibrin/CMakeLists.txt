
include_directories(${CUDA_INCLUDE_DIRS})
include_directories(${THRUST_INCLUDE_DIR})
include_directories(${SRC}) #include for higher source files include

MESSAGE("SRC in SRC_FIBRIN: " ${SRC})
MESSAGE("SRC in CUDA_INCLUDE_DIRS: " ${CUDA_INCLUDE_DIRS})

add_library(fibrinLib STATIC
	Advance_Positions_Fibrin.cu
	Advance_Positions_Fibrin.h
	functor_advance_pos.h
	Link_Nodes.cu
	Link_Nodes.h
	functor_link_nodes.h
	functor_de_link_nodes.h
	WLC_Force.cu
	WLC_Force.h
	functor_wlc.h
	Torsion_Force.cu
	Torsion_Force.h
	functor_torsion.h
	Params_Calc.cu
	Params_Calc.h
	functor_calc_strain_params.h)


target_compile_features(fibrinLib PUBLIC cxx_std_11)