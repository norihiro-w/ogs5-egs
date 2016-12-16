# Source files

set( HEADERS
	BoundaryCondition.h
	DistributionInfo.h
	DistributionTools.h
	eos.h
	fem_ele.h
	fem_ele_std.h
	fem_ele_vec.h
	FEMCondition.h
	FEMEnums.h
	femlib.h
	files0.h
	GeoInfo.h
	InitialCondition.h
	LinearFunctionData.h
	mathlib.h
	matrix_class.h
	mechanics_utils.h
	Output.h
	pcs_dm.h
	problem.h
	ProcessInfo.h
	prototyp.h
	rf_bc_new.h
	rf_fct.h
	rf_fluid_momentum.h
	rf_ic_new.h
	rf_kinreact.h
	rf_mfp_new.h
	rf_mmp_new.h
	rf_msp_new.h
	rf_node.h
	rf_num_new.h
	rf_out_new.h
	rf_pcs.h
	rf_pcs_TH.h
	rf_react.h
	rf_st_new.h
	rf_tim_new.h
	rfmat_cp.h
	solver.h
	SourceTerm.h
	SparseMatrixDOK.h
	Stiff_Bulirsch-Stoer.h
	tools.h
)

set( SOURCES
	BoundaryCondition.cpp
	DistributionInfo.cpp
	DistributionTools.cpp
	eos.cpp
	fem_ele.cpp
	fem_ele_std.cpp
	fem_ele_std1.cpp
	fem_ele_vec.cpp
	FEMCondition.cpp
	FEMEnums.cpp
	femlib.cpp
	files0.cpp
	GeoInfo.cpp
	InitialCondition.cpp
	LinearFunctionData.cpp
	mathlib.cpp
	matrix_class.cpp
	mechanics_utils.cpp
	Output.cpp
	pcs_dm.cpp
	problem.cpp
	ProcessInfo.cpp
	rf_bc_new.cpp
	rf_fct.cpp
	rf_fluid_momentum.cpp
	rf_ic_new.cpp
	rf_kinreact.cpp
	rf_mfp_new.cpp
	rf_mmp_new.cpp
	rf_msp_new.cpp
	rf_node.cpp
	rf_num_new.cpp
	rf_out_new.cpp
	rf_pcs.cpp
	rf_pcs_TH.cpp
	rf_react.cpp
	rf_st_new.cpp
	rf_tim_new.cpp
	rfmat_cp.cpp
	SourceTerm.cpp
	SparseMatrixDOK.cpp
	Stiff_Bulirsch-Stoer.cpp
	tools.cpp
)

IF(OGS_FEM_PETSC)
	SET( SOURCES ${SOURCES} rf_pcs1.cpp fct_mpi.h fct_mpi.cpp)
ELSE(OGS_FEM_PETSC)
	SET( SOURCES ${SOURCES} par_ddc.h par_ddc.cpp)
	## New sparse class. WW
	IF(OGS_FEM_MKL OR OGS_FEM_LIS)
		SET( SOURCES ${SOURCES} equation_class.h equation_class.cpp )
	ELSE()
		SET( SOURCES ${SOURCES} solver.h solver.cpp matrix_routines.h matrix_routines.cpp)
	ENDIF()
ENDIF()

IF(NOT OGS_ONLY_TH)
	SET( SOURCES ${SOURCES} rf_random_walk.h rf_random_walk.cpp)
ENDIF( )



include_directories(
	${CMAKE_SOURCE_DIR}/FEM
	${CMAKE_SOURCE_DIR}/Base
	${CMAKE_SOURCE_DIR}/MathLib
	${CMAKE_SOURCE_DIR}/GEO
	${CMAKE_SOURCE_DIR}/FileIO
	${CMAKE_SOURCE_DIR}/FileIO/MeshIO
	${CMAKE_SOURCE_DIR}/FileIO/FEMIO
	${CMAKE_SOURCE_DIR}/MSH
	${CMAKE_SOURCE_DIR}/ThirdParty
)

IF (LIS)
  INCLUDE_DIRECTORIES(${LIS_INCLUDE_DIR})
ENDIF (LIS)

IF (MKL)
  INCLUDE_DIRECTORIES(${MKL_INCLUDES})
ENDIF (MKL)

IF (PARALLEL_USE_MPI)
	INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
ENDIF (PARALLEL_USE_MPI)

# Create the library
add_library( FEM STATIC ${HEADERS} ${SOURCES} )

target_link_libraries(FEM
	Base
	FileIO
	GEO
	MSH
)


if (WIN32)
	if (LIS)
		target_link_libraries( FEM ${LIS_LIBRARIES} ) 
		set_target_properties( FEM PROPERTIES LINK_FLAGS "/NODEFAULTLIB:libcmt" )
	endif (LIS)
	if (MKL)
		target_link_libraries( FEM ${MKL_LIBRARIES} )
	endif (MKL)
	if (PARALLEL_USE_MPI)
		SET( CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} ${MPI_COMPILE_FLAGS} )
		target_link_libraries( FEM ${MPI_LIBRARIES} )
	endif (PARALLEL_USE_MPI)
endif (WIN32)

if (UNIX)
	if (LIS)
		target_link_libraries( FEM ${LIS_LIBRARIES} )
	endif (LIS)

	if (MKL)
		target_link_libraries( FEM ${MKL_LIBRARIES} )
	endif (MKL)
	if (PARALLEL_USE_MPI)
#		SET( CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} ${MPI_COMPILE_FLAGS} )
		target_link_libraries( FEM ${MPI_LIBRARIES} )
	endif (PARALLEL_USE_MPI)
endif (UNIX)

