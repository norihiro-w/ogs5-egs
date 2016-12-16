# Source files

APPEND_SOURCE_FILES(SOURCES)

IF(OGS_FEM_PETSC)
	LIST(REMOVE_ITEM SOURCES equation_class.h equation_class.cpp)
ELSE()
	LIST(REMOVE_ITEM SOURCES fct_mpi.h fct_mpi.cpp)
ENDIF()

IF(OGS_ONLY_TH)
	LIST(REMOVE_ITEM SOURCES rf_random_walk.h rf_random_walk.cpp)
ENDIF( )

include_directories(
	${CMAKE_SOURCE_DIR}/Base
	${CMAKE_SOURCE_DIR}/GEO
	${CMAKE_SOURCE_DIR}/FEM
	${CMAKE_SOURCE_DIR}/FileIO
	${CMAKE_SOURCE_DIR}/FileIO/FEMIO
	${CMAKE_SOURCE_DIR}/FileIO/MeshIO
	${CMAKE_SOURCE_DIR}/MathLib
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
add_library( FEM STATIC ${SOURCES} )

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

