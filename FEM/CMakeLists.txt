# Source files

APPEND_SOURCE_FILES(SOURCES)

IF(NOT OGS_FEM_PETSC)
	LIST(REMOVE_ITEM SOURCES rf_pcs1.cpp fct_mpi.h fct_mpi.cpp)
ENDIF()

IF(OGS_FEM_PETSC)
	LIST(REMOVE_ITEM SOURCES equation_class.h equation_class.cpp)
ENDIF()

IF(OGS_ONLY_TH)
	LIST(REMOVE_ITEM SOURCES rf_random_walk.h rf_random_walk.cpp)
ENDIF( )

APPEND_SOURCE_FILES(SOURCES FileIO)


include_directories(
	${CMAKE_SOURCE_DIR}/FEM
	${CMAKE_SOURCE_DIR}/FEM/FileIO
	${CMAKE_SOURCE_DIR}/Base
	${CMAKE_SOURCE_DIR}/Math
	${CMAKE_SOURCE_DIR}/GEO
	${CMAKE_SOURCE_DIR}/GEO/FileIO
	${CMAKE_SOURCE_DIR}/MSH
	${CMAKE_SOURCE_DIR}/MSH/FileIO
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
	GEO
	MSH
)

if (LIS)
	target_link_libraries( FEM ${LIS_LIBRARIES} )
	if (WIN32)
		set_target_properties( FEM PROPERTIES LINK_FLAGS "/NODEFAULTLIB:libcmt" )
	endif()
endif()

if (MKL)
	target_link_libraries( FEM ${MKL_LIBRARIES} )
endif ()

if (PARALLEL_USE_MPI)
	target_link_libraries( FEM ${MPI_LIBRARIES} )
endif ()
