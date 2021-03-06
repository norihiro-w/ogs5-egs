#### OGS-5 Project ####

# Specify minimum CMake version
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.3)

# Project name
PROJECT(OGS5-EGS)

###############
### Options ###
###############

# Configurations
OPTION(OGS_FEM_PETSC "Use PETSc parallel solver" OFF)
OPTION(OGS_FEM_LIS "Library of Iterative Solvers for Linear Systems configuration" OFF)
OPTION(OGS_FEM_MKL "Math kernel library configuration" OFF)
OPTION(OGS_FEM_Paralution "Paralution libarary for linear solve" OFF)

OPTION(OGS_ONLY_TH "Optimize OGS for TH simulations" OFF)

# Set default configuration when no other config is given
IF (NOT OGS_FEM_LIS AND NOT OGS_FEM_MKL AND NOT OGS_FEM_PETSC)
	MESSAGE (STATUS "No configuration specified. Default confuguration is used.")
        SET (OGS_FEM_LIS ON)
ENDIF ()

# Check if a valid OGS configuration is given (see Macros.cmake)
#CHECK_CONFIG()

# Additional options
OPTION(OGS_CMAKE_DEBUG "Show additional cmake debug information" OFF)
OPTION(OGS_BUILD_INFO "Should build information be generated" ON)
OPTION(CMAKE_CMD_ARGS "Build ogs5 storing the CMake command line args (hint: use cmake.ogs.sh" OFF )
IF(GCC)
	IF(GPROF_PATH)
		OPTION(OGS_PROFILE "Enables compiling with flags set for profiling with gprof." OFF)
		IF(OGS_PROFILE AND DOT_TOOL_PATH)
			OPTION(OGS_OUTPUT_PROFILE "Enables pdf output of profiling results." OFF)
		ENDIF() # OGS_PROFILE AND DOT_TOOL_PATH
	ENDIF() # GPROF_PATH

	OPTION_REQUIRES(OGS_COVERAGE "Enables code coverage measurements with gcov/lcov." OGS_BUILD_TESTS)
ENDIF() # GCC
OPTION (OGS_BUILD_TESTS "Enables building of tests." OFF)
OPTION (OGS_DELETE_EDGES_AFTER_INIT "Delete mesh edges after initialization if possible" OFF)

MARK_AS_ADVANCED(FORCE OGS_CMAKE_DEBUG OGS_BUILD_INFO CMAKE_CMD_ARGS)

IF(OGS_FEM_LIS)
	SET(LIS ON CACHE BOOL "Use Library of Iterative Solvers")
	SET(PARALLEL_USE_OPENMP ON CACHE BOOL "Use Library of Iterative Solvers")
	MARK_AS_ADVANCED(LIS, PARALLEL_USE_OPENMP)
ENDIF()

IF(OGS_FEM_MKL)
	SET(PARALLEL_USE_OPENMP ON CACHE BOOL "Use Library of Iterative Solvers")
	SET(LIS ON CACHE BOOL "Use Library of Iterative Solvers")
	SET(MKL ON CACHE BOOL "Use Math Kernel Library")
	MARK_AS_ADVANCED(PARALLEL_USE_OPENMP LIS MKL)
ENDIF(OGS_FEM_MKL)

IF(OGS_FEM_Paralution)
	SET(PARALLEL_USE_OPENMP ON CACHE BOOL "Use Library of Iterative Solvers")
	MARK_AS_ADVANCED(PARALLEL_USE_OPENMP)
ENDIF()


#########################
### Finding libraries ###
#########################

INCLUDE(scripts/cmake/Functions.cmake)
INCLUDE(scripts/CMakeConfiguration/CMakeSetup.cmake)
INCLUDE(scripts/CMakeConfiguration/CompilerSetup.cmake)
INCLUDE (scripts/CMakeConfiguration/Find.cmake)
INCLUDE(scripts/CMakeConfiguration/GeneralProjectSetup.cmake)



###########################################################################
### OGS version information. Adjust these if you release a new version. ###
###########################################################################
SET (OGS_VERSION_MAJOR 5)
SET (OGS_VERSION_MINOR 3)
SET (OGS_VERSION_PATCH 3-egs.1.10)
SET (OGS_RELEASE_PERSONS "NW")
SET (OGS_VERSION "${OGS_VERSION_MAJOR}.${OGS_VERSION_MINOR}.${OGS_VERSION_PATCH}(${OGS_RELEASE_PERSONS})")
SET (OGS_DATE "2017-01-05")

# Check for number of processors
INCLUDE(ProcessorCount)
ProcessorCount(PROCESSOR_COUNT)
IF(PROCESSOR_COUNT EQUAL 0)
	MESSAGE(WARNING "Processor count could not be detected. Setting to one processor.")
	SET(PROCESSOR_COUNT 1)
ELSE()
	MESSAGE(STATUS "Number of processors: ${PROCESSOR_COUNT}")
ENDIF() # PROCESSOR_COUNT

IF(OGS_ONLY_TH)
	ADD_DEFINITIONS(-DOGS_ONLY_TH)
ENDIF()

IF(OGS_FEM_PETSC)
	MESSAGE (STATUS "Configuring for PETSc" )

	SET(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/scripts/cmake/findPETSC)

	FIND_PACKAGE(PETSc REQUIRED QUIET)

	ADD_DEFINITIONS(-DUSE_PETSC)
	ADD_DEFINITIONS(-DPETSC_VERSION_NUMBER=PETSC_VERSION_MAJOR*1000+PETSC_VERSION_MINOR*10)


	INCLUDE_DIRECTORIES( ${PETSC_INCLUDES} )
	FIND_PACKAGE(MPI)
	IF(MPI_FOUND)
		SET(CMAKE_C_COMPILER ${MPI_COMPILER})
		SET(CMAKE_CXX_COMPILER ${MPI_COMPILER})
	ELSE(MPI_FOUND)
		MESSAGE (FATAL_ERROR "Aborting: MPI implementation is not found!")
	ENDIF(MPI_FOUND)
ENDIF(OGS_FEM_PETSC)

IF(OGS_FEM_Paralution)
	set(NEW_EQS ON)
	ADD_DEFINITIONS(-DNEW_EQS -DUSE_PARALUTION)
	INCLUDE_DIRECTORIES (${Paralution_INCLUDE})
ENDIF()

IF(LIS)
	MESSAGE (STATUS "Configuring for FEM command line with LIS" )
	INCLUDE_DIRECTORIES(${LIS_INCLUDE_DIR})
	ADD_DEFINITIONS(-DNEW_EQS -DIPMGEMPLUGIN -DLIS)
	IF(OGS_USE_LONG)
		ADD_DEFINITIONS(-D_LONGLONG)
	ENDIF()
ENDIF()

IF(MKL)
	MESSAGE (STATUS	"Configuring for FEM command line with MKL" )
	INCLUDE_DIRECTORIES(${MKL_INCLUDES})
	ADD_DEFINITIONS(-DNO_ERROR_CONTROL -DNEW_EQS -DMKL)
	IF(OGS_USE_LONG)
		ADD_DEFINITIONS(-DMKL_ILP64)
	ENDIF()
ENDIF()


IF(OGS_USE_LONG)
	MESSAGE (STATUS	"Use long long int" )
	ADD_DEFINITIONS(-DOGS_USE_LONG)
ENDIF()

IF(MSVC)
	ADD_DEFINITIONS(-DMSVC)
ENDIF() # MSVC

IF(OGS_DELETE_EDGES_AFTER_INIT)
	ADD_DEFINITIONS(-DOGS_DELETE_EDGES_AFTER_INIT)
ENDIF()

IF(OGS_SAVE_MEMORY)
	ADD_DEFINITIONS(-DOGS_DELETE_EDGES_AFTER_INIT)
	ADD_DEFINITIONS(-DOGS_SAVE_MEMORY)
ENDIF()

# Add subdirectories with the projects
#ADD_SUBDIRECTORY( ThirdParty )
INCLUDE_DIRECTORIES (SYSTEM ${CMAKE_CURRENT_SOURCE_DIR}/ThirdParty)
IF(NOT USE_EXTERNAL_EIGEN)
	INCLUDE_DIRECTORIES (${CMAKE_CURRENT_SOURCE_DIR}/ThirdParty/eigen3)
ENDIF()

ADD_SUBDIRECTORY( Base )
ADD_SUBDIRECTORY( Math )
ADD_SUBDIRECTORY( GEO )
ADD_SUBDIRECTORY( MSH )
ADD_SUBDIRECTORY( FEM )
ADD_SUBDIRECTORY( OGS )


# Create Configure.h
CONFIGURE_FILE (Base/Configure.h.in ${PROJECT_BINARY_DIR}/Base/Configure.h)
CONFIGURE_FILE (Base/BuildInfo.h.in ${PROJECT_BINARY_DIR}/Base/BuildInfo.h)

IF (BENCHMARK_DIR_FOUND OR GTEST_FOUND)
	MESSAGE(STATUS "ENABLE_TESTING")
	ENABLE_TESTING()
ENDIF (BENCHMARK_DIR_FOUND OR GTEST_FOUND)

IF (OGS_BUILD_TESTS)
	ADD_SUBDIRECTORY(ThirdParty/gtest)
	ADD_SUBDIRECTORY(tests)
ENDIF (OGS_BUILD_TESTS)

IF (OGS_CMAKE_DEBUG)
	INCLUDE(ListAllCMakeVariableValues)
	list_all_cmake_variable_values()
ENDIF (OGS_CMAKE_DEBUG)

IF (BENCHMARK_DIR_FOUND)
	SET(BENCHMARK_TIMEOUT 1800) # in s, 30 minutes timeout on normal benchmarks
	SET(EXCEEDING_BENCHMARK_TIMEOUT 86400) # 1 day timeout on exceeding benchmarks
	STRING (COMPARE EQUAL "$ENV{HOSTNAME}" "dev2.intern.ufz.de" HOST_IS_DEV2)
	IF (${HOST_IS_DEV2})
		MESSAGE (STATUS "On dev2: file comparing enabled")
	ENDIF (${HOST_IS_DEV2})

	MESSAGE(STATUS "Looking for benchmark cases")
	FILE (GLOB BENCHMARK_CONFIGS "${BENCHMARK_DIR_FOUND}/*.cmake")
#	FILE (GLOB BENCHMARK_CONFIGS "${PROJECT_SOURCE_DIR}/../benchmarks/*.cmake")
	FOREACH (BENCHMARK_CONFIG ${BENCHMARK_CONFIGS})
		MESSAGE(STATUS "INCLUDE ${BENCHMARK_CONFIG}")
		INCLUDE ("${BENCHMARK_CONFIG}")
	ENDFOREACH(BENCHMARK_CONFIG)

	UNSET (COPY_BENCHMARKS_TO_REF CACHE)
ENDIF (BENCHMARK_DIR_FOUND)

