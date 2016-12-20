############################
### Find OGS directories ###
############################

#UNSET(BENCHMARK_DIR_FOUND CACHE)
IF(NOT BENCHMARK_DIR_FOUND)
	FIND_PATH (BENCHMARK_DIR_FOUND copy.py HINTS ${BENCHMARK_DIR})
ELSE()
	FIND_PATH (BENCHMARK_DIR_FOUND copy.py ${PROJECT_SOURCE_DIR}/../benchmarks)
ENDIF()

IF(NOT BENCHMARK_REF_DIR_FOUND)
	SET(BENCHMARK_REF_DIR_FOUND  "${BENCHMARK_REF_DIR}/")
ELSE()
	SET(BENCHMARK_REF_DIR_FOUND  ${BENCHMARK_DIR_FOUND}/../benchmarks/)
ENDIF()

get_filename_component(AbsPathBenchRefDur ${BENCHMARK_REF_DIR_FOUND} ABSOLUTE)
SET(BENCHMARK_REF_DIR_FOUND ${AbsPathBenchRefDur})

MESSAGE(STATUS "Benchmarks directory")
MESSAGE(STATUS "- intput files dir: ${BENCHMARK_DIR_FOUND}")
MESSAGE(STATUS "- reference files dir: ${BENCHMARK_REF_DIR_FOUND}")


######################
### Find libraries ###
######################

## pthread ##
SET ( CMAKE_THREAD_PREFER_PTHREAD ON CACHE BOOL "" )
FIND_PACKAGE( Threads )
IF ( CMAKE_USE_PTHREADS_INIT AND NOT HAVE_PTHREADS)
	SET (HAVE_PTHREADS TRUE CACHE BOOL "Is PThreads found.")
	MESSAGE (STATUS "pthread library found." )
ENDIF ()
IF(HAVE_PTHREADS)
  ADD_DEFINITIONS(-DHAVE_PTHREADS)
ENDIF()
MARK_AS_ADVANCED(CMAKE_THREAD_PREFER_PTHREAD)


IF(MKL)
	# Find MKLlib
	FIND_PACKAGE( MKL REQUIRED )
	INCLUDE_DIRECTORIES (${MKL_INCLUDE_DIR})
ENDIF(MKL)

IF(LIS)
	# Find LISlib
	FIND_PACKAGE( LIS REQUIRED )
	set (NEW_EQS ON)
	add_definitions(
		-o3
		-DIPMGEMPLUGIN
	)
ENDIF(LIS)

# Find OpenMP
IF(PARALLEL_USE_OPENMP)
	FIND_PACKAGE( OpenMP )
	IF (OPENMP_FOUND)
		MESSAGE(STATUS "Found OpenMP")
		SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}" )
		SET( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}" )
		SET( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lgomp" )
	ELSE()
		MESSAGE(STATUS "Could not find OpenMP")
	ENDIF()
ENDIF(PARALLEL_USE_OPENMP)

IF(PARALLEL_USE_MPI)
	IF (WIN32)
#		MESSAGE (FATAL_ERROR "Aborting: MPI is only supported under UNIX/LINUX!")
		#ADD_DEFINITIONS(-DMPICH_IGNORE_CXX_SEEK)
		FIND_PACKAGE(MPI REQUIRED)
	ENDIF(WIN32)
	IF(UNIX)

# If there is an mpi compiler find it and interogate (farther below) it for the include
# and lib dirs otherwise we will continue to search from ${_MPI_BASE_DIR}.

		IF( ${CMAKE_MAJOR_VERSION} EQUAL 2 AND ${CMAKE_MINOR_VERSION} LESS 8)
			find_program(MPI_COMPILER
				NAMES mpic++ mpicxx mpiCC mpicc
				HINTS "${_MPI_BASE_DIR}"
				PATH_SUFFIXES bin
				DOC "MPI compiler. Used only to detect MPI compilation flags.")
			IF(MPI_COMPILER)

			MESSAGE (STATUS  "CMake version is less than 2.8, MPI compiler is set directly" )
			mark_as_advanced(MPI_COMPILER)
				SET(CMAKE_C_COMPILER ${MPI_COMPILER})
				SET(CMAKE_CXX_COMPILER ${MPI_COMPILER})
			ENDIF(MPI_COMPILER)
		ELSE( ${CMAKE_MAJOR_VERSION}  EQUAL 2 AND ${CMAKE_MINOR_VERSION} LESS 8)
			FIND_PACKAGE(MPI REQUIRED)
		ENDIF( ${CMAKE_MAJOR_VERSION} EQUAL 2 AND ${CMAKE_MINOR_VERSION} LESS 8)
	ENDIF(UNIX)
ENDIF(PARALLEL_USE_MPI)

######################
### Find tools     ###
######################

# Find Git
FIND_PACKAGE(Git)

# msysGit on Windows
IF(WIN32 AND GIT_FOUND)
	FIND_PACKAGE(MsysGit)
ENDIF() # WIN32 AND GIT_FOUND

# Find gnu profiler gprof
FIND_PROGRAM(GPROF_PATH gprof DOC "GNU profiler gprof")


## Unix tools ##
# Date
FIND_PROGRAM(DATE_TOOL_PATH date PATHS ${MSYSGIT_BIN_DIR})
# Grep
FIND_PROGRAM(GREP_TOOL_PATH grep PATHS ${MSYSGIT_BIN_DIR})
# Unzip
FIND_PROGRAM(UNZIP_TOOL_PATH unzip PATHS ${MSYSGIT_BIN_DIR})

# Hide these variables for the CMake user
MARK_AS_ADVANCED(GPROF_PATH UNZIP_TOOL_PATH)
