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

IF(DEFINED EXAMPLEDATA_DIR)
	FIND_PATH (EXAMPLEDATA_DIR_FOUND points.gli ${EXAMPLEDATA_DIR})
ELSE()
	FIND_PATH (EXAMPLEDATA_DIR_FOUND points.gli ${PROJECT_SOURCE_DIR}/../ExampleData)
ENDIF()

FIND_PATH (OGS_LIBS_DIR_FOUND geotiff.lib ${PROJECT_SOURCE_DIR}/../Libs/libgeotiff)

IF(DEFINED TESTDATA_DIR)
	FIND_PATH(TESTDATA_DIR_FOUND testdata.dummy ${TESTDATA_DIR})
ELSE()
	FIND_PATH(TESTDATA_DIR_FOUND testdata.dummy ${PROJECT_SOURCE_DIR}/../testdata)
ENDIF()

# Find precompiled libraries (for BRNS GEMS LIS)
FIND_PATH (OGS_PRECOMPILED_LIBS_DIR_FOUND GEMS3_rl.lib ${PROJECT_SOURCE_DIR}/../Libs/precompiled)
IF (OGS_PRECOMPILED_LIBS_DIR_FOUND)
#    MESSAGE(STATUS "OGS_PRECOMPILED_LIBS_DIR_FOUND=${OGS_PRECOMPILED_LIBS_DIR_FOUND}")
#	INCLUDE_DIRECTORIES (${PROJECT_SOURCE_DIR}/../Libs/precompiled)
#	LINK_DIRECTORIES (${PROJECT_SOURCE_DIR}/../Libs/precompiled)
	INCLUDE_DIRECTORIES (${OGS_PRECOMPILED_LIBS_DIR_FOUND})
	LINK_DIRECTORIES (${OGS_PRECOMPILED_LIBS_DIR_FOUND})
ELSE (OGS_PRECOMPILED_LIBS_DIR_FOUND)
	IF (WIN32)
		IF (OGS_FEM_BRNS OR OGS_FEM_GEMS OR OGS_FEM_CHEMAPP)
			MESSAGE (FATAL_ERROR "Precompiled libraries not found! Make sure to also check out the trunk/Libs directory beneath your sources directory.")
		ENDIF (OGS_FEM_BRNS OR OGS_FEM_GEMS OR OGS_FEM_CHEMAPP)
	ENDIF (WIN32)
ENDIF (OGS_PRECOMPILED_LIBS_DIR_FOUND)

######################
### Find libraries ###
######################
#FIND_PACKAGE (PythonInterp)
#FIND_PACKAGE( Shapelib )
#IF(Shapelib_FOUND)
#	ADD_DEFINITIONS(-DShapelib_FOUND)
#ENDIF() # Shapelib_FOUND

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

## boost (see FindBoost.cmake for more options) ##
IF (UNIX AND GCC AND OGS_FEM_GEMS)
       set(Boost_USE_STATIC_LIBS    OFF)
ELSE()
       set(Boost_USE_STATIC_LIBS    ON)
ENDIF()

set(Boost_USE_MULTITHREADED      ON)
set(Boost_USE_STATIC_RUNTIME    OFF)

IF(NOT OGS_FEM_GEMS)
	FIND_PACKAGE( Boost COMPONENTS filesystem system regex)
ELSE()
	# Boost with threads is required for GEMS
	FIND_PACKAGE( Boost COMPONENTS filesystem system regex thread )
ENDIF()


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

IF(OGS_FEM_Paralution)
	FIND_PACKAGE( Paralution REQUIRED )
ENDIF()

# Find OpenMP
IF(PARALLEL_USE_OPENMP)
	FIND_PACKAGE( OpenMP REQUIRED )
	SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}" )
	SET( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}" )
	SET( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lgomp" )
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

# Find Python interpreter
find_package (PythonInterp)
message(STATUS "PythonInterp - ${PYTHONINTERP_FOUND}")

# Find Subversion
FIND_PACKAGE(Subversion)

# Find Git
FIND_PACKAGE(Git)

# msysGit on Windows
IF(WIN32 AND GIT_FOUND)
	FIND_PACKAGE(MsysGit)
ENDIF() # WIN32 AND GIT_FOUND

# Find dot tool from graphviz
FIND_PROGRAM(DOT_TOOL_PATH dot DOC "Dot tool from graphviz")

# Find doxygen
FIND_PACKAGE(Doxygen)

# Find gnu profiler gprof
FIND_PROGRAM(GPROF_PATH gprof DOC "GNU profiler gprof")

FIND_PACKAGE(cppcheck)

# Find Exuberant ctags or BBEdit for code completion
FIND_PROGRAM(CTAGS_TOOL_PATH ctags DOC "Exuberant ctags")
FIND_PROGRAM(BBEDIT_TOOL_PATH bbedit DOC "BBEdit Editor")
IF(BBEDIT_TOOL_PATH)
	ADD_CUSTOM_TARGET(ctags
		bbedit --maketags
		WORKING_DIRECTORY ${CMAKE_SOURCES_DIR}
		COMMENT "Creating tags..." VERBATIM
	)
	ADD_CUSTOM_COMMAND(TARGET ctags POST_BUILD
		COMMAND mv -f tags ../tags
		WORKING_DIRECTORY ${CMAKE_SOURCES_DIR}
		COMMENT "Moving tags..." VERBATIM
	)
ELSE()
	IF(CTAGS_TOOL_PATH)
		ADD_CUSTOM_TARGET(ctags
			ctags -R --fields=+iamS -f ${CMAKE_SOURCES_DIR}/../tags
			WORKING_DIRECTORY ${CMAKE_SOURCES_DIR}
			COMMENT "Creating tags..." VERBATIM
		)
	ENDIF()
ENDIF()

## Unix tools ##
# Date
FIND_PROGRAM(DATE_TOOL_PATH date PATHS ${MSYSGIT_BIN_DIR})
# Grep
FIND_PROGRAM(GREP_TOOL_PATH grep PATHS ${MSYSGIT_BIN_DIR})
# Unzip
FIND_PROGRAM(UNZIP_TOOL_PATH unzip PATHS ${MSYSGIT_BIN_DIR})

# Hide these variables for the CMake user
MARK_AS_ADVANCED(DOT_TOOL_PATH GPROF_PATH CTAGS_TOOL_PATH BBEDIT_TOOL_PATH
	UNZIP_TOOL_PATH
)
########################
### Find other stuff ###
########################

# Check if on Jenkins
IF(NOT $ENV{JENKINS_URL} STREQUAL "")
	SET(JENKINS_URL $ENV{JENKINS_URL})
	SET(JENKINS_JOB_NAME $ENV{JOB_NAME})
ENDIF()
