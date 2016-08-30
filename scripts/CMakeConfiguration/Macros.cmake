# check 64 bit
if( CMAKE_SIZEOF_VOID_P EQUAL 4 )
	set( HAVE_64_BIT 0 )
	set( BITS 32 )
else( CMAKE_SIZEOF_VOID_P EQUAL 4 )
	set( HAVE_64_BIT 1 )
	add_definitions(-DHAVE_64_BIT)
	set( BITS 64)
endif( CMAKE_SIZEOF_VOID_P EQUAL 4 )


# Visual Studio detection
if (${CMAKE_GENERATOR} STREQUAL "Visual Studio 8 2005"
	OR ${CMAKE_GENERATOR} STREQUAL "Visual Studio 9 2008"
	OR ${CMAKE_GENERATOR} STREQUAL "Visual Studio 10"
	OR ${CMAKE_GENERATOR} STREQUAL "NMake Makefiles")

	set (VS32 TRUE)
	set (VS64 FALSE)
	message (STATUS "Generator: Visual Studio 32 Bit")

endif (${CMAKE_GENERATOR} STREQUAL "Visual Studio 8 2005"
	OR ${CMAKE_GENERATOR} STREQUAL "Visual Studio 9 2008"
	OR ${CMAKE_GENERATOR} STREQUAL "Visual Studio 10"
	OR ${CMAKE_GENERATOR} STREQUAL "NMake Makefiles")

if (${CMAKE_GENERATOR} STREQUAL "Visual Studio 8 2005 Win64"
	OR ${CMAKE_GENERATOR} STREQUAL "Visual Studio 9 2008 Win64"
	OR ${CMAKE_GENERATOR} STREQUAL "Visual Studio 10 Win64")

	set (VS32 FALSE)
	set (VS64 TRUE)
	message (STATUS "Generator: Visual Studio 64 Bit")

endif (${CMAKE_GENERATOR} STREQUAL "Visual Studio 8 2005 Win64"
	OR ${CMAKE_GENERATOR} STREQUAL "Visual Studio 9 2008 Win64"
	OR ${CMAKE_GENERATOR} STREQUAL "Visual Studio 10 Win64")

# Convert environment variables
if (NOT $ENV{LIBRARIES_DIR} STREQUAL "")
	STRING(REGEX REPLACE "\\\\" "/" LIBRARIES_DIR $ENV{LIBRARIES_DIR})
endif (NOT $ENV{LIBRARIES_DIR} STREQUAL "")

MACRO(COPY_FILE_IF_CHANGED in_file out_file target)
	IF(${in_file} IS_NEWER_THAN ${out_file})
		#message("Copying file: ${in_file} to: ${out_file}")
		ADD_CUSTOM_COMMAND (
				TARGET     ${target}
				POST_BUILD
				COMMAND    ${CMAKE_COMMAND}
				ARGS       -E copy ${in_file} ${out_file}
		)
		ENDIF(${in_file} IS_NEWER_THAN ${out_file})
ENDMACRO(COPY_FILE_IF_CHANGED)

MACRO(COPY_FILE_INTO_DIRECTORY_IF_CHANGED in_file out_dir target)
		GET_FILENAME_COMPONENT(file_name ${in_file} NAME)
		COPY_FILE_IF_CHANGED(${in_file} ${out_dir}/${file_name}
${target})
ENDMACRO(COPY_FILE_INTO_DIRECTORY_IF_CHANGED)

#Copies all the files from in_file_list into the out_dir.
# sub-trees are ignored (files are stored in same out_dir)
MACRO(COPY_FILES_INTO_DIRECTORY_IF_CHANGED in_file_list out_dir target)
	FOREACH(in_file ${in_file_list})
				COPY_FILE_INTO_DIRECTORY_IF_CHANGED(${in_file}
${out_dir} ${target})
		ENDFOREACH(in_file)
ENDMACRO(COPY_FILES_INTO_DIRECTORY_IF_CHANGED)

MACRO(COPY_FILE_INTO_EXECUTABLE_DIRECTORY in_file target)
	IF (WIN32)
		COPY_FILE_INTO_DIRECTORY_IF_CHANGED(
			${CMAKE_CURRENT_SOURCE_DIR}/${in_file}
			${EXECUTABLE_OUTPUT_PATH}/Debug
			${target}
		)
		COPY_FILE_INTO_DIRECTORY_IF_CHANGED(
			"${CMAKE_CURRENT_SOURCE_DIR}/${in_file}"
			${EXECUTABLE_OUTPUT_PATH}/Release
			${target}
		)
	ELSE (WIN32)
		COPY_FILE_INTO_DIRECTORY_IF_CHANGED(
			${CMAKE_CURRENT_SOURCE_DIR}/${in_file}
			${EXECUTABLE_OUTPUT_PATH}
			${target}
		)
	ENDIF (WIN32)
ENDMACRO(COPY_FILE_INTO_EXECUTABLE_DIRECTORY)

# Adds a benchmark run.
# authorName Your short name
# benchmarkName Relative path in benchmarks directory
# ogsConfiguration E.g. "OGS_FEM"
# numProcesses Number of processes for mpirun, Please set to 1 for non-MPI benchmarks
# Additional arguments add output files to compare
FUNCTION (ADD_BENCHMARK authorName benchmarkName ogsConfiguration numProcesses)

  SET (CONFIG_MATCH FALSE)
  IF (${ogsConfiguration} STREQUAL "OGS_FEM" AND OGS_FEM)
	SET (CONFIG_MATCH TRUE)
  ENDIF (${ogsConfiguration} STREQUAL "OGS_FEM" AND OGS_FEM)
  IF (${ogsConfiguration} STREQUAL "OGS_FEM_SP" AND OGS_FEM_SP)
	SET (CONFIG_MATCH TRUE)
  ENDIF (${ogsConfiguration} STREQUAL "OGS_FEM_SP" AND OGS_FEM_SP)
  IF (${ogsConfiguration} STREQUAL "OGS_FEM_GEMS" AND OGS_FEM_GEMS)
	SET (CONFIG_MATCH TRUE)
  ENDIF (${ogsConfiguration} STREQUAL "OGS_FEM_GEMS" AND OGS_FEM_GEMS)
  IF (${ogsConfiguration} STREQUAL "OGS_FEM_BRNS" AND OGS_FEM_BRNS)
	SET (CONFIG_MATCH TRUE)
  ENDIF (${ogsConfiguration} STREQUAL "OGS_FEM_BRNS" AND OGS_FEM_BRNS)
  IF (${ogsConfiguration} STREQUAL "OGS_FEM_PQC" AND OGS_FEM_PQC)
	SET (CONFIG_MATCH TRUE)
  ENDIF (${ogsConfiguration} STREQUAL "OGS_FEM_PQC" AND OGS_FEM_PQC)
  IF (UNIX) # Only supported on Linux
	IF (${ogsConfiguration} STREQUAL "OGS_FEM_LIS" AND OGS_FEM_LIS)
	  SET (CONFIG_MATCH TRUE)
	ENDIF (${ogsConfiguration} STREQUAL "OGS_FEM_LIS" AND OGS_FEM_LIS)
	IF (${ogsConfiguration} STREQUAL "OGS_FEM_MKL" AND OGS_FEM_MKL)
	  SET (CONFIG_MATCH TRUE)
	ENDIF (${ogsConfiguration} STREQUAL "OGS_FEM_MKL" AND OGS_FEM_MKL)
	IF (${ogsConfiguration} STREQUAL "OGS_FEM_MPI" AND OGS_FEM_MPI)
	  SET (CONFIG_MATCH TRUE)
	ENDIF (${ogsConfiguration} STREQUAL "OGS_FEM_MPI" AND OGS_FEM_MPI)
	IF (${ogsConfiguration} STREQUAL "OGS_FEM_PETSC" AND OGS_FEM_PETSC)
		SET (CONFIG_MATCH TRUE)
	ENDIF ()
  ENDIF (UNIX)

  IF (CONFIG_MATCH)
	IF (WIN32)
	  SET (ogsExe ${EXECUTABLE_OUTPUT_PATH}/Release/ogs)
	ELSE (WIN32)
	  SET (ogsExe ${EXECUTABLE_OUTPUT_PATH}/ogs)
	ENDIF(WIN32)

	# Set timeout
	STRING(REGEX MATCH "EXCEEDING" benchmarkExceeding ${benchmarkName})
	IF(benchmarkExceeding)
		SET(THIS_BENCHMARK_TIMEOUT ${EXCEEDING_BENCHMARK_TIMEOUT})
	ELSE()
		SET(THIS_BENCHMARK_TIMEOUT ${BENCHMARK_TIMEOUT})
	ENDIF()

	STRING (REGEX MATCH "[^/]+$" benchmarkStrippedName ${benchmarkName})
	STRING (LENGTH ${benchmarkName} benchmarkNameLength)
	STRING (LENGTH ${benchmarkStrippedName} benchmarkStrippedNameLength)
	MATH (EXPR substringLength ${benchmarkNameLength}-${benchmarkStrippedNameLength})
	STRING (SUBSTRING ${benchmarkName} 0 ${substringLength} benchmarkDir)
	STRING (REPLACE "/" "_" benchmarkNameUnderscore ${benchmarkName})
	STRING (REPLACE "_LONG_" "_" benchmarkNameUnderscore ${benchmarkNameUnderscore})
	STRING (REPLACE "_EXCEEDING_" "_" benchmarkNameUnderscore ${benchmarkNameUnderscore})

	# Delete output files on testing
	FOREACH (entry ${ARGN})
		SET (FILES_TO_DELETE "${FILES_TO_DELETE} ${entry}")
	ENDFOREACH (entry ${ARGN})
	IF(OGS_PROFILE)
		SET (FILES_TO_DELETE "${FILES_TO_DELETE} gmon.out")
	ENDIF()

	find_program(CSVDIFF_TOOL_PATH csvdiff)
	#find_program(NUMDIFF_TOOL_PATH numdiff)
	#message(STATUS "NUMDIFF_TOOL_PATH=${NUMDIFF_TOOL_PATH}")
	FOREACH (entry ${ARGN})
		LIST (APPEND OUTPUT_FILES ${entry})
	ENDFOREACH (entry ${ARGN})
	string (REPLACE ";" " " OUTPUT_FILES "${OUTPUT_FILES}")

	# Adds a benchmark run. This calls AddTest.cmake to execute several steps.
	ADD_TEST (
		${authorName}_BENCHMARK_${benchmarkName}
		${CMAKE_COMMAND}
		-DPROJECT_SOURCE_DIR=${PROJECT_SOURCE_DIR}
		-DEXECUTABLE_OUTPUT_PATH=${EXECUTABLE_OUTPUT_PATH}
		-DbenchmarkStrippedName=${benchmarkStrippedName}
		-DbenchmarkDir=${benchmarkDir}
		-DFILES_TO_DELETE=${FILES_TO_DELETE}
		-DOGS_PROFILE=${OGS_PROFILE}
		-DOGS_OUTPUT_PROFILE=${OGS_OUTPUT_PROFILE}
		-DGPROF_PATH=${GPROF_PATH}
		-DDOT_TOOL_PATH=${DOT_TOOL_PATH}
		-DBENCHMARK_TIMEOUT=${THIS_BENCHMARK_TIMEOUT}
		-DOGS_FEM_CONFIG=${ogsConfiguration}
		-DNUM_PROCESSES=${numProcesses}
		-DBENCHMARK_DIR_FOUND=${BENCHMARK_DIR_FOUND}
		-P ${PROJECT_SOURCE_DIR}/scripts/CMakeConfiguration/AddBenchmark.cmake
	)

	ADD_TEST (
		${authorName}_FILECOMPARE_${benchmarkName}
		${CMAKE_COMMAND}
		-DPROJECT_SOURCE_DIR=${PROJECT_SOURCE_DIR}
		-DEXECUTABLE_OUTPUT_PATH=${EXECUTABLE_OUTPUT_PATH}
		-DbenchmarkStrippedName=${benchmarkStrippedName}
		-DbenchmarkDir=${benchmarkDir}
		-DFILES_TO_DELETE=${FILES_TO_DELETE}
		-DOGS_PROFILE=${OGS_PROFILE}
		-DOGS_OUTPUT_PROFILE=${OGS_OUTPUT_PROFILE}
		-DGPROF_PATH=${GPROF_PATH}
		-DDOT_TOOL_PATH=${DOT_TOOL_PATH}
		-DBENCHMARK_TIMEOUT=${THIS_BENCHMARK_TIMEOUT}
		-DOGS_FEM_CONFIG=${ogsConfiguration}
		-DNUM_PROCESSES=${numProcesses}
		-DBENCHMARK_DIR_FOUND=${BENCHMARK_DIR_FOUND}
		-DNUMDIFF_TOOL_PATH=${CSVDIFF_TOOL_PATH}
		-DOUTPUT_FILES=${OUTPUT_FILES}
		-DBENCHMARK_REF_DIR=${BENCHMARK_REF_DIR_FOUND}
		-P ${PROJECT_SOURCE_DIR}/scripts/CMakeConfiguration/CompareBenchmark.cmake
	)

	# compare file differences with python script
#	IF (PYTHONINTERP_FOUND)
#		FILE (REMOVE ${BENCHMARK_DIR_FOUND}/results/temp/temp_${benchmarkNameUnderscore}.txt)
#		FOREACH (entry ${ARGN})
#			FILE (APPEND ${BENCHMARK_DIR_FOUND}/results/temp/temp_${benchmarkNameUnderscore}.txt "${entry}\n")
#		ENDFOREACH (entry ${ARGN})
#		ADD_TEST (
#			${authorName}_FILECOMPARE_${benchmarkName}
#			${CMAKE_COMMAND} -E chdir ${BENCHMARK_DIR_FOUND}/results
#			${PYTHON_EXECUTABLE}
#			${PROJECT_SOURCE_DIR}/scripts/compare.py
#			temp/temp_${benchmarkNameUnderscore}.txt
#			${BENCHMARK_REF_DIR_FOUND}
#			${authorName}_${benchmarkNameUnderscore}.html
#			../
#		)
#	ENDIF (PYTHONINTERP_FOUND)

  # copy benchmark output files to reference directory
  IF (COPY_BENCHMARKS_TO_REF)
	FOREACH (entry ${ARGN})
	  CONFIGURE_FILE( ${BENCHMARK_DIR_FOUND}/${entry} ${BENCHMARK_REF_DIR_FOUND}/${entry} COPYONLY)
	ENDFOREACH (entry ${ARGN})
  ENDIF (COPY_BENCHMARKS_TO_REF)

  ENDIF (CONFIG_MATCH)

ENDFUNCTION (ADD_BENCHMARK authorName benchmarkName ogsConfiguration filesToCompare numProcesses)

# Checks if a valid ogs configuration is given
FUNCTION(CHECK_CONFIG)

	SET(configs
		"${OGS_FEM_MKL}"
		"${OGS_FEM_LIS}"
		"${OGS_FEM_PETSC}")

	SET(counter 0)

	FOREACH(config ${configs})
		IF (config)
			MATH(EXPR counter "${counter} + 1")
		ENDIF ()
	ENDFOREACH()
	IF (counter EQUAL 0)
		MESSAGE(STATUS "No configuration specified. Assuming default configuration...")
		SET(OGS_FEM ON)
	ENDIF (counter EQUAL 0)

	IF (counter GREATER 1)
		MESSAGE(FATAL_ERROR "Error: More than one OGS configuration given. Please use only one of the following configurations:
			OGS_FEM_MKL
			OGS_FEM_LIS
			OGS_FEM_PETSC")
	ENDIF (counter GREATER 1)

ENDFUNCTION()

# Creates one ctest for each googletest found in source files passed as arguments
# number two onwards. Argument one specifies the testrunner executable.
MACRO(ADD_GOOGLE_TESTS executable)
	FOREACH ( source ${ARGN} )
		FILE(READ "${source}" contents)
		STRING(REGEX MATCHALL "TEST_?F?\\(([A-Za-z_0-9 ,]+)\\)" found_tests ${contents})
		FOREACH(hit ${found_tests})
			STRING(REGEX REPLACE ".*\\(([A-Za-z_0-9]+)[, ]*([A-Za-z_0-9]+)\\).*" "\\1.\\2" test_name ${hit})
			ADD_TEST(${test_name} ${executable}  --gtest_output=xml --gtest_filter=${test_name} ${MI3CTestingDir})
			# message ("Adding test: ${test_name}")
		ENDFOREACH(hit)
	ENDFOREACH()
ENDMACRO()

# copies the model files to the build dir and adds them as targets so that
# the build files are re-built when the source model files change
MACRO ( UPDATE_MODEL_FILES dirOUT fileLIST )
	GET_FILENAME_COMPONENT( _tdir ${CMAKE_CURRENT_SOURCE_DIR} NAME )
	#message (STATUS "Copying files to ${dirOUT} from ${fileLIST}.\n")
	FOREACH ( _file1 ${${fileLIST}} )
		SET( _file ${CMAKE_CURRENT_SOURCE_DIR}/${_file1} )
		GET_FILENAME_COMPONENT( _fdest ${_file} NAME )
		SET( dest ${dirOUT}/${_fdest} )
		#message( STATUS "Copying ${_file} to ${dest} \n" )

		#message (STATUS "Adding targets ${_tdir}.${_fdest}\n").\n")
		ADD_CUSTOM_TARGET( ${_tdir}.${_fdest}
			${CMAKE_COMMAND} -E copy_if_different
			${_file} ${dest})

		ADD_DEPENDENCIES( testrunner ${_tdir}.${_fdest} )
	ENDFOREACH(_file1 ${${fileLIST}})
ENDMACRO ( UPDATE_MODEL_FILES dirOUT fileLIST)
