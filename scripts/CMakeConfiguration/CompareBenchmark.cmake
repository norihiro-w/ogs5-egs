
# file comparison with numdiff or cmake
if(NUMDIFF_TOOL_PATH)
	set(TESTER_COMMAND ${NUMDIFF_TOOL_PATH})
	set(TESTER_ARGS --absolute-tolerance=1e-5 --relative-tolerance=1e-4 "-s'\n, '")
	message("Using numdiff")
else()
	set(TESTER_COMMAND ${CMAKE_COMMAND} -E compare_files)
	set(TESTER_ARGS "")
	message("Using cmake diff")
endif()
set(TESTER_COMMAND ${TESTER_COMMAND} )

set(SCRIPT_EXIT_CODE 0)
set(BENCHMARK_OUTPUT_DIRECTORY "${BENCHMARK_DIR_FOUND}/results")
set(NUMDIFF_OUTPUT_FILE ${BENCHMARK_OUTPUT_DIRECTORY}/${benchmarkStrippedName}.numdiff)
execute_process(COMMAND ${CMAKE_COMMAND} -E remove ${NUMDIFF_OUTPUT_FILE})
message(STATUS "OUTPUT_FILES: ${OUTPUT_FILES}")
separate_arguments(OUTPUT_FILES) # reconstruct list
foreach(OUTPUT_FILE ${OUTPUT_FILES})
	execute_process(COMMAND ${TESTER_COMMAND} ${TESTER_ARGS}
		${BENCHMARK_REF_DIR}/${OUTPUT_FILE}
		${BENCHMARK_DIR_FOUND}/${OUTPUT_FILE}
		RESULT_VARIABLE EXIT_CODE
		OUTPUT_QUIET
	)

	if(EXIT_CODE GREATER 0)
		if(NUMDIFF_TOOL_PATH)
			execute_process(COMMAND ${TESTER_COMMAND} ${TESTER_ARGS} -E -S
				${BENCHMARK_REF_DIR}/${OUTPUT_FILE}
				${BENCHMARK_DIR_FOUND}/${OUTPUT_FILE}
				ERROR_VARIABLE NUMDIFF_ERROR
				OUTPUT_VARIABLE NUMDIFF_OUT
			)
			file(APPEND ${NUMDIFF_OUTPUT_FILE} "### Numdiff error ###\n${NUMDIFF_ERROR}\n\n")
			file(APPEND ${NUMDIFF_OUTPUT_FILE} "### Numdiff output ###\n${NUMDIFF_OUT}\n\n")
		endif()
		set(SCRIPT_EXIT_CODE 1)
		message(STATUS "Benchmark file compare of ${OUTPUT_FILE} failed.")
	endif()
endforeach()

if(SCRIPT_EXIT_CODE GREATER 0)
	if(NUMDIFF_TOOL_PATH)
		message(FATAL_ERROR "Benchmark ${benchmarkDir}/${benchmarkStrippedName} failed.\n See ${NUMDIFF_OUTPUT_FILE} for details.")
	else()
		message(FATAL_ERROR "Benchmark ${benchmarkDir}/${benchmarkStrippedName} failed.")
	endif()
endif()
