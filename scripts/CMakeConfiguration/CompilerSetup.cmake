INCLUDE(ResetConfigurations)        # To Debug, Release, RelWithDbgInfo
INCLUDE(SetDefaultBuildType)
INCLUDE(DisableCompilerFlag)
SET_DEFAULT_BUILD_TYPE(Release)
INCLUDE(MSVCMultipleProcessCompile) # /MP Switch for VS

# Set compiler helper variables

if(${CMAKE_CXX_COMPILER_ID} MATCHES "Clang")
	set(COMPILER_IS_CLANG TRUE CACHE INTERNAL "")
elseif(${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
	set(COMPILER_IS_GCC TRUE CACHE INTERNAL "")
elseif(${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
	set(COMPILER_IS_INTEL TRUE CACHE INTERNAL "")
elseif(${CMAKE_CXX_COMPILER_ID} STREQUAL "MSVC")
	set(COMPILER_IS_MSVC TRUE CACHE INTERNAL "")
endif() # CMAKE_CXX_COMPILER_ID

# Better Clang warning suppression, see http://www.openwalnut.org/issues/230
if(NOT COMPILER_IS_MSVC)
	set( CMAKE_INCLUDE_SYSTEM_FLAG_CXX "-isystem" CACHE STRING "" FORCE )
endif()

find_program(CCACHE_FOUND ccache)
if(CCACHE_FOUND)
	message(STATUS "CCACHE found")
	set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ccache)
	set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK ccache)
endif(CCACHE_FOUND)

IF (WIN32)
	## For Visual Studio compiler
	IF (MSVC)
		ADD_DEFINITIONS(-D_CRT_SECURE_NO_WARNINGS -D_CRT_NONSTDC_NO_WARNINGS
			-D_CRT_XNONSTDC_NO_WARNINGS)
		# Sets warning level 3 and ignores some warnings
		SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W3 /wd4290 /wd4267")
		# Allow big object files generated for template heavy code
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /bigobj")
		SET(GCC OFF)

		DisableCompilerFlag(DEBUG /RTC1)

		# Set $PATH to Visual Studio bin directory. Needed for finding dumpbin.exe
		IF (MSVC80)
			SET(ENV{PATH} "$ENV{PATH};$ENV{VS80COMNTOOLS}..\\..\\VC\\bin")
		ENDIF ()
		IF (MSVC90)
			SET(ENV{PATH} "$ENV{PATH};$ENV{VS90COMNTOOLS}..\\..\\VC\\bin")
		ENDIF ()
		IF (MSVC10)
			SET(ENV{PATH} "$ENV{PATH};$ENV{VS100COMNTOOLS}..\\..\\VC\\bin")
		ENDIF ()

	ELSE (MSVC)
#FOR CYGWIN.  25.02.2010. WW
		MESSAGE (STATUS "Might be GCC under cygwin.")
		SET(GCC ON)
#		MESSAGE (FATAL_ERROR "Aborting: On Windows only the Visual Studio compiler is supported!")
	ENDIF (MSVC)
ENDIF (WIN32)

### For GNU C/CXX. WW
IF(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_COMPILER_IS_GNUCC)
	SET(GCC ON)
	IF( NOT CMAKE_BUILD_TYPE STREQUAL "Debug" )
		SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -march=native -DNDEBUG")
	ENDIF()
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
	#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pedantic-errors") # disable GCC extensions
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-nonansi-builtins") # Disable built-in declarations of functions that are not mandated by ANSI/ISO C
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wextra")
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wwrite-strings")
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unknown-pragmas -Wno-unused-result")
	#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wshadow")
	#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wconversion")
	#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated")
	IF (NOT (GCC_VERSION VERSION_LESS 4.8) ) 
	  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-local-typedefs") # suppress warnings in Eigen
	ENDIF()
  
	# would be cool: -Woverloaded-virtual, would be overkill: -Weffc++
	ADD_DEFINITIONS(-DGCC)

	IF (OGS_PROFILE)
		IF( NOT CMAKE_BUILD_TYPE STREQUAL "Release" )
			MESSAGE(STATUS "When using profiling you should set CMAKE_BUILD_TYPE to Release.")
		ENDIF()
		SET(PROFILE_FLAGS "-pg -fno-omit-frame-pointer -O2 -DNDEBUG -fno-inline-functions -fno-inline-functions-called-once -fno-optimize-sibling-calls")
		SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${PROFILE_FLAGS}")
	ENDIF (OGS_PROFILE)

	IF(OGS_BUILD_INFO)
		# Get compiler version info, only first line, remove linebreaks
		EXECUTE_PROCESS(COMMAND ${CMAKE_CXX_COMPILER} --version
			COMMAND sed -n 1p
			OUTPUT_VARIABLE GCC_VERSION)
		STRING(REPLACE "\n" "" GCC_VERSION ${GCC_VERSION})
	ENDIF() # OGS_BUILD_INFO
	
ENDIF() # CMAKE_COMPILER_IS_GNUCXX OR CMAKE_COMPILER_IS_GNUCC

IF(COMPILER_IS_CLANG)
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -Wno-self-assign -Wno-deprecated-register")
ENDIF()

IF(OGS_COVERAGE)
	INCLUDE(CodeCoverage)
ENDIF()
