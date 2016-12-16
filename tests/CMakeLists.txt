IF( MSVC ) # VS2012 doesn't support correctly the tuples yet
	ADD_DEFINITIONS( /D _VARIADIC_MAX=10 )
ENDIF()

# Add tests which do not need testdata here:
SET ( SOURCES
	testrunner.cpp
	testBase.cpp
)

INCLUDE_DIRECTORIES(
	${CMAKE_SOURCE_DIR}/ThirdParty/gtest
	${CMAKE_SOURCE_DIR}/Base
	${CMAKE_SOURCE_DIR}/FileIO
	${CMAKE_SOURCE_DIR}/GEO
	${CMAKE_SOURCE_DIR}/Math
	${CMAKE_SOURCE_DIR}/FEM
	${CMAKE_SOURCE_DIR}/MSH
)

ADD_EXECUTABLE (testrunner ${SOURCES})
TARGET_LINK_LIBRARIES(testrunner
	GTest
	Base
	FileIO
	GEO
	Math
	FEM
	MSH
)

IF (MSVC)
	SET(testrunnerExe ${EXECUTABLE_OUTPUT_PATH}/Release/testrunner)
ELSE ()
	SET(testrunnerExe ${EXECUTABLE_OUTPUT_PATH}/testrunner)
	IF("${CMAKE_SYSTEM_NAME}" MATCHES "Linux")
		TARGET_LINK_LIBRARIES( testrunner pthread )
	ENDIF() # "${CMAKE_SYSTEM_NAME}" MATCHES "Linux"
ENDIF() # MSVC

ADD_GOOGLE_TESTS ( ${testrunnerExe} ${SOURCES})

