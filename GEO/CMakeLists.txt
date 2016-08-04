# Source files
APPEND_SOURCE_FILES(SOURCES)
APPEND_SOURCE_FILES(SOURCES FileIO)

# Create the library
add_library(GEO ${SOURCES})

include_directories(
	${CMAKE_SOURCE_DIR}/GEO
	${CMAKE_SOURCE_DIR}/GEO/FileIO
	${CMAKE_SOURCE_DIR}/Base
	${CMAKE_SOURCE_DIR}/MathLib
	${CMAKE_SOURCE_DIR}/MathLib/LinAlg
)


target_link_libraries ( GEO
	Base
	MathLib
)

