# Source files
APPEND_SOURCE_FILES(SOURCES)
APPEND_SOURCE_FILES(SOURCES FileIO)

# Create the library
add_library(GEO ${SOURCES})

include_directories(
	${CMAKE_SOURCE_DIR}/Base
	${CMAKE_SOURCE_DIR}/GEO
	${CMAKE_SOURCE_DIR}/GEO/FileIO
	${CMAKE_SOURCE_DIR}/Math
	${CMAKE_SOURCE_DIR}/Math/LinAlg
)

target_link_libraries (GEO
	Base
	Math
)

set_target_properties(GEO PROPERTIES LINKER_LANGUAGE CXX)
