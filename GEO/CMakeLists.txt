# Source files
GET_SOURCE_FILES(SOURCES_GEO)
set(SOURCES ${SOURCES_GEO})

# Create the library
add_library(GEO ${SOURCES})

include_directories(
	${CMAKE_SOURCE_DIR}/GEO
	${CMAKE_SOURCE_DIR}/Base
	${CMAKE_SOURCE_DIR}/FileIO
	${CMAKE_SOURCE_DIR}/Math
	${CMAKE_SOURCE_DIR}/Math/LinAlg
	${CMAKE_SOURCE_DIR}/FEM
	${CMAKE_SOURCE_DIR}/MSH
)


target_link_libraries (
	GEO
	Base
	FEM
	Math
)

