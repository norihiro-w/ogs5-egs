APPEND_SOURCE_FILES(SOURCES)
APPEND_SOURCE_FILES(SOURCES FileIO)

if(NOT OGS_FEM_PETSC)
	LIST(REMOVE_ITEM SOURCES msh_lib_ext.cpp)
endif()

include_directories(
	${CMAKE_SOURCE_DIR}/Base
	${CMAKE_SOURCE_DIR}/GEO
	${CMAKE_SOURCE_DIR}/Math
	${CMAKE_SOURCE_DIR}/MSH
	${CMAKE_SOURCE_DIR}/MSH/FileIO
)

# Create the library
add_library(MSH ${SOURCES})

target_link_libraries( MSH
	GEO
)

set_target_properties(MSH PROPERTIES LINKER_LANGUAGE CXX)
