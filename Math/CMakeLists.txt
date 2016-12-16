# Source files
APPEND_SOURCE_FILES(SOURCES)
APPEND_SOURCE_FILES(SOURCES InterpolationAlgorithms)
APPEND_SOURCE_FILES(SOURCES LinAlg)

IF(OGS_FEM_PETSC)
	APPEND_SOURCE_FILES(SOURCES PETSC)
ENDIF()

# Create the library
add_library( Math ${SOURCES} )

set_target_properties(Math PROPERTIES LINKER_LANGUAGE CXX)

include_directories(
	${CMAKE_SOURCE_DIR}/Math
	${CMAKE_SOURCE_DIR}/Math/LinAlg
	${CMAKE_SOURCE_DIR}/Math/PETSC
	${CMAKE_SOURCE_DIR}/Base
)

target_link_libraries (Math
	Base
)

