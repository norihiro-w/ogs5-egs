# Source files
APPEND_SOURCE_FILES(SOURCES)
APPEND_SOURCE_FILES(SOURCES InterpolationAlgorithms)
APPEND_SOURCE_FILES(SOURCES LinAlg)
IF(OGS_FEM_PETSC)
    APPEND_SOURCE_FILES(SOURCES PETSC)
ENDIF()

# Create the library
add_library( MathLib ${SOURCES} )

set_target_properties(MathLib PROPERTIES LINKER_LANGUAGE CXX)

include_directories(
	${CMAKE_SOURCE_DIR}/MathLib
	${CMAKE_SOURCE_DIR}/MathLib/LinAlg
	${CMAKE_SOURCE_DIR}/MathLib/PETSC
	${CMAKE_SOURCE_DIR}/Base
)

target_link_libraries (	MathLib
	Base
)

