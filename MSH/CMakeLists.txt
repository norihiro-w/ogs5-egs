set(SOURCES
	msh_core.h
	msh_edge.h
	msh_elem.h
	msh_faces.h
	msh_lib.h
	msh_mesh.h
	msh_node.h
	MSHEnums.h
	MeshNodesAlongPolyline.h
	msh_core.cpp
	msh_edge.cpp
	msh_elem.cpp
	msh_faces.cpp
	msh_lib.cpp
	msh_mesh.cpp
	msh_node.cpp
	MSHEnums.cpp
	MeshNodesAlongPolyline.cpp
)


if(OGS_FEM_PETSC)
	set(SOURCES ${SOURCES} msh_lib_ext.cpp)
endif()

include_directories(
	${CMAKE_SOURCE_DIR}/MSH
	${CMAKE_SOURCE_DIR}/Base
	${CMAKE_SOURCE_DIR}/GEO
	${CMAKE_SOURCE_DIR}/MathLib
	${CMAKE_SOURCE_DIR}/FileIO
	${CMAKE_SOURCE_DIR}/FEM
)

# Create the library
add_library(MSH  ${SOURCES})

target_link_libraries( MSH
	FEM
	GEO
)
