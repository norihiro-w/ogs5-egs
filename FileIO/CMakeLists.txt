# Source files
SET( SOURCES
	Writer.cpp
	OGSIOVer4.cpp
	readNonBlankLineFromInputStream.cpp
	FEMIO/GeoIO.cpp
	FEMIO/ProcessIO.cpp
	FEMIO/BoundaryConditionIO.cpp
	MeshIO/LegacyVtkInterface.cpp
	MeshIO/vtk.cpp
	MeshIO/OGSMeshIO.cpp
)

# Header files
SET( HEADERS
	Writer.h
	OGSIOVer4.h
	readNonBlankLineFromInputStream.h
	FEMIO/GeoIO.h
	FEMIO/ProcessIO.h
	FEMIO/BoundaryConditionIO.h
	MathIO/CRSIO.h
	MeshIO/LegacyVtkInterface.h
	MeshIO/vtk.h
	MeshIO/OGSMeshIO.h
)

# Create the library
ADD_LIBRARY( FileIO STATIC
	${SOURCES}
	${HEADERS}
)

INCLUDE_DIRECTORIES(
	${CMAKE_SOURCE_DIR}/FileIO
	${CMAKE_SOURCE_DIR}/FileIO/FEMIO
	${CMAKE_SOURCE_DIR}/FileIO/MathIO
	${CMAKE_SOURCE_DIR}/FileIO/MeshIO
	${CMAKE_SOURCE_DIR}/Base
	${CMAKE_SOURCE_DIR}/FEM
	${CMAKE_SOURCE_DIR}/GEO
	${CMAKE_SOURCE_DIR}/MathLib
	${CMAKE_SOURCE_DIR}/MSH
)

TARGET_LINK_LIBRARIES( FileIO
	Base
	#	MathLib
	GEO
	FEM
)

