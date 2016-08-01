# Source files
GET_SOURCE_FILES(SOURCES_Base)
set(SOURCES ${SOURCES_Base})

# Create the library
add_library( Base ${SOURCES} )

set_target_properties(Base PROPERTIES LINKER_LANGUAGE CXX)

include_directories(
	${CMAKE_SOURCE_DIR}/Base
)

