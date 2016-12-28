# - Try to find Paralution
# Once done, this will define
#
#  Paralution_FOUND
#  Paralution_INCLUDE_DIR
#  Paralution_LIBRARY

SET(Paralution_DIR
    "${Paralution_DIR}"
    CACHE
    PATH
    "Directory to search for Paralution library")
IF (NOT Paralution_DIR)
    MESSAGE (FATAL_ERROR "Please set Paralution_DIR root diretory path" )
ENDIF()

find_path(Paralution_INCLUDE
	NAMES paralution.hpp
	HINTS ${Paralution_DIR} PATH_SUFFIXES include)

find_library(Paralution_LIBRARY
	NAMES libparalution
	HINTS ${Paralution_DIR} PATH_SUFFIXES lib)

if(NOT Paralution_LIBRARY STREQUAL "Paralution_LIBRARY-NOTFOUND" AND NOT Paralution_INCLUDE STREQUAL "Paralution_INCLUDE-NOTFOUND")
	set(Paralution_FOUND ON)
else()
	message (STATUS "Warning: Paralution not found!")
endif ()

MESSAGE (STATUS "Paralution lib: ${Paralution_LIBRARY}")

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Paralution DEFAULT_MSG Paralution_INCLUDE Paralution_LIBRARY)
mark_as_advanced(Paralution_INCLUDE Paralution_LIBRARY)

