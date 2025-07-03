# Find the FastJet includes and libraries.

include(FindPackageHandleStandardArgs)

find_library(FastJet_LIBRARY NAMES FastJet fastjet DOC "The FastJet library")

find_path(
    FastJet_INCLUDE_DIR
    fastjet/version.hh
    DOC "The FastJet include directory"
)

if(${FastJet_INCLUDE_DIR} STREQUAL "FastJet_INCLUDE_DIR-NOTFOUND")
    message(FATAL_ERROR "FastJet include directory not found")
endif()

file(READ "${FastJet_INCLUDE_DIR}/fastjet/config_auto.h" FastJet_VERSION_FILE)
string(
    REGEX MATCH
    "#define[ \t]+FASTJET_PACKAGE_VERSION[ \t]+\"([0-9]+\.[0-9]+\.[0-9]+)\""
    _
    ${FastJet_VERSION_FILE}
)
if(NOT CMAKE_MATCH_1)
    message(FATAL_ERROR "Failed to extract FastJet version from config_auto.h")
endif()

set(FastJet_VERSION ${CMAKE_MATCH_1})

find_package_handle_standard_args(
    FastJet
    REQUIRED_VARS FastJet_LIBRARY FastJet_INCLUDE_DIR
    VERSION_VAR FastJet_VERSION
)

add_library(FastJet SHARED IMPORTED)
add_library(FastJet::FastJet ALIAS FastJet)

set_property(TARGET FastJet PROPERTY IMPORTED_LOCATION ${FastJet_LIBRARY})
set_property(
    TARGET FastJet
    PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${FastJet_INCLUDE_DIR}
)

mark_as_advanced(FastJet_FOUND FastJet_INCLUDE_DIR FastJet_LIBRARY)

if(NOT FastJet_FIND_QUIETLY)
    if(FastJet_FOUND)
        message(STATUS "Found FastJet ${FastJet_VERSION} at ${FastJet_LIBRARY}")
    else()
        message(FATAL_ERROR "FastJet not found")
    endif()
endif()
