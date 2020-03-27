# Find the HepPID include directory and library.
#
# This module defines the `HepPID` imported target that encodes all
# necessary information in its target properties.

list(APPEND CMAKE_PREFIX_PATH $ENV{HepPID_DIR})

find_library(
  HepPID_LIBRARY
  NAMES HepPID libHepPI.so
  DOC "The HepPID library")
find_path(
  HepPID_INCLUDE_DIR
  NAMES HepPID/Version.hh
  DOC "The HepPID include directory")

find_package_handle_standard_args(
  HepPID
  REQUIRED_VARS HepPID_LIBRARY HepPID_INCLUDE_DIR)

add_library(HepPID SHARED IMPORTED)
set_property(TARGET HepPID PROPERTY IMPORTED_LOCATION ${HepPID_LIBRARY})
set_property(TARGET HepPID PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${HepPID_INCLUDE_DIR})

mark_as_advanced(HepPID_FOUND HepPID_INCLUDE_DIR HepPID_LIBRARY)
