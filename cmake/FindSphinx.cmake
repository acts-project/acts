# Find the Sphinx documentation generator.
#
# This module defines the following variables:
#   Sphinx_FOUND - whether Sphinx was found or not
#   Sphinx_EXECUTABLE - Sphinx executable if if was found
#

find_program(
  Sphinx_EXECUTABLE
  NAMES sphinx-build sphinx-build-3
  DOC "Sphinx documentation generator")

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  Sphinx
  REQUIRED_VARS Sphinx_EXECUTABLE)
