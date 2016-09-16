# - Config file for the ACTS package
#
# It include the cmake targets of all defined components. Additionally, the 
# following variables are set:
# - ACTS_INCLUDE_DIR ... directory containing the ACTS header files
# - ACTS_LIBRARY_DIR ... directory containing the ACTS libraries
# - compiler flags used during the compilation of ACTS:
#  + ACTS_CXX_FLAGS
#  + ACTS_CXX_FLAGS_DEBUG
#  + ACTS_CXX_FLAGS_MINSIZEREL
#  + ACTS_CXX_FLAGS_RELEASE
#  + ACTS_CXX_FLAGS_RELWITHDEBINFO
# - linker flags used for ACTS:
#  + ACTS_EXE_LINKER_FLAGS_DEBUG
#  + ACTS_SHARED_LINKER_FLAGS_DEBUG

# - Init CMakePackageConfigHelpers
@PACKAGE_INIT@

set(_supported_components @_supported_components@)

# print information about found version
if (NOT (${ACTS_FIND_QUIETLY}))
  message (STATUS "found ACTS version ${ACTS_VERSION}")
  message (STATUS "supported components are:")
  foreach (_scomp ${_supported_components})
    message (STATUS "   - ${_scomp}")
  endforeach ()
endif()

# check for all required components
foreach (_comp ${ACTS_FIND_COMPONENTS})
  # check if this component is supported
  if (NOT ";${_supported_components};" MATCHES ";${_comp};")
    # not supported, but required -> fail
    if (${ACTS_FIND_REQUIRED_${_comp}})
      set (ACTS_FOUND False)
      message (FATAL_ERROR "required component \"${_comp}\" not found")
    # not supported and optional -> skip
    else ()
      list (REMOVE_ITEM ACTS_FIND_COMPONENTS ${_comp})
      if (NOT (${ACTS_FIND_QUIETLY}))
      	message (STATUS "optional component \"${_comp}\" not found")
      endif ()
    endif ()
  endif ()
endforeach ()

if (NOT (${ACTS_FIND_QUIETLY}))
  message (STATUS "loading components:")
endif ()
foreach (_comp ${ACTS_FIND_COMPONENTS})
  if (NOT (${ACTS_FIND_QUIETLY}))
    message (STATUS "   - ${_comp}")
  endif ()
  # - Include the targets file to create the imported targets that a client can
  # link to (libraries) or execute (programs)
  include ("${CMAKE_CURRENT_LIST_DIR}/ACTS${_comp}Targets.cmake")
endforeach ()

# - Create relocatable paths
# NOTE: Do not strictly need paths as all usage requirements are encoded in
# the imported targets created later.
set_and_check (ACTS_INCLUDE_DIR "@PACKAGE_CMAKE_INSTALL_INCLUDEDIR@")
set_and_check (ACTS_LIBRARY_DIR "@PACKAGE_CMAKE_INSTALL_LIBDIR@")

# set ACTS compilation and linker flags
set (ACTS_CXX_FLAGS "@ACTS_CXX_FLAGS@")
set (ACTS_CXX_FLAGS_DEBUG "@ACTS_CXX_FLAGS_DEBUG@")
set (ACTS_CXX_FLAGS_MINSIZEREL "@ACTS_CXX_FLAGS_MINSIZEREL@")
set (ACTS_CXX_FLAGS_RELEASE "@ACTS_CXX_FLAGS_RELEASE@")
set (ACTS_CXX_FLAGS_RELWITHDEBINFO "@ACTS_CXX_FLAGS_RELWITHDEBINFO@")
set (ACTS_EXE_LINKER_FLAGS_DEBUG "@ACTS_EXE_LINKER_FLAGS_DEBUG@")
set (ACTS_SHARED_LINKER_FLAGS_DEBUG "@ACTS_SHARED_LINKER_FLAGS_DEBUG@")