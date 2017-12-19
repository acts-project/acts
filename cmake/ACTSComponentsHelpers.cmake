# Provide helper functions to simplify the handling of (optional) components
#
# Components must always be placed in a separate directory which can than
# be added using either of the two provided functions
#
#     add_component(<SUBDIR>)
#     add_component_if(<SUBDIR> ...) # for optional components
#
# In both cases the subdirectory names is used to identify the component.
# The following helper macros are also available
#
#     propagate_components_to_parent()
#     print_components()
#
# and the list of components is stored in the global `_supported_components`
# variable.

set(_supported_components)

# add an optional directory and register its name as a component
function(add_component_if path)
  get_filename_component(_name "${path}" NAME)
  file(RELATIVE_PATH _rel ${PROJECT_SOURCE_DIR} "${CMAKE_CURRENT_SOURCE_DIR}/${path}")
  if(${ARGN})
    add_subdirectory(${path})
    list(APPEND _supported_components "${_name}")
    set(_supported_components "${_supported_components}" PARENT_SCOPE)
    message(STATUS "Enable component '${_name}'")
  else()
    message(STATUS "Disable component '${_name}'")
  endif()
endfunction()

# add a directory and register its name as a component
macro(add_component path)
  add_component_if(${path} TRUE)
endmacro()

# propagate the list of components to the parent scope
macro(propagate_components_to_parent)
  set(_supported_components "${_supported_components}" PARENT_SCOPE)
endmacro()

# print a list of components that
macro(print_components)
  message(STATUS "The following components are being build:")
  foreach(comp ${_supported_components})
    message (STATUS "  ${comp}")
  endforeach()
endmacro()
