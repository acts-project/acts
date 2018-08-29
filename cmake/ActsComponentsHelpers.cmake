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
  set(_active ${ARGV1})
  set(_explicit_name ${ARGV2})
  get_filename_component(_dirname "${path}" NAME)
  if("${_explicit_name}" STREQUAL "")
    set(_name ${_dirname})
  else()
    set(_name ${_explicit_name})
  endif()
  file(RELATIVE_PATH _rel ${PROJECT_SOURCE_DIR} "${CMAKE_CURRENT_SOURCE_DIR}/${path}")
  if(${_active})
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
  set(_explicit_name ${ARGV1})
  if("${_explicit_name}" STREQUAL "")
    add_component_if(${path} TRUE)
  else()
    add_component_if(${path} TRUE ${_explicit_name})
  endif()
endmacro()

# propagate the list of components to the parent scope
macro(propagate_components_to_parent)
  set(_supported_components "${_supported_components}" PARENT_SCOPE)
endmacro()

# add an optional subdirectory that is **not** registered as a component
function(add_subdirectory_if path)
  file(RELATIVE_PATH _rel ${PROJECT_SOURCE_DIR} "${CMAKE_CURRENT_SOURCE_DIR}/${path}")
  if(${ARGN})
    add_subdirectory(${path})
  else()
    message(STATUS "Ignore '${_rel}' subdirectory")
  endif()
endfunction()
