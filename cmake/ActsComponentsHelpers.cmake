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
function(add_component_if path name)
  file(RELATIVE_PATH _rel ${PROJECT_SOURCE_DIR} "${CMAKE_CURRENT_SOURCE_DIR}/${path}")
  if(${ARGN})
    add_subdirectory(${path})
    list(APPEND _supported_components "${name}")
    set(_supported_components "${_supported_components}" PARENT_SCOPE)
    message(STATUS "Enable component '${name}'")
  else()
    message(STATUS "Disable component '${name}'")
  endif()
endfunction()

# add a directory and register its name as a component
macro(add_component path name)
    add_component_if(${path} ${name} TRUE)
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
