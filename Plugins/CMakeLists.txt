add_subdirectory(DD4hepPlugins)
add_subdirectory(TGeoPlugins)
# propagate list of supported components
set (_supported_components "${_supported_components}" PARENT_SCOPE)

