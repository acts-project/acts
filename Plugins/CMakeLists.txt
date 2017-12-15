add_subdirectory (DD4hepPlugin)
add_subdirectory (TGeoPlugin)
add_subdirectory (MaterialPlugin)
# propagate list of supported components
set (_supported_components "${_supported_components}" PARENT_SCOPE)
