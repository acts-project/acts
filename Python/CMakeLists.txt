# Add the utilities as ActsPythonUtilties
add_subdirectory(Utilities)

# Add the core module
add_subdirectory(Core)

# Add the plugins module, individual plugins
# are added in respecting the corresponding
# build flags
add_subdirectory(Plugins)

# Add the examples module
add_subdirectory(Examples)
