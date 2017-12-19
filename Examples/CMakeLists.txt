# define library target
add_library (ACTSCustomLogger SHARED src/CustomDefaultLogger.cpp)
target_include_directories (ACTSCustomLogger PRIVATE ${CMAKE_SOURCE_DIR}/Core/include)

# set installation directories
install (TARGETS ACTSCustomLogger
         EXPORT ACTSExamplesTargets
         LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
         COMPONENT Examples)
