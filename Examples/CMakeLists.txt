add_library (ActsCustomLogger SHARED src/CustomDefaultLogger.cpp)
target_include_directories (ActsCustomLogger PRIVATE ${PROJECT_SOURCE_DIR}/Core/include)

install (TARGETS ActsCustomLogger
         EXPORT ActsExamplesTargets
         LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
         COMPONENT Examples)
