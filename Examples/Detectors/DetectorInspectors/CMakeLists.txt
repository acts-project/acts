if (ACTS_BUILD_PLUGIN_ACTSVG AND ACTS_BUILD_PLUGIN_JSON)
  add_library(ActsExamplesDetectoInspectors SHARED src/SurfaceIndexing.cpp)
  target_include_directories(
    ActsExamplesDetectoInspectors
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
  target_link_libraries(
    ActsExamplesDetectoInspectors
    PUBLIC ActsCore ActsPluginJson ActsPluginActSVG)
  install(
    TARGETS ActsExamplesDetectoInspectors
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()

  