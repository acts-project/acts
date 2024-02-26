add_library(
  ActsPluginDetray SHARED#not sure about "shared" library type sos
  src/Detray.cpp)
  #include/*.hpp

target_include_directories(
  ActsPluginDetray
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/>
    $<INSTALL_INTERFACE:include>)
    
target_link_libraries(
  ActsPluginDetray 
  PUBLIC 
    ActsCore
    detray::core)

install(
  TARGETS ActsPluginDetray
  EXPORT ActsPluginDetrayTargets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
  
install(
  DIRECTORY include/Acts
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
