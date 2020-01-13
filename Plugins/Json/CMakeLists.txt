if (ACTS_USE_EXTERNAL_NLOHMANN_JSON)
  message(STATUS "Using external install of nlohmann::json")
  find_package(nlohmann_json 3.2.0 REQUIRED)

else()
  message(STATUS "Will build nlohmann::json as part of this build")

  set(_source "https://github.com/nlohmann/json/archive/v3.7.3.tar.gz")
  set(_md5 "846bbc611ce9ecd7d45d6554679245e1")

  include(FetchContent)
  FetchContent_Declare(
    nlohmann_json
    URL ${_source}
    URL_MD5 ${_md5}
  )

  set(JSON_BuildTests OFF CACHE INTERNAL "")
  FetchContent_MakeAvailable(nlohmann_json)

endif()

file(GLOB_RECURSE src_files "src/*.cpp" "include/*.hpp")

add_library(
  ActsJsonPlugin SHARED ${src_files})

target_include_directories(
  ActsJsonPlugin
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/>
    $<INSTALL_INTERFACE:include>)

target_link_libraries(
  ActsJsonPlugin PUBLIC ActsCore nlohmann_json::nlohmann_json)

install(
  TARGETS ActsJsonPlugin
  EXPORT ActsJsonPluginTargets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(
  DIRECTORY include/Acts
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
