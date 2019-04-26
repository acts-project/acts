add_library(ActsJsonPlugin INTERFACE)

# Reference for the json header file, obtained from
# https://github.com/nlohmann/json/releases
set(_json_dep_header_sha256_ref d2eeb25d2e95bffeb08ebb7704cdffd2e8fca7113eba9a0b38d60a5c391ea09a)

# Check if hash of file on disk matches the stated hash of the release
file(READ "include/Acts/Plugins/Json/lib/json.hpp" _json_dep_header)
string(SHA256 _json_dep_header_sha256 "${_json_dep_header}")

if (NOT ${_json_dep_header_sha256} STREQUAL ${_json_dep_header_sha256_ref})
  message(FATAL_ERROR "json.hpp hash does not match reference!")
endif()

target_include_directories(ActsJsonPlugin SYSTEM INTERFACE ${Boost_INCLUDE_DIRS} ${EIGEN_INCLUDE_DIRS})
target_include_directories(ActsJsonPlugin INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/> $<INSTALL_INTERFACE:include>)

install(
  TARGETS ActsJsonPlugin
  EXPORT ActsJsonPluginTargets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(
  DIRECTORY include/Acts
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
