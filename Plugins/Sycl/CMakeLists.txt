add_library(
  ActsPluginSycl SHARED
  # header files
  include/Acts/Plugins/Sycl/Seeding/detail/Types.hpp
  include/Acts/Plugins/Sycl/Seeding/CreateSeedsForGroupSycl.hpp
  include/Acts/Plugins/Sycl/Seeding/DeviceExperimentCuts.hpp
  include/Acts/Plugins/Sycl/Seeding/Seedfinder.hpp
  include/Acts/Plugins/Sycl/Seeding/Seedfinder.ipp
  include/Acts/Plugins/Sycl/Utilities/DeviceSelector.hpp
  include/Acts/Plugins/Sycl/Utilities/ListPlatforms.hpp
  include/Acts/Plugins/Sycl/Utilities/QueueWrapper.hpp
  # implementation files
  src/Seeding/CreateSeedsForGroupSycl.cpp
  src/Utilities/ListPlatforms.cpp
  src/Utilities/DeviceSelector.cpp
  src/Utilities/QueueWrapper.cpp
)

target_include_directories(
  ActsPluginSycl
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
  
target_link_libraries(
  ActsPluginSycl
  PUBLIC ActsCore)

acts_target_setup_sycl(ActsPluginSycl DEPENDENCY PRIVATE)

install(
  TARGETS ActsPluginSycl
  EXPORT ActsPluginSyclTargets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
  
install(
  DIRECTORY include/Acts
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
  