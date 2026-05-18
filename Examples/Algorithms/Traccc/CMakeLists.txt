acts_add_library(ExamplesTraccc SHARED
  src/TracccChain.cpp
  src/TracccSeqAlg.cpp
  src/ActsMeasToTracccAlg.cpp
  src/ActsSpToTracccAlg.cpp
)

target_include_directories(
  ActsExamplesTraccc
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

target_link_libraries(
  ActsExamplesTraccc
  PUBLIC
    Acts::Core
    Acts::ExamplesFramework
    traccc::core
    detray::core
    detray::io
    traccc::io
    vecmem::core
)

if(ACTS_ENABLE_CUDA)
  target_compile_definitions(
    ActsExamplesTraccc
    PRIVATE
      ACTS_ENABLE_CUDA
  )
  target_link_libraries(
    ActsExamplesTraccc
    PRIVATE
      traccc::cuda
      traccc::device_common
      vecmem::cuda
  )
endif()