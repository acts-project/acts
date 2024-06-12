add_library(
  ActsExamplesTracccCommon SHARED
  src/TracccChainAlgorithmBase.cpp
)

target_include_directories(
  ActsExamplesTracccCommon
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

target_link_libraries(
  ActsExamplesTracccCommon
  PUBLIC
    ActsPluginTraccc
    ActsExamplesFramework
    ActsExamplesDigitization
)

install(
  TARGETS ActsExamplesTracccCommon
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
