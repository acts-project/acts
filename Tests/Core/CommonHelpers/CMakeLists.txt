add_library(ActsTestsCommonHelpers INTERFACE)
target_include_directories(ActsTestsCommonHelpers INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/> $<INSTALL_INTERFACE:include>)
target_include_directories(ActsTestsCommonHelpers INTERFACE ActsCore)

install(TARGETS ActsTestsCommonHelpers LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(DIRECTORY include/Acts DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
