add_library(ActsTestsCommon INTERFACE)
target_include_directories(ActsTestsCommon INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/> $<INSTALL_INTERFACE:include>)
target_include_directories(ActsTestsCommon INTERFACE ActsCore)

install(TARGETS ActsTestsCommon LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(DIRECTORY include/Acts DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
